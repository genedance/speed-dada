//! Cross-sample DADA strategies: pooled (single joint denoise) and
//! pseudo-pooled (two-pass: per-sample, then per-sample with cross-sample
//! priors).
//!
//! Split out of `dada.rs` to keep that file under the 500-line cap.

use crate::{
    dada::{dada, dada_with_priors, Asv, DadaConfig},
    derep::UniqueSeq,
    error_model::ErrorModel,
    pool::PoolStore,
    Dada2Error,
};
use rayon::prelude::*;
use std::cmp::Reverse;

/// Run DADA on multiple samples with cross-sample pooling.
///
/// Sequences from all samples are pooled into a single [`PoolStore`],
/// denoised jointly, then re-split back to per-sample ASV lists using provenance
/// stored during accumulation.
///
/// Returns one `Vec<Asv>` per input sample, in the same order.
///
/// # Errors
/// Returns [`Dada2Error`] if the pool store fails or [`dada`] fails.
pub fn dada_pooled(
    samples: &[&[UniqueSeq]],
    error_model: &ErrorModel,
    cfg: &DadaConfig,
) -> Result<Vec<Vec<Asv>>, Dada2Error> {
    let n_samples = samples.len();

    let mut store = PoolStore::new(500_000)?;
    for (i, sample) in samples.iter().enumerate() {
        store.add_sample(i, sample)?;
    }

    let (pooled_uniques, pool_entries) = store.into_pooled_uniques()?;
    let pooled_asvs = dada(&pooled_uniques, error_model, cfg)?;

    // Build a lookup: sequence → pool entry index.
    let mut seq_to_entry: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::new();
    for (idx, u) in pooled_uniques.iter().enumerate() {
        seq_to_entry.insert(&u.seq, idx);
    }

    // Re-split ASVs back to per-sample by allocating abundance proportional
    // to each sample's contribution to the pooled entry's count.
    let mut per_sample: Vec<std::collections::HashMap<Vec<u8>, u32>> = (0..n_samples)
        .map(|_| std::collections::HashMap::new())
        .collect();

    for asv in &pooled_asvs {
        if let Some(&entry_idx) = seq_to_entry.get(asv.sequence.as_slice()) {
            let entry = &pool_entries[entry_idx];
            #[allow(clippy::cast_precision_loss)]
            let total = f64::from(entry.total_count);
            for &(sample_idx, sample_count) in &entry.per_sample {
                if sample_idx >= n_samples || total == 0.0 {
                    continue;
                }
                #[allow(clippy::cast_precision_loss)]
                let frac = f64::from(sample_count) / total;
                #[allow(
                    clippy::cast_possible_truncation,
                    clippy::cast_sign_loss,
                    clippy::cast_precision_loss
                )]
                let alloc = (f64::from(asv.abundance) * frac).round() as u32;
                if alloc > 0 {
                    *per_sample[sample_idx]
                        .entry(asv.sequence.clone())
                        .or_insert(0) += alloc;
                }
            }
        }
    }

    let result: Vec<Vec<Asv>> = per_sample
        .into_iter()
        .map(|map| {
            let mut v: Vec<Asv> = map
                .into_iter()
                .map(|(sequence, abundance)| Asv {
                    sequence,
                    abundance,
                })
                .collect();
            v.sort_unstable_by_key(|a| Reverse(a.abundance));
            v
        })
        .collect();

    Ok(result)
}

/// Run DADA with pseudo-pooling — Callahan's two-pass cross-sample scheme.
///
/// **Pass 1**: per-sample DADA in parallel (no priors). Collect every ASV
/// detected in any sample → `priors` set.
///
/// **Pass 2**: per-sample DADA in parallel again, but sequences in `priors`
/// auto-promote on first encounter (bypass the Poisson `omega_a` test). This
/// recovers rare-but-real ASVs that occur across multiple samples but would
/// fail the abundance test in any single sample.
///
/// Per-sample passes parallelise across the rayon thread pool — much faster
/// than the single-threaded greedy promotion of [`dada_pooled`].
///
/// # Errors
/// Returns [`Dada2Error`] if any per-sample [`dada`] call fails.
pub fn dada_pseudo(
    samples: &[&[UniqueSeq]],
    error_model: &ErrorModel,
    cfg: &DadaConfig,
) -> Result<Vec<Vec<Asv>>, Dada2Error> {
    // Pass 1 — per-sample, no priors, run in parallel.
    let pass1: Vec<Vec<Asv>> = samples
        .par_iter()
        .map(|s| dada(s, error_model, cfg))
        .collect::<Result<_, _>>()?;

    // Union of every ASV across all samples becomes the prior set.
    let prior_seqs: Vec<Vec<u8>> = {
        let mut set: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
        for asvs in &pass1 {
            for a in asvs {
                set.insert(a.sequence.clone());
            }
        }
        set.into_iter().collect()
    };
    let priors_set: std::collections::HashSet<&[u8]> =
        prior_seqs.iter().map(std::vec::Vec::as_slice).collect();

    // Pass 2 — per-sample with priors, run in parallel.
    let pass2: Vec<Vec<Asv>> = samples
        .par_iter()
        .map(|s| dada_with_priors(s, error_model, cfg, &priors_set))
        .collect::<Result<_, _>>()?;

    Ok(pass2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::derep::UniqueSeq;
    use crate::error_model::ErrorModel;

    fn make_unique(seq: &str, count: u32) -> UniqueSeq {
        let len = seq.len();
        UniqueSeq {
            seq: seq.bytes().collect(),
            count,
            qual_sum: vec![40.0 * f64::from(count); len],
        }
    }

    #[test]
    fn pooled_two_samples_same_sequence() {
        let seq = "ACGTACGTACGTACGT";
        let s1 = vec![make_unique(seq, 100)];
        let s2 = vec![make_unique(seq, 50)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig {
            omega_a: 0.5,
            ..Default::default()
        };

        let result = dada_pooled(&[&s1, &s2], &em, &cfg).unwrap();
        assert_eq!(result.len(), 2);
        let total_s1: u32 = result[0].iter().map(|a| a.abundance).sum();
        let total_s2: u32 = result[1].iter().map(|a| a.abundance).sum();
        assert!(total_s1 > 0, "sample 0 should have non-zero abundance");
        assert!(total_s2 > 0, "sample 1 should have non-zero abundance");
        let ratio = f64::from(total_s1) / f64::from(total_s2);
        assert!(
            ratio > 1.5 && ratio < 2.5,
            "expected ~2:1 ratio, got {ratio:.2}"
        );
    }
}
