//! Stage 5 — Core DADA denoising algorithm.
//!
//! Implements the greedy partition + EM refinement described in
//! Callahan et al. 2016 (Nature Methods, Suppl. Note 1).

use crate::{
    Dada2Error,
    derep::UniqueSeq,
    error_model::{ErrorModel, base_index},
    pool::PoolStore,
    Phred,
};
use rayon::prelude::*;
use std::cmp::Reverse;
use statrs::distribution::{Discrete, Poisson};

/// A denoised Amplicon Sequence Variant.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Asv {
    /// The inferred true sequence.
    pub sequence: Vec<u8>,
    /// Total read abundance assigned to this ASV.
    pub abundance: u32,
}

/// Configuration for the DADA algorithm.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DadaConfig {
    /// Abundance p-value threshold: sequences with p < `omega_a` are accepted as new ASVs.
    pub omega_a: f64,
    /// Pool reads from all samples before denoising (pseudo-pooling).
    pub pool: bool,
    /// Maximum EM iterations.
    pub max_iter: usize,
    /// Log-likelihood convergence tolerance.
    pub tol: f64,
    /// RNG seed (unused currently, reserved for future tie-breaking).
    pub seed: u64,
}

impl Default for DadaConfig {
    fn default() -> Self {
        Self {
            omega_a: 1e-40,
            pool: false,
            max_iter: 16,
            tol: 1e-6,
            seed: 42,
        }
    }
}

/// Run the DADA algorithm on a dereplicated sample.
///
/// # Errors
/// Returns [`Dada2Error::InvalidInput`] if `uniques` is empty.
pub fn dada(
    uniques: &[UniqueSeq],
    error_model: &ErrorModel,
    cfg: &DadaConfig,
) -> Result<Vec<Asv>, Dada2Error> {
    if uniques.is_empty() {
        return Err(Dada2Error::InvalidInput("no unique sequences supplied to dada".into()));
    }

    let total_reads: u64 = uniques.iter().map(|u| u64::from(u.count)).sum();

    // Initialise: the most-abundant unique seq is the first cluster centre
    let mut centres: Vec<Vec<u8>> = vec![uniques[0].seq.clone()];
    let mut assignments: Vec<usize> = vec![0usize; uniques.len()];

    let mut prev_ll = f64::NEG_INFINITY;

    for iter in 0..cfg.max_iter {
        // Greedy promotion: process uniques in decreasing count order (already
        // sorted by derep_fastq).  For each candidate we look up its CURRENT
        // best centre — including any centres promoted earlier in this same
        // pass — so that once a true ASV is promoted, its close neighbours
        // (error reads) are tested against that ASV and are not spuriously
        // promoted themselves.
        let mut live_centre_set: std::collections::HashSet<Vec<u8>> =
            centres.iter().cloned().collect();
        for u in uniques.iter() {
            if live_centre_set.contains(&u.seq) {
                continue;
            }
            let ci = best_centre(u, &centres, error_model);
            // Compute log_prob without holding a reference into centres,
            // so that centres.push below is not blocked by the borrow.
            let log_prob = seq_log_likelihood(u, &centres[ci], error_model);
            #[allow(clippy::cast_precision_loss)]
            let log_lambda = (total_reads as f64).ln() + log_prob;
            if is_significant_log(u64::from(u.count), log_lambda, cfg.omega_a) {
                live_centre_set.insert(u.seq.clone());
                centres.push(u.seq.clone());
            }
        }

        // Re-assign with updated centres
        let new_assignments2: Vec<usize> = uniques
            .par_iter()
            .map(|u| best_centre(u, &centres, error_model))
            .collect();

        // Compute total log-likelihood
        let ll: f64 = uniques
            .iter()
            .zip(new_assignments2.iter())
            .map(|(u, &ci)| {
                let centre = &centres[ci.min(centres.len() - 1)];
                let ll_per_read = seq_log_likelihood(u, centre, error_model);
                ll_per_read * f64::from(u.count)
            })
            .sum();

        assignments = new_assignments2;

        let delta = (ll - prev_ll).abs();
        let n_centres = centres.len();
        let max_iter = cfg.max_iter;
        log::info!(
            "dada: iter {iter}/{max_iter}, {n_centres} centres, ΔlogL = {delta:.2e}"
        );
        if delta < cfg.tol {
            break;
        }
        prev_ll = ll;
    }

    // Collect ASVs
    let mut asv_counts: std::collections::HashMap<Vec<u8>, u32> = std::collections::HashMap::new();
    for (u, &ci) in uniques.iter().zip(assignments.iter()) {
        let centre = centres[ci.min(centres.len() - 1)].clone();
        *asv_counts.entry(centre).or_insert(0) += u.count;
    }

    let mut asvs: Vec<Asv> = asv_counts
        .into_iter()
        .map(|(sequence, abundance)| Asv { sequence, abundance })
        .collect();
    asvs.sort_unstable_by_key(|a| std::cmp::Reverse(a.abundance));

    Ok(asvs)
}

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

    // 1. Build a PoolStore and accumulate all samples
    let mut store = PoolStore::new(500_000)?;
    for (i, sample) in samples.iter().enumerate() {
        store.add_sample(i, sample)?;
    }

    // 2. Merge into pooled uniques + provenance entries
    let (pooled_uniques, pool_entries) = store.into_pooled_uniques()?;

    // 3. Run DADA on the merged pool
    let pooled_asvs = dada(&pooled_uniques, error_model, cfg)?;

    // 4. Build a lookup: sequence → pool entry index
    let mut seq_to_entry: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::new();
    for (idx, u) in pooled_uniques.iter().enumerate() {
        seq_to_entry.insert(&u.seq, idx);
    }

    // 5. Re-split ASVs back to per-sample
    let mut per_sample: Vec<std::collections::HashMap<Vec<u8>, u32>> =
        (0..n_samples).map(|_| std::collections::HashMap::new()).collect();

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
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_precision_loss)]
                let alloc = (f64::from(asv.abundance) * frac).round() as u32;
                if alloc > 0 {
                    *per_sample[sample_idx]
                        .entry(asv.sequence.clone())
                        .or_insert(0) += alloc;
                }
            }
        }
    }

    // Convert hashmaps to sorted Vec<Asv>
    let result: Vec<Vec<Asv>> = per_sample
        .into_iter()
        .map(|map| {
            let mut v: Vec<Asv> = map
                .into_iter()
                .map(|(sequence, abundance)| Asv { sequence, abundance })
                .collect();
            v.sort_unstable_by_key(|a| Reverse(a.abundance));
            v
        })
        .collect();

    Ok(result)
}

/// Return the index of the centre with the highest likelihood for `u`.
fn best_centre(u: &UniqueSeq, centres: &[Vec<u8>], em: &ErrorModel) -> usize {
    centres
        .iter()
        .enumerate()
        .map(|(i, c)| (i, seq_log_likelihood(u, c, em)))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map_or(0, |(i, _)| i)
}

/// Log-likelihood of `u.seq` given `centre` under the error model.
fn seq_log_likelihood(u: &UniqueSeq, centre: &[u8], em: &ErrorModel) -> f64 {
    let len = u.seq.len().min(centre.len());
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    (0..len)
        .map(|i| {
            let tb = base_index(centre[i]);
            let ob = base_index(u.seq[i]);
            let q = Phred(u.mean_qual(i) as u8);
            em.p_error(tb, ob, q).max(1e-300).ln()
        })
        .sum()
}

/// Poisson abundance significance test working entirely in log-space.
///
/// Returns `true` if P(X ≥ count | Poisson(exp(log_lambda))) < `omega_a`.
///
/// When `log_lambda` is very negative (sequence far from centre), `exp(log_lambda)`
/// underflows to 0.0 in f64.  We fall back to the log-space approximation:
///   log P(X = count) ≈ count · log_lambda − log(count!)
/// which stays representable even for log_lambda ≈ −1000.
fn is_significant_log(count: u64, log_lambda: f64, omega_a: f64) -> bool {
    if count == 0 {
        return false;
    }
    let lambda = log_lambda.exp();
    if lambda >= 1e-14 {
        let Ok(dist) = Poisson::new(lambda) else {
            return false;
        };
        // P(X >= count) = 1 - CDF(count-1)
        let p_val: f64 = 1.0 - (0..count).map(|k| dist.pmf(k)).sum::<f64>();
        return p_val < omega_a;
    }
    // Log-space path: P(X >= count) ≈ P(X = count) when lambda is tiny.
    // log P(X = count) = count * log_lambda - log(count!)
    #[allow(clippy::cast_precision_loss)]
    let log_p = (count as f64) * log_lambda - log_factorial(count);
    log_p < omega_a.ln()
}

fn log_factorial(n: u64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    if n <= 20 {
        #[allow(clippy::cast_precision_loss)]
        return (1..=n).map(|k| (k as f64).ln()).sum();
    }
    // Stirling approximation
    #[allow(clippy::cast_precision_loss)]
    let nf = n as f64;
    nf * nf.ln() - nf + 0.5 * (2.0 * std::f64::consts::PI * nf).ln()
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
        // 2 samples each with the same true sequence at different counts.
        // Pool → should recover the ASV, and both samples get abundance back.
        let seq = "ACGTACGTACGTACGT";
        let s1 = vec![make_unique(seq, 100)];
        let s2 = vec![make_unique(seq, 50)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig { omega_a: 0.5, ..Default::default() };

        let result = dada_pooled(&[&s1, &s2], &em, &cfg).unwrap();
        assert_eq!(result.len(), 2);
        // Both samples should have at least one ASV
        let total_s1: u32 = result[0].iter().map(|a| a.abundance).sum();
        let total_s2: u32 = result[1].iter().map(|a| a.abundance).sum();
        assert!(total_s1 > 0, "sample 0 should have non-zero abundance");
        assert!(total_s2 > 0, "sample 1 should have non-zero abundance");
        // Ratio should be roughly 2:1
        let ratio = f64::from(total_s1) / f64::from(total_s2);
        assert!(ratio > 1.5 && ratio < 2.5, "expected ~2:1 ratio, got {ratio:.2}");
    }

    #[test]
    fn single_cluster_identical_reads() {
        // All reads are the same → should converge to exactly 1 ASV
        let uniques = vec![make_unique("ACGTACGTACGT", 1000)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig::default();
        let asvs = dada(&uniques, &em, &cfg).unwrap();
        assert_eq!(asvs.len(), 1);
        assert_eq!(asvs[0].sequence, b"ACGTACGTACGT");
        assert_eq!(asvs[0].abundance, 1000);
    }
}
