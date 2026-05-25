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

    // Precompute per-position log-probability lookup tables for every unique sequence.
    // logp_table[u][i][tb] = log P(u.seq[i] | true_base=tb, mean_qual_at_i)
    // Stored as f32 (half the RAM of f64); seq_ll upcasts to f64 for accumulation.
    let logp_table: Vec<Vec<[f32; 4]>> = uniques
        .par_iter()
        .map(|u| precompute_logp(u, error_model))
        .collect();

    // Initialise: most-abundant unique seq is the first cluster centre
    let mut centres: Vec<Vec<u8>> = vec![uniques[0].seq.clone()];
    let mut assignments: Vec<usize> = vec![0usize; uniques.len()];

    // Track which sequences are already centres; maintained incrementally across
    // all iterations (centres only grow, never shrink).
    let mut centre_set: std::collections::HashSet<&[u8]> =
        std::collections::HashSet::from([uniques[0].seq.as_slice()]);

    let mut prev_ll = f64::NEG_INFINITY;

    for iter in 0..cfg.max_iter {
        // Greedy promotion: process uniques in decreasing count order (already
        // sorted by derep_fastq).  For each candidate we look up its CURRENT
        // best centre — including any centres promoted earlier in this same
        // pass — so that once a true ASV is promoted, its close neighbours
        // (error reads) are tested against that ASV and are not spuriously
        // promoted themselves.
        for (u_idx, u) in uniques.iter().enumerate() {
            if centre_set.contains(u.seq.as_slice()) {
                continue;
            }
            let ci = best_centre(&logp_table[u_idx], &centres);
            let log_prob = seq_ll(&logp_table[u_idx], &centres[ci]);
            #[allow(clippy::cast_precision_loss)]
            let log_lambda = (total_reads as f64).ln() + log_prob;
            if is_significant_log(u64::from(u.count), log_lambda, cfg.omega_a) {
                centre_set.insert(u.seq.as_slice());
                centres.push(u.seq.clone());
            }
        }

        // Re-assign all uniques to their nearest centre (parallel).
        let new_assignments: Vec<usize> = logp_table
            .par_iter()
            .map(|logp| best_centre(logp, &centres))
            .collect();

        // Total log-likelihood (parallel sum).
        let ll: f64 = logp_table
            .par_iter()
            .zip(new_assignments.par_iter())
            .zip(uniques.par_iter())
            .map(|((logp, &ci), u)| {
                seq_ll(logp, &centres[ci.min(centres.len() - 1)]) * f64::from(u.count)
            })
            .sum();

        assignments = new_assignments;

        let delta = (ll - prev_ll).abs();
        let n_centres = centres.len();
        let max_iter = cfg.max_iter;
        log::info!("dada: iter {iter}/{max_iter}, {n_centres} centres, ΔlogL = {delta:.2e}");
        if delta < cfg.tol {
            break;
        }
        prev_ll = ll;
    }

    // Collect ASVs using centre indices to avoid cloning sequences more than once.
    let mut asv_counts: std::collections::HashMap<usize, u32> =
        std::collections::HashMap::new();
    for (u, &ci) in uniques.iter().zip(assignments.iter()) {
        let ci = ci.min(centres.len() - 1);
        *asv_counts.entry(ci).or_insert(0) += u.count;
    }

    let mut asvs: Vec<Asv> = asv_counts
        .into_iter()
        .map(|(ci, abundance)| Asv { sequence: centres[ci].clone(), abundance })
        .collect();
    asvs.sort_unstable_by_key(|a| Reverse(a.abundance));

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

    let mut store = PoolStore::new(500_000)?;
    for (i, sample) in samples.iter().enumerate() {
        store.add_sample(i, sample)?;
    }

    let (pooled_uniques, pool_entries) = store.into_pooled_uniques()?;
    let pooled_asvs = dada(&pooled_uniques, error_model, cfg)?;

    // Build a lookup: sequence → pool entry index
    let mut seq_to_entry: std::collections::HashMap<&[u8], usize> =
        std::collections::HashMap::new();
    for (idx, u) in pooled_uniques.iter().enumerate() {
        seq_to_entry.insert(&u.seq, idx);
    }

    // Re-split ASVs back to per-sample
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

/// Precompute per-position log-probability table for a single unique sequence.
///
/// `result[i][tb]` = log P(u.seq[i] | `true_base=tb`, `mean_qual_at_i`).
/// Stored as `f32` — halves the logp-table RAM cost vs `f64` with no loss of
/// clustering accuracy (scores compared relatively, not absolutely).
fn precompute_logp(u: &UniqueSeq, em: &ErrorModel) -> Vec<[f32; 4]> {
    (0..u.seq.len())
        .map(|i| {
            let ob = base_index(u.seq[i]);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let q = Phred(u.mean_qual(i) as u8);
            #[allow(clippy::cast_possible_truncation)]
            std::array::from_fn(|tb| em.log_p_error(tb as u8, ob, q) as f32)
        })
        .collect()
}

/// Log-likelihood of unique `logp` given `centre` — pure table lookups, no transcendentals.
/// Upcasts f32 entries to f64 before accumulating to preserve sum precision.
fn seq_ll(logp: &[[f32; 4]], centre: &[u8]) -> f64 {
    logp.iter()
        .zip(centre.iter())
        .map(|(lp, &cb)| f64::from(lp[base_index(cb) as usize]))
        .sum()
}

/// Return the index of the centre with the highest log-likelihood for `logp`.
fn best_centre(logp: &[[f32; 4]], centres: &[Vec<u8>]) -> usize {
    centres
        .iter()
        .enumerate()
        .map(|(i, c)| (i, seq_ll(logp, c)))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map_or(0, |(i, _)| i)
}

/// Poisson abundance significance test working entirely in log-space.
///
/// Returns `true` if P(X ≥ count | `Poisson(exp(log_lambda))`) < `omega_a`.
fn is_significant_log(count: u64, log_lambda: f64, omega_a: f64) -> bool {
    if count == 0 {
        return false;
    }
    let lambda = log_lambda.exp();
    if lambda >= 1e-14 {
        let Ok(dist) = Poisson::new(lambda) else {
            return false;
        };
        let p_val: f64 = 1.0 - (0..count).map(|k| dist.pmf(k)).sum::<f64>();
        return p_val < omega_a;
    }
    // Log-space path: P(X >= count) ≈ P(X = count) when lambda is tiny.
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
        let seq = "ACGTACGTACGTACGT";
        let s1 = vec![make_unique(seq, 100)];
        let s2 = vec![make_unique(seq, 50)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig { omega_a: 0.5, ..Default::default() };

        let result = dada_pooled(&[&s1, &s2], &em, &cfg).unwrap();
        assert_eq!(result.len(), 2);
        let total_s1: u32 = result[0].iter().map(|a| a.abundance).sum();
        let total_s2: u32 = result[1].iter().map(|a| a.abundance).sum();
        assert!(total_s1 > 0, "sample 0 should have non-zero abundance");
        assert!(total_s2 > 0, "sample 1 should have non-zero abundance");
        let ratio = f64::from(total_s1) / f64::from(total_s2);
        assert!(ratio > 1.5 && ratio < 2.5, "expected ~2:1 ratio, got {ratio:.2}");
    }

    #[test]
    fn single_cluster_identical_reads() {
        let uniques = vec![make_unique("ACGTACGTACGT", 1000)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig::default();
        let asvs = dada(&uniques, &em, &cfg).unwrap();
        assert_eq!(asvs.len(), 1);
        assert_eq!(asvs[0].sequence, b"ACGTACGTACGT");
        assert_eq!(asvs[0].abundance, 1000);
    }
}
