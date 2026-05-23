//! Stage 5 — Core DADA denoising algorithm.
//!
//! Implements the greedy partition + EM refinement described in
//! Callahan et al. 2016 (Nature Methods, Suppl. Note 1).

use crate::{
    Dada2Error,
    derep::UniqueSeq,
    error_model::{ErrorModel, base_index},
    Phred,
};
use rayon::prelude::*;
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
    /// Abundance p-value threshold: sequences with p < omega_a are accepted as new ASVs.
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
pub fn run_dada(
    uniques: &[UniqueSeq],
    error_model: &ErrorModel,
    cfg: &DadaConfig,
) -> Result<Vec<Asv>, Dada2Error> {
    if uniques.is_empty() {
        return Err(Dada2Error::InvalidInput("no unique sequences supplied to run_dada".into()));
    }

    let total_reads: u64 = uniques.iter().map(|u| u64::from(u.count)).sum();

    // Initialise: the most-abundant unique seq is the first cluster centre
    let mut centres: Vec<Vec<u8>> = vec![uniques[0].seq.clone()];
    let mut assignments: Vec<usize> = vec![0usize; uniques.len()];

    let mut prev_ll = f64::NEG_INFINITY;

    for _iter in 0..cfg.max_iter {
        // E-step: assign each unique to its best centre
        let new_assignments: Vec<usize> = uniques
            .par_iter()
            .map(|u| best_centre(u, &centres, error_model))
            .collect();

        // Promote unique sequences that are significantly over-abundant
        for (i, u) in uniques.iter().enumerate() {
            if new_assignments[i] < centres.len() {
                continue;
            }
            // Already assigned to an existing centre — check if it should become one
            let centre = &centres[new_assignments[i].min(centres.len() - 1)];
            let err_rate = error_rate(u, centre, error_model);
            let lambda = (total_reads as f64) * err_rate;
            if is_significant(u64::from(u.count), lambda, cfg.omega_a) {
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

        if (ll - prev_ll).abs() < cfg.tol {
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
    asvs.sort_unstable_by(|a, b| b.abundance.cmp(&a.abundance));

    Ok(asvs)
}

/// Return the index of the centre with the highest likelihood for `u`.
fn best_centre(u: &UniqueSeq, centres: &[Vec<u8>], em: &ErrorModel) -> usize {
    centres
        .iter()
        .enumerate()
        .map(|(i, c)| (i, seq_log_likelihood(u, c, em)))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Compute the mean per-base error rate between `u` and `centre`.
fn error_rate(u: &UniqueSeq, centre: &[u8], em: &ErrorModel) -> f64 {
    let len = u.seq.len().min(centre.len());
    if len == 0 {
        return 1.0;
    }
    let rate: f64 = (0..len)
        .map(|i| {
            let tb = base_index(centre[i]);
            let ob = base_index(u.seq[i]);
            let q = Phred(u.mean_qual(i) as u8);
            em.p_error(tb, ob, q)
        })
        .sum::<f64>()
        / len as f64;
    rate.max(1e-300)
}

/// Log-likelihood of `u.seq` given `centre` under the error model.
fn seq_log_likelihood(u: &UniqueSeq, centre: &[u8], em: &ErrorModel) -> f64 {
    let len = u.seq.len().min(centre.len());
    (0..len)
        .map(|i| {
            let tb = base_index(centre[i]);
            let ob = base_index(u.seq[i]);
            let q = Phred(u.mean_qual(i) as u8);
            em.p_error(tb, ob, q).max(1e-300).ln()
        })
        .sum()
}

/// Poisson abundance significance test.
///
/// Returns `true` if P(X ≥ count | Poisson(lambda)) < omega_a.
fn is_significant(count: u64, lambda: f64, omega_a: f64) -> bool {
    if lambda <= 0.0 || count == 0 {
        return false;
    }
    let dist = match Poisson::new(lambda) {
        Ok(d) => d,
        Err(_) => return false,
    };
    // P(X >= count) = 1 - CDF(count-1)
    let p_val: f64 = 1.0
        - (0..count)
            .map(|k| dist.pmf(k))
            .sum::<f64>();
    p_val < omega_a
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
    fn single_cluster_identical_reads() {
        // All reads are the same → should converge to exactly 1 ASV
        let uniques = vec![make_unique("ACGTACGTACGT", 1000)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig::default();
        let asvs = run_dada(&uniques, &em, &cfg).unwrap();
        assert_eq!(asvs.len(), 1);
        assert_eq!(asvs[0].sequence, b"ACGTACGTACGT");
        assert_eq!(asvs[0].abundance, 1000);
    }
}
