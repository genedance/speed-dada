//! Stage 5 — Core DADA denoising algorithm (public API + EM loop).
//!
//! Implements the greedy partition + EM refinement described in
//! Callahan et al. 2016 (Nature Methods, Suppl. Note 1). The scoring
//! primitives live in `dada_scoring.rs`; the cross-sample strategies
//! `dada_pooled` and `dada_pseudo` live in `dada_pool.rs`.

use crate::{
    dada_scoring::{
        best_centre_for_promotion, best_centre_serial_packed, gap_corrected_ll, is_significant_log,
        pack_seq, precompute_logp, seq_ll_packed,
    },
    derep::UniqueSeq,
    error_model::ErrorModel,
    Dada2Error,
};
use rayon::prelude::*;
use std::cmp::Reverse;

// Re-export the multi-sample strategies so callers can still write
// `use speeddada_core::dada::{dada_pooled, dada_pseudo, ...}` — the split
// is an internal refactor, not a public API change.
pub use crate::dada_pool::{dada_pooled, dada_pseudo};

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
    /// Log-probability of a single insertion or deletion in `gap_corrected_ll`.
    ///
    /// The substitution-only scorer mis-scores reads with a sequencer-introduced
    /// indel as having `L - p` substitution errors after the indel position p,
    /// which spuriously promotes the indel artefact as a new ASV. The gap-
    /// corrected scorer adds a single-indel-alignment model with this constant
    /// penalty; an indel-artefact read scores ≈ centre's log-likelihood + this
    /// constant instead of accumulating substitution penalties.
    ///
    /// Default: `ln(1e-3) ≈ -6.91`, matching R dada2's empirical indel rate on
    /// Illumina paired-end data.
    pub gap_log_p: f64,
}

impl Default for DadaConfig {
    fn default() -> Self {
        Self {
            omega_a: 1e-40,
            pool: false,
            max_iter: 16,
            tol: 1e-6,
            seed: 42,
            gap_log_p: -6.907_755_278_982_137, // (1e-3_f64).ln()
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
    let empty = std::collections::HashSet::new();
    dada_with_priors(uniques, error_model, cfg, &empty)
}

/// Variant of [`dada`] that auto-promotes sequences listed in `priors`,
/// bypassing the Poisson abundance test (only requires presence in `uniques`).
///
/// This is the inner used by [`dada_pseudo`] for the second pass; sequences
/// that were detected as ASVs in at least one sample become priors for every
/// sample, allowing rare-but-real ASVs to be recovered consistently across
/// samples even when their per-sample abundance is below `omega_a`.
///
/// # Errors
/// Returns [`Dada2Error::InvalidInput`] if `uniques` is empty.
pub fn dada_with_priors(
    uniques: &[UniqueSeq],
    error_model: &ErrorModel,
    cfg: &DadaConfig,
    priors: &std::collections::HashSet<&[u8]>,
) -> Result<Vec<Asv>, Dada2Error> {
    if uniques.is_empty() {
        return Err(Dada2Error::InvalidInput(
            "no unique sequences supplied to dada".into(),
        ));
    }

    let total_reads: u64 = uniques.iter().map(|u| u64::from(u.count)).sum();

    // Precompute per-position log-probability lookup tables for every unique sequence.
    // logp_table[u][i][tb] = log P(u.seq[i] | true_base=tb, mean_qual_at_i)
    // Stored as f32 (half the RAM of f64); seq_ll upcasts to f64 for accumulation.
    let logp_table: Vec<Vec<[f32; 4]>> = uniques
        .par_iter()
        .map(|u| precompute_logp(u, error_model))
        .collect();

    // Centres are stored as indices into `uniques` rather than cloned byte
    // sequences — the i-th cluster centre's sequence is `uniques[centre_idx[i]].seq`.
    // We also maintain `centre_packed[i]` — the centre's bases as 2-bit
    // indices — so the scoring hot path can skip the per-position
    // `base_index` match.
    let mut centre_idx: Vec<usize> = vec![0];
    let mut centre_packed: Vec<Vec<u8>> = vec![pack_seq(&uniques[0].seq)];
    let mut assignments: Vec<usize> = vec![0usize; uniques.len()];

    // Track which sequences are already centres; maintained incrementally.
    let mut centre_set: std::collections::HashSet<&[u8]> =
        std::collections::HashSet::from([uniques[0].seq.as_slice()]);

    let mut prev_ll = f64::NEG_INFINITY;

    for iter in 0..cfg.max_iter {
        let n_centres_before = centre_idx.len();
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
            // Prior bypass: cross-sample priors auto-promote without the
            // Poisson abundance test (pseudo-pooling pass 2).
            if !priors.is_empty() && priors.contains(u.seq.as_slice()) {
                centre_set.insert(u.seq.as_slice());
                centre_idx.push(u_idx);
                centre_packed.push(pack_seq(&u.seq));
                continue;
            }
            // best_centre_for_promotion returns both the winning index AND
            // its log-likelihood — avoiding the redundant seq_ll call that
            // the previous version did to feed the significance test.
            let (best_i, log_prob_sub) =
                best_centre_for_promotion(&logp_table[u_idx], &centre_packed);
            // Gap-corrected likelihood: if the candidate is an indel artefact
            // of `centre_packed[best_i]`, the single-indel alignment scores
            // far better than the substitution-only path. Taking the max
            // suppresses spurious promotion of indel artefacts as new ASVs
            // and aligns the algorithm with R dada2's banded-aligner output.
            let log_prob_gap =
                gap_corrected_ll(&logp_table[u_idx], &centre_packed[best_i], cfg.gap_log_p);
            let log_prob = log_prob_sub.max(log_prob_gap);
            #[allow(clippy::cast_precision_loss)]
            let log_lambda = (total_reads as f64).ln() + log_prob;
            if is_significant_log(u64::from(u.count), log_lambda, cfg.omega_a) {
                centre_set.insert(u.seq.as_slice());
                centre_idx.push(u_idx);
                centre_packed.push(pack_seq(&u.seq));
            }
        }
        let n_centres_added = centre_idx.len() - n_centres_before;

        // Re-assign all uniques to their nearest centre. Outer par_iter over
        // uniques already saturates the thread pool, so each per-unique
        // best_centre call uses the SERIAL bound-pruned path — nested
        // rayon parallelism would cost more in fork overhead than it saves.
        let new_assignments: Vec<usize> = logp_table
            .par_iter()
            .map(|logp| best_centre_serial_packed(logp, &centre_packed).0)
            .collect();

        // Total log-likelihood (parallel sum).
        let ll: f64 = logp_table
            .par_iter()
            .zip(new_assignments.par_iter())
            .zip(uniques.par_iter())
            .map(|((logp, &ci), u)| {
                let safe_ci = ci.min(centre_packed.len() - 1);
                seq_ll_packed(logp, &centre_packed[safe_ci]) * f64::from(u.count)
            })
            .sum();

        let assignments_stable = new_assignments == assignments;
        assignments = new_assignments;

        let delta = (ll - prev_ll).abs();
        let n_centres = centre_idx.len();
        let max_iter = cfg.max_iter;
        log::info!("dada: iter {iter}/{max_iter}, {n_centres} centres, ΔlogL = {delta:.2e}");
        // Early-exit when (a) no centres were added this iteration AND
        // (b) assignments did not change. The next iteration would re-test
        // the same non-centre uniques against the same centres and produce
        // identical results.
        if iter > 0 && n_centres_added == 0 && assignments_stable {
            break;
        }
        if delta < cfg.tol {
            break;
        }
        prev_ll = ll;
    }

    // Collect ASVs by aggregating counts per centre.
    let mut asv_counts: std::collections::HashMap<usize, u32> = std::collections::HashMap::new();
    for (u, &ci) in uniques.iter().zip(assignments.iter()) {
        let ci = ci.min(centre_idx.len() - 1);
        *asv_counts.entry(ci).or_insert(0) += u.count;
    }

    let mut asvs: Vec<Asv> = asv_counts
        .into_iter()
        .map(|(ci, abundance)| Asv {
            sequence: uniques[centre_idx[ci]].seq.clone(),
            abundance,
        })
        .collect();
    asvs.sort_unstable_by_key(|a| Reverse(a.abundance));

    Ok(asvs)
}

/// Run DADA independently per sample, parallelised across samples via Rayon.
///
/// Identical to calling [`dada`] once per sample, but the per-sample calls
/// run concurrently across the Rayon thread pool. This is what should back
/// `dada(list_of_dereps, pool=FALSE)` in language bindings — the alternative
/// (host-language for-loop) leaves cross-sample parallelism unused.
///
/// # Errors
/// Returns [`Dada2Error`] if any per-sample [`dada`] call fails.
pub fn dada_many(
    samples: &[&[UniqueSeq]],
    error_model: &ErrorModel,
    cfg: &DadaConfig,
) -> Result<Vec<Vec<Asv>>, Dada2Error> {
    samples
        .par_iter()
        .map(|s| dada(s, error_model, cfg))
        .collect()
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
        let uniques = vec![make_unique("ACGTACGTACGT", 1000)];
        let em = ErrorModel::illumina_default();
        let cfg = DadaConfig::default();
        let asvs = dada(&uniques, &em, &cfg).unwrap();
        assert_eq!(asvs.len(), 1);
        assert_eq!(asvs[0].sequence, b"ACGTACGTACGT");
        assert_eq!(asvs[0].abundance, 1000);
    }
}
