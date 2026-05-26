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
    let mut asv_counts: std::collections::HashMap<usize, u32> =
        std::collections::HashMap::new();
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

/// Pack a sequence's bytes to 2-bit base indices (0=A, 1=C, 2=G, 3=T, N→0).
///
/// Done once per centre when it's promoted; subsequent scoring iterates the
/// packed array directly, skipping the per-position `base_index` match.
fn pack_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| base_index(b)).collect()
}

/// Log-likelihood of `logp` against a pre-packed centre (values 0-3 per pos).
///
/// Pure table indexing — no `base_index` match in the inner loop. Upcasts
/// f32 entries to f64 before accumulating to preserve sum precision.
fn seq_ll_packed(logp: &[[f32; 4]], packed: &[u8]) -> f64 {
    logp.iter()
        .zip(packed.iter())
        .map(|(lp, &p)| f64::from(lp[p as usize]))
        .sum()
}

/// Best log-likelihood under the model "one indel anywhere in the query".
///
/// Without this, a sequencer-introduced indel at position p in the query
/// makes the substitution-only scorer count ~(L-p) phantom substitutions,
/// spuriously promoting the indel artefact as a new ASV. R dada2 handles
/// indels natively through a banded Needleman-Wunsch aligner; this function
/// closes most of that gap analytically without an aligner.
///
/// Computes:
///   for each breakpoint p in `[0, L]`:
///     - **insertion** model: score `query[0..p]` against `centre[0..p]`,
///       then add `gap_log_p` for the inserted base, then score
///       `query[p+1..L]` against `centre[p..L-1]` (centre shifted left).
///     - **deletion** model: score `query[0..p]` against `centre[0..p]`,
///       add `gap_log_p` for the missing base, then score `query[p..L-1]`
///       against `centre[p+1..L]` (centre shifted right).
///   returns the max over all `p` and both models.
///
/// Implemented in O(L) via prefix sums of the four possible per-position
/// scores: `(query=q[i], true=c[i])`, `(query=q[i], true=c[i-1])`,
/// `(query=q[i], true=c[i+1])`, etc. — done with one pass each.
#[allow(clippy::cast_possible_wrap)]
fn gap_corrected_ll(logp: &[[f32; 4]], centre: &[u8], gap_log_p: f64) -> f64 {
    let n = logp.len().min(centre.len());
    if n < 2 {
        return f64::NEG_INFINITY;
    }
    // Prefix sums of score(query[i] | true=centre[i]) — same as seq_ll_packed
    // truncated to length p.
    let mut prefix = vec![0.0f64; n + 1];
    for i in 0..n {
        prefix[i + 1] = prefix[i] + f64::from(logp[i][centre[i] as usize]);
    }
    // Suffix sums under the INSERTION model: score(query[i+1] | true=centre[i])
    // for i in [p, n-1] — i.e. centre shifted LEFT by one relative to query.
    // After the gap, query[p+1..n] is scored against centre[p..n-1].
    let mut ins_suffix = vec![0.0f64; n + 1];
    for i in (1..n).rev() {
        ins_suffix[i] = ins_suffix[i + 1] + f64::from(logp[i][centre[i - 1] as usize]);
    }
    // Suffix sums under the DELETION model: score(query[i] | true=centre[i+1])
    // for i in [p, n-2] — centre shifted RIGHT by one.
    // After the gap, query[p..n-1] is scored against centre[p+1..n].
    let mut del_suffix = vec![0.0f64; n + 1];
    for i in (0..n - 1).rev() {
        del_suffix[i] = del_suffix[i + 1] + f64::from(logp[i][centre[i + 1] as usize]);
    }
    let mut best = f64::NEG_INFINITY;
    for p in 0..n {
        // Insertion at position p in the query: prefix[p] + gap + ins_suffix[p+1]
        let ins_ll = prefix[p] + gap_log_p + ins_suffix[p + 1];
        if ins_ll > best {
            best = ins_ll;
        }
        // Deletion at position p in the query: prefix[p] + gap + del_suffix[p]
        let del_ll = prefix[p] + gap_log_p + del_suffix[p];
        if del_ll > best {
            best = del_ll;
        }
    }
    best
}

/// Centre count above which `best_centre_for_promotion` fans out across
/// rayon threads. Higher than the previous threshold (64) because the
/// per-call work in promotion at K=64 is too small to amortise rayon fork
/// overhead. At K=256, per-call work is ~256 × ~50 prune-avg positions
/// ≈ 13 k ops, enough to make fan-out worthwhile.
const BEST_CENTRE_PAR_THRESHOLD: usize = 256;

/// Return `(argmax_centre_idx, log_likelihood)` for the candidate `logp`
/// against the current centre set, when called from the SERIAL promotion
/// loop. All 13 idle cores can fan out across the K centres.
///
/// Dispatches to a serial bound-pruning path for small K, or a parallel
/// reduce for large K. Both preserve the argmax exactly.
fn best_centre_for_promotion(
    logp: &[[f32; 4]],
    centre_packed: &[Vec<u8>],
) -> (usize, f64) {
    debug_assert!(!centre_packed.is_empty());
    if centre_packed.len() < BEST_CENTRE_PAR_THRESHOLD {
        best_centre_serial_packed(logp, centre_packed)
    } else {
        best_centre_parallel_packed(logp, centre_packed)
    }
}

/// Bound-pruned serial scan. Used inside re-assignment (where outer
/// `par_iter` over uniques already saturates the thread pool) and as the
/// promotion path for small centre sets.
///
/// Log-probabilities are `<= 0`, so a partial sum is a non-increasing upper
/// bound on the final sum. The branch is checked **every position** — modern
/// branch prediction makes the per-position check essentially free, and
/// non-matching centres typically exit within ~5–10 positions instead of
/// running to the next 16-position boundary.
fn best_centre_serial_packed(
    logp: &[[f32; 4]],
    centre_packed: &[Vec<u8>],
) -> (usize, f64) {
    let mut best_ll = seq_ll_packed(logp, &centre_packed[0]);
    let mut best_i = 0usize;
    for (i, packed) in centre_packed.iter().enumerate().skip(1) {
        let n = logp.len().min(packed.len());
        let mut acc = 0.0f64;
        let mut pruned = false;
        for pos in 0..n {
            acc += f64::from(logp[pos][packed[pos] as usize]);
            if acc < best_ll {
                pruned = true;
                break;
            }
        }
        // `>=` (not `>`) keeps the LAST centre on exact ties — matching
        // the original `centres.iter().max_by(partial_cmp)` semantics.
        if !pruned && acc >= best_ll {
            best_ll = acc;
            best_i = i;
        }
    }
    (best_i, best_ll)
}

/// Parallel exhaustive scan across centres — used only from the serial
/// promotion loop when K exceeds [`BEST_CENTRE_PAR_THRESHOLD`].
///
/// Each rayon worker computes `seq_ll` for its slice of centres; the
/// reduction picks the global argmax with the same `>=` tie-break.
///
/// MUST NOT be called from inside another `par_iter` (e.g. re-assignment) —
/// nested rayon parallelism with this granularity costs more in fork
/// overhead than it saves.
fn best_centre_parallel_packed(
    logp: &[[f32; 4]],
    centre_packed: &[Vec<u8>],
) -> (usize, f64) {
    centre_packed
        .par_iter()
        .enumerate()
        .map(|(i, packed)| (i, seq_ll_packed(logp, packed)))
        .reduce(
            || (0usize, f64::NEG_INFINITY),
            |a, b| {
                // `>=` keeps later-indexed centre on ties, matching
                // best_centre_serial_packed. Rayon reduces left-to-right
                // within a slice, so index ordering is deterministic.
                if b.1 >= a.1 {
                    b
                } else {
                    a
                }
            },
        )
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

    /// Reference O(centres × seq_len) implementation — no pruning. Used to
    /// verify that the pruned `best_centre_serial_packed` returns identical
    /// argmax + log-likelihood.
    fn best_centre_exhaustive_packed(
        logp: &[[f32; 4]],
        centre_packed: &[Vec<u8>],
    ) -> (usize, f64) {
        centre_packed
            .iter()
            .enumerate()
            .map(|(i, packed)| (i, seq_ll_packed(logp, packed)))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, f64::NEG_INFINITY))
    }

    fn make_unique_seq(seed: u64, len: usize) -> UniqueSeq {
        // SplitMix64 — full 64-bit mixing so distinct seeds give distinct streams.
        let bases = b"ACGT";
        let mut s = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut v = Vec::with_capacity(len);
        for _ in 0..len {
            s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = s;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            v.push(bases[(z as usize) % 4]);
        }
        UniqueSeq { seq: v, count: 1, qual_sum: vec![30.0; len] }
    }

    #[test]
    fn best_centre_pruning_matches_exhaustive() {
        // Build a small fixture: 8 candidate centres of length 64 nt + queries.
        let em = ErrorModel::illumina_default();
        let uniques: Vec<UniqueSeq> =
            (0..8u64).map(|i| make_unique_seq(i * 100 + 1, 64)).collect();
        let centre_packed: Vec<Vec<u8>> =
            uniques.iter().map(|u| pack_seq(&u.seq)).collect();

        for qseed in 1..=20u64 {
            let q = make_unique_seq(qseed * 7919, 64);
            let logp = precompute_logp(&q, &em);
            let pruned = best_centre_serial_packed(&logp, &centre_packed);
            let exhaustive = best_centre_exhaustive_packed(&logp, &centre_packed);
            assert_eq!(pruned, exhaustive, "argmax/ll mismatch for qseed={qseed}");
        }
    }

    #[test]
    fn best_centre_parallel_matches_serial() {
        // Fixture with K = 300 centres > BEST_CENTRE_PAR_THRESHOLD (256), so
        // best_centre_for_promotion dispatches to the parallel path.
        let em = ErrorModel::illumina_default();
        let uniques: Vec<UniqueSeq> =
            (0..300u64).map(|i| make_unique_seq(i * 991 + 13, 64)).collect();
        let centre_packed: Vec<Vec<u8>> =
            uniques.iter().map(|u| pack_seq(&u.seq)).collect();
        assert!(centre_packed.len() >= BEST_CENTRE_PAR_THRESHOLD);

        for qseed in 1..=30u64 {
            let q = make_unique_seq(qseed * 7919, 64);
            let logp = precompute_logp(&q, &em);
            let parallel = best_centre_for_promotion(&logp, &centre_packed);
            let serial = best_centre_serial_packed(&logp, &centre_packed);
            assert_eq!(parallel, serial, "argmax mismatch at qseed={qseed}");
        }
    }

    #[test]
    fn gap_corrected_ll_no_indel_returns_substitution_score() {
        // If the query and centre are identical, both substitution-only and
        // gap-corrected paths should produce comparable likelihoods — the gap
        // model adds gap_log_p penalty so it should be strictly worse for the
        // no-indel case.
        let em = ErrorModel::illumina_default();
        let q = make_unique("AAAACCCCGGGGTTTT", 1);
        let logp = precompute_logp(&q, &em);
        let centre_packed = pack_seq(&q.seq);
        let sub_ll = seq_ll_packed(&logp, &centre_packed);
        let gap_ll = gap_corrected_ll(&logp, &centre_packed, -6.9);
        // gap_ll is best single-indel score; with identical seqs, any indel
        // alignment introduces mismatches, so gap_ll must be ≤ sub_ll.
        assert!(
            gap_ll <= sub_ll,
            "gap_ll ({gap_ll}) should not exceed sub_ll ({sub_ll}) for identical seqs"
        );
    }

    #[test]
    fn gap_corrected_ll_detects_single_insertion() {
        // Query is centre with one base inserted at the start.
        // The substitution-only scorer treats every position as a mismatch
        // (catastrophic LL); the gap-corrected scorer should recognise the
        // insertion and score close to the all-match LL plus one gap penalty.
        let em = ErrorModel::illumina_default();
        let centre_seq = b"AAAACCCCGGGGTTTT".to_vec();
        let mut query_seq = vec![b'A']; // extra A at start
        query_seq.extend_from_slice(&centre_seq[..centre_seq.len() - 1]); // drop last so lengths match
        let query = UniqueSeq {
            seq: query_seq,
            count: 1,
            qual_sum: vec![40.0; centre_seq.len()],
        };
        let logp = precompute_logp(&query, &em);
        let centre_packed = pack_seq(&centre_seq);
        let sub_ll = seq_ll_packed(&logp, &centre_packed);
        let gap_ll = gap_corrected_ll(&logp, &centre_packed, -6.9);
        // gap_ll should be much higher (closer to 0) than sub_ll since the
        // insertion model captures most positions as matches.
        assert!(
            gap_ll > sub_ll + 10.0,
            "gap_ll ({gap_ll}) should be much higher than sub_ll ({sub_ll}) when query has an insertion"
        );
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
