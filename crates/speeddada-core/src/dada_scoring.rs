//! Internal scoring primitives for the DADA denoising algorithm.
//!
//! Split out of `dada.rs` to keep that file under the 500-line cap. All items
//! are `pub(crate)` so the sibling modules (`dada`, `dada_pool`) can use them
//! without exposing them in the public API.

use crate::{
    derep::UniqueSeq,
    error_model::{base_index, ErrorModel},
    Phred,
};
use rayon::prelude::*;
use statrs::distribution::{Discrete, Poisson};

/// Centre count above which `best_centre_for_promotion` fans out across
/// rayon threads. Higher than the previous threshold (64) because the
/// per-call work in promotion at K=64 is too small to amortise rayon fork
/// overhead. At K=256, per-call work is ~256 × ~50 prune-avg positions
/// ≈ 13 k ops, enough to make fan-out worthwhile.
pub(crate) const BEST_CENTRE_PAR_THRESHOLD: usize = 256;

/// Precompute per-position log-probability table for a single unique sequence.
///
/// `result[i][tb]` = log P(u.seq[i] | `true_base=tb`, `mean_qual_at_i`).
/// Stored as `f32` — halves the logp-table RAM cost vs `f64` with no loss of
/// clustering accuracy (scores compared relatively, not absolutely).
pub(crate) fn precompute_logp(u: &UniqueSeq, em: &ErrorModel) -> Vec<[f32; 4]> {
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
pub(crate) fn pack_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| base_index(b)).collect()
}

/// Log-likelihood of `logp` against a pre-packed centre (values 0-3 per pos).
///
/// Pure table indexing — no `base_index` match in the inner loop. Upcasts
/// f32 entries to f64 before accumulating to preserve sum precision.
pub(crate) fn seq_ll_packed(logp: &[[f32; 4]], packed: &[u8]) -> f64 {
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
pub(crate) fn gap_corrected_ll(logp: &[[f32; 4]], centre: &[u8], gap_log_p: f64) -> f64 {
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

/// Return `(argmax_centre_idx, log_likelihood)` for the candidate `logp`
/// against the current centre set, when called from the SERIAL promotion
/// loop. All idle cores can fan out across the K centres.
///
/// Dispatches to a serial bound-pruning path for small K, or a parallel
/// reduce for large K. Both preserve the argmax exactly.
pub(crate) fn best_centre_for_promotion(
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
pub(crate) fn best_centre_serial_packed(
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
pub(crate) fn best_centre_parallel_packed(
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
pub(crate) fn is_significant_log(count: u64, log_lambda: f64, omega_a: f64) -> bool {
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
            // `& 3` keeps only the low two bits, so the `as usize` cast is
            // exact on any pointer width (avoids clippy::cast_possible_truncation).
            v.push(bases[(z & 3) as usize]);
        }
        UniqueSeq {
            seq: v,
            count: 1,
            qual_sum: vec![30.0; len],
        }
    }

    /// Reference `O(centres × seq_len)` implementation — no pruning. Used to
    /// verify that the pruned `best_centre_serial_packed` returns identical
    /// argmax + log-likelihood.
    fn best_centre_exhaustive_packed(logp: &[[f32; 4]], centre_packed: &[Vec<u8>]) -> (usize, f64) {
        centre_packed
            .iter()
            .enumerate()
            .map(|(i, packed)| (i, seq_ll_packed(logp, packed)))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, f64::NEG_INFINITY))
    }

    #[test]
    fn best_centre_pruning_matches_exhaustive() {
        let em = ErrorModel::illumina_default();
        let uniques: Vec<UniqueSeq> = (0..8u64)
            .map(|i| make_unique_seq(i * 100 + 1, 64))
            .collect();
        let centre_packed: Vec<Vec<u8>> = uniques.iter().map(|u| pack_seq(&u.seq)).collect();

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
        let em = ErrorModel::illumina_default();
        let uniques: Vec<UniqueSeq> = (0..300u64)
            .map(|i| make_unique_seq(i * 991 + 13, 64))
            .collect();
        let centre_packed: Vec<Vec<u8>> = uniques.iter().map(|u| pack_seq(&u.seq)).collect();
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
        let em = ErrorModel::illumina_default();
        let q = make_unique("AAAACCCCGGGGTTTT", 1);
        let logp = precompute_logp(&q, &em);
        let centre_packed = pack_seq(&q.seq);
        let sub_ll = seq_ll_packed(&logp, &centre_packed);
        let gap_ll = gap_corrected_ll(&logp, &centre_packed, -6.9);
        assert!(
            gap_ll <= sub_ll,
            "gap_ll ({gap_ll}) should not exceed sub_ll ({sub_ll}) for identical seqs"
        );
    }

    #[test]
    fn gap_corrected_ll_detects_single_insertion() {
        let em = ErrorModel::illumina_default();
        let centre_seq = b"AAAACCCCGGGGTTTT".to_vec();
        let mut query_seq = vec![b'A'];
        query_seq.extend_from_slice(&centre_seq[..centre_seq.len() - 1]);
        let query = UniqueSeq {
            seq: query_seq,
            count: 1,
            qual_sum: vec![40.0; centre_seq.len()],
        };
        let logp = precompute_logp(&query, &em);
        let centre_packed = pack_seq(&centre_seq);
        let sub_ll = seq_ll_packed(&logp, &centre_packed);
        let gap_ll = gap_corrected_ll(&logp, &centre_packed, -6.9);
        assert!(
            gap_ll > sub_ll + 10.0,
            "gap_ll ({gap_ll}) should be much higher than sub_ll ({sub_ll}) when query has an insertion"
        );
    }
}
