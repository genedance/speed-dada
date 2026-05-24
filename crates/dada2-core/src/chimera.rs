//! Stage 7 — Bimera (chimera) detection and removal.
//!
//! For each candidate sequence, tests all pairs of more-abundant sequences
//! as potential parents. A sequence is a bimera if there exist parents P1, P2
//! such that the candidate's left arm matches P1 and the right arm matches P2,
//! with the crossover at the first mismatch position.

use crate::{
    align::{first_mismatch, range_equal},
    Dada2Error,
};
use rayon::prelude::*;

/// Minimum arm length in bases for a valid bimera call.
const MIN_ARM_LEN: usize = 8;

/// Remove bimeric sequences from a list of (sequence, abundance) pairs.
///
/// Sequences are processed in descending abundance order.  A sequence is
/// retained if no bimera parents can be found among all more-abundant sequences.
///
/// # Errors
/// Always succeeds; returns `Ok` for API uniformity.
pub fn remove_bimera_denovo(
    seqs: &[(Vec<u8>, u32)],
) -> Result<Vec<(Vec<u8>, u32)>, Dada2Error> {
    let n_in = seqs.len();
    // Sort descending by abundance (stable for determinism)
    let mut sorted = seqs.to_vec();
    sorted.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    // Compute keep/reject in parallel; collect as booleans to preserve sorted order
    // without a second sort pass.
    let keep: Vec<bool> = sorted
        .par_iter()
        .enumerate()
        .map(|(i, (seq, abund))| {
            if i == 0 {
                return true; // most abundant is never a chimera
            }
            let parents: Vec<&[u8]> = sorted[..i]
                .iter()
                .filter(|(_, pa)| *pa > *abund)
                .map(|(s, _)| s.as_slice())
                .collect();
            !is_bimera(seq, &parents)
        })
        .collect();

    // Consume `sorted` directly — no clone per item, no second sort needed.
    let out: Vec<(Vec<u8>, u32)> = sorted
        .into_iter()
        .zip(keep)
        .filter_map(|(item, keep)| if keep { Some(item) } else { None })
        .collect();
    let n_out = out.len();
    log::info!("remove_bimeras: {n_in} sequences → {n_out} after chimera removal");
    Ok(out)
}

/// Return `true` if `candidate` is a bimera of any pair in `parents`.
fn is_bimera(candidate: &[u8], parents: &[&[u8]]) -> bool {
    for i in 0..parents.len() {
        for j in 0..parents.len() {
            if i == j {
                continue;
            }
            if check_bimera_pair(candidate, parents[i], parents[j]) {
                return true;
            }
        }
    }
    false
}

/// Test whether `candidate` is formed by left-arm from `p1` and right-arm from `p2`.
fn check_bimera_pair(candidate: &[u8], p1: &[u8], p2: &[u8]) -> bool {
    let len = candidate.len().min(p1.len()).min(p2.len());
    if len < MIN_ARM_LEN * 2 {
        return false;
    }

    // Find first position where candidate diverges from p1
    let Some(crossover) = first_mismatch(&candidate[..len], &p1[..len]) else {
        return false; // candidate identical to p1 — not a bimera
    };

    if crossover < MIN_ARM_LEN {
        return false;
    }

    let right_len = len - crossover;
    if right_len < MIN_ARM_LEN {
        return false;
    }

    // Right arm must match p2 from crossover onwards
    range_equal(candidate, p2, crossover, len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_bimera_detected() {
        // p1 = AAAAAAAAAACCCCCCCCCC (10xA + 10xC)
        // p2 = AAAAAAAAAATTTTTTTTTT (10xA + 10xT)
        // chimera = AAAAAAAAAATTTTTTTTTT formed from p1[0..10] + p2[10..20]
        let p1: Vec<u8> = b"AAAAAAAAAACCCCCCCCCC".to_vec();
        let p2: Vec<u8> = b"AAAAAAAAAATTTTTTTTTT".to_vec();
        // Bimera: first 10 from p1, last 10 from p2 = p2 in this case
        let bimera: Vec<u8> = [&p1[..10], &p2[10..]].concat();

        let seqs = vec![
            (p1, 1000),
            (p2, 900),
            (bimera, 50),
        ];
        let result = remove_bimera_denovo(&seqs).unwrap();
        // Bimera should be removed
        assert_eq!(result.len(), 2, "bimera should have been removed");
    }

    #[test]
    fn non_chimera_retained() {
        let s1: Vec<u8> = b"AAAAAAAAAACCCCCCCCCC".to_vec();
        let s2: Vec<u8> = b"GGGGGGGGGGTTTTTTTTTT".to_vec();
        let seqs = vec![(s1, 100), (s2, 80)];
        let result = remove_bimera_denovo(&seqs).unwrap();
        assert_eq!(result.len(), 2);
    }
}
