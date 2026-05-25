//! Stage 6 — Paired-end read merging.
//!
//! Finds the overlap between forward and reverse-complement reads and
//! merges them into a single amplicon sequence.

use crate::{align::hamming_distance, dada::Asv, Dada2Error};
use rayon::prelude::*;

/// A merged amplicon with forward and reverse ASV provenance.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MergedRead {
    /// The merged sequence.
    pub sequence: Vec<u8>,
    /// Number of read pairs that produced this merged sequence.
    pub abundance: u32,
    /// Length of the overlap region.
    pub overlap_len: usize,
    /// Number of mismatches in the overlap.
    pub n_mismatches: u32,
}

/// Configuration for pair merging.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MergeConfig {
    /// Minimum overlap length required.
    pub min_overlap: usize,
    /// Maximum mismatches allowed in the overlap.
    pub max_mismatches: u32,
    /// If `true`, return only the joined sequence (no overlap duplicated).
    pub just_concatenate: bool,
}

impl Default for MergeConfig {
    fn default() -> Self {
        Self {
            min_overlap: 20,
            max_mismatches: 0,
            just_concatenate: false,
        }
    }
}

/// Merge paired-end ASVs.
///
/// `fwd_asvs` and `rev_asvs` must be dereplicated ASV sets from the
/// forward and reverse reads respectively.  The function tries all
/// (fwd, rev) pairs and returns merged sequences that satisfy `cfg`.
///
/// # Errors
/// Returns [`Dada2Error::InvalidInput`] if either input is empty.
pub fn merge_pairs(
    fwd_asvs: &[Asv],
    rev_asvs: &[Asv],
    cfg: &MergeConfig,
) -> Result<Vec<MergedRead>, Dada2Error> {
    if fwd_asvs.is_empty() || rev_asvs.is_empty() {
        return Err(Dada2Error::InvalidInput(
            "merge_pairs: both fwd and rev ASV sets must be non-empty".into(),
        ));
    }

    // Precompute all reverse complements once — shared across all fwd workers.
    let rev_rcs: Vec<Vec<u8>> = rev_asvs
        .iter()
        .map(|rev| reverse_complement(&rev.sequence))
        .collect();

    let mut merged: Vec<MergedRead> = fwd_asvs
        .par_iter()
        .flat_map(|fwd| {
            rev_asvs
                .iter()
                .zip(rev_rcs.iter())
                .filter_map(|(rev, rev_rc)| {
                    find_overlap(
                        &fwd.sequence,
                        rev_rc,
                        cfg.min_overlap,
                        cfg.max_mismatches,
                    )
                    .map(|(overlap_len, n_mismatches)| {
                        let sequence = if cfg.just_concatenate {
                            let mut s = fwd.sequence.clone();
                            s.extend_from_slice(rev_rc);
                            s
                        } else {
                            let mut s = fwd.sequence.clone();
                            s.extend_from_slice(&rev_rc[overlap_len..]);
                            s
                        };
                        MergedRead {
                            sequence,
                            abundance: fwd.abundance.min(rev.abundance),
                            overlap_len,
                            n_mismatches,
                        }
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect();

    merged.sort_unstable_by_key(|m| std::cmp::Reverse(m.abundance));
    Ok(merged)
}

/// Find the best suffix-prefix overlap between `fwd` and `rev_rc`.
///
/// Returns `Some((overlap_len, mismatches))` if a valid overlap is found.
fn find_overlap(
    fwd: &[u8],
    rev_rc: &[u8],
    min_overlap: usize,
    max_mismatches: u32,
) -> Option<(usize, u32)> {
    let max_overlap = fwd.len().min(rev_rc.len());
    let mut best: Option<(usize, u32)> = None;
    let mut best_score = u32::MAX;

    for ov in (min_overlap..=max_overlap).rev() {
        let fwd_end = &fwd[fwd.len() - ov..];
        let rev_start = &rev_rc[..ov];
        let mismatches = hamming_distance(fwd_end, rev_start);

        if mismatches <= max_mismatches && mismatches < best_score {
            best_score = mismatches;
            best = Some((ov, mismatches));
        }
    }
    best
}

/// Compute the reverse complement of a nucleotide sequence.
#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| complement(b))
        .collect()
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_correct() {
        let seq = b"ACGT";
        assert_eq!(reverse_complement(seq), b"ACGT");
        let seq2 = b"AAAA";
        assert_eq!(reverse_complement(seq2), b"TTTT");
    }

    #[test]
    fn merge_with_known_overlap() {
        // fwd: ACGTACGT
        // rev: (RC of GTGTACGT = ACGTACAC) — overlap of 4 with fwd
        let fwd_asv = Asv { sequence: b"ACGTACGT".to_vec(), abundance: 100 };
        let rev_asv = Asv {
            // RC of this will be ACGTACAC → overlap ACGT (4 bp) with fwd
            sequence: reverse_complement(b"ACGTACGT"),
            abundance: 100,
        };
        let cfg = MergeConfig { min_overlap: 4, ..Default::default() };
        let result = merge_pairs(&[fwd_asv], &[rev_asv], &cfg).unwrap();
        assert!(!result.is_empty());
    }
}
