//! SIMD-friendly sequence alignment primitives.
//!
//! The core functions are written as tight scalar loops that LLVM reliably
//! auto-vectorises to SSE4/AVX2/NEON when compiled with `target-cpu=native`
//! (see `.cargo/config.toml`).  No unsafe blocks are required.

/// Count the number of positions where `a[i] != b[i]`.
///
/// Processes up to `min(a.len(), b.len())` positions.
/// LLVM auto-vectorises this to SIMD on all major architectures.
#[inline]
#[must_use]
pub fn hamming_distance(a: &[u8], b: &[u8]) -> u32 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| u32::from(x != y))
        .sum()
}

/// Find the first position where `a[i] != b[i]`.
///
/// Returns `None` if the sequences are identical up to `min(a.len(), b.len())`.
#[inline]
#[must_use]
pub fn first_mismatch(a: &[u8], b: &[u8]) -> Option<usize> {
    a.iter().zip(b.iter()).position(|(x, y)| x != y)
}

/// Check whether `a[start..end] == b[start..end]`.
#[inline]
#[must_use]
pub fn range_equal(a: &[u8], b: &[u8], start: usize, end: usize) -> bool {
    let end = end.min(a.len()).min(b.len());
    if start >= end {
        return true;
    }
    a[start..end] == b[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hamming_identical() {
        assert_eq!(hamming_distance(b"ACGT", b"ACGT"), 0);
    }

    #[test]
    fn hamming_all_different() {
        assert_eq!(hamming_distance(b"AAAA", b"TTTT"), 4);
    }

    #[test]
    fn hamming_partial_overlap() {
        assert_eq!(hamming_distance(b"ACGT", b"ACGG"), 1);
    }

    #[test]
    fn first_mismatch_none() {
        assert!(first_mismatch(b"ACGT", b"ACGT").is_none());
    }

    #[test]
    fn first_mismatch_at_pos_2() {
        assert_eq!(first_mismatch(b"AAAT", b"AACT"), Some(2));
    }
}
