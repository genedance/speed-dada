//! Optional primer/adapter trimming before quality filtering.

use crate::{align::hamming_distance, filter::FilterStats, Dada2Error};
use std::path::Path;

/// Configuration for primer trimming.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PrimerConfig {
    /// Forward primer sequence (5′→3′).
    pub fwd_primer: Vec<u8>,
    /// Reverse primer sequence (5′→3′, as it appears in the read).
    pub rev_primer: Vec<u8>,
    /// Maximum mismatches allowed when locating a primer.
    pub max_mismatches: u32,
    /// Minimum bases of primer that must be present at the read start.
    pub min_overlap: usize,
}

/// Trim primers from all reads in `input`, write trimmed reads to `output`.
///
/// Reads where the primer cannot be located within `max_mismatches` are discarded.
///
/// # Errors
/// Returns [`Dada2Error`] on I/O or parse failure.
pub fn trim_primers(
    cfg: &PrimerConfig,
    input: &Path,
    output: &Path,
) -> Result<FilterStats, Dada2Error> {
    use needletail::parse_fastx_file;
    use std::io::Write;

    let mut reader = parse_fastx_file(input)
        .map_err(|e| Dada2Error::Parse(format!("cannot open {}: {e}", input.display())))?;
    let mut writer = std::io::BufWriter::new(std::fs::File::create(output)?);
    let mut reads_in = 0u64;
    let mut reads_out = 0u64;

    while let Some(rec_result) = reader.next() {
        let rec = rec_result.map_err(|e| Dada2Error::Parse(e.to_string()))?;
        reads_in += 1;

        let id = std::str::from_utf8(rec.id())
            .map_err(|e| Dada2Error::Parse(e.to_string()))?
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_owned();
        let seq = rec.seq().to_vec();
        let qual = rec.qual().map(<[u8]>::to_vec).unwrap_or_default();

        if let Some((trimmed_seq, trimmed_qual)) = apply_primer_trim(&seq, &qual, cfg) {
            writeln!(writer, "@{id}")?;
            writer.write_all(&trimmed_seq)?;
            write!(writer, "\n+\n")?;
            writer.write_all(&trimmed_qual)?;
            writeln!(writer)?;
            reads_out += 1;
        }
    }

    Ok(FilterStats {
        reads_in,
        reads_out,
    })
}

/// Locate and trim primers from a single read.
///
/// Returns `Some((trimmed_seq, trimmed_qual))` if both primers are found
/// within tolerance, `None` if the read should be discarded.
fn apply_primer_trim(seq: &[u8], qual: &[u8], cfg: &PrimerConfig) -> Option<(Vec<u8>, Vec<u8>)> {
    let fwd_len = cfg.fwd_primer.len();
    let rev_len = cfg.rev_primer.len();

    // --- Forward primer search at the 5' end ---
    let fwd_search_end = (fwd_len + 5).min(seq.len());
    let fwd_trim_start = if fwd_len == 0 {
        0
    } else {
        find_primer_position(
            seq,
            &cfg.fwd_primer,
            0,
            fwd_search_end,
            cfg.max_mismatches,
            cfg.min_overlap,
        )?
    };

    // --- Reverse primer search at the 3' end ---
    let rev_trim_end = if rev_len == 0 {
        seq.len()
    } else {
        // The reverse primer appears near the 3' end of the read.
        // Search window: last (rev_len + 5) bases.
        let rev_search_start = seq.len().saturating_sub(rev_len + 5);
        find_primer_end(
            seq,
            &cfg.rev_primer,
            rev_search_start,
            cfg.max_mismatches,
            cfg.min_overlap,
        )?
    };

    if fwd_trim_start >= rev_trim_end {
        return None;
    }

    Some((
        seq[fwd_trim_start..rev_trim_end].to_vec(),
        if qual.len() >= rev_trim_end {
            qual[fwd_trim_start..rev_trim_end].to_vec()
        } else {
            Vec::new()
        },
    ))
}

/// Scan `seq[0..search_end]` for `primer` starting at offset 0.
///
/// Returns the position immediately after the primer (i.e., where the insert begins),
/// or `None` if not found within `max_mismatches`.
fn find_primer_position(
    seq: &[u8],
    primer: &[u8],
    search_start: usize,
    search_end: usize,
    max_mismatches: u32,
    min_overlap: usize,
) -> Option<usize> {
    let primer_len = primer.len();
    if primer_len == 0 || seq.len() < min_overlap {
        return None;
    }

    // Try each starting offset within the search window
    for offset in search_start..search_end.saturating_sub(min_overlap).saturating_add(1) {
        let available = search_end.min(seq.len()).saturating_sub(offset);
        let compare_len = primer_len.min(available);
        if compare_len < min_overlap {
            break;
        }
        let dist = hamming_distance(&seq[offset..offset + compare_len], &primer[..compare_len]);
        if dist <= max_mismatches {
            return Some(offset + compare_len);
        }
    }
    None
}

/// Find the end of `primer` in `seq` near the 3′ end.
///
/// Returns the index *before* the primer begins (i.e., where the insert ends),
/// or `None` if not found within `max_mismatches`.
fn find_primer_end(
    seq: &[u8],
    primer: &[u8],
    search_start: usize,
    max_mismatches: u32,
    min_overlap: usize,
) -> Option<usize> {
    let primer_len = primer.len();
    if primer_len == 0 {
        return Some(seq.len());
    }
    if seq.len() < min_overlap {
        return None;
    }

    for offset in search_start..=seq.len().saturating_sub(min_overlap) {
        let available = seq.len().saturating_sub(offset);
        let compare_len = primer_len.min(available);
        if compare_len < min_overlap {
            break;
        }
        let dist = hamming_distance(&seq[offset..offset + compare_len], &primer[..compare_len]);
        if dist <= max_mismatches {
            return Some(offset);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::write_test_fastq as write_fastq;
    use tempfile::NamedTempFile;

    #[test]
    fn primers_trimmed_correctly() {
        // FWD primer = "AAAA", REV primer = "TTTT"
        // Insert = "CCCCCCCCCC" (10 C's)
        let fwd = b"AAAA".to_vec();
        let rev = b"TTTT".to_vec();
        let read = "AAAACCCCCCCCCCTTTT";
        let qual = "I".repeat(read.len());

        let inp = write_fastq(&[("r1", read, &qual)]);
        let out = NamedTempFile::new().unwrap();

        let cfg = PrimerConfig {
            fwd_primer: fwd,
            rev_primer: rev,
            max_mismatches: 0,
            min_overlap: 4,
        };
        let stats = trim_primers(&cfg, inp.path(), out.path()).unwrap();
        assert_eq!(stats.reads_in, 1);
        assert_eq!(stats.reads_out, 1);

        let content = std::fs::read_to_string(out.path()).unwrap();
        let seq_line = content.lines().nth(1).unwrap();
        assert_eq!(seq_line, "CCCCCCCCCC");
    }

    #[test]
    fn wrong_primer_discarded() {
        // Read starts with "GGGG" instead of "AAAA" primer → discard
        let read = "GGGGCCCCCCCCCCTTTT";
        let qual = "I".repeat(read.len());

        let inp = write_fastq(&[("r1", read, &qual)]);
        let out = NamedTempFile::new().unwrap();

        let cfg = PrimerConfig {
            fwd_primer: b"AAAA".to_vec(),
            rev_primer: b"TTTT".to_vec(),
            max_mismatches: 0,
            min_overlap: 4,
        };
        let stats = trim_primers(&cfg, inp.path(), out.path()).unwrap();
        assert_eq!(stats.reads_in, 1);
        assert_eq!(stats.reads_out, 0);
    }
}
