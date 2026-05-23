//! Stage 2 — Quality filtering and adapter trimming.
//!
//! Mirrors dada2's `filterAndTrim()` function.

use crate::{Dada2Error, Phred};
use std::path::{Path, PathBuf};

/// Configuration for the filter-and-trim stage.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FilterConfig {
    /// Truncate reads to this length (0 = no truncation).
    pub trunc_len: usize,
    /// Discard reads shorter than this after truncation.
    pub min_len: usize,
    /// Maximum number of expected errors allowed per read.
    pub max_ee: f64,
    /// Minimum Phred quality score; bases below this truncate the read.
    pub trunc_q: u8,
    /// Remove `n_left` bases from the left end before any other processing.
    pub trim_left: usize,
    /// Remove `n_right` bases from the right end before any other processing.
    pub trim_right: usize,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            trunc_len: 0,
            min_len: 20,
            max_ee: 2.0,
            trunc_q: 2,
            trim_left: 0,
            trim_right: 0,
        }
    }
}

/// Statistics returned by [`filter_and_trim`].
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FilterStats {
    /// Total reads supplied as input.
    pub reads_in: u64,
    /// Reads that passed all filters.
    pub reads_out: u64,
}

/// Filter and trim reads from `input` and write passing reads to `output`.
///
/// Processes one record at a time using needletail streaming — constant RAM.
///
/// # Errors
/// Returns [`Dada2Error`] on I/O or parse failure.
pub fn filter_and_trim(
    cfg: &FilterConfig,
    input: &Path,
    output: &Path,
) -> Result<FilterStats, Dada2Error> {
    use needletail::parse_fastx_file;
    use std::io::Write;

    if cfg.max_ee < 0.0 {
        return Err(Dada2Error::InvalidInput("max_ee must be >= 0".into()));
    }

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
        let seq_bytes = rec.seq().to_vec();
        let qual_bytes = rec.qual().map(|q| q.to_vec()).unwrap_or_default();

        if let Some((seq, qual)) = apply_filters_owned(seq_bytes, qual_bytes, cfg) {
            write!(writer, "@{id}\n")?;
            writer.write_all(&seq)?;
            write!(writer, "\n+\n")?;
            writer.write_all(&qual)?;
            writeln!(writer)?;
            reads_out += 1;
        }
    }
    Ok(FilterStats { reads_in, reads_out })
}

/// Filter and trim multiple FASTQ files in parallel.
///
/// Processes each `(input, output)` pair on a separate Rayon thread.
/// Returns one [`FilterStats`] per pair in the same order as input.
///
/// # Errors
/// Returns the first error encountered.
pub fn filter_and_trim_many(
    cfg: &FilterConfig,
    pairs: &[(PathBuf, PathBuf)],
) -> Result<Vec<FilterStats>, Dada2Error> {
    use rayon::prelude::*;
    pairs
        .par_iter()
        .map(|(inp, out)| filter_and_trim(cfg, inp, out))
        .collect()
}

/// Apply all filter steps to owned seq/qual vecs.
///
/// Returns `Some((seq, qual))` if the read passes, or `None` if it should be dropped.
fn apply_filters_owned(
    mut seq: Vec<u8>,
    mut qual: Vec<u8>,
    cfg: &FilterConfig,
) -> Option<(Vec<u8>, Vec<u8>)> {
    // Trim ends first
    if cfg.trim_left > 0 {
        if seq.len() <= cfg.trim_left {
            return None;
        }
        seq = seq[cfg.trim_left..].to_vec();
        qual = qual[cfg.trim_left..].to_vec();
    }
    if cfg.trim_right > 0 {
        let len = seq.len();
        if len <= cfg.trim_right {
            return None;
        }
        seq.truncate(len - cfg.trim_right);
        qual.truncate(len - cfg.trim_right);
    }

    // Truncate at low-quality base
    if let Some(pos) = qual
        .iter()
        .position(|&q| q.saturating_sub(33) < cfg.trunc_q)
    {
        seq.truncate(pos);
        qual.truncate(pos);
    }

    // Fixed-length truncation
    if cfg.trunc_len > 0 && seq.len() > cfg.trunc_len {
        seq.truncate(cfg.trunc_len);
        qual.truncate(cfg.trunc_len);
    }

    // Minimum length check
    if seq.len() < cfg.min_len {
        return None;
    }

    // Expected-error filter
    let ee: f64 = qual
        .iter()
        .map(|&q| Phred::from_ascii(q).error_prob())
        .sum();
    if ee > cfg.max_ee {
        return None;
    }

    Some((seq, qual))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::read_fastq;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp_fastq(records: &[(&str, &str, &str)]) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        for (id, seq, qual) in records {
            writeln!(f, "@{id}\n{seq}\n+\n{qual}").unwrap();
        }
        f
    }

    #[test]
    fn test_trunc_len_exact() {
        // Quality string: all 'I' = Phred 40 — well above any threshold
        let inp = write_temp_fastq(&[("r1", "ACGTACGTACGT", "IIIIIIIIIIII")]);
        let out = NamedTempFile::new().unwrap();
        let cfg = FilterConfig {
            trunc_len: 6,
            min_len: 1,
            max_ee: 100.0,
            trunc_q: 0,
            trim_left: 0,
            trim_right: 0,
        };
        let stats = filter_and_trim(&cfg, inp.path(), out.path()).unwrap();
        assert_eq!(stats.reads_in, 1);
        assert_eq!(stats.reads_out, 1);
        let recs = read_fastq(out.path()).unwrap();
        assert_eq!(recs[0].seq.len(), 6);
    }

    #[test]
    fn test_max_ee_filter() {
        // Low-quality read: '!' = Phred 0 → error_prob = 1.0 per base
        let inp = write_temp_fastq(&[("r1", "ACGT", "!!!!")]);
        let out = NamedTempFile::new().unwrap();
        let cfg = FilterConfig { max_ee: 0.1, min_len: 1, ..Default::default() };
        let stats = filter_and_trim(&cfg, inp.path(), out.path()).unwrap();
        assert_eq!(stats.reads_out, 0);
    }

    #[test]
    fn test_filter_many_parallel() {
        let inp1 = write_temp_fastq(&[("r1", "ACGTACGT", "IIIIIIII")]);
        let inp2 = write_temp_fastq(&[("r2", "TTTTTTTT", "IIIIIIII")]);
        let out1 = NamedTempFile::new().unwrap();
        let out2 = NamedTempFile::new().unwrap();
        let cfg = FilterConfig { min_len: 4, ..Default::default() };
        let pairs = vec![
            (inp1.path().to_path_buf(), out1.path().to_path_buf()),
            (inp2.path().to_path_buf(), out2.path().to_path_buf()),
        ];
        let stats = filter_and_trim_many(&cfg, &pairs).unwrap();
        assert_eq!(stats.len(), 2);
        assert_eq!(stats[0].reads_out, 1);
        assert_eq!(stats[1].reads_out, 1);
    }
}
