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
        let qual_bytes = rec.qual().map(<[u8]>::to_vec).unwrap_or_default();

        if let Some((seq, qual)) = apply_filters_owned(seq_bytes, qual_bytes, cfg) {
            writeln!(writer, "@{id}")?;
            writer.write_all(&seq)?;
            write!(writer, "\n+\n")?;
            writer.write_all(&qual)?;
            writeln!(writer)?;
            reads_out += 1;
        }
    }
    log::info!("filter_and_trim: {reads_in} reads in, {reads_out} passed");
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

/// Statistics for paired-end `filter_and_trim`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FilterStatsPaired {
    /// Total pairs supplied as input.
    pub reads_in: u64,
    /// Pairs where both reads passed all filters.
    pub pairs_out: u64,
    /// Pairs discarded because the forward read failed.
    pub fwd_failed: u64,
    /// Pairs discarded because the reverse read failed.
    pub rev_failed: u64,
    /// Pairs discarded because both reads failed.
    pub both_failed: u64,
}

/// Filter and trim paired-end FASTQ files in lock-step.
///
/// A pair is written to output only if both forward and reverse reads pass
/// all filters. Discards both if either fails.
///
/// # Errors
/// Returns [`Dada2Error`] on I/O or parse failure.
pub fn filter_and_trim_paired(
    cfg_fwd: &FilterConfig,
    cfg_rev: &FilterConfig,
    r1_in: &Path,
    r2_in: &Path,
    r1_out: &Path,
    r2_out: &Path,
) -> Result<FilterStatsPaired, Dada2Error> {
    use needletail::parse_fastx_file;
    use std::io::Write;

    let mut reader1 = parse_fastx_file(r1_in)
        .map_err(|e| Dada2Error::Parse(format!("cannot open {}: {e}", r1_in.display())))?;
    let mut reader2 = parse_fastx_file(r2_in)
        .map_err(|e| Dada2Error::Parse(format!("cannot open {}: {e}", r2_in.display())))?;

    let mut writer1 = std::io::BufWriter::new(std::fs::File::create(r1_out)?);
    let mut writer2 = std::io::BufWriter::new(std::fs::File::create(r2_out)?);

    let mut reads_in = 0u64;
    let mut pairs_out = 0u64;
    let mut fwd_failed = 0u64;
    let mut rev_failed = 0u64;
    let mut both_failed = 0u64;

    loop {
        let rec1_opt = reader1.next();
        let rec2_opt = reader2.next();

        match (rec1_opt, rec2_opt) {
            (None, None) => break,
            (Some(r1), Some(r2)) => {
                reads_in += 1;

                let r1 = r1.map_err(|e| Dada2Error::Parse(e.to_string()))?;
                let r2 = r2.map_err(|e| Dada2Error::Parse(e.to_string()))?;

                let id1 = std::str::from_utf8(r1.id())
                    .map_err(|e| Dada2Error::Parse(e.to_string()))?
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_owned();
                let id2 = std::str::from_utf8(r2.id())
                    .map_err(|e| Dada2Error::Parse(e.to_string()))?
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_owned();

                let fwd_pass =
                    apply_filters_owned(r1.seq().to_vec(), r1.qual().map(<[u8]>::to_vec).unwrap_or_default(), cfg_fwd);
                let rev_pass =
                    apply_filters_owned(r2.seq().to_vec(), r2.qual().map(<[u8]>::to_vec).unwrap_or_default(), cfg_rev);

                match (fwd_pass, rev_pass) {
                    (Some((seq1, qual1)), Some((seq2, qual2))) => {
                        writeln!(writer1, "@{id1}")?;
                        writer1.write_all(&seq1)?;
                        write!(writer1, "\n+\n")?;
                        writer1.write_all(&qual1)?;
                        writeln!(writer1)?;

                        writeln!(writer2, "@{id2}")?;
                        writer2.write_all(&seq2)?;
                        write!(writer2, "\n+\n")?;
                        writer2.write_all(&qual2)?;
                        writeln!(writer2)?;

                        pairs_out += 1;
                    }
                    (None, None) => both_failed += 1,
                    (None, Some(_)) => fwd_failed += 1,
                    (Some(_), None) => rev_failed += 1,
                }
            }
            _ => {
                return Err(Dada2Error::Parse(
                    "paired FASTQ files have different numbers of records".into(),
                ));
            }
        }
    }

    log::info!(
        "filter_and_trim_paired: {reads_in} pairs in, {pairs_out} passed, \
         fwd_failed={fwd_failed} rev_failed={rev_failed} both_failed={both_failed}"
    );

    Ok(FilterStatsPaired { reads_in, pairs_out, fwd_failed, rev_failed, both_failed })
}

/// Filter and trim multiple paired-end FASTQ file pairs in parallel.
///
/// Each tuple is `(r1_in, r2_in, r1_out, r2_out)`.
/// Processes all pairs in parallel using Rayon.
///
/// # Errors
/// Returns the first error encountered.
pub fn filter_and_trim_paired_many(
    cfg_fwd: &FilterConfig,
    cfg_rev: &FilterConfig,
    pairs: &[(PathBuf, PathBuf, PathBuf, PathBuf)],
) -> Result<Vec<FilterStatsPaired>, Dada2Error> {
    use rayon::prelude::*;
    pairs
        .par_iter()
        .map(|(r1_in, r2_in, r1_out, r2_out)| {
            filter_and_trim_paired(cfg_fwd, cfg_rev, r1_in, r2_in, r1_out, r2_out)
        })
        .collect()
}

/// Apply all filter steps to owned seq/qual vecs.
///
/// Returns `Some((seq, qual))` if the read passes, or `None` if it should be dropped.
pub(crate) fn apply_filters_owned(
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
    fn test_filter_and_trim_paired_basic() {
        // 5 pairs; read 3 (index 2) of R1 has bad quality → fwd_failed == 1, pairs_out == 4
        let good_seq = "ACGTACGTACGT";
        let good_qual = "IIIIIIIIIIII"; // Phred 40

        let mut r1_lines = Vec::new();
        let mut r2_lines = Vec::new();
        for i in 0..5 {
            let (_seq, qual) = if i == 2 {
                // Force EE failure: '!' = Phred 0 → EE = 12 per base
                (good_seq, "!!!!!!!!!!!!".into())
            } else {
                (good_seq, good_qual.to_owned())
            };
            r1_lines.push(format!("@r{i}\n{good_seq}\n+\n{qual}"));
            r2_lines.push(format!("@r{i}\n{good_seq}\n+\n{good_qual}"));
        }
        let r1_in = write_temp_fastq(
            &r1_lines
                .iter()
                .map(|s| {
                    let parts: Vec<&str> = s.lines().collect();
                    (parts[0].trim_start_matches('@'), parts[1], parts[3])
                })
                .collect::<Vec<_>>(),
        );
        let r2_in = write_temp_fastq(
            &r2_lines
                .iter()
                .map(|s| {
                    let parts: Vec<&str> = s.lines().collect();
                    (parts[0].trim_start_matches('@'), parts[1], parts[3])
                })
                .collect::<Vec<_>>(),
        );
        let r1_out = NamedTempFile::new().unwrap();
        let r2_out = NamedTempFile::new().unwrap();

        let cfg_fwd = FilterConfig { max_ee: 0.5, min_len: 1, ..Default::default() };
        let cfg_rev = FilterConfig { max_ee: 100.0, min_len: 1, ..Default::default() };

        let stats = filter_and_trim_paired(
            &cfg_fwd,
            &cfg_rev,
            r1_in.path(),
            r2_in.path(),
            r1_out.path(),
            r2_out.path(),
        )
        .unwrap();

        assert_eq!(stats.reads_in, 5);
        assert_eq!(stats.pairs_out, 4);
        assert_eq!(stats.fwd_failed, 1);
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
