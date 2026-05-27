//! Stage 1 — Per-cycle quality statistics from FASTQ files.

use crate::Dada2Error;
use std::path::Path;

/// Per-cycle quality summary statistics for a FASTQ file.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct QualityProfile {
    /// Number of reads sampled.
    pub n_reads: u64,
    /// Mean Phred score per cycle.
    pub cycle_mean: Vec<f64>,
    /// 25th-percentile Phred per cycle.
    pub cycle_p25: Vec<f64>,
    /// Median Phred per cycle.
    pub cycle_p50: Vec<f64>,
    /// 75th-percentile Phred per cycle.
    pub cycle_p75: Vec<f64>,
    /// Number of reads reaching each cycle (decreases for variable-length reads).
    pub cycle_count: Vec<u64>,
}

/// Maximum Phred quality score we track in histograms.
const MAX_PHRED: usize = 42; // 0..=41

/// Compute per-cycle quality statistics from a FASTQ file.
///
/// Reads at most `n_reads` records. Uses O(`max_read_len` × 42) memory.
///
/// # Errors
/// Returns [`Dada2Error`] on I/O or parse failure.
pub fn quality_profile(path: &Path, n_reads: usize) -> Result<QualityProfile, Dada2Error> {
    use needletail::parse_fastx_file;

    let mut reader = parse_fastx_file(path)
        .map_err(|e| Dada2Error::Parse(format!("cannot open {}: {e}", path.display())))?;

    // hist[cycle][phred] = count of reads with that Phred at that cycle
    let mut hist: Vec<Vec<u64>> = Vec::new();
    // cycle_count[cycle] = number of reads reaching this cycle
    let mut cycle_count: Vec<u64> = Vec::new();
    let mut n_read = 0u64;

    while let Some(rec_result) = reader.next() {
        if n_reads > 0 && n_read >= n_reads as u64 {
            break;
        }
        let rec = rec_result.map_err(|e| Dada2Error::Parse(e.to_string()))?;
        let Some(qual) = rec.qual() else { continue };
        n_read += 1;

        // Extend histogram if this read is longer than previous ones
        if qual.len() > hist.len() {
            hist.resize_with(qual.len(), || vec![0u64; MAX_PHRED]);
            cycle_count.resize(qual.len(), 0u64);
        }

        for (i, &qc) in qual.iter().enumerate() {
            let phred = qc.saturating_sub(33) as usize;
            let phred_clamped = phred.min(MAX_PHRED - 1);
            hist[i][phred_clamped] += 1;
            cycle_count[i] += 1;
        }
    }

    let n_cycles = hist.len();
    let mut cycle_mean = vec![0.0f64; n_cycles];
    let mut cycle_p25 = vec![0.0f64; n_cycles];
    let mut cycle_p50 = vec![0.0f64; n_cycles];
    let mut cycle_p75 = vec![0.0f64; n_cycles];

    for (i, counts) in hist.iter().enumerate() {
        let total = cycle_count[i];
        if total == 0 {
            continue;
        }

        // Mean
        #[allow(clippy::cast_precision_loss)]
        let sum: f64 = counts
            .iter()
            .enumerate()
            .map(|(q, &c)| q as f64 * c as f64)
            .sum();
        #[allow(clippy::cast_precision_loss)]
        let total_f = total as f64;
        cycle_mean[i] = sum / total_f;

        // Percentiles via cumulative sum
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let p25_target = (total_f * 0.25).ceil() as u64;
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let p50_target = (total_f * 0.50).ceil() as u64;
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let p75_target = (total_f * 0.75).ceil() as u64;

        let mut cum = 0u64;
        let mut p25_set = false;
        let mut p50_set = false;
        for (q, &c) in counts.iter().enumerate() {
            cum += c;
            #[allow(clippy::cast_precision_loss)]
            if !p25_set && cum >= p25_target {
                cycle_p25[i] = q as f64;
                p25_set = true;
            }
            #[allow(clippy::cast_precision_loss)]
            if !p50_set && cum >= p50_target {
                cycle_p50[i] = q as f64;
                p50_set = true;
            }
            #[allow(clippy::cast_precision_loss)]
            if cum >= p75_target {
                cycle_p75[i] = q as f64;
                break;
            }
        }
    }

    Ok(QualityProfile {
        n_reads: n_read,
        cycle_mean,
        cycle_p25,
        cycle_p50,
        cycle_p75,
        cycle_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn uniform_quality_profile() {
        // 10 reads of length 20, all with quality 'I' = Phred 40
        let mut f = NamedTempFile::new().unwrap();
        let seq = "ACGTACGTACGTACGTACGT";
        let qual = "IIIIIIIIIIIIIIIIIIII";
        for i in 0..10 {
            writeln!(f, "@r{i}\n{seq}\n+\n{qual}").unwrap();
        }
        let profile = quality_profile(f.path(), 0).unwrap();
        assert_eq!(profile.cycle_mean.len(), 20, "expected 20 cycles");
        for &mean in &profile.cycle_mean {
            assert!((mean - 40.0).abs() < 0.01, "expected mean 40.0, got {mean}");
        }
        assert_eq!(profile.n_reads, 10);
    }
}
