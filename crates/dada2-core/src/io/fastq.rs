//! Streaming FASTQ parser built on top of `needletail`.

use crate::Dada2Error;
use std::path::Path;

/// A single FASTQ record held in memory.
#[derive(Debug, Clone)]
pub struct FastqRecord {
    /// Sequence identifier (without leading `@`).
    pub id: String,
    /// Raw nucleotide bytes (uppercase A/C/G/T/N).
    pub seq: Vec<u8>,
    /// Phred+33 encoded quality string (same length as `seq`).
    pub qual: Vec<u8>,
}

impl FastqRecord {
    /// Return the mean Phred quality score across all bases.
    pub fn mean_quality(&self) -> f64 {
        if self.qual.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.qual.iter().map(|&q| f64::from(q.saturating_sub(33))).sum();
        sum / self.qual.len() as f64
    }

    /// Truncate sequence and quality to `len` bases.
    pub fn truncate(&mut self, len: usize) {
        self.seq.truncate(len);
        self.qual.truncate(len);
    }

    /// Return `true` if any window of `width` bases has mean quality < `min_q`.
    pub fn has_low_quality_window(&self, min_q: u8, width: usize) -> bool {
        if width == 0 || self.qual.len() < width {
            return false;
        }
        self.qual.windows(width).any(|w| {
            let mean: f64 = w.iter().map(|&q| f64::from(q.saturating_sub(33))).sum::<f64>()
                / w.len() as f64;
            mean < f64::from(min_q)
        })
    }
}

/// Read all records from a FASTQ file into a `Vec`.
///
/// # Errors
/// Returns [`Dada2Error::Io`] or [`Dada2Error::Parse`] on failure.
pub fn read_fastq(path: &Path) -> Result<Vec<FastqRecord>, Dada2Error> {
    use needletail::parse_fastx_file;

    let mut reader = parse_fastx_file(path)
        .map_err(|e| Dada2Error::Parse(format!("cannot open {}: {e}", path.display())))?;

    let mut records = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| Dada2Error::Parse(e.to_string()))?;
        let id = std::str::from_utf8(rec.id())
            .map_err(|e| Dada2Error::Parse(e.to_string()))?
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_owned();
        let seq = rec.seq().to_vec();
        let qual = rec.qual().map(|q| q.to_vec()).unwrap_or_default();
        records.push(FastqRecord { id, seq, qual });
    }
    Ok(records)
}

/// Write records to a FASTQ file.
///
/// # Errors
/// Returns [`Dada2Error::Io`] on failure.
pub fn write_fastq(path: &Path, records: &[FastqRecord]) -> Result<(), Dada2Error> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    for r in records {
        write!(f, "@{}\n", r.id)?;
        f.write_all(&r.seq)?;
        write!(f, "\n+\n")?;
        f.write_all(&r.qual)?;
        writeln!(f)?;
    }
    Ok(())
}
