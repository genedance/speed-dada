//! Sequence table: sample × ASV count matrix.

use crate::{dada::Asv, Dada2Error};
use std::collections::HashMap;
use std::path::Path;

/// Sample × ASV abundance matrix.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SequenceTable {
    /// Sample names in row order.
    pub samples: Vec<String>,
    /// ASV sequences in column order (hex-encoded for JSON compatibility).
    pub sequences: Vec<String>,
    /// `counts[sample_idx][asv_idx]` = abundance.
    pub counts: Vec<Vec<u32>>,
}

impl SequenceTable {
    /// Build a sequence table from per-sample DADA results.
    ///
    /// `sample_names` and `results` must be the same length.
    #[must_use]
    pub fn new(sample_names: &[&str], results: &[Vec<Asv>]) -> Self {
        // Collect all unique ASV sequences in stable order
        let mut seq_index: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut sequences: Vec<Vec<u8>> = Vec::new();

        for sample_asvs in results {
            for asv in sample_asvs {
                if !seq_index.contains_key(&asv.sequence) {
                    let idx = sequences.len();
                    seq_index.insert(asv.sequence.clone(), idx);
                    sequences.push(asv.sequence.clone());
                }
            }
        }

        let n_asvs = sequences.len();
        let mut counts: Vec<Vec<u32>> = results
            .iter()
            .map(|sample_asvs| {
                let mut row = vec![0u32; n_asvs];
                for asv in sample_asvs {
                    if let Some(&col) = seq_index.get(&asv.sequence) {
                        row[col] = asv.abundance;
                    }
                }
                row
            })
            .collect();

        // Ensure lengths match even if results is shorter than sample_names
        while counts.len() < sample_names.len() {
            counts.push(vec![0u32; n_asvs]);
        }

        let hex_seqs: Vec<String> = sequences.iter().map(|s| crate::bytes_to_hex(s)).collect();

        Self {
            samples: sample_names.iter().map(|s| (*s).to_owned()).collect(),
            sequences: hex_seqs,
            counts,
        }
    }

    /// Write tab-separated: rows = samples, columns = ASV sequences.
    ///
    /// # Errors
    /// Returns [`Dada2Error::Io`] on write failure.
    pub fn to_tsv(&self, path: &Path) -> Result<(), Dada2Error> {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);

        // Header row
        write!(f, "sample")?;
        for seq in &self.sequences {
            write!(f, "\t{seq}")?;
        }
        writeln!(f)?;

        // Data rows
        for (i, name) in self.samples.iter().enumerate() {
            write!(f, "{name}")?;
            if let Some(row) = self.counts.get(i) {
                for &cnt in row {
                    write!(f, "\t{cnt}")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }

    /// Serialise to JSON string.
    ///
    /// # Errors
    /// Returns [`Dada2Error::Parse`] if serialisation fails.
    pub fn to_json(&self) -> Result<String, Dada2Error> {
        serde_json::to_string(self).map_err(|e| Dada2Error::Parse(e.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dada::Asv;

    fn asv(seq: &str, abundance: u32) -> Asv {
        Asv { sequence: seq.bytes().collect(), abundance }
    }

    #[test]
    fn sequence_table_shape_and_values() {
        // 2 samples, 3 ASVs between them
        let s1 = vec![asv("AAAA", 100), asv("CCCC", 50)];
        let s2 = vec![asv("CCCC", 30), asv("TTTT", 20)];

        let table = SequenceTable::new(&["s1", "s2"], &[s1, s2]);

        assert_eq!(table.samples.len(), 2);
        assert_eq!(table.sequences.len(), 3, "expected 3 unique ASVs");
        assert_eq!(table.counts.len(), 2);
        assert_eq!(table.counts[0].len(), 3);
        assert_eq!(table.counts[1].len(), 3);

        // Check that total counts are preserved
        let total_s1: u32 = table.counts[0].iter().sum();
        let total_s2: u32 = table.counts[1].iter().sum();
        assert_eq!(total_s1, 150);
        assert_eq!(total_s2, 50);
    }
}
