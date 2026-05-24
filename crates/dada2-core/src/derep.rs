//! Stage 4 — Dereplication: collapse identical sequences and count abundances.

use crate::{Dada2Error, io::fastq::FastqRecord};
use std::collections::HashMap;

/// A unique sequence with its observed count and representative quality scores.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct UniqueSeq {
    /// The canonical sequence.
    pub seq: Vec<u8>,
    /// Number of reads with this exact sequence.
    pub count: u32,
    /// Per-base quality sums (used to compute mean quality per position).
    pub qual_sum: Vec<f64>,
}

impl UniqueSeq {
    /// Mean Phred quality at position `i` (returns 0 if out of bounds).
    #[must_use]
    pub fn mean_qual(&self, i: usize) -> f64 {
        if i >= self.qual_sum.len() || self.count == 0 {
            return 0.0;
        }
        self.qual_sum[i] / f64::from(self.count)
    }
}

/// Dereplicate a slice of FASTQ records.
///
/// Returns a list of [`UniqueSeq`] sorted by descending abundance.
///
/// # Errors
/// Always succeeds; returns `Ok` for API consistency.
pub fn derep_fastq(records: &[FastqRecord]) -> Result<Vec<UniqueSeq>, Dada2Error> {
    let mut map: HashMap<Vec<u8>, (u32, Vec<f64>)> = HashMap::new();

    for rec in records {
        let entry = map.entry(rec.seq.clone()).or_insert_with(|| {
            let zeros = vec![0.0f64; rec.seq.len()];
            (0, zeros)
        });
        entry.0 += 1;
        for (i, &qc) in rec.qual.iter().enumerate() {
            if i < entry.1.len() {
                entry.1[i] += f64::from(qc.saturating_sub(33));
            }
        }
    }

    let mut uniques: Vec<UniqueSeq> = map
        .into_iter()
        .map(|(seq, (count, qual_sum))| UniqueSeq { seq, count, qual_sum })
        .collect();

    // Deterministic sort: descending count, then lexicographic sequence as tiebreak
    uniques.sort_unstable_by(|a, b| b.count.cmp(&a.count).then_with(|| a.seq.cmp(&b.seq)));

    let total: u64 = uniques.iter().map(|u| u64::from(u.count)).sum();
    let n_uniq = uniques.len();
    log::info!("dereplicate: {total} reads → {n_uniq} unique sequences");

    Ok(uniques)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::FastqRecord;

    fn make(seq: &str) -> FastqRecord {
        FastqRecord {
            id: "x".into(),
            seq: seq.bytes().collect(),
            qual: vec![b'I'; seq.len()],
        }
    }

    #[test]
    fn dedup_count_correctness() {
        let records = vec![
            make("AAAA"),
            make("CCCC"),
            make("AAAA"),
            make("AAAA"),
            make("TTTT"),
        ];
        let uniq = derep_fastq(&records).unwrap();
        assert_eq!(uniq.len(), 3, "expected 3 unique sequences");

        // Most abundant first
        assert_eq!(uniq[0].seq, b"AAAA");
        assert_eq!(uniq[0].count, 3);

        let cccc = uniq.iter().find(|u| u.seq == b"CCCC").unwrap();
        assert_eq!(cccc.count, 1);
    }
}
