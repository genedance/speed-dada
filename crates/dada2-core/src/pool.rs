//! Disk-backed pooled dereplication for pool=TRUE DADA processing.
//!
//! Instead of holding every sample's reads in RAM simultaneously, each
//! sample's dereplicated sequences are written to a temp file.  The pooled
//! DADA step pages sequences in on demand, keeping RAM proportional to the
//! number of UNIQUE sequences (much smaller than raw read count).

use crate::{derep::UniqueSeq, Dada2Error};
use std::{
    collections::BTreeMap,
    fs,
    io::{BufRead, BufReader, Write},
    path::PathBuf,
};

/// One entry in the pooled store: the combined count and per-sample breakdown.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PoolEntry {
    /// Total count across all samples.
    pub total_count: u32,
    /// (`sample_index`, count) pairs for reassignment after DADA.
    pub per_sample: Vec<(usize, u32)>,
    /// Per-position quality sum accumulated from all samples.
    pub qual_sum: Vec<f64>,
}

/// A disk-backed store that accumulates unique sequences from multiple samples.
pub struct PoolStore {
    /// In-memory index: seq → entry.  Flushed to disk when > `flush_threshold` entries.
    mem: BTreeMap<Vec<u8>, PoolEntry>,
    /// Path to the backing temp directory.
    dir: tempfile::TempDir,
    /// Paths to flushed chunk files.
    chunks: Vec<PathBuf>,
    /// Flush when in-memory entries exceed this count.
    flush_threshold: usize,
}

impl PoolStore {
    /// Create a new empty pool store.
    ///
    /// `flush_threshold` controls how many unique sequences are held in RAM
    /// before a chunk is flushed to disk.  A value of `500_000` uses roughly
    /// 200 MB for 200-bp sequences.
    ///
    /// # Errors
    /// Returns [`Dada2Error::Io`] if the temp directory cannot be created.
    pub fn new(flush_threshold: usize) -> Result<Self, Dada2Error> {
        Ok(Self {
            mem: BTreeMap::new(),
            dir: tempfile::TempDir::new()?,
            chunks: Vec::new(),
            flush_threshold,
        })
    }

    /// Add all unique sequences from `sample_idx` into the pool.
    ///
    /// # Errors
    /// Returns [`Dada2Error::Io`] on flush failure.
    pub fn add_sample(&mut self, sample_idx: usize, uniques: &[UniqueSeq]) -> Result<(), Dada2Error> {
        for u in uniques {
            let entry = self.mem.entry(u.seq.clone()).or_insert_with(|| PoolEntry {
                total_count: 0,
                per_sample: Vec::new(),
                qual_sum: vec![0.0; u.seq.len()],
            });
            entry.total_count += u.count;
            entry.per_sample.push((sample_idx, u.count));
            for (i, &q) in u.qual_sum.iter().enumerate() {
                if i < entry.qual_sum.len() {
                    entry.qual_sum[i] += q;
                }
            }
        }

        if self.mem.len() >= self.flush_threshold {
            self.flush_chunk()?;
        }
        Ok(())
    }

    /// Flush the current in-memory index to a chunk file.
    fn flush_chunk(&mut self) -> Result<(), Dada2Error> {
        if self.mem.is_empty() {
            return Ok(());
        }
        let chunk_path = self.dir.path().join(format!("chunk_{}.jsonl", self.chunks.len()));
        let mut f = std::io::BufWriter::new(fs::File::create(&chunk_path)?);
        for (seq, entry) in &self.mem {
            let line = serde_json::to_string(&(seq, entry))
                .map_err(|e| Dada2Error::Parse(e.to_string()))?;
            writeln!(f, "{line}")?;
        }
        self.chunks.push(chunk_path);
        self.mem.clear();
        Ok(())
    }

    /// Merge all chunks and the current in-memory state into a final `Vec<UniqueSeq>`
    /// suitable for passing to [`crate::dada::dada`].
    ///
    /// # Errors
    /// Returns [`Dada2Error::Io`] or [`Dada2Error::Parse`] on read failure.
    pub fn into_pooled_uniques(mut self) -> Result<(Vec<UniqueSeq>, Vec<PoolEntry>), Dada2Error> {
        // Flush any remaining in-memory data
        self.flush_chunk()?;

        // Re-merge all chunks
        let mut merged: BTreeMap<Vec<u8>, PoolEntry> = BTreeMap::new();

        for chunk_path in &self.chunks {
            let f = BufReader::new(fs::File::open(chunk_path)?);
            for line in f.lines() {
                let line = line?;
                let (seq, entry): (Vec<u8>, PoolEntry) = serde_json::from_str(&line)
                    .map_err(|e| Dada2Error::Parse(e.to_string()))?;
                let m = merged.entry(seq).or_insert_with(|| PoolEntry {
                    total_count: 0,
                    per_sample: Vec::new(),
                    qual_sum: entry.qual_sum.clone(),
                });
                m.total_count += entry.total_count;
                m.per_sample.extend_from_slice(&entry.per_sample);
                for (i, &q) in entry.qual_sum.iter().enumerate() {
                    if i < m.qual_sum.len() {
                        m.qual_sum[i] += q;
                    }
                }
            }
        }

        let mut sorted: Vec<(Vec<u8>, PoolEntry)> = merged.into_iter().collect();
        sorted.sort_unstable_by(|a, b| {
            b.1.total_count.cmp(&a.1.total_count).then_with(|| a.0.cmp(&b.0))
        });

        let mut uniques = Vec::with_capacity(sorted.len());
        let mut entries = Vec::with_capacity(sorted.len());

        for (seq, entry) in sorted {
            let qual_sum = entry.qual_sum.clone();
            uniques.push(UniqueSeq { seq, count: entry.total_count, qual_sum });
            entries.push(entry);
        }

        Ok((uniques, entries))
    }
}
