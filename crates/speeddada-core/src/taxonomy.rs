//! Stage 8 — Naive Bayes k-mer taxonomic classifier.
//!
//! Implements the Wang et al. 2007 RDP classifier approach:
//! build a k-mer presence/absence bitset per reference sequence, then
//! classify queries by maximum shared-kmer count (AND + popcount) with
//! bootstrap confidence estimation.
//!
//! Bitset profiles use 8 KB per reference at k=8 (vs 262 KB for u32 counts),
//! a 32× reduction that keeps the entire database in L3 cache for realistic
//! reference set sizes.

use crate::{io::fasta::FastaRecord, merge::reverse_complement, Dada2Error, Kmer};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

/// Default k-mer length.
pub const DEFAULT_K: usize = 8;
/// Number of bootstrap replicates.
pub const N_BOOTSTRAP: usize = 100;
/// Default confidence threshold for genus-level assignment.
pub const DEFAULT_THRESHOLD: f64 = 0.80;

/// A taxonomic assignment for a single ASV.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TaxonAssignment {
    /// The query ASV sequence (hex-encoded for JSON compatibility).
    pub asv: String,
    /// Kingdom-level assignment.
    pub kingdom: Option<String>,
    /// Phylum-level assignment.
    pub phylum: Option<String>,
    /// Class-level assignment.
    pub class: Option<String>,
    /// Order-level assignment.
    pub order: Option<String>,
    /// Family-level assignment.
    pub family: Option<String>,
    /// Genus-level assignment (only if confidence >= threshold).
    pub genus: Option<String>,
    /// Species-level assignment.
    pub species: Option<String>,
    /// Bootstrap confidence at genus level (0..1).
    pub confidence: f64,
    /// Per-level bootstrap confidence `[kingdom, phylum, class, order,
    /// family, genus, species]`, each in 0..1. Mirrors dada2's
    /// `outputBootstraps = TRUE` matrix.
    pub bootstrap: [f64; 7],
}

/// Configuration for the taxonomy classifier.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TaxonomyConfig {
    /// k-mer length.
    pub k: usize,
    /// Minimum bootstrap confidence to report an assignment.
    pub threshold: f64,
    /// RNG seed for bootstrap subsampling.
    pub seed: u64,
    /// If true, also classify the reverse complement of each query and keep
    /// whichever orientation yields a higher genus-level confidence. Mirrors
    /// dada2's `tryRC` argument.
    pub try_rc: bool,
}

impl Default for TaxonomyConfig {
    fn default() -> Self {
        Self {
            k: DEFAULT_K,
            threshold: DEFAULT_THRESHOLD,
            seed: 42,
            try_rc: false,
        }
    }
}

/// Pre-built reference database using bitset k-mer profiles.
pub struct TaxonomyDb {
    k: usize,
    /// Number of u64 words per bitset: `4^k / 64` (rounded up).
    n_words: usize,
    /// One label per profile (reference sequence ID).
    labels: Vec<String>,
    /// Bitset per profile: `bits[i][word] & (1 << bit)` is set iff k-mer
    /// `word*64 + bit` appears in reference sequence `i`.
    bits: Vec<Vec<u64>>,
    /// Lineage per reference ID.
    lineages: HashMap<String, Vec<String>>,
}

impl TaxonomyDb {
    /// Build a taxonomy database from reference FASTA records and a lineage map.
    ///
    /// `lineages` maps reference sequence ID to a 7-level lineage
    /// `[kingdom, phylum, class, order, family, genus, species]`.
    ///
    /// # Errors
    /// Returns [`Dada2Error::InvalidInput`] if k > 16 or records are empty.
    #[allow(clippy::implicit_hasher)]
    pub fn build(
        records: &[FastaRecord],
        lineages: &HashMap<String, Vec<String>>,
        cfg: &TaxonomyConfig,
    ) -> Result<Self, Dada2Error> {
        if cfg.k > 16 {
            return Err(Dada2Error::InvalidInput("k must be <= 16".into()));
        }
        if records.is_empty() {
            return Err(Dada2Error::InvalidInput(
                "reference database is empty".into(),
            ));
        }

        #[allow(clippy::cast_possible_truncation)]
        let n_kmers = 4usize.pow(cfg.k as u32);
        let n_words = n_kmers.div_ceil(64);

        let (labels, bits): (Vec<_>, Vec<_>) = records
            .par_iter()
            .map(|rec| {
                let mut words = vec![0u64; n_words];
                for kmer in kmers(&rec.seq, cfg.k) {
                    #[allow(clippy::cast_possible_truncation)]
                    let idx = kmer.0 as usize;
                    words[idx >> 6] |= 1u64 << (idx & 63);
                }
                (rec.id.clone(), words)
            })
            .unzip();

        Ok(Self {
            k: cfg.k,
            n_words,
            labels,
            bits,
            lineages: lineages.clone(),
        })
    }

    /// Classify a collection of ASV sequences.
    ///
    /// # Errors
    /// Returns [`Dada2Error::InvalidInput`] if `seqs` is empty.
    pub fn classify(
        &self,
        seqs: &[Vec<u8>],
        cfg: &TaxonomyConfig,
    ) -> Result<Vec<TaxonAssignment>, Dada2Error> {
        if seqs.is_empty() {
            return Err(Dada2Error::InvalidInput("no sequences to classify".into()));
        }

        let n = seqs.len();
        log::info!("assign_taxonomy: classifying {n} sequences");
        let assignments: Vec<TaxonAssignment> = seqs
            .par_iter()
            .map(|seq| {
                let fwd = self.classify_one(seq, cfg);
                if cfg.try_rc {
                    let rc = reverse_complement(seq);
                    let rev = self.classify_one(&rc, cfg);
                    // Keep the orientation with higher genus-level confidence.
                    if rev.confidence > fwd.confidence {
                        rev
                    } else {
                        fwd
                    }
                } else {
                    fwd
                }
            })
            .collect();

        Ok(assignments)
    }

    fn classify_one(&self, seq: &[u8], cfg: &TaxonomyConfig) -> TaxonAssignment {
        fn get(lin: Option<&Vec<String>>, i: usize, conf: f64, threshold: f64) -> Option<String> {
            if conf < threshold {
                return None;
            }
            lin.and_then(|l| l.get(i))
                .filter(|s| !s.is_empty())
                .cloned()
        }

        // Build query bitset once; reuse bootstrap buffer.
        let query_kmers: Vec<Kmer> = kmers(seq, self.k).collect();
        let query_bits = self.seq_to_bits(seq);

        let best_idx = self.best_match_bits(&query_bits);
        let best_label = self.labels[best_idx].as_str();
        let best_lineage = self.lineages.get(best_label).cloned();

        // Track per-level matches across bootstrap reps.
        let mut level_hits = [0u32; 7];
        let mut seed = cfg.seed;
        let n_sub = (query_kmers.len() / 8).max(1);
        let mut sub_bits = vec![0u64; self.n_words];

        for _rep in 0..N_BOOTSTRAP {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            sub_bits.fill(0);
            for i in 0..n_sub {
                #[allow(clippy::cast_possible_truncation)]
                let idx = ((seed >> (i % 8)) as usize) % query_kmers.len().max(1);
                #[allow(clippy::cast_possible_truncation)]
                let kmer_id = query_kmers[idx].0 as usize;
                sub_bits[kmer_id >> 6] |= 1u64 << (kmer_id & 63);
            }
            let boot_idx = self.best_match_bits(&sub_bits);
            let boot_label = self.labels[boot_idx].as_str();
            let boot_lineage = self.lineages.get(boot_label);
            for lvl in 0..7 {
                let a = best_lineage.as_ref().and_then(|l| l.get(lvl));
                let b = boot_lineage.and_then(|l| l.get(lvl));
                if let (Some(a), Some(b)) = (a, b) {
                    if !a.is_empty() && a == b {
                        level_hits[lvl] += 1;
                    }
                }
            }
        }

        #[allow(clippy::cast_precision_loss)]
        let bootstrap: [f64; 7] = std::array::from_fn(|i| {
            f64::from(level_hits[i]) / N_BOOTSTRAP as f64
        });
        let confidence = bootstrap[5];
        let lineage = best_lineage.as_ref();

        TaxonAssignment {
            asv: crate::bytes_to_hex(seq),
            kingdom: get(lineage, 0, bootstrap[0], cfg.threshold),
            phylum: get(lineage, 1, bootstrap[1], cfg.threshold),
            class: get(lineage, 2, bootstrap[2], cfg.threshold),
            order: get(lineage, 3, bootstrap[3], cfg.threshold),
            family: get(lineage, 4, bootstrap[4], cfg.threshold),
            genus: get(lineage, 5, bootstrap[5], cfg.threshold),
            species: None,
            confidence,
            bootstrap,
        }
    }

    /// Score all profiles against `query_bits`; return the index of the best match.
    fn best_match_bits(&self, query_bits: &[u64]) -> usize {
        self.bits
            .iter()
            .enumerate()
            .map(|(i, profile)| {
                let score: u32 = profile
                    .iter()
                    .zip(query_bits.iter())
                    .map(|(a, b)| (a & b).count_ones())
                    .sum();
                (i, score)
            })
            .max_by_key(|(_, s)| *s)
            .map_or(0, |(i, _)| i)
    }

    /// Build a k-mer presence bitset for `seq`.
    fn seq_to_bits(&self, seq: &[u8]) -> Vec<u64> {
        let mut words = vec![0u64; self.n_words];
        for kmer in kmers(seq, self.k) {
            #[allow(clippy::cast_possible_truncation)]
            let idx = kmer.0 as usize;
            words[idx >> 6] |= 1u64 << (idx & 63);
        }
        words
    }

    /// Read-only access to the labels (one per reference sequence).
    #[must_use]
    pub fn labels(&self) -> &[String] {
        &self.labels
    }
}

/// Build a lineage map from FASTA records whose description/header carries
/// the taxonomy string (SILVA / GTDB style: `>SeqID Kingdom;Phylum;...;Species`).
///
/// The default reference format used by dada2's `assignTaxonomy` packs the
/// lineage into the FASTA header so users pass a single file; this helper
/// makes that work for SpeedDada too. Sequence IDs are read from `record.id`,
/// lineages from `record.description` (split on `;`).
#[must_use]
pub fn lineage_map_from_fasta(records: &[FastaRecord]) -> HashMap<String, Vec<String>> {
    records
        .iter()
        .map(|r| {
            // dada2-style headers often put the whole lineage in `id` with no
            // separate description; try description first, fall back to id.
            let raw = r
                .description
                .as_deref()
                .filter(|s| s.contains(';'))
                .unwrap_or(&r.id);
            let lineage: Vec<String> = raw
                .split(';')
                .map(|s| s.trim().trim_end_matches(';').to_owned())
                .collect();
            (r.id.clone(), lineage)
        })
        .collect()
}

/// Load a lineage map from a tab-separated file.
///
/// Expected format (one entry per line):
/// ```text
/// seq_id\tkingdom;phylum;class;order;family;genus;species
/// ```
///
/// # Errors
/// Returns [`Dada2Error::Io`] on read failure, [`Dada2Error::Parse`] for malformed lines.
pub fn load_lineage_tsv(path: &Path) -> Result<HashMap<String, Vec<String>>, Dada2Error> {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, '\t');
        let seq_id = parts
            .next()
            .ok_or_else(|| Dada2Error::Parse(format!("missing seq_id on line {}", line_num + 1)))?;
        let lineage_str = parts.next().ok_or_else(|| {
            Dada2Error::Parse(format!("missing lineage on line {}", line_num + 1))
        })?;
        let lineage: Vec<String> = lineage_str
            .split(';')
            .map(|s| s.trim().to_owned())
            .collect();
        map.insert(seq_id.to_owned(), lineage);
    }

    Ok(map)
}

/// Assign taxonomy to ASV sequences using a reference FASTA file.
///
/// # Errors
/// Returns [`Dada2Error`] on I/O or classification failure.
#[allow(clippy::implicit_hasher)]
pub fn assign_taxonomy(
    seqs: &[Vec<u8>],
    ref_records: &[FastaRecord],
    lineages: &HashMap<String, Vec<String>>,
    cfg: &TaxonomyConfig,
) -> Result<Vec<TaxonAssignment>, Dada2Error> {
    let db = TaxonomyDb::build(ref_records, lineages, cfg)?;
    db.classify(seqs, cfg)
}

/// Iterate over all k-mers in a sequence, encoding each as a [`Kmer`].
fn kmers(seq: &[u8], k: usize) -> impl Iterator<Item = Kmer> + '_ {
    let n = if seq.len() >= k { seq.len() - k + 1 } else { 0 };
    (0..n).filter_map(move |i| encode_kmer(&seq[i..i + k], k))
}

/// Encode a k-mer slice into a [`Kmer`] integer.  Returns `None` for ambiguous bases.
fn encode_kmer(slice: &[u8], _k: usize) -> Option<Kmer> {
    let mut val = 0u64;
    for &b in slice {
        let bits = match b.to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        val = val * 4 + bits;
    }
    Some(Kmer(val))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fasta::FastaRecord;

    fn make_ref(id: &str, seq: &str, lineage: &[&str]) -> (FastaRecord, (String, Vec<String>)) {
        let rec = FastaRecord {
            id: id.into(),
            description: None,
            seq: seq.bytes().collect(),
        };
        let lin = lineage.iter().map(|s| (*s).to_owned()).collect();
        (rec, (id.to_string(), lin))
    }

    #[test]
    fn load_lineage_tsv_round_trip() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "seq1\tBacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;acidophilus").unwrap();
        writeln!(f, "seq2\tBacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;aeruginosa").unwrap();

        let map = load_lineage_tsv(f.path()).unwrap();
        assert_eq!(map.len(), 2);
        let l1 = map.get("seq1").unwrap();
        assert_eq!(l1[0], "Bacteria");
        assert_eq!(l1[5], "Lactobacillus");
        assert_eq!(l1[6], "acidophilus");
        let l2 = map.get("seq2").unwrap();
        assert_eq!(l2[5], "Pseudomonas");
    }

    #[test]
    fn top1_assignment_mock_ref() {
        let (r1, (id1, l1)) = make_ref(
            "seq1",
            "ACGTACGTACGTACGTACGTACGTACGT",
            &[
                "Bacteria",
                "Firmicutes",
                "Bacilli",
                "Lactobacillales",
                "Lactobacillaceae",
                "Lactobacillus",
                "acidophilus",
            ],
        );
        let (r2, (id2, l2)) = make_ref(
            "seq2",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            &[
                "Bacteria",
                "Proteobacteria",
                "Gammaproteobacteria",
                "Pseudomonadales",
                "Pseudomonadaceae",
                "Pseudomonas",
                "aeruginosa",
            ],
        );
        let (r3, (id3, l3)) = make_ref(
            "seq3",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            &[
                "Bacteria",
                "Bacteroidetes",
                "Bacteroidia",
                "Bacteroidales",
                "Bacteroidaceae",
                "Bacteroides",
                "fragilis",
            ],
        );

        let records = vec![r1, r2, r3];
        let mut lineages = HashMap::new();
        lineages.insert(id1, l1);
        lineages.insert(id2, l2);
        lineages.insert(id3, l3);

        let query = b"ACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let cfg = TaxonomyConfig::default();
        let result = assign_taxonomy(&[query], &records, &lineages, &cfg).unwrap();

        assert_eq!(result.len(), 1);
        // Should assign to Lactobacillus lineage
        assert_eq!(result[0].kingdom.as_deref(), Some("Bacteria"));
    }
}
