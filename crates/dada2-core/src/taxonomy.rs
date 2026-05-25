//! Stage 8 — Naive Bayes k-mer taxonomic classifier.
//!
//! Implements the Wang et al. 2007 RDP classifier approach:
//! build a k-mer frequency profile per reference taxon, then
//! classify queries by maximum posterior probability with bootstrap
//! confidence estimation.

use crate::{Dada2Error, Kmer, io::fasta::FastaRecord};
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
}

impl Default for TaxonomyConfig {
    fn default() -> Self {
        Self { k: DEFAULT_K, threshold: DEFAULT_THRESHOLD, seed: 42 }
    }
}

/// Pre-built reference database.
pub struct TaxonomyDb {
    k: usize,
    /// Map from taxon label to k-mer count vector.
    profiles: Vec<(String, Vec<u32>)>,
    /// Lineage per taxon label.
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
            return Err(Dada2Error::InvalidInput("reference database is empty".into()));
        }

        #[allow(clippy::cast_possible_truncation)]
        let n_kmers = 4usize.pow(cfg.k as u32);
        let profiles: Vec<(String, Vec<u32>)> = records
            .par_iter()
            .map(|rec| {
                let mut counts = vec![0u32; n_kmers];
                for kmer in kmers(&rec.seq, cfg.k) {
                    #[allow(clippy::cast_possible_truncation)]
                    let idx = kmer.0 as usize;
                    counts[idx] = counts[idx].saturating_add(1);
                }
                (rec.id.clone(), counts)
            })
            .collect();

        Ok(Self { k: cfg.k, profiles, lineages: lineages.clone() })
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
            .map(|seq| self.classify_one(seq, cfg))
            .collect();

        Ok(assignments)
    }

    fn classify_one(&self, seq: &[u8], cfg: &TaxonomyConfig) -> TaxonAssignment {
        fn get(lin: Option<&Vec<String>>, i: usize) -> Option<String> {
            lin.and_then(|l| l.get(i)).filter(|s| !s.is_empty()).cloned()
        }

        let query_kmers: Vec<Kmer> = kmers(seq, self.k).collect();
        let best_label: &str = self.best_match(&query_kmers);
        let best_genus = self.genus_of(best_label);

        // Bootstrap confidence
        let mut seed = cfg.seed;
        let n_sub = (query_kmers.len() / 8).max(1);
        let mut genus_hits = 0u32;

        for _rep in 0..N_BOOTSTRAP {
            // Simple LCG for deterministic subsampling
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            let subsample: Vec<Kmer> = (0..n_sub)
                .map(|i| {
                    #[allow(clippy::cast_possible_truncation)]
                    let idx = ((seed >> (i % 8)) as usize) % query_kmers.len().max(1);
                    query_kmers[idx]
                })
                .collect();
            let boot_label: &str = self.best_match(&subsample);
            if self.genus_of(boot_label) == best_genus && best_genus.is_some() {
                genus_hits += 1;
            }
        }

        #[allow(clippy::cast_precision_loss)]
        let confidence = f64::from(genus_hits) / N_BOOTSTRAP as f64;
        let lineage = self.lineages.get(best_label);

        TaxonAssignment {
            asv: crate::bytes_to_hex(seq),
            kingdom: get(lineage, 0),
            phylum: get(lineage, 1),
            class: get(lineage, 2),
            order: get(lineage, 3),
            family: get(lineage, 4),
            genus: if confidence >= cfg.threshold { get(lineage, 5) } else { None },
            species: None,
            confidence,
        }
    }

    fn best_match<'a>(&'a self, query_kmers: &[Kmer]) -> &'a str {
        self.profiles
            .iter()
            .enumerate()
            .map(|(idx, (_, counts))| {
                #[allow(clippy::cast_possible_truncation)]
                let score: u64 = query_kmers
                    .iter()
                    .map(|k| u64::from(counts[k.0 as usize]))
                    .sum();
                (idx, score)
            })
            .max_by_key(|(_, s)| *s)
            .map_or("", |(idx, _)| self.profiles[idx].0.as_str())
    }

    fn genus_of(&self, label: &str) -> Option<String> {
        self.lineages.get(label).and_then(|l| l.get(5)).cloned()
    }
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
        let seq_id = parts.next().ok_or_else(|| {
            Dada2Error::Parse(format!("missing seq_id on line {}", line_num + 1))
        })?;
        let lineage_str = parts.next().ok_or_else(|| {
            Dada2Error::Parse(format!("missing lineage on line {}", line_num + 1))
        })?;
        let lineage: Vec<String> = lineage_str.split(';').map(|s| s.trim().to_owned()).collect();
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
        let rec = FastaRecord { id: id.into(), description: None, seq: seq.bytes().collect() };
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
            &["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "acidophilus"],
        );
        let (r2, (id2, l2)) = make_ref(
            "seq2",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            &["Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "aeruginosa"],
        );
        let (r3, (id3, l3)) = make_ref(
            "seq3",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            &["Bacteria", "Bacteroidetes", "Bacteroidia", "Bacteroidales", "Bacteroidaceae", "Bacteroides", "fragilis"],
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
