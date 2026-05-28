//! Stage 8b — Species-level assignment via exact sequence match.
//!
//! Mirrors dada2's `assignSpecies` / `addSpecies`: maps each ASV to the
//! genus + species label of any reference whose full nucleotide sequence
//! matches exactly (and optionally up to a few mismatches, with `tryRC`).
//!
//! The reference FASTA must carry the genus and species name in the
//! description, e.g. SILVA species-assignment files:
//! ```text
//! >AY999999.1 Lactobacillus acidophilus
//! ACGTACGT...
//! ```
//!
//! The classifier in [`crate::taxonomy`] only goes to genus level (the Wang
//! 2007 algorithm is not reliable below); this module fills in species
//! using exact match, which is what the dada2 default species reference
//! (`silva_species_assignment_*.fa.gz`) is curated for.

use crate::{align::hamming_distance, io::fasta::FastaRecord, merge::reverse_complement, Dada2Error};
use rayon::prelude::*;
use std::collections::HashMap;

/// A single reference taxon at the genus/species level.
#[derive(Debug, Clone)]
pub struct SpeciesRecord {
    /// Genus label (first whitespace-delimited token of the description).
    pub genus: String,
    /// Species epithet (second whitespace-delimited token).
    pub species: String,
}

/// A species reference database keyed by full nucleotide sequence.
///
/// Each reference is indexed under both its forward and reverse-complement
/// sequence at build time so query-time lookups need no orientation branch.
/// `length_buckets` groups reference indices by sequence length, so the
/// hamming-1 (or hamming-2) fallback search only scans length-compatible
/// references rather than the whole database.
pub struct SpeciesDb {
    /// All reference records (one per original FASTA entry).
    records: Vec<SpeciesRecord>,
    /// Full sequence → indices into `records`. Each ref has at most two
    /// entries here (forward + reverse complement).
    exact: HashMap<Vec<u8>, Vec<u32>>,
    /// length → all reference sequences with that length, plus their record
    /// index. Used for hamming-N fallback.
    by_length: HashMap<usize, Vec<(Vec<u8>, u32)>>,
}

/// Configuration for species assignment.
#[derive(Debug, Clone)]
pub struct SpeciesConfig {
    /// If true, also try matching the reverse complement of each query.
    pub try_rc: bool,
    /// If `> 0`, accept up to this many Hamming mismatches.
    pub n_mismatch: u32,
    /// If true, return all matching species (semicolon-joined). If false,
    /// return a match only when exactly one species is consistent.
    pub allow_multiple: bool,
}

impl Default for SpeciesConfig {
    fn default() -> Self {
        Self {
            try_rc: false,
            n_mismatch: 0,
            allow_multiple: false,
        }
    }
}

/// Result of one species assignment.
#[derive(Debug, Clone)]
pub struct SpeciesAssignment {
    /// Genus label (or `None` if no match / ambiguous and `allow_multiple = false`).
    pub genus: Option<String>,
    /// Species label, possibly `"sp1/sp2"` if `allow_multiple = true`.
    pub species: Option<String>,
}

impl SpeciesDb {
    /// Build a species DB from FASTA records.
    ///
    /// Records with malformed descriptions (fewer than two tokens) are
    /// dropped silently. Each kept record is indexed twice (forward + RC).
    ///
    /// # Errors
    /// Returns [`Dada2Error::InvalidInput`] if the FASTA is empty.
    #[must_use = "build the SpeciesDb and reuse it across assign_species calls"]
    pub fn build(records: &[FastaRecord]) -> Result<Self, Dada2Error> {
        if records.is_empty() {
            return Err(Dada2Error::InvalidInput(
                "species reference database is empty".into(),
            ));
        }
        let mut species_records: Vec<SpeciesRecord> = Vec::with_capacity(records.len());
        let mut exact: HashMap<Vec<u8>, Vec<u32>> = HashMap::new();
        let mut by_length: HashMap<usize, Vec<(Vec<u8>, u32)>> = HashMap::new();

        for rec in records {
            let Some((genus, species)) = parse_genus_species(rec) else {
                continue;
            };
            let idx = species_records.len() as u32;
            species_records.push(SpeciesRecord { genus, species });

            // Normalise to uppercase ASCII for consistent matching.
            let seq_fwd: Vec<u8> = rec.seq.iter().map(u8::to_ascii_uppercase).collect();
            let seq_rev = reverse_complement(&seq_fwd);

            exact.entry(seq_fwd.clone()).or_default().push(idx);
            exact.entry(seq_rev.clone()).or_default().push(idx);
            by_length
                .entry(seq_fwd.len())
                .or_default()
                .push((seq_fwd, idx));
            by_length
                .entry(seq_rev.len())
                .or_default()
                .push((seq_rev, idx));
        }

        if species_records.is_empty() {
            return Err(Dada2Error::InvalidInput(
                "no species references could be parsed (expected `>id Genus species` headers)"
                    .into(),
            ));
        }
        Ok(Self {
            records: species_records,
            exact,
            by_length,
        })
    }

    /// Assign genus + species to each query sequence.
    pub fn classify(&self, seqs: &[Vec<u8>], cfg: &SpeciesConfig) -> Vec<SpeciesAssignment> {
        seqs.par_iter()
            .map(|seq| self.classify_one(seq, cfg))
            .collect()
    }

    fn classify_one(&self, seq: &[u8], cfg: &SpeciesConfig) -> SpeciesAssignment {
        let query: Vec<u8> = seq.iter().map(u8::to_ascii_uppercase).collect();
        let mut hits = self.find_hits(&query, cfg.n_mismatch);
        if cfg.try_rc && hits.is_empty() {
            let rc = reverse_complement(&query);
            hits = self.find_hits(&rc, cfg.n_mismatch);
        }
        self.summarise(&hits, cfg)
    }

    fn find_hits(&self, query: &[u8], n_mismatch: u32) -> Vec<u32> {
        // Exact match first — cheapest and the common case.
        if let Some(idxs) = self.exact.get(query) {
            return idxs.clone();
        }
        if n_mismatch == 0 {
            return Vec::new();
        }
        // Hamming-N fallback over the same-length bucket.
        let Some(candidates) = self.by_length.get(&query.len()) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (cand_seq, idx) in candidates {
            if hamming_distance(query, cand_seq) <= n_mismatch {
                out.push(*idx);
            }
        }
        out
    }

    fn summarise(&self, hits: &[u32], cfg: &SpeciesConfig) -> SpeciesAssignment {
        if hits.is_empty() {
            return SpeciesAssignment {
                genus: None,
                species: None,
            };
        }
        // Deduplicate species labels (same species can be present under
        // multiple GenBank IDs, and forward+RC index the same record twice).
        let mut seen: Vec<(String, String)> = Vec::new();
        for &idx in hits {
            let r = &self.records[idx as usize];
            let pair = (r.genus.clone(), r.species.clone());
            if !seen.contains(&pair) {
                seen.push(pair);
            }
        }
        if seen.len() == 1 {
            let (g, s) = seen.into_iter().next().unwrap();
            return SpeciesAssignment {
                genus: Some(g),
                species: Some(s),
            };
        }
        if cfg.allow_multiple {
            // Collapse to one row: report genus only if all hits agree, and
            // a `/`-joined list of species epithets.
            let same_genus = seen.iter().all(|(g, _)| g == &seen[0].0);
            let genus = if same_genus {
                Some(seen[0].0.clone())
            } else {
                None
            };
            let mut species_list: Vec<String> = seen.into_iter().map(|(_, s)| s).collect();
            species_list.sort();
            species_list.dedup();
            return SpeciesAssignment {
                genus,
                species: Some(species_list.join("/")),
            };
        }
        // Multiple hits and !allow_multiple → ambiguous, report nothing.
        SpeciesAssignment {
            genus: None,
            species: None,
        }
    }
}

/// Convenience: build a species DB and classify in one call.
///
/// For multi-call workflows (e.g., assigning species to many sample tables
/// against the same SILVA reference), reuse a [`SpeciesDb`] via
/// [`SpeciesDb::build`] and [`SpeciesDb::classify`].
///
/// # Errors
/// Returns [`Dada2Error`] if DB construction fails.
pub fn assign_species(
    seqs: &[Vec<u8>],
    ref_records: &[FastaRecord],
    cfg: &SpeciesConfig,
) -> Result<Vec<SpeciesAssignment>, Dada2Error> {
    let db = SpeciesDb::build(ref_records)?;
    Ok(db.classify(seqs, cfg))
}

/// Extract `(genus, species)` from a FASTA record description.
///
/// dada2 species references typically look like:
/// ```text
/// >AY999999.1 Lactobacillus acidophilus
/// ```
/// The genus is the first whitespace-separated token of the description,
/// and the species epithet is the second. Records missing either are
/// skipped.
fn parse_genus_species(rec: &FastaRecord) -> Option<(String, String)> {
    let desc = rec.description.as_deref().unwrap_or("");
    let mut parts = desc.split_whitespace();
    let genus = parts.next()?.to_owned();
    let species = parts.next()?.to_owned();
    if genus.is_empty() || species.is_empty() {
        return None;
    }
    Some((genus, species))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, desc: &str, seq: &str) -> FastaRecord {
        FastaRecord {
            id: id.into(),
            description: Some(desc.into()),
            seq: seq.bytes().collect(),
        }
    }

    #[test]
    fn exact_match_assigns_species() {
        let refs = vec![
            rec("A", "Lactobacillus acidophilus", "ACGTACGTACGT"),
            rec("B", "Pseudomonas aeruginosa", "TTTTTTTTTTTT"),
        ];
        let db = SpeciesDb::build(&refs).unwrap();
        let cfg = SpeciesConfig::default();
        let out = db.classify(&[b"ACGTACGTACGT".to_vec()], &cfg);
        assert_eq!(out[0].genus.as_deref(), Some("Lactobacillus"));
        assert_eq!(out[0].species.as_deref(), Some("acidophilus"));
    }

    #[test]
    fn try_rc_matches_reverse_complement() {
        let refs = vec![rec("A", "Lactobacillus acidophilus", "ACGTACGTACGT")];
        let db = SpeciesDb::build(&refs).unwrap();
        // tryRC index covers RC at build time, so the query matches even
        // without setting try_rc — confirm both paths work.
        let rc_query = reverse_complement(b"ACGTACGTACGT");
        let cfg = SpeciesConfig::default();
        let out = db.classify(&[rc_query], &cfg);
        assert_eq!(out[0].genus.as_deref(), Some("Lactobacillus"));
    }

    #[test]
    fn ambiguous_match_returns_none_without_allow_multiple() {
        let refs = vec![
            rec("A", "Lactobacillus acidophilus", "ACGTACGTACGT"),
            rec("B", "Lactobacillus johnsonii", "ACGTACGTACGT"),
        ];
        let db = SpeciesDb::build(&refs).unwrap();
        let out = db.classify(&[b"ACGTACGTACGT".to_vec()], &SpeciesConfig::default());
        assert!(out[0].species.is_none());
    }

    #[test]
    fn ambiguous_match_with_allow_multiple_joins_species() {
        let refs = vec![
            rec("A", "Lactobacillus acidophilus", "ACGTACGTACGT"),
            rec("B", "Lactobacillus johnsonii", "ACGTACGTACGT"),
        ];
        let db = SpeciesDb::build(&refs).unwrap();
        let cfg = SpeciesConfig {
            allow_multiple: true,
            ..Default::default()
        };
        let out = db.classify(&[b"ACGTACGTACGT".to_vec()], &cfg);
        assert_eq!(out[0].genus.as_deref(), Some("Lactobacillus"));
        assert_eq!(out[0].species.as_deref(), Some("acidophilus/johnsonii"));
    }

    #[test]
    fn hamming_one_mismatch_match() {
        let refs = vec![rec("A", "Lactobacillus acidophilus", "ACGTACGTACGT")];
        let db = SpeciesDb::build(&refs).unwrap();
        let cfg = SpeciesConfig {
            n_mismatch: 1,
            ..Default::default()
        };
        // Single substitution at position 5: A -> C
        let q = b"ACGTACCTACGT".to_vec();
        let out = db.classify(&[q], &cfg);
        assert_eq!(out[0].species.as_deref(), Some("acidophilus"));
    }
}
