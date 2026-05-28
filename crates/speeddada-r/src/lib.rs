//! `extendr` bindings for speeddada-core — exposes the DADA2 pipeline as an R package.
#![allow(non_snake_case)]
//!
//! Each function mirrors a step in the R dada2 API so the R wrapper layer can
//! call it via `.Call("wrap__<fn>", ...)`.  Errors are surfaced by panicking;
//! extendr wraps every call in `std::panic::catch_unwind` and converts panics
//! into R `stop()` errors automatically.

use extendr_api::prelude::*;

use speeddada_core::{
    chimera::remove_bimera_denovo,
    dada::{
        dada as dada_core, dada_many as dada_many_core, dada_pooled as dada_pooled_core,
        dada_pseudo as dada_pseudo_core, Asv, DadaConfig,
    },
    derep::{derep_fastq, UniqueSeq},
    error_model::{learn_errors, ErrorLearningConfig, ErrorModel},
    filter::{filter_and_trim_many, filter_and_trim_paired_many, FilterConfig},
    io::{fasta::read_fasta, fastq::read_fastq},
    merge::{merge_pairs, reverse_complement, MergeConfig},
    primer::{trim_primers, PrimerConfig},
    quality_profile::quality_profile,
    species::{assign_species, SpeciesConfig},
    taxonomy::{lineage_map_from_fasta, load_lineage_tsv, TaxonomyConfig, TaxonomyDb},
};
use std::path::{Path, PathBuf};

// ── ErrorModel external pointer ───────────────────────────────────────────────

/// Opaque wrapper so [`ErrorModel`] can live in an R `externalptr`.
pub struct RErrorModel(pub ErrorModel);

// Safety: ErrorModel contains only Array2<f64> + primitives, all Send.
#[allow(unsafe_code)]
unsafe impl Send for RErrorModel {}

/// Empty impl block required so extendr registers RErrorModel as an R type.
#[extendr]
impl RErrorModel {}

// ── Derep external pointer ────────────────────────────────────────────────────

/// Opaque wrapper so a per-sample dereplicated `Vec<UniqueSeq>` can live in
/// an R `externalptr`.
///
/// This carries the FULL `UniqueSeq` (including `qual_sum`) across the FFI
/// boundary, so the dada algorithm uses *real* per-position quality from
/// the FASTQ data — not a hardcoded Phred-30 placeholder that the previous
/// (seq, count)-only representation forced us to fabricate.
pub struct RDereped(pub Vec<UniqueSeq>);

#[allow(unsafe_code)]
unsafe impl Send for RDereped {}

#[extendr]
impl RDereped {}

// ── TaxonomyDb external pointer ───────────────────────────────────────────────

/// Opaque wrapper so a built [`TaxonomyDb`] (bitset k-mer profiles + lineage
/// map) can live in an R `externalptr` across many `assignTaxonomy` calls.
pub struct RTaxonomyDb(pub TaxonomyDb);

#[allow(unsafe_code)]
unsafe impl Send for RTaxonomyDb {}

#[extendr]
impl RTaxonomyDb {}

/// Reconstruct [`Asv`] vec from parallel `seqs` / `counts` vectors.
fn asvs_from_r(seqs: &[String], counts: &[i32]) -> Vec<Asv> {
    seqs.iter()
        .zip(counts.iter())
        .map(|(s, &c)| Asv {
            sequence: s.as_bytes().to_vec(),
            abundance: c as u32,
        })
        .collect()
}

// ── 1. filterAndTrim ─────────────────────────────────────────────────────────

/// Filter and trim FASTQ files (single- or paired-end).
///
/// `rev` / `rev_out` are `NULL` for single-end runs.  Returns a list
/// `(reads_in, reads_out, rownames)` — the R wrapper reshapes it into the
/// `[n × 2]` matrix that dada2's `filterAndTrim` returns.
#[allow(clippy::too_many_arguments)]
#[extendr]
fn filterAndTrim(
    fwd: Vec<String>,
    fwd_out: Vec<String>,
    rev: Nullable<Vec<String>>,
    rev_out: Nullable<Vec<String>>,
    trunc_len_fwd: i32,
    trunc_len_rev: i32,
    trim_left_fwd: i32,
    trim_left_rev: i32,
    max_ee_fwd: f64,
    max_ee_rev: f64,
    trunc_q: i32,
    min_len: i32,
) -> List {
    let cfg_fwd = FilterConfig {
        trunc_len: trunc_len_fwd as usize,
        min_len: min_len as usize,
        max_ee: max_ee_fwd,
        trunc_q: trunc_q as u8,
        trim_left: trim_left_fwd as usize,
        trim_right: 0,
    };

    let rownames: Vec<String> = fwd
        .iter()
        .map(|p| {
            Path::new(p)
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| p.clone())
        })
        .collect();

    let (reads_in, reads_out): (Vec<i32>, Vec<i32>) = match rev {
        Nullable::NotNull(rev_paths) => {
            let rev_out_paths = match rev_out {
                Nullable::NotNull(ro) => ro,
                Nullable::Null => fwd_out.clone(),
            };
            let cfg_rev = FilterConfig {
                trunc_len: trunc_len_rev as usize,
                min_len: min_len as usize,
                max_ee: max_ee_rev,
                trunc_q: trunc_q as u8,
                trim_left: trim_left_rev as usize,
                trim_right: 0,
            };
            let pairs: Vec<(PathBuf, PathBuf, PathBuf, PathBuf)> = fwd
                .iter()
                .zip(rev_paths.iter())
                .zip(fwd_out.iter())
                .zip(rev_out_paths.iter())
                .map(|(((f, r), fo), ro)| {
                    (
                        PathBuf::from(f),
                        PathBuf::from(r),
                        PathBuf::from(fo),
                        PathBuf::from(ro),
                    )
                })
                .collect();
            let stats = filter_and_trim_paired_many(&cfg_fwd, &cfg_rev, &pairs)
                .unwrap_or_else(|e| panic!("filterAndTrim: {e}"));
            stats
                .iter()
                .map(|s| (s.reads_in as i32, s.pairs_out as i32))
                .unzip()
        }
        Nullable::Null => {
            let pairs: Vec<(PathBuf, PathBuf)> = fwd
                .iter()
                .zip(fwd_out.iter())
                .map(|(f, fo)| (PathBuf::from(f), PathBuf::from(fo)))
                .collect();
            let stats = filter_and_trim_many(&cfg_fwd, &pairs)
                .unwrap_or_else(|e| panic!("filterAndTrim: {e}"));
            stats
                .iter()
                .map(|s| (s.reads_in as i32, s.reads_out as i32))
                .unzip()
        }
    };

    list!(
        reads_in = reads_in,
        reads_out = reads_out,
        rownames = rownames
    )
}

// ── 2. learnErrors ───────────────────────────────────────────────────────────

/// Learn error rates from FASTQ files.
///
/// `nbases` is converted to a read count via `nbases / 150.0`.
/// Returns an opaque [`RErrorModel`] external pointer.
#[extendr]
fn learnErrors(fls: Vec<String>, nbases: f64) -> ExternalPtr<RErrorModel> {
    let n_reads = ((nbases / 150.0).round() as usize).clamp(10_000, 10_000_000);
    let mut all_records = Vec::new();
    for path in &fls {
        let recs = read_fastq(Path::new(path)).unwrap_or_else(|e| panic!("learnErrors: {e}"));
        all_records.extend(recs);
        if all_records.len() >= n_reads {
            break;
        }
    }
    let cfg = ErrorLearningConfig {
        n_reads,
        ..Default::default()
    };
    let model = learn_errors(&all_records, &cfg).unwrap_or_else(|_| ErrorModel::illumina_default());
    ExternalPtr::new(RErrorModel(model))
}

// ── 3. derepFastq ────────────────────────────────────────────────────────────

/// Dereplicate one FASTQ file.
///
/// Returns a list `(seqs, counts, ptr)` — the R wrapper builds the named
/// integer `$uniques` vector (back-compat introspection), and stores `ptr`
/// as `.rust_ptr` on the derep object. Downstream `dada()` calls extract
/// the full `Vec<UniqueSeq>` (including per-position quality) through that
/// pointer, instead of round-tripping through `(seq, count)` tuples that
/// dropped quality information.
#[extendr]
fn derepFastq(path: &str) -> List {
    let records = read_fastq(Path::new(path)).unwrap_or_else(|e| panic!("derepFastq: {e}"));
    let uniques = derep_fastq(&records).unwrap_or_else(|e| panic!("derepFastq: {e}"));
    let seqs: Vec<String> = uniques
        .iter()
        .map(|u| String::from_utf8_lossy(&u.seq).into_owned())
        .collect();
    let counts: Vec<i32> = uniques.iter().map(|u| u.count as i32).collect();
    let ptr = ExternalPtr::new(RDereped(uniques));
    list!(seqs = seqs, counts = counts, ptr = ptr)
}

// ── 4. dada ──────────────────────────────────────────────────────────────────

/// Run DADA denoising on a single sample.
///
/// Takes the opaque dereplicated sample handle from `wrap__derepFastq`
/// (carrying per-position quality), the error model, and DADA parameters.
/// Returns `(seqs = character, counts = integer)` for the denoised ASVs.
#[extendr]
fn dada(
    derep: ExternalPtr<RDereped>,
    err: ExternalPtr<RErrorModel>,
    omega_a: f64,
    pool: bool,
) -> List {
    let cfg = DadaConfig {
        omega_a,
        pool,
        ..Default::default()
    };
    let asvs = dada_core(&derep.0, &err.0, &cfg).unwrap_or_else(|e| panic!("dada: {e}"));
    let out_seqs: Vec<String> = asvs
        .iter()
        .map(|a| String::from_utf8_lossy(&a.sequence).into_owned())
        .collect();
    let out_counts: Vec<i32> = asvs.iter().map(|a| a.abundance as i32).collect();
    list!(seqs = out_seqs, counts = out_counts)
}

// ── 4a. dada_many ───────────────────────────────────────────────────────────

/// Take a list of derep externalptrs and a wrapper closure that runs a
/// core multi-sample dada call (dada_many/dada_pooled/dada_pseudo), then
/// flattens the result back to parallel sample_idx / seqs / counts vectors.
fn dispatch_multi_sample<F>(
    dereps: List,
    err: &speeddada_core::error_model::ErrorModel,
    omega_a: f64,
    kernel: F,
) -> List
where
    F: FnOnce(
        &[&[UniqueSeq]],
        &speeddada_core::error_model::ErrorModel,
        &DadaConfig,
    ) -> Result<Vec<Vec<Asv>>, speeddada_core::Dada2Error>,
{
    // Robj for each list element; downcast each to ExternalPtr<RDereped>.
    let owned: Vec<ExternalPtr<RDereped>> = dereps
        .iter()
        .map(|(_, r)| {
            ExternalPtr::<RDereped>::try_from(r)
                .unwrap_or_else(|e| panic!("expected list of derep externalptrs: {e:?}"))
        })
        .collect();
    let refs: Vec<&[UniqueSeq]> = owned.iter().map(|p| p.0.as_slice()).collect();
    let cfg = DadaConfig {
        omega_a,
        ..Default::default()
    };
    let result = kernel(&refs, err, &cfg).unwrap_or_else(|e| panic!("dada multi-sample: {e}"));

    let mut out_sample_idx: Vec<i32> = Vec::new();
    let mut out_seqs: Vec<String> = Vec::new();
    let mut out_counts: Vec<i32> = Vec::new();
    for (si, asvs) in result.iter().enumerate() {
        for a in asvs {
            out_sample_idx.push(si as i32);
            out_seqs.push(String::from_utf8_lossy(&a.sequence).into_owned());
            out_counts.push(a.abundance as i32);
        }
    }
    list!(
        sample_idx = out_sample_idx,
        seqs = out_seqs,
        counts = out_counts
    )
}

/// Run DADA per-sample across multiple samples, parallelised via Rayon.
///
/// Takes a list of derep externalptrs (one per sample) carrying per-position
/// quality. Each sample is denoised independently; no cross-sample priors.
#[extendr]
fn dada_many(dereps: List, err: ExternalPtr<RErrorModel>, omega_a: f64) -> List {
    dispatch_multi_sample(dereps, &err.0, omega_a, dada_many_core)
}

// ── 4b. dada_pooled ─────────────────────────────────────────────────────────

/// Run DADA denoising on multiple samples with cross-sample pooling.
///
/// Takes a list of derep externalptrs (one per sample) carrying per-position
/// quality. Returns parallel vectors `(sample_idx, seqs, counts)` of
/// per-sample ASVs after pooled denoising.
#[extendr]
fn dada_pooled(dereps: List, err: ExternalPtr<RErrorModel>, omega_a: f64) -> List {
    dispatch_multi_sample(dereps, &err.0, omega_a, |refs, err, cfg| {
        let mut cfg = cfg.clone();
        cfg.pool = true;
        dada_pooled_core(refs, err, &cfg)
    })
}

// ── 4c. dada_pseudo ─────────────────────────────────────────────────────────

/// Run DADA with two-pass pseudo-pooling across samples.
///
/// Takes a list of derep externalptrs (one per sample) carrying per-position
/// quality. Returns parallel vectors `(sample_idx, seqs, counts)` of
/// per-sample ASVs after pseudo-pool denoising.
#[extendr]
fn dada_pseudo(dereps: List, err: ExternalPtr<RErrorModel>, omega_a: f64) -> List {
    dispatch_multi_sample(dereps, &err.0, omega_a, dada_pseudo_core)
}

// ── 5. mergePairs ────────────────────────────────────────────────────────────

/// Merge paired-end ASVs.
///
/// Returns a list-formatted data.frame with columns `sequence`, `abundance`,
/// `accept`, `nmatch`, `nmismatch`, `nindel` — matches dada2's `mergePairs`.
#[extendr]
fn mergePairs(
    fwd_seqs: Vec<String>,
    fwd_counts: Vec<i32>,
    rev_seqs: Vec<String>,
    rev_counts: Vec<i32>,
    min_overlap: i32,
    max_mismatch: i32,
    just_concatenate: bool,
) -> List {
    let fwd_asvs = asvs_from_r(&fwd_seqs, &fwd_counts);
    let rev_asvs = asvs_from_r(&rev_seqs, &rev_counts);
    let cfg = MergeConfig {
        min_overlap: min_overlap as usize,
        max_mismatches: max_mismatch as u32,
        just_concatenate,
    };
    let merged =
        merge_pairs(&fwd_asvs, &rev_asvs, &cfg).unwrap_or_else(|e| panic!("mergePairs: {e}"));

    let mut sequences = Vec::with_capacity(merged.len());
    let mut abundances = Vec::with_capacity(merged.len());
    let mut nmatch = Vec::with_capacity(merged.len());
    let mut nmismatch = Vec::with_capacity(merged.len());
    for m in &merged {
        sequences.push(String::from_utf8_lossy(&m.sequence).into_owned());
        abundances.push(m.abundance as i32);
        nmatch.push(m.overlap_len as i32);
        nmismatch.push(m.n_mismatches as i32);
    }
    let n = merged.len();
    list!(
        sequence = sequences,
        abundance = abundances,
        accept = vec![true; n],
        nmatch = nmatch,
        nmismatch = nmismatch,
        nindel = vec![0i32; n]
    )
}

// ── 6. makeSequenceTable ─────────────────────────────────────────────────────

/// Build a sample × ASV count matrix.
///
/// Inputs are flat parallel vectors covering all samples.  Returns a list
/// `(data, seqs, samples)` that the R wrapper reshapes into an integer matrix.
#[extendr]
fn makeSequenceTable(
    sample_names: Vec<String>,
    all_seqs: Vec<String>,
    all_counts: Vec<i32>,
    all_sample_idx: Vec<i32>,
    order_by_abundance: bool,
) -> List {
    use std::collections::HashMap;

    let mut asv_index: HashMap<&str, usize> = HashMap::new();
    let mut asvs: Vec<&str> = Vec::new();
    for s in &all_seqs {
        let next = asv_index.len();
        asv_index.entry(s.as_str()).or_insert_with(|| {
            asvs.push(s.as_str());
            next
        });
    }

    let n_samples = sample_names.len();
    let n_asvs = asvs.len();
    let mut mat = vec![0i32; n_samples * n_asvs];
    for ((seq, &cnt), &si) in all_seqs.iter().zip(&all_counts).zip(&all_sample_idx) {
        let ai = asv_index[seq.as_str()];
        mat[si as usize + n_samples * ai] += cnt;
    }

    let (final_mat, asv_strings): (Vec<i32>, Vec<String>) = if order_by_abundance {
        let col_sums: Vec<i32> = (0..n_asvs)
            .map(|ai| (0..n_samples).map(|si| mat[si + n_samples * ai]).sum())
            .collect();
        let mut order: Vec<usize> = (0..n_asvs).collect();
        order.sort_by(|&a, &b| col_sums[b].cmp(&col_sums[a]));
        let mut reordered = vec![0i32; n_samples * n_asvs];
        for (new_ai, &old_ai) in order.iter().enumerate() {
            for si in 0..n_samples {
                reordered[si + n_samples * new_ai] = mat[si + n_samples * old_ai];
            }
        }
        let new_asvs: Vec<String> = order.iter().map(|&i| asvs[i].to_owned()).collect();
        (reordered, new_asvs)
    } else {
        (mat, asvs.iter().map(|s| s.to_string()).collect())
    };

    list!(data = final_mat, seqs = asv_strings, samples = sample_names)
}

// ── 7. removeBimeraDenovo ────────────────────────────────────────────────────

/// Identify chimeric ASVs.
///
/// Returns an integer vector (1 = keep, 0 = chimera) parallel to `seqs`.
#[extendr]
fn removeBimeraDenovo(seqs: Vec<String>, counts: Vec<i32>) -> Vec<i32> {
    let pairs: Vec<(Vec<u8>, u32)> = seqs
        .iter()
        .zip(counts.iter())
        .map(|(s, &c)| (s.as_bytes().to_vec(), c as u32))
        .collect();

    let clean = remove_bimera_denovo(&pairs).unwrap_or_else(|e| panic!("removeBimeraDenovo: {e}"));
    let kept: std::collections::HashSet<&[u8]> = clean.iter().map(|(s, _)| s.as_slice()).collect();

    seqs.iter()
        .map(|s| {
            if kept.contains(s.as_bytes()) {
                1i32
            } else {
                0i32
            }
        })
        .collect()
}

// ── 8. assignTaxonomy ────────────────────────────────────────────────────────

/// Build a taxonomy reference database from a FASTA file.
///
/// If `lineage_tsv` is the empty string, lineages are parsed from the FASTA
/// header descriptions (SILVA / GTDB style `;`-separated). Returns an opaque
/// `externalptr` that downstream `assignTaxonomy` calls reuse so the costly
/// bitset profile build runs once per session.
#[extendr]
fn buildTaxonomyDb(
    ref_fasta: &str,
    lineage_tsv: &str,
    k: i32,
    threshold: f64,
    seed: f64,
) -> ExternalPtr<RTaxonomyDb> {
    let records =
        read_fasta(Path::new(ref_fasta)).unwrap_or_else(|e| panic!("buildTaxonomyDb: {e}"));
    let lineages = if lineage_tsv.is_empty() {
        lineage_map_from_fasta(&records)
    } else {
        load_lineage_tsv(Path::new(lineage_tsv))
            .unwrap_or_else(|e| panic!("buildTaxonomyDb: {e}"))
    };
    let cfg = TaxonomyConfig {
        k: k.max(2) as usize,
        threshold,
        seed: seed as u64,
        try_rc: false,
    };
    let db = TaxonomyDb::build(&records, &lineages, &cfg)
        .unwrap_or_else(|e| panic!("buildTaxonomyDb: {e}"));
    ExternalPtr::new(RTaxonomyDb(db))
}

/// Assign taxonomy to a set of ASV sequences against a pre-built database.
///
/// Returns parallel character vectors for each of the seven taxonomic levels
/// (`Kingdom..Species`), a `confidence` numeric vector (genus-level
/// bootstrap), and a flat `bootstrap` numeric vector of length `7 * n_seqs`
/// (row-major: `[asv_i * 7 + level]`). The R wrapper reshapes the bootstrap
/// vector into the `[n × 7]` matrix dada2 returns when `outputBootstraps =
/// TRUE`.
#[extendr]
fn assignTaxonomy(
    seqs: Vec<String>,
    db: ExternalPtr<RTaxonomyDb>,
    min_boot: f64,
    try_rc: bool,
) -> List {
    let query: Vec<Vec<u8>> = seqs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let cfg = TaxonomyConfig {
        k: 0, // unused at classify-time; DB was built with its own k
        threshold: (min_boot / 100.0).clamp(0.0, 1.0),
        seed: 42,
        try_rc,
    };
    let assignments = db
        .0
        .classify(&query, &cfg)
        .unwrap_or_else(|e| panic!("assignTaxonomy: {e}"));

    let n = assignments.len();
    let mut kingdom = Vec::with_capacity(n);
    let mut phylum = Vec::with_capacity(n);
    let mut class_ = Vec::with_capacity(n);
    let mut order = Vec::with_capacity(n);
    let mut family = Vec::with_capacity(n);
    let mut genus = Vec::with_capacity(n);
    let mut species = Vec::with_capacity(n);
    let mut confidence = Vec::with_capacity(n);
    let mut bootstrap = Vec::with_capacity(n * 7);
    let na = String::from("NA");
    for a in &assignments {
        kingdom.push(a.kingdom.clone().unwrap_or_else(|| na.clone()));
        phylum.push(a.phylum.clone().unwrap_or_else(|| na.clone()));
        class_.push(a.class.clone().unwrap_or_else(|| na.clone()));
        order.push(a.order.clone().unwrap_or_else(|| na.clone()));
        family.push(a.family.clone().unwrap_or_else(|| na.clone()));
        genus.push(a.genus.clone().unwrap_or_else(|| na.clone()));
        species.push(a.species.clone().unwrap_or_else(|| na.clone()));
        confidence.push(a.confidence);
        bootstrap.extend(a.bootstrap.iter().map(|x| x * 100.0));
    }

    list!(
        Kingdom = kingdom,
        Phylum = phylum,
        Class = class_,
        Order = order,
        Family = family,
        Genus = genus,
        Species = species,
        confidence = confidence,
        bootstrap = bootstrap
    )
}

// ── 8b. assignSpecies ────────────────────────────────────────────────────────

/// Exact-match (with optional Hamming-tolerance) species assignment.
///
/// Returns parallel `(genus, species)` character vectors over the input
/// sequences. Empty strings mean "no match" (or "ambiguous and
/// `allow_multiple = FALSE`").
#[extendr]
fn assignSpecies(
    seqs: Vec<String>,
    ref_fasta: &str,
    allow_multiple: bool,
    try_rc: bool,
    n_mismatch: i32,
) -> List {
    let ref_records =
        read_fasta(Path::new(ref_fasta)).unwrap_or_else(|e| panic!("assignSpecies: {e}"));
    let cfg = SpeciesConfig {
        try_rc,
        n_mismatch: n_mismatch.max(0) as u32,
        allow_multiple,
    };
    let query: Vec<Vec<u8>> = seqs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let out = assign_species(&query, &ref_records, &cfg)
        .unwrap_or_else(|e| panic!("assignSpecies: {e}"));

    let mut genus = Vec::with_capacity(out.len());
    let mut species = Vec::with_capacity(out.len());
    for a in &out {
        genus.push(a.genus.clone().unwrap_or_default());
        species.push(a.species.clone().unwrap_or_default());
    }
    list!(Genus = genus, Species = species)
}

// ── 9. removePrimers ─────────────────────────────────────────────────────────

/// Trim PCR primers from FASTQ reads.
///
/// `fwd_primer` / `rev_primer` are empty strings when the corresponding
/// side should be skipped. Returns parallel `(reads_in, reads_out, rownames)`
/// vectors over the input files — the R wrapper reshapes them into the
/// `[n × 2]` matrix that dada2's `removePrimers` returns.
#[allow(clippy::too_many_arguments)]
#[extendr]
fn removePrimers(
    fn_: Vec<String>,
    fout: Vec<String>,
    primer_fwd: String,
    primer_rev: String,
    max_mismatch: i32,
    min_overlap: i32,
    orient: bool,
) -> List {
    let _ = orient; // Currently only single-orientation matching is implemented.
    let cfg = PrimerConfig {
        fwd_primer: primer_fwd.as_bytes().to_vec(),
        rev_primer: primer_rev.as_bytes().to_vec(),
        max_mismatches: max_mismatch.max(0) as u32,
        min_overlap: min_overlap.max(1) as usize,
    };

    if fn_.len() != fout.len() {
        panic!("removePrimers: length(fn) must equal length(fout)");
    }

    let mut reads_in = Vec::with_capacity(fn_.len());
    let mut reads_out = Vec::with_capacity(fn_.len());
    let rownames: Vec<String> = fn_
        .iter()
        .map(|p| {
            Path::new(p)
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| p.clone())
        })
        .collect();

    for (input, output) in fn_.iter().zip(fout.iter()) {
        let stats = trim_primers(&cfg, Path::new(input), Path::new(output))
            .unwrap_or_else(|e| panic!("removePrimers: {e}"));
        reads_in.push(stats.reads_in as i32);
        reads_out.push(stats.reads_out as i32);
    }

    list!(
        reads_in = reads_in,
        reads_out = reads_out,
        rownames = rownames
    )
}

// ── 9b. qualityProfile ───────────────────────────────────────────────────────

/// Compute per-cycle Phred quality statistics for one FASTQ file.
///
/// Returns parallel numeric vectors of length `n_cycles`:
///   - `position` (1-based cycle index)
///   - `mean`, `q25`, `q50`, `q75` (Phred scores)
///   - `count` (number of reads reaching that cycle)
/// plus a scalar `n_reads`.
#[extendr]
fn qualityProfile(path: &str, n_reads: i32) -> List {
    let limit = if n_reads <= 0 { 0 } else { n_reads as usize };
    let profile = quality_profile(Path::new(path), limit)
        .unwrap_or_else(|e| panic!("qualityProfile: {e}"));
    let n = profile.cycle_mean.len();
    let position: Vec<i32> = (1..=n as i32).collect();
    let count: Vec<i32> = profile
        .cycle_count
        .iter()
        .map(|&c| c.min(i64::from(i32::MAX) as u64) as i32)
        .collect();
    list!(
        position = position,
        mean = profile.cycle_mean,
        q25 = profile.cycle_p25,
        q50 = profile.cycle_p50,
        q75 = profile.cycle_p75,
        count = count,
        n_reads = profile.n_reads as i32
    )
}

// ── 10. rc ───────────────────────────────────────────────────────────────────

/// Reverse-complement a vector of DNA strings.
///
/// Non-ACGT bases pass through as `N` (matches the original dada2 behaviour
/// of leaving ambiguous bases alone after RC). Vectorised so a single FFI
/// crossing handles all ASVs from a sequence table.
#[extendr]
fn rc(seqs: Vec<String>) -> Vec<String> {
    seqs.iter()
        .map(|s| String::from_utf8_lossy(&reverse_complement(s.as_bytes())).into_owned())
        .collect()
}

// ── Module registration ───────────────────────────────────────────────────────

extendr_module! {
    mod SpeedDada;
    impl RErrorModel;
    impl RDereped;
    impl RTaxonomyDb;
    fn filterAndTrim;
    fn learnErrors;
    fn derepFastq;
    fn dada;
    fn dada_many;
    fn dada_pooled;
    fn dada_pseudo;
    fn mergePairs;
    fn makeSequenceTable;
    fn removeBimeraDenovo;
    fn removePrimers;
    fn buildTaxonomyDb;
    fn assignTaxonomy;
    fn assignSpecies;
    fn qualityProfile;
    fn rc;
}
