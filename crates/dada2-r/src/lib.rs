//! `extendr` bindings for dada2-core — exposes the DADA2 pipeline as an R package.
#![allow(non_snake_case)]
//!
//! Each function mirrors a step in the R dada2 API so the R wrapper layer can
//! call it via `.Call("wrap__<fn>", ...)`.  Errors are surfaced by panicking;
//! extendr wraps every call in `std::panic::catch_unwind` and converts panics
//! into R `stop()` errors automatically.

use extendr_api::prelude::*;

use dada2_core::{
    chimera::remove_bimera_denovo,
    dada::{dada as dada_core, Asv, DadaConfig},
    derep::{derep_fastq, UniqueSeq},
    error_model::{learn_errors, ErrorLearningConfig, ErrorModel},
    filter::{filter_and_trim_many, filter_and_trim_paired_many, FilterConfig},
    io::fastq::read_fastq,
    merge::{merge_pairs, MergeConfig},
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

// ── Shared helpers ────────────────────────────────────────────────────────────

/// Reconstruct dereplicated sequences from parallel `seqs` / `counts` vectors.
fn uniques_from_r(seqs: &[String], counts: &[i32]) -> Vec<UniqueSeq> {
    seqs.iter()
        .zip(counts.iter())
        .map(|(s, &c)| {
            let seq = s.as_bytes().to_vec();
            let len = seq.len();
            UniqueSeq { seq, count: c as u32, qual_sum: vec![30.0 * f64::from(c); len] }
        })
        .collect()
}

/// Reconstruct [`Asv`] vec from parallel `seqs` / `counts` vectors.
fn asvs_from_r(seqs: &[String], counts: &[i32]) -> Vec<Asv> {
    seqs.iter()
        .zip(counts.iter())
        .map(|(s, &c)| Asv { sequence: s.as_bytes().to_vec(), abundance: c as u32 })
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
                    (PathBuf::from(f), PathBuf::from(r), PathBuf::from(fo), PathBuf::from(ro))
                })
                .collect();
            let stats = filter_and_trim_paired_many(&cfg_fwd, &cfg_rev, &pairs)
                .unwrap_or_else(|e| panic!("filterAndTrim: {e}"));
            stats.iter().map(|s| (s.reads_in as i32, s.pairs_out as i32)).unzip()
        }
        Nullable::Null => {
            let pairs: Vec<(PathBuf, PathBuf)> = fwd
                .iter()
                .zip(fwd_out.iter())
                .map(|(f, fo)| (PathBuf::from(f), PathBuf::from(fo)))
                .collect();
            let stats = filter_and_trim_many(&cfg_fwd, &pairs)
                .unwrap_or_else(|e| panic!("filterAndTrim: {e}"));
            stats.iter().map(|s| (s.reads_in as i32, s.reads_out as i32)).unzip()
        }
    };

    list!(reads_in = reads_in, reads_out = reads_out, rownames = rownames)
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
        let recs =
            read_fastq(Path::new(path)).unwrap_or_else(|e| panic!("learnErrors: {e}"));
        all_records.extend(recs);
        if all_records.len() >= n_reads {
            break;
        }
    }
    let cfg = ErrorLearningConfig { n_reads, ..Default::default() };
    let model =
        learn_errors(&all_records, &cfg).unwrap_or_else(|_| ErrorModel::illumina_default());
    ExternalPtr::new(RErrorModel(model))
}

// ── 3. derepFastq ────────────────────────────────────────────────────────────

/// Dereplicate one FASTQ file.
///
/// Returns a list `(seqs = character, counts = integer)` — the R wrapper
/// builds the named integer `$uniques` vector and attaches class `"derep"`.
#[extendr]
fn derepFastq(path: &str) -> List {
    let records =
        read_fastq(Path::new(path)).unwrap_or_else(|e| panic!("derepFastq: {e}"));
    let uniques =
        derep_fastq(&records).unwrap_or_else(|e| panic!("derepFastq: {e}"));
    let seqs: Vec<String> =
        uniques.iter().map(|u| String::from_utf8_lossy(&u.seq).into_owned()).collect();
    let counts: Vec<i32> = uniques.iter().map(|u| u.count as i32).collect();
    list!(seqs = seqs, counts = counts)
}

// ── 4. dada ──────────────────────────────────────────────────────────────────

/// Run DADA denoising on a single sample.
///
/// Accepts a dereplicated sample as parallel `seqs` / `counts` vectors.
/// Returns `(seqs = character, counts = integer)` for the denoised ASVs.
#[extendr]
fn dada(
    seqs: Vec<String>,
    counts: Vec<i32>,
    err: ExternalPtr<RErrorModel>,
    omega_a: f64,
    pool: bool,
) -> List {
    let uniques = uniques_from_r(&seqs, &counts);
    let cfg = DadaConfig { omega_a, pool, ..Default::default() };
    let asvs =
        dada_core(&uniques, &err.0, &cfg).unwrap_or_else(|e| panic!("dada: {e}"));
    let out_seqs: Vec<String> =
        asvs.iter().map(|a| String::from_utf8_lossy(&a.sequence).into_owned()).collect();
    let out_counts: Vec<i32> = asvs.iter().map(|a| a.abundance as i32).collect();
    list!(seqs = out_seqs, counts = out_counts)
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
    let merged = merge_pairs(&fwd_asvs, &rev_asvs, &cfg)
        .unwrap_or_else(|e| panic!("mergePairs: {e}"));

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
        sequence  = sequences,
        abundance = abundances,
        accept    = vec![true; n],
        nmatch    = nmatch,
        nmismatch = nmismatch,
        nindel    = vec![0i32; n]
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

    let clean = remove_bimera_denovo(&pairs)
        .unwrap_or_else(|e| panic!("removeBimeraDenovo: {e}"));
    let kept: std::collections::HashSet<&[u8]> =
        clean.iter().map(|(s, _)| s.as_slice()).collect();

    seqs.iter().map(|s| if kept.contains(s.as_bytes()) { 1i32 } else { 0i32 }).collect()
}

// ── Module registration ───────────────────────────────────────────────────────

extendr_module! {
    mod dada2rs;
    impl RErrorModel;
    fn filterAndTrim;
    fn learnErrors;
    fn derepFastq;
    fn dada;
    fn mergePairs;
    fn makeSequenceTable;
    fn removeBimeraDenovo;
}
