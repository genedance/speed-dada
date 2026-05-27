//! PyO3 pipeline function bindings for speeddada-core.
// PyO3 functions returning PyResult document errors through Python exceptions,
// not Rust doc-sections; and PyO3 macros consume return values automatically.
#![allow(clippy::missing_errors_doc, clippy::must_use_candidate)]

use crate::types::{
    PyDadaResult, PyDerepResult, PyErrorModel, PyFilterConfig, PyFilterStats, PyFilterStatsPaired,
    PyMergedRead, PyQualityProfile, PySequenceTable, PyTaxonAssignment,
};
use pyo3::prelude::*;
use rayon::prelude::*;
use speeddada_core::{
    chimera::remove_bimera_denovo,
    dada::{dada, dada_many, dada_pooled, dada_pseudo, DadaConfig},
    derep::{derep_fastq, derep_fastq_path},
    error_model::{learn_errors, ErrorLearningConfig},
    filter::{filter_and_trim, filter_and_trim_paired, FilterConfig},
    io::{
        fasta::read_fasta,
        fastq::{read_fastq, read_fastq_n},
    },
    merge::{merge_pairs, MergeConfig},
    primer::{trim_primers, PrimerConfig},
    quality_profile::quality_profile,
    runtime::RuntimeConfig,
    sequence_table::SequenceTable,
    taxonomy::{assign_taxonomy, load_lineage_tsv, TaxonomyConfig},
};
use std::{collections::HashMap, path::Path, path::PathBuf};

/// Return the crate version string.
#[pyfunction]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Configure the rayon thread pool based on available hardware.
///
/// Parameters
/// ----------
/// n_threads : int, optional
///     Override the number of worker threads. 0 means auto-detect.
/// mb_per_thread : int, optional
///     Assumed RAM cost per worker thread in MiB (default 512).
///     Use 800 for DADA-heavy runs; 64 for filter/taxonomy-only runs.
///
/// Returns
/// -------
/// tuple[int, int | None]
///     ``(n_threads_configured, available_ram_mb)``
#[pyfunction]
#[pyo3(name = "configure_runtime")]
#[pyo3(signature = (n_threads = 0, mb_per_thread = 512))]
pub fn configure_runtime_py(n_threads: usize, mb_per_thread: u64) -> (usize, Option<u64>) {
    let cfg = if n_threads == 0 {
        RuntimeConfig::detect_with(mb_per_thread)
    } else {
        RuntimeConfig::detect_with(mb_per_thread).with_threads(n_threads)
    };
    // Ignore errors — pool may already be initialised (idempotent).
    cfg.apply().ok();
    (cfg.n_threads, cfg.mem_available_mb)
}

/// Build a sequence table from per-sample DADA results.
///
/// Parameters
/// ----------
/// sample_names : list[str]
/// results : list[DadaResult]
///
/// Returns
/// -------
/// SequenceTable
#[pyfunction]
#[pyo3(name = "make_sequence_table")]
#[allow(clippy::needless_pass_by_value, clippy::unnecessary_wraps)]
pub fn make_sequence_table_py(
    sample_names: Vec<String>,
    results: Vec<PyRef<'_, PyDadaResult>>,
) -> PySequenceTable {
    let name_refs: Vec<&str> = sample_names.iter().map(String::as_str).collect();
    let asvs: Vec<_> = results.iter().map(|r| r.0.clone()).collect();
    PySequenceTable(SequenceTable::new(&name_refs, &asvs))
}

/// Trim primers from FASTQ reads.
///
/// Parameters
/// ----------
/// fwd_primer : str | bytes
/// rev_primer : str | bytes
/// input_path : str
/// output_path : str
/// max_mismatches : int, optional
/// min_overlap : int, optional
///
/// Returns
/// -------
/// FilterStats
#[pyfunction]
#[pyo3(name = "trim_primers")]
#[pyo3(signature = (fwd_primer, rev_primer, input_path, output_path, max_mismatches = 0, min_overlap = 10))]
pub fn trim_primers_py(
    py: Python<'_>,
    fwd_primer: Vec<u8>,
    rev_primer: Vec<u8>,
    input_path: &str,
    output_path: &str,
    max_mismatches: u32,
    min_overlap: usize,
) -> PyResult<PyFilterStats> {
    let inp = input_path.to_owned();
    let out = output_path.to_owned();
    let stats = py
        .allow_threads(move || {
            let cfg = PrimerConfig {
                fwd_primer,
                rev_primer,
                max_mismatches,
                min_overlap,
            };
            trim_primers(&cfg, Path::new(&inp), Path::new(&out))
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyFilterStats {
        reads_in: stats.reads_in,
        reads_out: stats.reads_out,
    })
}

/// Compute per-cycle quality statistics from a FASTQ file.
///
/// Parameters
/// ----------
/// fastq_path : str
/// n_reads : int, optional
///     Maximum number of reads to sample (0 = all reads).
///
/// Returns
/// -------
/// QualityProfile
#[pyfunction]
#[pyo3(name = "quality_profile")]
#[pyo3(signature = (fastq_path, n_reads = 500_000))]
pub fn quality_profile_py(
    py: Python<'_>,
    fastq_path: &str,
    n_reads: usize,
) -> PyResult<PyQualityProfile> {
    let path = fastq_path.to_owned();
    let profile = py
        .allow_threads(move || quality_profile(Path::new(&path), n_reads))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyQualityProfile {
        n_reads: profile.n_reads,
        cycle_mean: profile.cycle_mean,
        cycle_p25: profile.cycle_p25,
        cycle_p50: profile.cycle_p50,
        cycle_p75: profile.cycle_p75,
        cycle_count: profile.cycle_count,
    })
}

/// Filter and trim FASTQ reads.
///
/// Parameters
/// ----------
/// config : FilterConfig
/// input_path : str
/// output_path : str
///
/// Returns
/// -------
/// FilterStats
#[pyfunction]
#[pyo3(name = "filter_and_trim")]
pub fn filter_and_trim_py(
    py: Python<'_>,
    config: &PyFilterConfig,
    input_path: &str,
    output_path: &str,
) -> PyResult<PyFilterStats> {
    let cfg = config.0.clone();
    let inp = input_path.to_owned();
    let out = output_path.to_owned();
    let stats = py
        .allow_threads(move || filter_and_trim(&cfg, Path::new(&inp), Path::new(&out)))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyFilterStats {
        reads_in: stats.reads_in,
        reads_out: stats.reads_out,
    })
}

/// Filter and trim paired-end FASTQ files in lock-step.
///
/// Parameters
/// ----------
/// config_fwd : FilterConfig
/// config_rev : FilterConfig
/// r1_in : str
/// r2_in : str
/// r1_out : str
/// r2_out : str
///
/// Returns
/// -------
/// FilterStatsPaired
#[pyfunction]
#[pyo3(name = "filter_and_trim_paired")]
pub fn filter_and_trim_paired_py(
    py: Python<'_>,
    config_fwd: &PyFilterConfig,
    config_rev: &PyFilterConfig,
    r1_in: &str,
    r2_in: &str,
    r1_out: &str,
    r2_out: &str,
) -> PyResult<PyFilterStatsPaired> {
    let cfg_fwd = config_fwd.0.clone();
    let cfg_rev = config_rev.0.clone();
    let r1_in = r1_in.to_owned();
    let r2_in = r2_in.to_owned();
    let r1_out = r1_out.to_owned();
    let r2_out = r2_out.to_owned();
    let stats = py
        .allow_threads(move || {
            filter_and_trim_paired(
                &cfg_fwd,
                &cfg_rev,
                Path::new(&r1_in),
                Path::new(&r2_in),
                Path::new(&r1_out),
                Path::new(&r2_out),
            )
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyFilterStatsPaired {
        reads_in: stats.reads_in,
        pairs_out: stats.pairs_out,
        fwd_failed: stats.fwd_failed,
        rev_failed: stats.rev_failed,
        both_failed: stats.both_failed,
    })
}

/// Learn error rates from FASTQ files.
///
/// Parameters
/// ----------
/// fastq_paths : list[str]
/// n_reads : int, optional
///
/// Returns
/// -------
/// ErrorModel
#[pyfunction]
#[pyo3(name = "learn_errors")]
#[pyo3(signature = (fastq_paths, n_reads = 1_000_000))]
pub fn learn_errors_py(
    py: Python<'_>,
    fastq_paths: Vec<String>,
    n_reads: usize,
) -> PyResult<PyErrorModel> {
    let model = py
        .allow_threads(move || {
            let mut all_records = Vec::new();
            for path in &fastq_paths {
                let recs = read_fastq(Path::new(path))?;
                all_records.extend(recs);
                if all_records.len() >= n_reads {
                    break;
                }
            }
            let cfg = ErrorLearningConfig {
                n_reads,
                ..Default::default()
            };
            learn_errors(&all_records, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyErrorModel(model))
}

/// Dereplicate a FASTQ file.
///
/// Returns
/// -------
/// DerepResult
///     Opaque handle carrying per-position quality across the FFI boundary
///     into the dada functions. Indexable as (sequence: bytes, count: int)
///     tuples for back-compat introspection.
#[pyfunction]
#[pyo3(name = "derep_fastq")]
pub fn derep_fastq_py(py: Python<'_>, fastq_path: &str) -> PyResult<PyDerepResult> {
    let path = fastq_path.to_owned();
    let uniques = py
        .allow_threads(move || {
            let records = read_fastq(Path::new(&path))?;
            derep_fastq(&records)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyDerepResult(uniques))
}

/// Run the DADA denoising algorithm.
///
/// Parameters
/// ----------
/// derep : DerepResult
///     Output of :py:func:`derep_fastq`.
/// error_model : ErrorModel
/// omega_a : float, optional
/// pool : bool, optional
///
/// Returns
/// -------
/// DadaResult
#[pyfunction]
#[pyo3(name = "dada")]
#[pyo3(signature = (derep, error_model, omega_a = 1e-40, pool = false))]
pub fn dada_py(
    py: Python<'_>,
    derep: PyRef<PyDerepResult>,
    error_model: &PyErrorModel,
    omega_a: f64,
    pool: bool,
) -> PyResult<PyDadaResult> {
    let inner_em = error_model.0.clone();
    let uniques = derep.0.clone();
    let asvs = py
        .allow_threads(move || {
            let cfg = DadaConfig {
                omega_a,
                pool,
                ..Default::default()
            };
            dada(&uniques, &inner_em, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyDadaResult(asvs))
}

/// Run DADA on multiple samples with cross-sample pooling.
///
/// Parameters
/// ----------
/// samples : list[DerepResult]
///     One per sample (output of :py:func:`derep_fastq`).
/// error_model : ErrorModel
/// omega_a : float, optional
///
/// Returns
/// -------
/// list[DadaResult]
#[pyfunction]
#[pyo3(name = "dada_pooled")]
#[pyo3(signature = (samples, error_model, omega_a = 1e-40))]
pub fn dada_pooled_py(
    py: Python<'_>,
    samples: Vec<PyRef<PyDerepResult>>,
    error_model: &PyErrorModel,
    omega_a: f64,
) -> PyResult<Vec<PyDadaResult>> {
    let inner_em = error_model.0.clone();
    let converted: Vec<Vec<speeddada_core::derep::UniqueSeq>> =
        samples.iter().map(|d| d.0.clone()).collect();
    let results = py
        .allow_threads(move || {
            let refs: Vec<&[speeddada_core::derep::UniqueSeq]> =
                converted.iter().map(Vec::as_slice).collect();
            let cfg = DadaConfig {
                omega_a,
                ..Default::default()
            };
            dada_pooled(&refs, &inner_em, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(results.into_iter().map(PyDadaResult).collect())
}

/// Run DADA per-sample for multiple samples, parallelised across Rayon.
///
/// Parameters
/// ----------
/// samples : list[DerepResult]
/// error_model : ErrorModel
/// omega_a : float, optional
///
/// Returns
/// -------
/// list[DadaResult]
#[pyfunction]
#[pyo3(name = "dada_many")]
#[pyo3(signature = (samples, error_model, omega_a = 1e-40))]
pub fn dada_many_py(
    py: Python<'_>,
    samples: Vec<PyRef<PyDerepResult>>,
    error_model: &PyErrorModel,
    omega_a: f64,
) -> PyResult<Vec<PyDadaResult>> {
    let inner_em = error_model.0.clone();
    let converted: Vec<Vec<speeddada_core::derep::UniqueSeq>> =
        samples.iter().map(|d| d.0.clone()).collect();
    let results = py
        .allow_threads(move || {
            let refs: Vec<&[speeddada_core::derep::UniqueSeq]> =
                converted.iter().map(Vec::as_slice).collect();
            let cfg = DadaConfig {
                omega_a,
                ..Default::default()
            };
            dada_many(&refs, &inner_em, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(results.into_iter().map(PyDadaResult).collect())
}

/// Run DADA with pseudo-pooling (Callahan two-pass cross-sample scheme).
///
/// Parameters
/// ----------
/// samples : list[DerepResult]
/// error_model : ErrorModel
/// omega_a : float, optional
///
/// Returns
/// -------
/// list[DadaResult]
#[pyfunction]
#[pyo3(name = "dada_pseudo")]
#[pyo3(signature = (samples, error_model, omega_a = 1e-40))]
pub fn dada_pseudo_py(
    py: Python<'_>,
    samples: Vec<PyRef<PyDerepResult>>,
    error_model: &PyErrorModel,
    omega_a: f64,
) -> PyResult<Vec<PyDadaResult>> {
    let inner_em = error_model.0.clone();
    let converted: Vec<Vec<speeddada_core::derep::UniqueSeq>> =
        samples.iter().map(|d| d.0.clone()).collect();
    let results = py
        .allow_threads(move || {
            let refs: Vec<&[speeddada_core::derep::UniqueSeq]> =
                converted.iter().map(Vec::as_slice).collect();
            let cfg = DadaConfig {
                omega_a,
                ..Default::default()
            };
            dada_pseudo(&refs, &inner_em, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(results.into_iter().map(PyDadaResult).collect())
}

/// Merge forward and reverse DADA results.
///
/// Parameters
/// ----------
/// fwd_dada : DadaResult
/// rev_dada : DadaResult
/// min_overlap : int, optional
///
/// Returns
/// -------
/// list[MergedRead]
///     Each element has ``.sequence``, ``.abundance``, ``.accept``,
///     ``.nmatch``, ``.nmismatch``, ``.nindel`` attributes (compatible
///     with R dada2 ``mergePairs`` column names).
#[pyfunction]
#[pyo3(name = "merge_pairs")]
#[pyo3(signature = (fwd_dada, rev_dada, min_overlap = 20))]
pub fn merge_pairs_py(
    py: Python<'_>,
    fwd_dada: &PyDadaResult,
    rev_dada: &PyDadaResult,
    min_overlap: usize,
) -> PyResult<Vec<PyMergedRead>> {
    let inner_fwd = fwd_dada.0.clone();
    let inner_rev = rev_dada.0.clone();
    let merged = py
        .allow_threads(move || {
            let cfg = MergeConfig {
                min_overlap,
                ..Default::default()
            };
            merge_pairs(&inner_fwd, &inner_rev, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(merged.into_iter().map(PyMergedRead::from).collect())
}

/// Remove bimeric (chimeric) sequences.
///
/// Parameters
/// ----------
/// seqs : list[tuple[bytes, int]]
///     (sequence, abundance) pairs.
///
/// Returns
/// -------
/// list[tuple[bytes, int]]
///     Non-chimeric (sequence, abundance) pairs.
#[pyfunction]
#[pyo3(name = "remove_bimera_denovo")]
pub fn remove_bimera_denovo_py(
    py: Python<'_>,
    seqs: Vec<(Vec<u8>, u32)>,
) -> PyResult<Vec<(Vec<u8>, u32)>> {
    py.allow_threads(move || remove_bimera_denovo(&seqs))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

/// Assign taxonomy to ASV sequences using a reference FASTA.
///
/// Parameters
/// ----------
/// seqs : list[bytes]
/// ref_fasta : str
///     Path to reference FASTA (headers must include lineage separated by `;`).
/// k : int, optional
/// lineage_tsv : str | None, optional
///     Path to a TSV file with ``seq_id\\tlineage`` format.
///
/// Returns
/// -------
/// list[TaxonAssignment]
#[pyfunction]
#[pyo3(name = "assign_taxonomy")]
#[pyo3(signature = (seqs, ref_fasta, k = 8, lineage_tsv = None))]
pub fn assign_taxonomy_py(
    py: Python<'_>,
    seqs: Vec<Vec<u8>>,
    ref_fasta: &str,
    k: usize,
    lineage_tsv: Option<&str>,
) -> PyResult<Vec<PyTaxonAssignment>> {
    let ref_path = ref_fasta.to_owned();
    let tsv_path = lineage_tsv.map(str::to_owned);
    let assignments = py
        .allow_threads(move || {
            let ref_records = read_fasta(Path::new(&ref_path))?;
            let lineages: HashMap<String, Vec<String>> = if let Some(tsv) = tsv_path {
                load_lineage_tsv(Path::new(&tsv))?
            } else {
                ref_records
                    .iter()
                    .map(|r| {
                        let lin: Vec<String> = r
                            .description
                            .as_deref()
                            .unwrap_or("")
                            .split(';')
                            .map(|s| s.trim().to_owned())
                            .collect();
                        (r.id.clone(), lin)
                    })
                    .collect()
            };
            let cfg = TaxonomyConfig {
                k,
                ..Default::default()
            };
            assign_taxonomy(&seqs, &ref_records, &lineages, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    Ok(assignments
        .into_iter()
        .map(|a| PyTaxonAssignment {
            asv: a.asv,
            kingdom: a.kingdom,
            phylum: a.phylum,
            class: a.class,
            order: a.order,
            family: a.family,
            genus: a.genus,
            species: a.species,
            confidence: a.confidence,
        })
        .collect())
}

/// Run the full DADA2 pipeline on a list of input FASTQ files.
///
/// Filters, learns errors, dereplicates, denoises, and removes bimeras,
/// entirely inside a GIL-free thread.
///
/// Parameters
/// ----------
/// input_paths : list[str]
/// output_dir : str
/// trunc_len : int, optional
/// max_ee : float, optional
/// omega_a : float, optional
///
/// Returns
/// -------
/// dict[str, int]
///     Mapping from hex-encoded ASV sequence to abundance.
#[pyfunction]
#[pyo3(name = "run_pipeline")]
#[pyo3(signature = (input_paths, output_dir, trunc_len=0, max_ee=2.0, omega_a=1e-40))]
pub fn run_pipeline_py(
    py: Python<'_>,
    input_paths: Vec<String>,
    output_dir: &str,
    trunc_len: usize,
    max_ee: f64,
    omega_a: f64,
) -> PyResult<HashMap<String, u32>> {
    let out_dir = output_dir.to_owned();
    let result = py
        .allow_threads(
            move || -> Result<HashMap<String, u32>, speeddada_core::Dada2Error> {
                let out_dir_path = PathBuf::from(&out_dir);
                std::fs::create_dir_all(&out_dir_path)?;

                let filter_cfg = FilterConfig {
                    trunc_len,
                    max_ee,
                    ..Default::default()
                };

                // Filter samples in parallel.
                let filtered_paths: Vec<PathBuf> = input_paths
                    .par_iter()
                    .enumerate()
                    .map(|(i, inp)| -> Result<PathBuf, speeddada_core::Dada2Error> {
                        let out = out_dir_path.join(format!("filtered_{i}.fastq"));
                        filter_and_trim(&filter_cfg, Path::new(inp), &out)?;
                        Ok(out)
                    })
                    .collect::<Result<_, _>>()?;

                // Sample reads for error learning — cap total at n_reads regardless of
                // sample count so peak RAM stays bounded for 1000-sample workflows.
                let em_cfg = ErrorLearningConfig::default();
                let per_file =
                    RuntimeConfig::reads_per_file(em_cfg.n_reads, filtered_paths.len(), 5_000);
                let mut all_records = Vec::with_capacity(em_cfg.n_reads);
                for path in &filtered_paths {
                    all_records.extend(read_fastq_n(path, per_file)?);
                    if all_records.len() >= em_cfg.n_reads {
                        break;
                    }
                }
                let error_model = learn_errors(&all_records, &em_cfg)?;
                drop(all_records); // release before parallel phase

                // Dereplicate and denoise each sample in parallel using streaming
                // derep — no intermediate Vec<FastqRecord>, no double file-read.
                let dada_cfg = DadaConfig {
                    omega_a,
                    ..Default::default()
                };
                let all_asvs: Vec<(Vec<u8>, u32)> = filtered_paths
                    .par_iter()
                    .map(
                        |path| -> Result<Vec<(Vec<u8>, u32)>, speeddada_core::Dada2Error> {
                            let uniques = derep_fastq_path(path)?;
                            let asvs = dada(&uniques, &error_model, &dada_cfg)?;
                            Ok(asvs
                                .into_iter()
                                .map(|a| (a.sequence, a.abundance))
                                .collect())
                        },
                    )
                    .collect::<Result<Vec<_>, _>>()?
                    .into_iter()
                    .flatten()
                    .collect();

                let clean = remove_bimera_denovo(&all_asvs)?;

                let mut out_map = HashMap::new();
                for (seq, abund) in clean {
                    let hex = speeddada_core::bytes_to_hex(&seq);
                    *out_map.entry(hex).or_insert(0) += abund;
                }
                Ok(out_map)
            },
        )
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(result)
}

/// Run the full DADA2 pipeline on multiple paired-end or single-end samples.
///
/// Designed for 100s–1000s of samples. Error learning is capped at
/// `n_reads_learn` total reads (spread evenly across samples), so peak
/// RAM is bounded regardless of sample count. Each sample is processed
/// independently through filter → derep (streaming) → denoise → merge
/// (if paired) → chimera removal.
///
/// Parameters
/// ----------
/// fwd_paths : list[str]
///     Paths to forward (R1) FASTQ files, one per sample.
/// rev_paths : list[str] | None
///     Paths to reverse (R2) FASTQ files (same order). ``None`` for single-end.
/// output_dir : str
/// trunc_len_fwd : int, optional
/// trunc_len_rev : int, optional
/// max_ee_fwd : float, optional
/// max_ee_rev : float, optional
/// min_overlap : int, optional
/// omega_a : float, optional
/// n_reads_learn : int, optional
///     Total reads to use for error learning across all samples.
///
/// Returns
/// -------
/// dict[str, dict[str, int]]
///     ``{sample_name: {asv_hex: abundance}}`` for every sample.
#[pyfunction]
#[pyo3(name = "run_pipeline_samples")]
#[pyo3(signature = (
    fwd_paths, rev_paths = None, output_dir = None,
    trunc_len_fwd = 0, trunc_len_rev = 0,
    max_ee_fwd = 2.0, max_ee_rev = 5.0,
    min_overlap = 20, omega_a = 1e-40, n_reads_learn = 1_000_000
))]
#[allow(clippy::too_many_arguments)]
pub fn run_pipeline_samples_py(
    py: Python<'_>,
    fwd_paths: Vec<String>,
    rev_paths: Option<Vec<String>>,
    output_dir: Option<String>,
    trunc_len_fwd: usize,
    trunc_len_rev: usize,
    max_ee_fwd: f64,
    max_ee_rev: f64,
    min_overlap: usize,
    omega_a: f64,
    n_reads_learn: usize,
) -> PyResult<HashMap<String, HashMap<String, u32>>> {
    // Cross-platform default: <OS temp>/speeddada_out instead of a Unix-only
    // hard-coded /tmp path (would crash on Windows).
    let out_dir = output_dir.unwrap_or_else(|| {
        std::env::temp_dir()
            .join("speeddada_out")
            .to_string_lossy()
            .into_owned()
    });
    let result = py
        .allow_threads(move || {
            run_samples_inner(
                fwd_paths,
                rev_paths,
                out_dir,
                trunc_len_fwd,
                trunc_len_rev,
                max_ee_fwd,
                max_ee_rev,
                min_overlap,
                omega_a,
                n_reads_learn,
            )
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(result)
}

#[allow(clippy::too_many_arguments, clippy::needless_pass_by_value)]
fn run_samples_inner(
    fwd_paths: Vec<String>,
    rev_paths: Option<Vec<String>>,
    out_dir: String,
    trunc_len_fwd: usize,
    trunc_len_rev: usize,
    max_ee_fwd: f64,
    max_ee_rev: f64,
    min_overlap: usize,
    omega_a: f64,
    n_reads_learn: usize,
) -> Result<HashMap<String, HashMap<String, u32>>, speeddada_core::Dada2Error> {
    use std::path::PathBuf;

    let out_dir_path = PathBuf::from(&out_dir);
    std::fs::create_dir_all(&out_dir_path)?;

    let n_samples = fwd_paths.len();
    let paired = rev_paths.is_some();
    let cfg_fwd = FilterConfig {
        trunc_len: trunc_len_fwd,
        max_ee: max_ee_fwd,
        ..Default::default()
    };
    let cfg_rev = FilterConfig {
        trunc_len: trunc_len_rev,
        max_ee: max_ee_rev,
        ..Default::default()
    };

    // ── 1. Filter all samples in parallel ────────────────────────────────────
    let filtered: Vec<(String, PathBuf, Option<PathBuf>)> = fwd_paths
        .par_iter()
        .enumerate()
        .map(|(i, fwd)| -> Result<_, speeddada_core::Dada2Error> {
            let stem = Path::new(fwd)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("sample");
            let name = format!("{stem}_{i}");
            let filt_fwd = out_dir_path.join(format!("filt_fwd_{i}.fastq"));
            if paired {
                let rev = rev_paths.as_ref().unwrap()[i].as_str();
                let filt_rev = out_dir_path.join(format!("filt_rev_{i}.fastq"));
                filter_and_trim_paired(
                    &cfg_fwd,
                    &cfg_rev,
                    Path::new(fwd),
                    Path::new(rev),
                    &filt_fwd,
                    &filt_rev,
                )?;
                Ok((name, filt_fwd, Some(filt_rev)))
            } else {
                filter_and_trim(&cfg_fwd, Path::new(fwd), &filt_fwd)?;
                Ok((name, filt_fwd, None))
            }
        })
        .collect::<Result<_, _>>()?;

    // ── 2. Learn errors from a capped sample across all filtered fwd files ───
    let em_cfg = ErrorLearningConfig {
        n_reads: n_reads_learn,
        ..Default::default()
    };
    let per_file = RuntimeConfig::reads_per_file(n_reads_learn, n_samples, 5_000);
    let mut learn_records = Vec::with_capacity(n_reads_learn);
    for (_, filt_fwd, _) in &filtered {
        learn_records.extend(read_fastq_n(filt_fwd, per_file)?);
        if learn_records.len() >= n_reads_learn {
            break;
        }
    }
    let error_model = learn_errors(&learn_records, &em_cfg)?;
    drop(learn_records);

    // ── 3. Denoise each sample in parallel (streaming derep, no double-read) ─
    let merge_cfg = MergeConfig {
        min_overlap,
        ..Default::default()
    };
    let dada_cfg = DadaConfig {
        omega_a,
        ..Default::default()
    };

    let results: HashMap<String, HashMap<String, u32>> = filtered
        .par_iter()
        .map(
            |(name, filt_fwd, filt_rev)| -> Result<_, speeddada_core::Dada2Error> {
                let fwd_uniq = derep_fastq_path(filt_fwd)?;
                let fwd_asvs = dada(&fwd_uniq, &error_model, &dada_cfg)?;

                let pairs: Vec<(Vec<u8>, u32)> = if let Some(rev_path) = filt_rev {
                    let rev_uniq = derep_fastq_path(rev_path)?;
                    let rev_asvs = dada(&rev_uniq, &error_model, &dada_cfg)?;
                    merge_pairs(&fwd_asvs, &rev_asvs, &merge_cfg)?
                        .into_iter()
                        .map(|m| (m.sequence, m.abundance))
                        .collect()
                } else {
                    fwd_asvs
                        .into_iter()
                        .map(|a| (a.sequence, a.abundance))
                        .collect()
                };

                let clean = remove_bimera_denovo(&pairs)?;
                let mut sample_map = HashMap::new();
                for (seq, abund) in clean {
                    let hex = speeddada_core::bytes_to_hex(&seq);
                    *sample_map.entry(hex).or_insert(0u32) += abund;
                }
                Ok((name.clone(), sample_map))
            },
        )
        .collect::<Result<_, _>>()?;

    Ok(results)
}

/// Initialise Rust-level logging. Call once at startup.
///
/// Parameters
/// ----------
/// level : str, optional
///     Log level: ``"error"``, ``"warn"``, ``"info"`` (default), ``"debug"``, ``"trace"``.
#[pyfunction]
#[pyo3(name = "init_logging")]
#[pyo3(signature = (level = "info"))]
pub fn init_logging_py(level: &str) {
    // SAFETY: single-threaded at import time; env var is read by env_logger::try_init.
    #[allow(unsafe_code)]
    unsafe {
        std::env::set_var("RUST_LOG", level);
    }
    env_logger::try_init().ok();
}
