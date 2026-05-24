//! `PyO3` bindings for dada2-core.
//!
//! Exposes the full pipeline as a Python module `dada2`.

// Python docstrings use snake_case parameter names that trigger doc_markdown.
// The remaining lints are suppressed where intentionally used.
#![allow(clippy::doc_markdown)]

use dada2_core::{
    chimera::remove_bimera_denovo,
    dada::{dada, dada_pooled, Asv, DadaConfig},
    derep::derep_fastq,
    error_model::{learn_errors, ErrorLearningConfig, ErrorModel},
    filter::{filter_and_trim, filter_and_trim_paired, FilterConfig},
    io::{fasta::read_fasta, fastq::read_fastq},
    merge::{merge_pairs, MergeConfig},
    primer::{trim_primers, PrimerConfig},
    quality_profile::quality_profile,
    sequence_table::SequenceTable,
    taxonomy::{assign_taxonomy, load_lineage_tsv, TaxonomyConfig},
};
use pyo3::prelude::*;
use std::{collections::HashMap, path::Path, path::PathBuf};

/// Return the crate version string.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

// ── FilterConfig ────────────────────────────────────────────────────────────

/// Configuration for the filter-and-trim stage.
///
/// Parameters
/// ----------
/// trunc_len : int
///     Truncate reads to this length (0 = no truncation).
/// min_len : int
///     Discard reads shorter than this after truncation.
/// max_ee : float
///     Maximum expected errors per read.
/// trunc_q : int
///     Truncate reads at the first base with quality below this threshold.
/// trim_left : int
///     Remove this many bases from the 5' end.
/// trim_right : int
///     Remove this many bases from the 3' end.
#[pyclass(name = "FilterConfig")]
#[derive(Clone)]
pub struct PyFilterConfig(FilterConfig);

#[pymethods]
impl PyFilterConfig {
    #[new]
    #[pyo3(signature = (trunc_len=0, min_len=20, max_ee=2.0, trunc_q=2, trim_left=0, trim_right=0))]
    fn new(
        trunc_len: usize,
        min_len: usize,
        max_ee: f64,
        trunc_q: u8,
        trim_left: usize,
        trim_right: usize,
    ) -> Self {
        Self(FilterConfig { trunc_len, min_len, max_ee, trunc_q, trim_left, trim_right })
    }

    fn __repr__(&self) -> String {
        format!(
            "FilterConfig(trunc_len={}, min_len={}, max_ee={}, trunc_q={})",
            self.0.trunc_len, self.0.min_len, self.0.max_ee, self.0.trunc_q
        )
    }
}

// ── FilterStats ─────────────────────────────────────────────────────────────

/// Statistics returned by filter_and_trim.
///
/// Attributes
/// ----------
/// reads_in : int
/// reads_out : int
#[pyclass(name = "FilterStats")]
pub struct PyFilterStats {
    #[pyo3(get)]
    pub reads_in: u64,
    #[pyo3(get)]
    pub reads_out: u64,
}

// ── FilterStatsPaired ───────────────────────────────────────────────────────

/// Statistics for paired-end filtering.
///
/// Attributes
/// ----------
/// reads_in : int
/// pairs_out : int
/// fwd_failed : int
/// rev_failed : int
/// both_failed : int
#[pyclass(name = "FilterStatsPaired")]
pub struct PyFilterStatsPaired {
    #[pyo3(get)]
    pub reads_in: u64,
    #[pyo3(get)]
    pub pairs_out: u64,
    #[pyo3(get)]
    pub fwd_failed: u64,
    #[pyo3(get)]
    pub rev_failed: u64,
    #[pyo3(get)]
    pub both_failed: u64,
}

// ── ErrorModel ──────────────────────────────────────────────────────────────

/// Learned parametric error model.
///
/// Methods
/// -------
/// plot_errors() -> dict
///     Return a dict with keys ``"quality"`` and ``"error_rates"`` suitable
///     for plotting with matplotlib.
#[pyclass(name = "ErrorModel")]
pub struct PyErrorModel(ErrorModel);

#[pymethods]
impl PyErrorModel {
    /// Return error rate data for plotting.
    fn plot_errors(&self) -> HashMap<String, Vec<f64>> {
        let mut quality = Vec::new();
        let mut error_rates = Vec::new();
        for q in 0..dada2_core::error_model::MAX_QUAL {
            #[allow(clippy::cast_precision_loss)]
            quality.push(q as f64);
            // Mean mismatch rate across all 12 substitution classes
            let rate: f64 = (0..16)
                .filter(|&r| r % 5 != 0)
                .map(|r| self.0.matrix[[r, q]])
                .sum::<f64>()
                / 12.0;
            error_rates.push(rate);
        }
        let mut m = HashMap::new();
        m.insert("quality".into(), quality);
        m.insert("error_rates".into(), error_rates);
        m
    }

    fn __repr__(&self) -> String {
        format!("ErrorModel(n_reads_used={})", self.0.n_reads_used)
    }
}

// ── DadaResult ──────────────────────────────────────────────────────────────

/// Result of the DADA denoising step.
///
/// Behaves like a list of (sequence: bytes, abundance: int) tuples.
#[pyclass(name = "DadaResult")]
pub struct PyDadaResult(Vec<Asv>);

#[pymethods]
impl PyDadaResult {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __getitem__(&self, idx: usize) -> PyResult<(Vec<u8>, u32)> {
        self.0
            .get(idx)
            .map(|a| (a.sequence.clone(), a.abundance))
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("index out of range"))
    }

    fn __repr__(&self) -> String {
        format!("DadaResult(n_asvs={})", self.0.len())
    }
}

// ── TaxonAssignment ──────────────────────────────────────────────────────────

/// Taxonomic assignment for a single ASV.
///
/// Attributes
/// ----------
/// asv : str, kingdom : str | None, phylum : str | None, ..., confidence : float
#[pyclass(name = "TaxonAssignment")]
pub struct PyTaxonAssignment {
    #[pyo3(get)] pub asv: String,
    #[pyo3(get)] pub kingdom: Option<String>,
    #[pyo3(get)] pub phylum: Option<String>,
    #[pyo3(get)] pub class: Option<String>,
    #[pyo3(get)] pub order: Option<String>,
    #[pyo3(get)] pub family: Option<String>,
    #[pyo3(get)] pub genus: Option<String>,
    #[pyo3(get)] pub species: Option<String>,
    #[pyo3(get)] pub confidence: f64,
}

#[pymethods]
impl PyTaxonAssignment {
    fn __repr__(&self) -> String {
        format!(
            "TaxonAssignment(genus={:?}, confidence={:.2})",
            self.genus, self.confidence
        )
    }
}

// ── SequenceTable ───────────────────────────────────────────────────────────

/// Sample × ASV abundance matrix.
///
/// Attributes
/// ----------
/// samples : list[str]
/// sequences : list[str]
///     Hex-encoded ASV sequences.
/// counts : list[list[int]]
///     ``counts[sample_idx][asv_idx]`` = abundance.
#[pyclass(name = "SequenceTable")]
pub struct PySequenceTable(SequenceTable);

#[pymethods]
impl PySequenceTable {
    /// Sample names.
    #[getter]
    fn samples(&self) -> Vec<String> {
        self.0.samples.clone()
    }

    /// Hex-encoded ASV sequences.
    #[getter]
    fn sequences(&self) -> Vec<String> {
        self.0.sequences.clone()
    }

    /// Abundance matrix.
    #[getter]
    fn counts(&self) -> Vec<Vec<u32>> {
        self.0.counts.clone()
    }

    /// Write as TSV to *path*.
    fn to_tsv(&self, path: &str) -> PyResult<()> {
        self.0
            .to_tsv(std::path::Path::new(path))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    /// Serialise to JSON string.
    fn to_json(&self) -> PyResult<String> {
        self.0
            .to_json()
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "SequenceTable(samples={}, sequences={})",
            self.0.samples.len(),
            self.0.sequences.len()
        )
    }
}

// ── QualityProfile ──────────────────────────────────────────────────────────

/// Per-cycle quality summary statistics for a FASTQ file.
///
/// Attributes
/// ----------
/// n_reads : int
/// cycle_mean : list[float]
/// cycle_p25 : list[float]
/// cycle_p50 : list[float]
/// cycle_p75 : list[float]
/// cycle_count : list[int]
#[pyclass(name = "QualityProfile")]
pub struct PyQualityProfile {
    #[pyo3(get)]
    pub n_reads: u64,
    #[pyo3(get)]
    pub cycle_mean: Vec<f64>,
    #[pyo3(get)]
    pub cycle_p25: Vec<f64>,
    #[pyo3(get)]
    pub cycle_p50: Vec<f64>,
    #[pyo3(get)]
    pub cycle_p75: Vec<f64>,
    #[pyo3(get)]
    pub cycle_count: Vec<u64>,
}

#[pymethods]
impl PyQualityProfile {
    fn __repr__(&self) -> String {
        format!(
            "QualityProfile(n_reads={}, n_cycles={})",
            self.n_reads,
            self.cycle_mean.len()
        )
    }
}

// ── Public pipeline functions ────────────────────────────────────────────────

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
fn make_sequence_table_py(
    sample_names: Vec<String>,
    results: Vec<PyRef<'_, PyDadaResult>>,
) -> PySequenceTable {
    let name_refs: Vec<&str> = sample_names.iter().map(String::as_str).collect();
    let asvs: Vec<Vec<Asv>> = results.iter().map(|r| r.0.clone()).collect();
    let table = SequenceTable::new(&name_refs, &asvs);
    PySequenceTable(table)
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
fn trim_primers_py(
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
            let cfg = PrimerConfig { fwd_primer, rev_primer, max_mismatches, min_overlap };
            trim_primers(&cfg, std::path::Path::new(&inp), std::path::Path::new(&out))
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyFilterStats { reads_in: stats.reads_in, reads_out: stats.reads_out })
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
fn quality_profile_py(
    py: Python<'_>,
    fastq_path: &str,
    n_reads: usize,
) -> PyResult<PyQualityProfile> {
    let path = fastq_path.to_owned();
    let profile = py
        .allow_threads(move || quality_profile(std::path::Path::new(&path), n_reads))
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
fn filter_and_trim_py(
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
    Ok(PyFilterStats { reads_in: stats.reads_in, reads_out: stats.reads_out })
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
fn filter_and_trim_paired_py(
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
                std::path::Path::new(&r1_in),
                std::path::Path::new(&r2_in),
                std::path::Path::new(&r1_out),
                std::path::Path::new(&r2_out),
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
fn learn_errors_py(
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
            let cfg = ErrorLearningConfig { n_reads, ..Default::default() };
            learn_errors(&all_records, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyErrorModel(model))
}

/// Dereplicate a FASTQ file.
///
/// Returns
/// -------
/// list[tuple[bytes, int]]
#[pyfunction]
#[pyo3(name = "derep_fastq")]
fn derep_fastq_py(py: Python<'_>, fastq_path: &str) -> PyResult<Vec<(Vec<u8>, u32)>> {
    let path = fastq_path.to_owned();
    let uniques = py
        .allow_threads(move || {
            let records = read_fastq(Path::new(&path))?;
            derep_fastq(&records)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(uniques.into_iter().map(|u| (u.seq, u.count)).collect())
}

/// Run the DADA denoising algorithm.
///
/// Parameters
/// ----------
/// derep : list[tuple[bytes, int]]
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
fn dada_py(
    py: Python<'_>,
    derep: Vec<(Vec<u8>, u32)>,
    error_model: &PyErrorModel,
    omega_a: f64,
    pool: bool,
) -> PyResult<PyDadaResult> {
    let inner_em = error_model.0.clone();
    let asvs = py
        .allow_threads(move || {
            let uniques: Vec<dada2_core::derep::UniqueSeq> = derep
                .into_iter()
                .map(|(seq, count)| {
                    let qual_sum = vec![30.0 * f64::from(count); seq.len()];
                    dada2_core::derep::UniqueSeq { seq, count, qual_sum }
                })
                .collect();
            let cfg = DadaConfig { omega_a, pool, ..Default::default() };
            dada(&uniques, &inner_em, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyDadaResult(asvs))
}

/// Run DADA on multiple samples with cross-sample pooling.
///
/// Parameters
/// ----------
/// samples : list[list[tuple[bytes, int]]]
///     One dereplicated sample per element (output of ``dereplicate``).
/// error_model : ErrorModel
/// omega_a : float, optional
///
/// Returns
/// -------
/// list[DadaResult]
///     One result per input sample.
#[pyfunction]
#[pyo3(name = "dada_pooled")]
#[pyo3(signature = (samples, error_model, omega_a = 1e-40))]
fn dada_pooled_py(
    py: Python<'_>,
    samples: Vec<Vec<(Vec<u8>, u32)>>,
    error_model: &PyErrorModel,
    omega_a: f64,
) -> PyResult<Vec<PyDadaResult>> {
    let inner_em = error_model.0.clone();
    let results = py
        .allow_threads(move || {
            let converted: Vec<Vec<dada2_core::derep::UniqueSeq>> = samples
                .into_iter()
                .map(|s| {
                    s.into_iter()
                        .map(|(seq, count)| {
                            let qual_sum = vec![30.0 * f64::from(count); seq.len()];
                            dada2_core::derep::UniqueSeq { seq, count, qual_sum }
                        })
                        .collect()
                })
                .collect();
            let refs: Vec<&[dada2_core::derep::UniqueSeq]> =
                converted.iter().map(Vec::as_slice).collect();
            let cfg = DadaConfig { omega_a, ..Default::default() };
            dada_pooled(&refs, &inner_em, &cfg)
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
/// list[tuple[bytes, int]]
#[pyfunction]
#[pyo3(name = "merge_pairs")]
#[pyo3(signature = (fwd_dada, rev_dada, min_overlap = 20))]
fn merge_pairs_py(
    py: Python<'_>,
    fwd_dada: &PyDadaResult,
    rev_dada: &PyDadaResult,
    min_overlap: usize,
) -> PyResult<Vec<(Vec<u8>, u32)>> {
    let inner_fwd = fwd_dada.0.clone();
    let inner_rev = rev_dada.0.clone();
    let merged = py
        .allow_threads(move || {
            let cfg = MergeConfig { min_overlap, ..Default::default() };
            merge_pairs(&inner_fwd, &inner_rev, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(merged.into_iter().map(|m| (m.sequence, m.abundance)).collect())
}

/// Remove bimeric (chimeric) sequences.
///
/// Parameters
/// ----------
/// seqs : list[tuple[bytes, int]]
///
/// Returns
/// -------
/// list[bytes]
#[pyfunction]
#[pyo3(name = "remove_bimera_denovo")]
fn remove_bimera_denovo_py(py: Python<'_>, seqs: Vec<(Vec<u8>, u32)>) -> PyResult<Vec<Vec<u8>>> {
    let result = py
        .allow_threads(move || remove_bimera_denovo(&seqs))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(result.into_iter().map(|(s, _)| s).collect())
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
///     Path to a TSV file with ``seq_id\\tlineage`` format. When provided,
///     lineages are loaded from this file instead of the FASTA headers.
///
/// Returns
/// -------
/// list[TaxonAssignment]
#[pyfunction]
#[pyo3(name = "assign_taxonomy")]
#[pyo3(signature = (seqs, ref_fasta, k = 8, lineage_tsv = None))]
fn assign_taxonomy_py(
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
            let cfg = TaxonomyConfig { k, ..Default::default() };
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
fn run_pipeline_py(
    py: Python<'_>,
    input_paths: Vec<String>,
    output_dir: &str,
    trunc_len: usize,
    max_ee: f64,
    omega_a: f64,
) -> PyResult<HashMap<String, u32>> {
    let out_dir = output_dir.to_owned();
    let result = py
        .allow_threads(move || -> Result<HashMap<String, u32>, dada2_core::Dada2Error> {
            let out_dir_path = PathBuf::from(&out_dir);
            std::fs::create_dir_all(&out_dir_path)?;

            let filter_cfg = FilterConfig {
                trunc_len,
                max_ee,
                ..Default::default()
            };

            // Filter all files and collect output paths
            let mut filtered_paths: Vec<PathBuf> = Vec::with_capacity(input_paths.len());
            for (i, inp) in input_paths.iter().enumerate() {
                let out = out_dir_path.join(format!("filtered_{i}.fastq"));
                filter_and_trim(&filter_cfg, Path::new(inp), &out)?;
                filtered_paths.push(out);
            }

            // Learn error model from filtered files
            let mut all_records = Vec::new();
            for path in &filtered_paths {
                let recs = read_fastq(path)?;
                all_records.extend(recs);
            }
            let em_cfg = ErrorLearningConfig::default();
            let error_model = learn_errors(&all_records, &em_cfg)?;

            // Dereplicate, denoise, collect all ASVs
            let dada_cfg = DadaConfig { omega_a, ..Default::default() };
            let mut all_asvs: Vec<(Vec<u8>, u32)> = Vec::new();
            for path in &filtered_paths {
                let records = read_fastq(path)?;
                let uniques = derep_fastq(&records)?;
                let asvs = dada(&uniques, &error_model, &dada_cfg)?;
                for asv in asvs {
                    all_asvs.push((asv.sequence, asv.abundance));
                }
            }

            // Remove bimeras from pooled ASVs
            let clean = remove_bimera_denovo(&all_asvs)?;

            // Return as hex-encoded sequence → abundance map
            let mut out_map = HashMap::new();
            for (seq, abund) in clean {
                use std::fmt::Write as _;
                let hex = seq.iter().fold(String::with_capacity(seq.len() * 2), |mut acc, b| {
                    let _ = write!(acc, "{b:02x}");
                    acc
                });
                *out_map.entry(hex).or_insert(0) += abund;
            }
            Ok(out_map)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(result)
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
fn init_logging_py(level: &str) {
    // SAFETY: single-threaded at import time; env var is read by env_logger::try_init.
    #[allow(unsafe_code)]
    unsafe {
        std::env::set_var("RUST_LOG", level);
    }
    env_logger::try_init().ok(); // ok() ignores already-initialised error
}

/// dada2 — high-performance ASV pipeline (Rust core).
#[pymodule]
fn dada2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(init_logging_py, m)?)?;
    m.add_function(wrap_pyfunction!(make_sequence_table_py, m)?)?;
    m.add_function(wrap_pyfunction!(quality_profile_py, m)?)?;
    m.add_function(wrap_pyfunction!(trim_primers_py, m)?)?;
    m.add_function(wrap_pyfunction!(filter_and_trim_py, m)?)?;
    m.add_function(wrap_pyfunction!(filter_and_trim_paired_py, m)?)?;
    m.add_function(wrap_pyfunction!(learn_errors_py, m)?)?;
    m.add_function(wrap_pyfunction!(derep_fastq_py, m)?)?;
    m.add_function(wrap_pyfunction!(dada_py, m)?)?;
    m.add_function(wrap_pyfunction!(dada_pooled_py, m)?)?;
    m.add_function(wrap_pyfunction!(merge_pairs_py, m)?)?;
    m.add_function(wrap_pyfunction!(remove_bimera_denovo_py, m)?)?;
    m.add_function(wrap_pyfunction!(assign_taxonomy_py, m)?)?;
    m.add_function(wrap_pyfunction!(run_pipeline_py, m)?)?;

    m.add_class::<PyFilterConfig>()?;
    m.add_class::<PyFilterStats>()?;
    m.add_class::<PyFilterStatsPaired>()?;
    m.add_class::<PyQualityProfile>()?;
    m.add_class::<PySequenceTable>()?;
    m.add_class::<PyErrorModel>()?;
    m.add_class::<PyDadaResult>()?;
    m.add_class::<PyTaxonAssignment>()?;

    m.add("__version__", version())?;
    Ok(())
}
