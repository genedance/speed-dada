//! PyO3 bindings for dada2-core.
//!
//! Exposes the full pipeline as a Python module `dada2`.

use dada2_core::{
    chimera::remove_bimeras,
    dada::{run_dada, Asv, DadaConfig},
    derep::dereplicate,
    error_model::{learn_errors, ErrorLearningConfig, ErrorModel},
    filter::{filter_and_trim, FilterConfig},
    io::{fasta::read_fasta, fastq::read_fastq},
    merge::{merge_pairs, MergeConfig},
    taxonomy::{assign_taxonomy, TaxonomyConfig},
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

// ── Public pipeline functions ────────────────────────────────────────────────

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
#[pyo3(name = "dereplicate")]
fn dereplicate_py(py: Python<'_>, fastq_path: &str) -> PyResult<Vec<(Vec<u8>, u32)>> {
    let path = fastq_path.to_owned();
    let uniques = py
        .allow_threads(move || {
            let records = read_fastq(Path::new(&path))?;
            dereplicate(&records)
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
#[pyo3(name = "run_dada")]
#[pyo3(signature = (derep, error_model, omega_a = 1e-40, pool = false))]
fn run_dada_py(
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
            run_dada(&uniques, &inner_em, &cfg)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(PyDadaResult(asvs))
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
#[pyo3(name = "remove_bimeras")]
fn remove_bimeras_py(py: Python<'_>, seqs: Vec<(Vec<u8>, u32)>) -> PyResult<Vec<Vec<u8>>> {
    let result = py
        .allow_threads(move || remove_bimeras(&seqs))
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
///
/// Returns
/// -------
/// list[TaxonAssignment]
#[pyfunction]
#[pyo3(name = "assign_taxonomy")]
#[pyo3(signature = (seqs, ref_fasta, k = 8))]
fn assign_taxonomy_py(
    py: Python<'_>,
    seqs: Vec<Vec<u8>>,
    ref_fasta: &str,
    k: usize,
) -> PyResult<Vec<PyTaxonAssignment>> {
    let ref_path = ref_fasta.to_owned();
    let assignments = py
        .allow_threads(move || {
            let ref_records = read_fasta(Path::new(&ref_path))?;
            let lineages: HashMap<String, Vec<String>> = ref_records
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
                .collect();
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
                let uniques = dereplicate(&records)?;
                let asvs = run_dada(&uniques, &error_model, &dada_cfg)?;
                for asv in asvs {
                    all_asvs.push((asv.sequence, asv.abundance));
                }
            }

            // Remove bimeras from pooled ASVs
            let clean = remove_bimeras(&all_asvs)?;

            // Return as hex-encoded sequence → abundance map
            let mut out_map = HashMap::new();
            for (seq, abund) in clean {
                let hex = seq.iter().map(|b| format!("{b:02x}")).collect::<String>();
                *out_map.entry(hex).or_insert(0) += abund;
            }
            Ok(out_map)
        })
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(result)
}

/// dada2 — high-performance ASV pipeline (Rust core).
#[pymodule]
fn dada2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(filter_and_trim_py, m)?)?;
    m.add_function(wrap_pyfunction!(learn_errors_py, m)?)?;
    m.add_function(wrap_pyfunction!(dereplicate_py, m)?)?;
    m.add_function(wrap_pyfunction!(run_dada_py, m)?)?;
    m.add_function(wrap_pyfunction!(merge_pairs_py, m)?)?;
    m.add_function(wrap_pyfunction!(remove_bimeras_py, m)?)?;
    m.add_function(wrap_pyfunction!(assign_taxonomy_py, m)?)?;
    m.add_function(wrap_pyfunction!(run_pipeline_py, m)?)?;

    m.add_class::<PyFilterConfig>()?;
    m.add_class::<PyFilterStats>()?;
    m.add_class::<PyErrorModel>()?;
    m.add_class::<PyDadaResult>()?;
    m.add_class::<PyTaxonAssignment>()?;

    m.add("__version__", version())?;
    Ok(())
}
