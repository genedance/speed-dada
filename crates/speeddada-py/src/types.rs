//! PyO3 class wrappers for dada2-core types.

use dada2_core::{
    dada::Asv,
    derep::UniqueSeq,
    error_model::ErrorModel,
    filter::FilterConfig,
    merge::MergedRead,
    sequence_table::SequenceTable,
};
use pyo3::prelude::*;
use std::collections::HashMap;

// ── FilterConfig ─────────────────────────────────────────────────────────────

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
pub struct PyFilterConfig(pub FilterConfig);

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

// ── FilterStats ──────────────────────────────────────────────────────────────

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

// ── FilterStatsPaired ────────────────────────────────────────────────────────

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

// ── ErrorModel ───────────────────────────────────────────────────────────────

/// Learned parametric error model.
///
/// Methods
/// -------
/// plot_errors() -> dict
///     Return a dict with keys ``"quality"`` and ``"error_rates"`` suitable
///     for plotting with matplotlib.
#[pyclass(name = "ErrorModel")]
pub struct PyErrorModel(pub ErrorModel);

#[pymethods]
impl PyErrorModel {
    /// Return error rate data for plotting.
    fn plot_errors(&self) -> HashMap<String, Vec<f64>> {
        let mut quality = Vec::new();
        let mut error_rates = Vec::new();
        for q in 0..dada2_core::error_model::MAX_QUAL {
            #[allow(clippy::cast_precision_loss)]
            quality.push(q as f64);
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

// ── DadaResult ───────────────────────────────────────────────────────────────

/// Result of the DADA denoising step.
///
/// Behaves like a list of (sequence: bytes, abundance: int) tuples.
#[pyclass(name = "DadaResult")]
pub struct PyDadaResult(pub Vec<Asv>);

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

// ── MergedRead ───────────────────────────────────────────────────────────────

/// A single merged paired-end read.
///
/// Attributes
/// ----------
/// sequence : bytes
/// abundance : int
/// accept : bool
///     Always True (rejected reads are not returned).
/// nmatch : int
///     Number of matching bases in the overlap region.
/// nmismatch : int
///     Number of mismatching bases in the overlap region.
/// nindel : int
///     Number of indels in the overlap (always 0 for this aligner).
#[pyclass(name = "MergedRead")]
pub struct PyMergedRead {
    #[pyo3(get)] pub sequence: Vec<u8>,
    #[pyo3(get)] pub abundance: u32,
    #[pyo3(get)] pub accept: bool,
    #[pyo3(get)] pub nmatch: usize,
    #[pyo3(get)] pub nmismatch: u32,
    #[pyo3(get)] pub nindel: u32,
}

impl From<MergedRead> for PyMergedRead {
    fn from(m: MergedRead) -> Self {
        Self {
            sequence: m.sequence,
            abundance: m.abundance,
            accept: true,
            nmatch: m.overlap_len,
            nmismatch: m.n_mismatches,
            nindel: 0,
        }
    }
}

#[pymethods]
impl PyMergedRead {
    fn __repr__(&self) -> String {
        format!("MergedRead(abundance={}, nmatch={}, nmismatch={})", self.abundance, self.nmatch, self.nmismatch)
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

// ── SequenceTable ────────────────────────────────────────────────────────────

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
pub struct PySequenceTable(pub SequenceTable);

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

// ── QualityProfile ───────────────────────────────────────────────────────────

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

// ── DerepResult ──────────────────────────────────────────────────────────────

/// Dereplicated FASTQ result — an opaque, Rust-owned container.
///
/// The previous Python representation (`list[tuple[bytes, int]]`) lost the
/// per-position quality summary that the DADA algorithm uses to estimate
/// per-base error rates, forcing the dada* functions to fabricate a flat
/// Phred-30 quality when reconstructing `UniqueSeq` on the Rust side.
/// `DerepResult` carries the full `Vec<UniqueSeq>` (including `qual_sum`)
/// across the FFI boundary so denoising uses real quality.
///
/// For back-compat, the class is still iterable / indexable as
/// `(sequence: bytes, abundance: int)` tuples.
#[pyclass(name = "DerepResult")]
pub struct PyDerepResult(pub Vec<UniqueSeq>);

#[pymethods]
impl PyDerepResult {
    fn __len__(&self) -> usize {
        self.0.len()
    }
    fn __getitem__(&self, i: usize) -> PyResult<(Vec<u8>, u32)> {
        let u = self
            .0
            .get(i)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("index out of range"))?;
        Ok((u.seq.clone(), u.count))
    }
    fn __repr__(&self) -> String {
        format!("DerepResult(n_uniques={})", self.0.len())
    }
}
