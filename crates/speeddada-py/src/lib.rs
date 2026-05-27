//! `PyO3` bindings for dada2-core.
//!
//! Exposes the full pipeline as a Python module `dada2`.

// Python docstrings use snake_case parameter names that trigger doc_markdown.
#![allow(clippy::doc_markdown)]

pub mod functions;
pub mod types;

use functions::{
    assign_taxonomy_py, configure_runtime_py, dada_many_py, dada_pooled_py,
    dada_pseudo_py, dada_py, derep_fastq_py, filter_and_trim_paired_py,
    filter_and_trim_py, init_logging_py, learn_errors_py, make_sequence_table_py,
    merge_pairs_py, quality_profile_py, remove_bimera_denovo_py, run_pipeline_py,
    run_pipeline_samples_py, trim_primers_py, version,
};
use types::{
    PyDadaResult, PyDerepResult, PyErrorModel, PyFilterConfig, PyFilterStats,
    PyFilterStatsPaired, PyMergedRead, PyQualityProfile, PySequenceTable,
    PyTaxonAssignment,
};
use pyo3::prelude::*;

/// dada2 — high-performance ASV pipeline (Rust core).
#[pymodule]
fn dada2(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(configure_runtime_py, m)?)?;
    m.add_function(wrap_pyfunction!(init_logging_py, m)?)?;
    m.add_function(wrap_pyfunction!(make_sequence_table_py, m)?)?;
    m.add_function(wrap_pyfunction!(quality_profile_py, m)?)?;
    m.add_function(wrap_pyfunction!(trim_primers_py, m)?)?;
    m.add_function(wrap_pyfunction!(filter_and_trim_py, m)?)?;
    m.add_function(wrap_pyfunction!(filter_and_trim_paired_py, m)?)?;
    m.add_function(wrap_pyfunction!(learn_errors_py, m)?)?;
    m.add_function(wrap_pyfunction!(derep_fastq_py, m)?)?;
    m.add_function(wrap_pyfunction!(dada_py, m)?)?;
    m.add_function(wrap_pyfunction!(dada_many_py, m)?)?;
    m.add_function(wrap_pyfunction!(dada_pooled_py, m)?)?;
    m.add_function(wrap_pyfunction!(dada_pseudo_py, m)?)?;
    m.add_function(wrap_pyfunction!(merge_pairs_py, m)?)?;
    m.add_function(wrap_pyfunction!(remove_bimera_denovo_py, m)?)?;
    m.add_function(wrap_pyfunction!(assign_taxonomy_py, m)?)?;
    m.add_function(wrap_pyfunction!(run_pipeline_py, m)?)?;
    m.add_function(wrap_pyfunction!(run_pipeline_samples_py, m)?)?;

    m.add_class::<PyFilterConfig>()?;
    m.add_class::<PyFilterStats>()?;
    m.add_class::<PyFilterStatsPaired>()?;
    m.add_class::<PyQualityProfile>()?;
    m.add_class::<PySequenceTable>()?;
    m.add_class::<PyErrorModel>()?;
    m.add_class::<PyDadaResult>()?;
    m.add_class::<PyDerepResult>()?;
    m.add_class::<PyTaxonAssignment>()?;
    m.add_class::<PyMergedRead>()?;

    m.add("__version__", version())?;
    Ok(())
}
