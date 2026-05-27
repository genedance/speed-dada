"""speeddada — high-performance DADA2 amplicon sequence variant pipeline.

A drop-in Python interface to the Rust-backed DADA2 pipeline. Re-exports the
compiled :mod:`speeddada._native` PyO3 module so end users write::

    import speeddada
    cfg = speeddada.FilterConfig(trunc_len=240, max_ee=2.0)
    speeddada.filter_and_trim(cfg, "R1.fastq.gz", "R1.filt.fastq.gz")

See the README and the documentation site for full examples.
"""
from __future__ import annotations

from ._native import (
    # Configuration / stats classes
    FilterConfig,
    FilterStats,
    FilterStatsPaired,
    QualityProfile,
    SequenceTable,
    ErrorModel,
    DadaResult,
    DerepResult,
    TaxonAssignment,
    MergedRead,
    # Pipeline functions
    configure_runtime,
    init_logging,
    quality_profile,
    trim_primers,
    filter_and_trim,
    filter_and_trim_paired,
    learn_errors,
    derep_fastq,
    dada,
    dada_many,
    dada_pooled,
    dada_pseudo,
    merge_pairs,
    remove_bimera_denovo,
    assign_taxonomy,
    make_sequence_table,
    run_pipeline,
    run_pipeline_samples,
    # Version
    __version__,
)

__all__ = [
    "FilterConfig",
    "FilterStats",
    "FilterStatsPaired",
    "QualityProfile",
    "SequenceTable",
    "ErrorModel",
    "DadaResult",
    "DerepResult",
    "TaxonAssignment",
    "MergedRead",
    "configure_runtime",
    "init_logging",
    "quality_profile",
    "trim_primers",
    "filter_and_trim",
    "filter_and_trim_paired",
    "learn_errors",
    "derep_fastq",
    "dada",
    "dada_many",
    "dada_pooled",
    "dada_pseudo",
    "merge_pairs",
    "remove_bimera_denovo",
    "assign_taxonomy",
    "make_sequence_table",
    "run_pipeline",
    "run_pipeline_samples",
    "__version__",
]
