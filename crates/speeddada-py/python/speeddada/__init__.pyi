"""Type stubs for speeddada — see speeddada/__init__.py for runtime re-exports."""
from __future__ import annotations

from typing import Any, Sequence

__version__: str

# ── Configuration / data classes ─────────────────────────────────────────────

class FilterConfig:
    def __init__(
        self,
        trunc_len: int = ...,
        min_len: int = ...,
        max_ee: float = ...,
        trunc_q: int = ...,
        trim_left: int = ...,
        trim_right: int = ...,
    ) -> None: ...

class FilterStats:
    reads_in: int
    reads_out: int

class FilterStatsPaired:
    reads_in: int
    pairs_out: int
    fwd_failed: int
    rev_failed: int
    both_failed: int

class QualityProfile:
    def __len__(self) -> int: ...
    def mean(self, position: int) -> float: ...

class SequenceTable:
    n_samples: int
    n_asvs: int
    samples: list[str]
    asvs: list[str]
    def counts(self) -> list[list[int]]: ...

class ErrorModel:
    """Opaque substitution-error model learned from FASTQ data."""

class DadaResult:
    """Per-sample DADA denoising output. Iterable as (sequence, abundance) tuples."""
    def __len__(self) -> int: ...
    def __getitem__(self, idx: int) -> tuple[str, int]: ...

class DerepResult:
    """Dereplicated sample (sequence, count) pairs."""
    def __len__(self) -> int: ...
    def __getitem__(self, idx: int) -> tuple[str, int]: ...

class TaxonAssignment:
    sequence: str
    lineage: list[str]
    bootstrap: list[float]

class MergedRead:
    sequence: str
    abundance: int
    overlap_len: int
    n_mismatches: int

# ── Pipeline functions ───────────────────────────────────────────────────────

def configure_runtime(n_threads: int = ...) -> None: ...
def init_logging(level: str = ...) -> None: ...
def quality_profile(path: str, n_reads: int = ...) -> QualityProfile: ...
def trim_primers(
    in_path: str,
    out_path: str,
    fwd_primer: bytes,
    rev_primer: bytes,
    max_mismatches: int = ...,
    min_overlap: int = ...,
) -> FilterStats: ...
def filter_and_trim(
    cfg: FilterConfig,
    in_path: str,
    out_path: str,
) -> FilterStats: ...
def filter_and_trim_paired(
    cfg_fwd: FilterConfig,
    cfg_rev: FilterConfig,
    fwd_in: str,
    rev_in: str,
    fwd_out: str,
    rev_out: str,
) -> FilterStatsPaired: ...
def learn_errors(paths: Sequence[str], n_reads: int = ...) -> ErrorModel: ...
def derep_fastq(path: str) -> DerepResult: ...
def dada(
    derep: DerepResult,
    err: ErrorModel,
    omega_a: float = ...,
    pool: bool = ...,
) -> DadaResult: ...
def dada_many(
    dereps: Sequence[DerepResult],
    err: ErrorModel,
    omega_a: float = ...,
) -> list[DadaResult]: ...
def dada_pooled(
    dereps: Sequence[DerepResult],
    err: ErrorModel,
    omega_a: float = ...,
) -> list[DadaResult]: ...
def dada_pseudo(
    dereps: Sequence[DerepResult],
    err: ErrorModel,
    omega_a: float = ...,
) -> list[DadaResult]: ...
def merge_pairs(
    fwd: DadaResult,
    rev: DadaResult,
    min_overlap: int = ...,
    max_mismatches: int = ...,
    just_concatenate: bool = ...,
) -> list[MergedRead]: ...
def remove_bimera_denovo(
    seqs_with_abundance: Sequence[tuple[str, int]],
) -> list[tuple[str, int]]: ...
def assign_taxonomy(
    queries: Sequence[str],
    reference_fasta: str,
    lineage_tsv: str,
    minBoot: int = ...,
    n_kmer: int = ...,
) -> list[TaxonAssignment]: ...
def make_sequence_table(
    sample_names: Sequence[str],
    per_sample: Sequence[Sequence[tuple[str, int]]],
) -> SequenceTable: ...
def run_pipeline(
    fwd_path: str,
    rev_path: str | None = ...,
    **kwargs: Any,
) -> dict[str, int]: ...
def run_pipeline_samples(
    fwd_paths: Sequence[str],
    rev_paths: Sequence[str] | None = ...,
    output_dir: str | None = ...,
    **kwargs: Any,
) -> dict[str, dict[str, int]]: ...
