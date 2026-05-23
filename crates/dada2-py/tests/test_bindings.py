"""Smoke tests for the dada2 Python bindings."""
from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path

import pytest

import dada2


FIXTURE_DIR = Path(__file__).parent.parent.parent.parent / "tests" / "integration" / "fixtures"
SAMPLE_FASTQ = FIXTURE_DIR / "sample_R1.fastq"


def test_version_semver():
    ver = dada2.__version__
    assert re.match(r"^\d+\.\d+\.\d+", ver), f"not semver: {ver}"


def test_filter_and_trim():
    with tempfile.NamedTemporaryFile(suffix=".fastq", delete=False) as tmp:
        out_path = tmp.name

    try:
        cfg = dada2.FilterConfig(trunc_len=30, min_len=10, max_ee=10.0)
        stats = dada2.filter_and_trim(cfg, str(SAMPLE_FASTQ), out_path)
        assert Path(out_path).exists()
        assert stats.reads_out > 0
    finally:
        os.unlink(out_path)


def test_learn_errors():
    model = dada2.learn_errors([str(SAMPLE_FASTQ)], n_reads=100)
    assert model is not None


def test_run_dada_remove_bimeras_roundtrip():
    derep = dada2.dereplicate(str(SAMPLE_FASTQ))
    assert len(derep) > 0

    model = dada2.learn_errors([str(SAMPLE_FASTQ)], n_reads=100)
    result = dada2.run_dada(derep, model, omega_a=1e-5)
    assert len(result) > 0

    seqs_with_abund = [(result[i][0], result[i][1]) for i in range(len(result))]
    clean = dada2.remove_bimeras(seqs_with_abund)
    assert len(clean) > 0
