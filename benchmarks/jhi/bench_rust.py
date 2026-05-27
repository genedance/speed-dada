#!/usr/bin/env python3
"""Benchmark: Python dada2 (Rust core via PyO3) on 3 paired JHI samples."""
import argparse
import json
import time
from pathlib import Path
from collections import defaultdict

import speeddada as dada2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", default="/tmp/bench_jhi_out/python")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    nt, _ = dada2.configure_runtime(n_threads=args.threads)
    stems = ["JHI-2025-Q1-A-004", "JHI-2025-Q1-A-009", "JHI-2025-Q1-A-010"]
    in_dir = Path(args.in_dir)
    fwd_in  = [str(in_dir / f"raw.{s}.1.fq.gz") for s in stems]
    rev_in  = [str(in_dir / f"raw.{s}.2.fq.gz") for s in stems]
    fwd_filt = [str(out / f"{s}_R1_filt.fastq.gz") for s in stems]
    rev_filt = [str(out / f"{s}_R2_filt.fastq.gz") for s in stems]

    print(f"[Python dada2] n_threads={nt}  samples={len(stems)}")
    t_total = time.perf_counter()

    print("[filter_and_trim_paired]")
    cfg_fwd = dada2.FilterConfig(trunc_len=240, min_len=50, max_ee=2.0, trunc_q=2)
    cfg_rev = dada2.FilterConfig(trunc_len=180, min_len=50, max_ee=4.0, trunc_q=2)
    t1 = time.perf_counter()
    per_sample = []
    for s, fi, ri, fo, ro in zip(stems, fwd_in, rev_in, fwd_filt, rev_filt):
        stats = dada2.filter_and_trim_paired(cfg_fwd, cfg_rev, fi, ri, fo, ro)
        per_sample.append({"sample": s, "reads_in": stats.reads_in,
                           "reads_out": stats.pairs_out})
        print(f"  {s}: in={stats.reads_in}  out={stats.pairs_out}")
    t_filter = time.perf_counter() - t1
    print(f"  total filter: {t_filter*1000:.1f} ms")

    print("[learn_errors]")
    t1 = time.perf_counter()
    errF = dada2.learn_errors(fwd_filt)
    errR = dada2.learn_errors(rev_filt)
    t_errors = time.perf_counter() - t1
    print(f"  done  ({t_errors*1000:.1f} ms)")

    print("[derep_fastq]")
    t1 = time.perf_counter()
    derep_fwd = [dada2.derep_fastq(p) for p in fwd_filt]
    derep_rev = [dada2.derep_fastq(p) for p in rev_filt]
    t_derep = time.perf_counter() - t1
    print(f"  done  ({t_derep*1000:.1f} ms)")

    print("[dada_pseudo]")
    t1 = time.perf_counter()
    dada_fwd = dada2.dada_pseudo(derep_fwd, errF)
    dada_rev = dada2.dada_pseudo(derep_rev, errR)
    t_dada = time.perf_counter() - t1
    n_asvF = sum(len(r) for r in dada_fwd)
    n_asvR = sum(len(r) for r in dada_rev)
    print(f"  fwd_asvs(total)={n_asvF}  rev_asvs(total)={n_asvR}  "
          f"({t_dada*1000:.1f} ms)")

    print("[merge_pairs]")
    t1 = time.perf_counter()
    merged = [dada2.merge_pairs(f, r, min_overlap=12)
              for f, r in zip(dada_fwd, dada_rev)]
    t_merge = time.perf_counter() - t1
    n_merged = sum(len(m) for m in merged)
    print(f"  total_merged_asvs={n_merged}  ({t_merge*1000:.1f} ms)")

    print("[chimera + sequence table]")
    t1 = time.perf_counter()
    counts = defaultdict(int)
    for sample_merged in merged:
        for m in sample_merged:
            if m.accept:
                seq = m.sequence.decode() if isinstance(m.sequence, bytes) else m.sequence
                counts[seq] += m.abundance
    agg = [(s.encode(), a) for s, a in counts.items()]
    clean = dada2.remove_bimera_denovo(agg)
    clean_d = {(s.decode() if isinstance(s, bytes) else s): a for s, a in clean}
    t_chimera = time.perf_counter() - t1
    print(f"  asvs_in={len(counts)}  asvs_out={len(clean_d)}  "
          f"({t_chimera*1000:.1f} ms)")

    t_total_ms = (time.perf_counter() - t_total) * 1000
    print(f"\nTotal Python dada2 time: {t_total_ms:.1f} ms")

    asvs = sorted(
        ({"sequence": s, "abundance": clean_d[s]} for s in clean_d),
        key=lambda x: -x["abundance"],
    )

    result = {
        "tool": "Python dada2",
        "n_threads": nt,
        "total_ms": round(t_total_ms, 1),
        "stages": {
            "filter_ms":       round(t_filter  * 1000, 1),
            "learn_errors_ms": round(t_errors  * 1000, 1),
            "derep_ms":        round(t_derep   * 1000, 1),
            "dada_ms":         round(t_dada    * 1000, 1),
            "merge_ms":        round(t_merge   * 1000, 1),
            "chimera_ms":      round(t_chimera * 1000, 1),
        },
        "samples": per_sample,
        "n_asvs_before_chimera": len(counts),
        "n_asvs_after_chimera":  len(clean_d),
        "total_abundance":       sum(clean_d.values()),
        "asvs": asvs,
    }
    (out / "rust_output.json").write_text(json.dumps(result, indent=2))
    print(f"\nWrote {out / 'rust_output.json'}")


if __name__ == "__main__":
    main()
