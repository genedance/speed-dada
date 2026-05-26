#!/usr/bin/env python3
"""
Benchmark: Python dada2 (Rust core via PyO3) on 6 paired manure samples.

Multi-sample workflow:
  filter_and_trim_paired (per-sample)
   → learn_errors (fwd pool + rev pool)
   → derep_fastq (per-sample)
   → dada (per-sample, fwd + rev)
   → merge_pairs (per-sample)
   → make_sequence_table → remove_bimera_denovo
"""
import argparse
import json
import time
from pathlib import Path

import dada2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--in-dir", default="/Users/alex/Downloads/raw_data_FIELD")
    ap.add_argument("--out-dir", default="/tmp/bench_field_out/python")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    nt, ram = dada2.configure_runtime(n_threads=args.threads)
    print(f"[Python dada2] n_threads={nt}  ram_mb={ram}  samples=6")

    samples = [f"T0-Manure-rep{i}" for i in range(1, 7)]
    in_dir  = Path(args.in_dir)
    fwd_in  = [str(in_dir / f"{s}_R1.fastq.gz") for s in samples]
    rev_in  = [str(in_dir / f"{s}_R2.fastq.gz") for s in samples]
    fwd_filt = [str(out / f"{s}_R1_filt.fastq.gz") for s in samples]
    rev_filt = [str(out / f"{s}_R2_filt.fastq.gz") for s in samples]

    t_total = time.perf_counter()

    # 1. Filter — Python bindings expose one-sample filter_and_trim_paired
    print("[filter_and_trim_paired]")
    cfg_fwd = dada2.FilterConfig(trunc_len=240, min_len=50, max_ee=2.0, trunc_q=2)
    cfg_rev = dada2.FilterConfig(trunc_len=180, min_len=50, max_ee=4.0, trunc_q=2)
    t1 = time.perf_counter()
    per_sample = []
    for s, fi, ri, fo, ro in zip(samples, fwd_in, rev_in, fwd_filt, rev_filt):
        stats = dada2.filter_and_trim_paired(cfg_fwd, cfg_rev, fi, ri, fo, ro)
        per_sample.append({"sample": s,
                           "reads_in":  stats.reads_in,
                           "reads_out": stats.pairs_out})
        print(f"  {s}: in={stats.reads_in}  out={stats.pairs_out}")
    t_filter = time.perf_counter() - t1
    print(f"  total filter: {t_filter*1000:.1f} ms")

    # 2. Learn errors — pool fwd files, then rev files
    print("[learn_errors]")
    t1 = time.perf_counter()
    errF = dada2.learn_errors(fwd_filt)
    errR = dada2.learn_errors(rev_filt)
    t_errors = time.perf_counter() - t1
    print(f"  done  ({t_errors*1000:.1f} ms)")

    # 3. Dereplicate per sample
    print("[derep_fastq]")
    t1 = time.perf_counter()
    derep_fwd = [dada2.derep_fastq(p) for p in fwd_filt]
    derep_rev = [dada2.derep_fastq(p) for p in rev_filt]
    t_derep = time.perf_counter() - t1
    print(f"  done  ({t_derep*1000:.1f} ms)")

    # 4. DADA — pseudo-pool (standard dada2 cross-sample practice)
    # Per-sample pass1 → collect ASV priors → per-sample pass2 with priors
    # auto-promoting. Both passes parallelise across samples via Rayon.
    print("[dada_pseudo]")
    t1 = time.perf_counter()
    dada_fwd = dada2.dada_pseudo(derep_fwd, errF, omega_a=1e-40)
    dada_rev = dada2.dada_pseudo(derep_rev, errR, omega_a=1e-40)
    t_dada = time.perf_counter() - t1
    n_asvF = sum(len(r) for r in dada_fwd)
    n_asvR = sum(len(r) for r in dada_rev)
    print(f"  fwd_asvs(total)={n_asvF}  rev_asvs(total)={n_asvR}  "
          f"({t_dada*1000:.1f} ms)")

    # 5. Merge per sample
    print("[merge_pairs]")
    t1 = time.perf_counter()
    merged = [dada2.merge_pairs(f, r, min_overlap=12) for f, r in zip(dada_fwd, dada_rev)]
    t_merge = time.perf_counter() - t1
    n_merged = sum(len(m) for m in merged)
    print(f"  total_merged_asvs={n_merged}  ({t_merge*1000:.1f} ms)")

    # 6. Sequence table + chimera removal
    # Convert merged → per-sample (seq, count) lists for make_sequence_table,
    # which in this binding takes a sample_names list and a list of DadaResults.
    # Since we have MergedRead objects, we'll build the count matrix ourselves
    # then run remove_bimera_denovo on aggregate counts.
    print("[chimera + sequence table]")
    t1 = time.perf_counter()
    from collections import defaultdict
    counts = defaultdict(lambda: [0] * len(samples))
    for i, sample_merged in enumerate(merged):
        for m in sample_merged:
            if m.accept:
                seq = m.sequence.decode() if isinstance(m.sequence, bytes) else m.sequence
                counts[seq][i] += m.abundance
    # Aggregate counts for chimera detection (pass (bytes, int) pairs)
    agg = [(seq.encode(), sum(per_sample_counts))
           for seq, per_sample_counts in counts.items()]
    clean = dada2.remove_bimera_denovo(agg)
    clean_seqs = {(s.decode() if isinstance(s, bytes) else s): a
                  for s, a in clean}
    t_chimera = time.perf_counter() - t1
    n_asvs_before = len(counts)
    n_asvs_after  = len(clean_seqs)
    print(f"  asvs_in={n_asvs_before}  asvs_out={n_asvs_after}  "
          f"({t_chimera*1000:.1f} ms)")

    t_total_ms = (time.perf_counter() - t_total) * 1000
    print(f"\nTotal Python dada2 time: {t_total_ms:.1f} ms")

    asvs = [{"sequence": seq, "abundance": clean_seqs[seq]}
            for seq in clean_seqs]
    asvs.sort(key=lambda x: -x["abundance"])

    result = {
        "tool":      "Python dada2",
        "n_threads": nt,
        "total_ms":  round(t_total_ms, 1),
        "stages": {
            "filter_ms":       round(t_filter  * 1000, 1),
            "learn_errors_ms": round(t_errors  * 1000, 1),
            "derep_ms":        round(t_derep   * 1000, 1),
            "dada_ms":         round(t_dada    * 1000, 1),
            "merge_ms":        round(t_merge   * 1000, 1),
            "chimera_ms":      round(t_chimera * 1000, 1),
        },
        "samples": per_sample,
        "n_asvs_before_chimera": n_asvs_before,
        "n_asvs_after_chimera":  n_asvs_after,
        "total_abundance":       sum(clean_seqs.values()),
        "asvs": asvs,
    }
    (out / "rust_output.json").write_text(json.dumps(result, indent=2))
    print(f"\nWrote {out / 'rust_output.json'}")


if __name__ == "__main__":
    main()
