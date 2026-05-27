"""
Benchmark: Python dada2 bindings (Rust core) on simulated 16S V3-V4 paired FASTQs.
Pipeline: filter_and_trim_paired → learn_errors → derep_fastq → dada → merge_pairs
          → remove_bimera_denovo
"""
import json, sys, time
from pathlib import Path

import speeddada as dada2

R1  = "/tmp/bench_fastq/R1.fastq"
R2  = "/tmp/bench_fastq/R2.fastq"
OUT = Path("/tmp/bench_out")
OUT.mkdir(exist_ok=True)

t_total = time.perf_counter()

# 1. Paired filter
cfg_fwd  = dada2.FilterConfig(trunc_len=230, min_len=150, max_ee=3.0)
cfg_rev  = dada2.FilterConfig(trunc_len=210, min_len=150, max_ee=5.0)
filt_r1  = str(OUT / "py_filt_R1.fastq")
filt_r2  = str(OUT / "py_filt_R2.fastq")
t1 = time.perf_counter()
fstats = dada2.filter_and_trim_paired(cfg_fwd, cfg_rev, R1, R2, filt_r1, filt_r2)
t_filter = time.perf_counter() - t1
print(f"[filter_and_trim_paired]  pairs_in=10000  pairs_out={fstats.pairs_out}"
      f"  fwd_fail={fstats.fwd_failed}  rev_fail={fstats.rev_failed}"
      f"  ({t_filter*1000:.1f} ms)")

# 2. Learn errors
t1 = time.perf_counter()
model = dada2.learn_errors([filt_r1, filt_r2], n_reads=20_000)
t_errors = time.perf_counter() - t1
print(f"[learn_errors]            ({t_errors*1000:.1f} ms)")

# 3. Dereplicate
t1 = time.perf_counter()
derep_fwd = dada2.derep_fastq(filt_r1)
derep_rev = dada2.derep_fastq(filt_r2)
t_derep = time.perf_counter() - t1
print(f"[derep_fastq]             fwd_uniq={len(derep_fwd)}  rev_uniq={len(derep_rev)}"
      f"  ({t_derep*1000:.1f} ms)")

# 4. DADA
t1 = time.perf_counter()
res_fwd = dada2.dada(derep_fwd, model, omega_a=1e-40)
res_rev = dada2.dada(derep_rev, model, omega_a=1e-40)
t_dada = time.perf_counter() - t1
print(f"[dada]                    fwd_asvs={len(res_fwd)}  rev_asvs={len(res_rev)}"
      f"  ({t_dada*1000:.1f} ms)")

# 5. Merge
t1 = time.perf_counter()
merged = dada2.merge_pairs(res_fwd, res_rev, min_overlap=12)
t_merge = time.perf_counter() - t1
print(f"[merge_pairs]             merged={len(merged)}  ({t_merge*1000:.1f} ms)")

# 6. Remove bimeras — pass (sequence, abundance) pairs
t1 = time.perf_counter()
clean = dada2.remove_bimera_denovo([(m.sequence, m.abundance) for m in merged])
t_chimera = time.perf_counter() - t1
print(f"[remove_bimera_denovo]    asvs_out={len(clean)}  ({t_chimera*1000:.1f} ms)")

t_total_ms = (time.perf_counter() - t_total) * 1000
print(f"\nTotal Python binding time: {t_total_ms:.1f} ms")

# Build result — clean is list[tuple[bytes, int]]
asvs = []
for seq, abund in clean:
    s = seq.decode() if isinstance(seq, bytes) else seq
    asvs.append({"sequence": s, "abundance": abund})
asvs.sort(key=lambda x: -x["abundance"])

result = {
    "tool":     "Python dada2",
    "total_ms": round(t_total_ms, 1),
    "stages": {
        "filter_ms":       round(t_filter  * 1000, 1),
        "learn_errors_ms": round(t_errors  * 1000, 1),
        "derep_ms":        round(t_derep   * 1000, 1),
        "dada_ms":         round(t_dada    * 1000, 1),
        "merge_ms":        round(t_merge   * 1000, 1),
        "chimera_ms":      round(t_chimera * 1000, 1),
    },
    "asvs": asvs,
}
(OUT / "rust_output.json").write_text(json.dumps(result, indent=2))

print(f"\nTop ASVs (Python dada2):")
for r in asvs[:5]:
    print(f"  {r['sequence'][:40]}...  abundance={r['abundance']}")
