"""
Benchmark the Rust dada2 pipeline on a pair of simulated FASTQs.
Pipeline: filter_and_trim_paired → learn_errors → derep_fastq → dada → merge_pairs → remove_bimera_denovo
Taxonomy is skipped.
"""
import sys, time, json
from pathlib import Path

sys.path.insert(0, str(Path("/home/fromage/dada2_test/.venv/lib/python3.11/site-packages")))
import dada2

R1 = "/tmp/bench_fastq/R1.fastq"
R2 = "/tmp/bench_fastq/R2.fastq"
OUT = Path("/tmp/bench_out")
OUT.mkdir(exist_ok=True)

t0 = time.perf_counter()

# 1. Paired filter
cfg_fwd = dada2.FilterConfig(trunc_len=110, min_len=50, max_ee=3.0)
cfg_rev = dada2.FilterConfig(trunc_len=110, min_len=50, max_ee=3.0)
filt_r1 = str(OUT / "filt_R1.fastq")
filt_r2 = str(OUT / "filt_R2.fastq")
fstats = dada2.filter_and_trim_paired(cfg_fwd, cfg_rev, R1, R2, filt_r1, filt_r2)
t_filter = time.perf_counter() - t0
print(f"[filter_and_trim_paired]  pairs_in=2000  pairs_out={fstats.pairs_out}  "
      f"fwd_fail={fstats.fwd_failed}  rev_fail={fstats.rev_failed}  "
      f"both_fail={fstats.both_failed}  ({t_filter*1000:.1f} ms)")

# 2. Learn errors (from both filtered reads)
t1 = time.perf_counter()
model = dada2.learn_errors([filt_r1, filt_r2], n_reads=5000)
t_err = time.perf_counter() - t1
print(f"[learn_errors]            ({t_err*1000:.1f} ms)")

# 3. Derep
t1 = time.perf_counter()
derep_fwd = dada2.derep_fastq(filt_r1)
derep_rev = dada2.derep_fastq(filt_r2)
t_derep = time.perf_counter() - t1
print(f"[derep_fastq]             fwd_uniq={len(derep_fwd)}  rev_uniq={len(derep_rev)}  ({t_derep*1000:.1f} ms)")

# 4. DADA
t1 = time.perf_counter()
res_fwd = dada2.dada(derep_fwd, model, omega_a=1e-40)
res_rev = dada2.dada(derep_rev, model, omega_a=1e-40)
t_dada = time.perf_counter() - t1
print(f"[dada]                    fwd_asvs={len(res_fwd)}  rev_asvs={len(res_rev)}  ({t_dada*1000:.1f} ms)")

# 5. Merge paired
t1 = time.perf_counter()
merged = dada2.merge_pairs(res_fwd, res_rev)
t_merge = time.perf_counter() - t1
print(f"[merge_pairs]             merged={len(merged)}  ({t_merge*1000:.1f} ms)")

# 6. Remove bimeras
t1 = time.perf_counter()
clean = dada2.remove_bimera_denovo(merged)
t_chim = time.perf_counter() - t1
total = time.perf_counter() - t0
print(f"[remove_bimera_denovo]    asvs_out={len(clean)}  ({t_chim*1000:.1f} ms)")
print(f"\nTotal Rust pipeline time: {total*1000:.1f} ms")

# Build merged dict for abundance lookup
merged_dict = {(m[0].decode() if isinstance(m[0], bytes) else m[0]): m[1] for m in merged}

results = []
for seq in clean:
    seq_str = seq.decode() if isinstance(seq, bytes) else seq
    results.append({"sequence": seq_str, "abundance": merged_dict.get(seq_str, 0)})
results.sort(key=lambda x: -x["abundance"])

(OUT / "rust_output.json").write_text(json.dumps(results, indent=2))
print(f"\nTop ASVs (Rust):")
for r in results[:5]:
    print(f"  {r['sequence'][:40]}...  abundance={r['abundance']}")
