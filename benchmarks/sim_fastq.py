"""Generate paired-end FASTQs with realistic per-cycle quality decay for dada2 benchmarking."""
import random, math, sys
from pathlib import Path

TRUE_SEQ = ("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")  # 120 bp

def rc(s):
    t = str.maketrans("ACGT", "TGCA")
    return s.translate(t)[::-1]

def cycle_error_rate(pos, read_len, base_err=0.001, tail_err=0.05):
    """Illumina-like: low error at 5' end, rising toward 3' end."""
    frac = pos / (read_len - 1)
    return base_err + (tail_err - base_err) * frac ** 2

def phred_char(err):
    q = min(40, max(2, int(-10 * math.log10(max(err, 1e-10)))))
    return chr(q + 33)

def make_reads(seq, n, seed):
    rng = random.Random(seed)
    L = len(seq)
    err_rates = [cycle_error_rate(i, L) for i in range(L)]
    qual_str = "".join(phred_char(e) for e in err_rates)
    lines = []
    for i in range(n):
        bases = list(seq)
        for j, b in enumerate(bases):
            if rng.random() < err_rates[j]:
                bases[j] = rng.choice([x for x in "ACGT" if x != b])
        s = "".join(bases)
        lines += [f"@read_{i}", s, "+", qual_str]
    return "\n".join(lines) + "\n"

n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/tmp/bench_fastq")
out.mkdir(exist_ok=True)

(out / "R1.fastq").write_text(make_reads(TRUE_SEQ, n, 42))
(out / "R2.fastq").write_text(make_reads(rc(TRUE_SEQ), n, 43))
print(f"Wrote {n} read pairs to {out}/  (len={len(TRUE_SEQ)}bp, err Q40->Q13 across read)")
