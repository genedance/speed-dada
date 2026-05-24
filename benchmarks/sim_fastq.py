#!/usr/bin/env python3
"""
Generate paired-end FASTQs simulating a 16S V3-V4 Illumina MiSeq run.

10 true ASVs, log-normal abundances, 10 000 read pairs.
Forward reads: 250 bp, Q35 → Q28 quadratic quality decay (MiSeq-like).
Reverse reads: 250 bp, Q28 → Q18 quadratic quality decay (R2 degrades more).
Amplicon:      420 bp (V3-V4 region between forward/reverse primers).
After truncation fwd=230 / rev=210:  230+210-420 = 20 bp overlap for merging.
~3 % of read pairs are bimeric (chimera of top-2 ASVs, crossover at pos 140).
"""
import math, random, sys
from pathlib import Path

SEED         = 42
N_PAIRS      = 10_000
AMPLICON_LEN = 420     # bp between primers
READ_LEN     = 250     # bp raw read length from sequencer
N_ASVS       = 10
CHIMERA_FRAC = 0.03

# ── Per-cycle error profiles ──────────────────────────────────────────────────

def _efwd(pos):
    """Q35 → Q28, quadratic (err 3.2e-4 → 1.6e-3)."""
    frac = pos / (READ_LEN - 1)
    return 3.2e-4 + (1.6e-3 - 3.2e-4) * frac ** 2

def _erev(pos):
    """Q28 → Q18, quadratic (err 1.6e-3 → 1.6e-2)."""
    frac = pos / (READ_LEN - 1)
    return 1.6e-3 + (1.6e-2 - 1.6e-3) * frac ** 2

def _phred(err):
    q = min(40, max(2, round(-10 * math.log10(max(err, 1e-10)))))
    return chr(q + 33)

ERR_FWD  = [_efwd(i) for i in range(READ_LEN)]
ERR_REV  = [_erev(i) for i in range(READ_LEN)]
QUAL_FWD = "".join(_phred(e) for e in ERR_FWD)
QUAL_REV = "".join(_phred(e) for e in ERR_REV)

# ── Helpers ───────────────────────────────────────────────────────────────────

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def _randseq(length, gc, rng):
    at = 1 - gc
    out = []
    for _ in range(length):
        r = rng.random()
        if   r < gc / 2:       out.append("G")
        elif r < gc:            out.append("C")
        elif r < gc + at / 2:  out.append("A")
        else:                   out.append("T")
    return "".join(out)

def _mutate(seq, err_rates, rng):
    bases = list(seq)
    for j, e in enumerate(err_rates):
        if rng.random() < e:
            bases[j] = rng.choice([x for x in "ACGT" if x != bases[j]])
    return "".join(bases)

def _pair(amp, rid, rng):
    fwd = _mutate(amp[:READ_LEN],                     ERR_FWD, rng)
    rev = _mutate(rc(amp[AMPLICON_LEN - READ_LEN :]), ERR_REV, rng)
    r1 = f"@{rid}\n{fwd}\n+\n{QUAL_FWD}"
    r2 = f"@{rid}\n{rev}\n+\n{QUAL_REV}"
    return r1, r2

# ── True sequences & abundances ───────────────────────────────────────────────

rng = random.Random(SEED)

true_seqs = [_randseq(AMPLICON_LEN, 0.50 + rng.uniform(-0.06, 0.06), rng)
             for _ in range(N_ASVS)]

raw    = sorted([rng.lognormvariate(0, 0.8) for _ in range(N_ASVS)], reverse=True)
n_true = round(N_PAIRS * (1 - CHIMERA_FRAC))
scale  = n_true / sum(raw)
counts = [max(1, round(x * scale)) for x in raw]
counts[0] += n_true - sum(counts)          # exact total

n_chimera   = N_PAIRS - sum(counts)
crossover   = AMPLICON_LEN // 3            # crossover at ~1/3 of amplicon
chimera_seq = true_seqs[0][:crossover] + true_seqs[1][crossover:]

# ── Build read pairs ──────────────────────────────────────────────────────────

pairs = []
idx = 0
for ai, (seq, cnt) in enumerate(zip(true_seqs, counts)):
    for _ in range(cnt):
        pairs.append(_pair(seq, f"r{idx}a{ai}", rng))
        idx += 1
for _ in range(n_chimera):
    pairs.append(_pair(chimera_seq, f"r{idx}ch", rng))
    idx += 1

rng.shuffle(pairs)
r1_out, r2_out = zip(*pairs)

# ── Write ─────────────────────────────────────────────────────────────────────

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/bench_fastq")
out.mkdir(parents=True, exist_ok=True)

(out / "R1.fastq").write_text("\n".join(r1_out) + "\n")
(out / "R2.fastq").write_text("\n".join(r2_out) + "\n")

q0f = round(-10 * math.log10(ERR_FWD[0]));  q1f = round(-10 * math.log10(ERR_FWD[-1]))
q0r = round(-10 * math.log10(ERR_REV[0]));  q1r = round(-10 * math.log10(ERR_REV[-1]))

print(f"Wrote {N_PAIRS:,} read pairs → {out}/")
print(f"  {N_ASVS} true ASVs  |  {n_chimera} chimeric reads (~{CHIMERA_FRAC*100:.0f}%)")
print(f"  Abundances: {counts}")
print(f"  Amplicon: {AMPLICON_LEN} bp  |  raw read: {READ_LEN} bp")
print(f"  Fwd quality: Q{q0f}→Q{q1f}  |  Rev quality: Q{q0r}→Q{q1r}")
print(f"  Recommended truncLen: fwd=230, rev=210  →  {230+210-AMPLICON_LEN} bp overlap")
