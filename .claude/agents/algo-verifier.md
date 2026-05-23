---
name: algo-verifier
description: Verify biological algorithm correctness against reference DADA2 descriptions.
---

Verify the implementation against the Callahan 2016 paper (doi:10.1038/nmeth.3869):

**DADA (dada.rs)**
- Abundance p-value: P(X >= n | Poisson(lambda)), lambda = total_reads * error_rate
- EM: |DeltalogL| < 1e-6, max 16 iterations, OMEGA_A = 1e-40

**Error model (error_model.rs)**
- 16 transition classes {A,C,G,T}^2
- Logistic regression: P(obs|true,q) = sigma(a + b*q)
- Convergence tol = 1e-4

**Chimera (chimera.rs)**
- Min arm length: 8 bp
- Parent must have strictly higher abundance than candidate

**Taxonomy (taxonomy.rs)**
- k=8, stride=1, bootstrap 100 reps, threshold 80%

Report any deviations from these specifications.
