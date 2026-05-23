# Skill: bio-algorithms

## Invoke with
"use the bio-algorithms skill"

## DADA core algorithm (Callahan 2016, Suppl. Note 1)
- Abundance p-value: p = P(X >= n_reads | Poisson(lambda))
  where lambda = total_reads * error_rate(seq_i, seq_j)
- Partition: greedily assign each unique seq to the centre
  that maximises P(reads | centre) * abundance
- EM convergence: |DeltalogL| < 1e-6 or maxIter = 16 (configurable)
- OMEGA_A default: 1e-40 (make this a function parameter)

## Error model (16-class logistic regression)
- 16 transitions: {A,C,G,T} x {A,C,G,T}
- Feature: [Phred score, is_transition (bool)]
- Fit with gradient descent; convergence tol = 1e-4

## Bimera detection rules
- Left arm: match from position 0 until first mismatch
- Right arm: match from last mismatch to end
- Both arms must be >= 8 bp (DADA2 default)
- Parent must have higher abundance than candidate

## Naive Bayes taxonomy (Wang et al. 2007)
- k = 8, stride = 1 (all k-mers)
- Bootstrap: subsample 1/8 of k-mers x 100 reps
- Report genus if bootstrap >= 80% (configurable threshold)
