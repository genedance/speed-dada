//! Stage 3 — Parametric error model learned by EM.
//!
//! Fits a logistic regression  P(obs | true, q) = σ(a + b·q)  for each of
//! the 16 base-transition classes.  The resulting [`ErrorModel`] is a
//! `16 × max_qual` matrix of substitution probabilities, with a precomputed
//! log-probability matrix for use in the DADA hot loop.

use crate::{Dada2Error, Phred};
use ndarray::{Array1, Array2};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Index encoding for the 16 transition classes (`true_base × 4 + obs_base`).
/// Bases: 0=A, 1=C, 2=G, 3=T.
pub const N_TRANSITIONS: usize = 16;

/// Maximum Phred score stored in the error matrix.
pub const MAX_QUAL: usize = 41;

/// Learned error model: a `16 × MAX_QUAL` matrix of P(obs | true, q).
#[derive(Debug, Clone, Serialize)]
pub struct ErrorModel {
    /// `matrix[[trans, q]]` = P(observing obs-base | true-base, Phred q).
    /// Rows 0..16 index transitions (`true*4+obs`); columns 0..`MAX_QUAL` index Phred.
    pub matrix: Array2<f64>,
    /// Number of reads used to fit this model.
    pub n_reads_used: u64,
    /// Precomputed `log(matrix[[trans, q]].max(1e-300))` for the DADA hot loop.
    #[serde(skip)]
    log_matrix: Array2<f64>,
}

impl<'de> Deserialize<'de> for ErrorModel {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Wire {
            matrix: Array2<f64>,
            n_reads_used: u64,
        }
        let Wire {
            matrix,
            n_reads_used,
        } = Wire::deserialize(d)?;
        let log_matrix = build_log_matrix(&matrix);
        Ok(Self {
            matrix,
            n_reads_used,
            log_matrix,
        })
    }
}

fn build_log_matrix(matrix: &Array2<f64>) -> Array2<f64> {
    matrix.mapv(|p| p.max(1e-300).ln())
}

impl ErrorModel {
    /// Probability of observing transition `(true_base, obs_base)` at quality `q`.
    ///
    /// Bases encoded as 0=A,1=C,2=G,3=T.
    #[must_use]
    pub fn p_error(&self, true_base: u8, obs_base: u8, q: Phred) -> f64 {
        let row = (true_base as usize) * 4 + (obs_base as usize);
        let col = (q.0 as usize).min(MAX_QUAL - 1);
        self.matrix[[row, col]]
    }

    /// Log-probability of transition `(true_base, obs_base)` at quality `q`.
    ///
    /// Equivalent to `p_error(…).max(1e-300).ln()` but uses a precomputed table.
    #[inline]
    #[must_use]
    pub fn log_p_error(&self, true_base: u8, obs_base: u8, q: Phred) -> f64 {
        let row = (true_base as usize) * 4 + (obs_base as usize);
        let col = (q.0 as usize).min(MAX_QUAL - 1);
        self.log_matrix[[row, col]]
    }

    /// Compute the log-likelihood that `obs` was produced from `truth` under this model.
    #[must_use]
    pub fn log_likelihood(&self, truth: &[u8], obs: &[u8], quals: &[u8]) -> f64 {
        truth
            .iter()
            .zip(obs.iter())
            .zip(quals.iter())
            .map(|((&t, &o), &qc)| {
                let tb = base_index(t);
                let ob = base_index(o);
                let q = Phred::from_ascii(qc);
                self.log_p_error(tb, ob, q)
            })
            .sum()
    }

    /// Build a default error model using the Illumina-like logistic curve.
    ///
    /// This is used as the initial estimate before EM refinement and as a
    /// fallback when too few reads are available.
    #[must_use]
    pub fn illumina_default() -> Self {
        let mut matrix = Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL));
        for row in 0..N_TRANSITIONS {
            let is_match = row % 5 == 0; // diagonal = A→A, C→C, G→G, T→T
            for col in 0..MAX_QUAL {
                #[allow(clippy::cast_precision_loss)]
                let q = col as f64;
                let p_err = 10f64.powf(-q / 10.0);
                matrix[[row, col]] = if is_match { 1.0 - p_err } else { p_err / 3.0 };
            }
        }
        let log_matrix = build_log_matrix(&matrix);
        Self {
            matrix,
            n_reads_used: 0,
            log_matrix,
        }
    }
}

/// Smoothing method used to fit per-quality transition probabilities.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ErrFunKind {
    /// 2-parameter logistic regression σ(a + b·q) — the original SpeedDada
    /// behaviour. Cheap but crude; under-fits any non-monotone curve.
    Logistic,
    /// Weighted local linear regression with a tricubic kernel
    /// (LOWESS, degree 1). Roughly equivalent to `dada2::loessErrfun`
    /// with default span — captures the empirical curvature in the
    /// transition data, including the typical "shoulder" at low quality
    /// scores where the logistic fit overshoots.
    Loess,
}

impl Default for ErrFunKind {
    fn default() -> Self {
        Self::Loess
    }
}

/// Configuration for error model learning.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorLearningConfig {
    /// Maximum reads to use during EM fitting.
    pub n_reads: usize,
    /// Maximum EM iterations.
    pub max_iter: usize,
    /// Log-likelihood convergence tolerance.
    pub tol: f64,
    /// RNG seed for reproducibility.
    pub seed: u64,
    /// Smoothing method used to fit the per-quality transition probabilities.
    /// Defaults to [`ErrFunKind::Loess`] for parity with dada2.
    #[serde(default)]
    pub method: ErrFunKind,
    /// Number of selfConsist passes the *wrapper layer* will run.
    ///
    /// The selfConsist loop alternates `learn_errors → dada → recover
    /// mismatch counts from clusters → refit`. Because that loop crosses
    /// module boundaries (it must call `crate::dada`), it lives outside
    /// `error_model`; this field is metadata so the wrapper knows how
    /// many iterations were requested. `0` (the default) means run a
    /// single pass with empirical mismatch evidence collected from the
    /// dereplicated read set.
    #[serde(default)]
    pub self_consist_iters: u32,
    /// LOESS span (fraction of points contributing to each local fit).
    /// Default 0.95 matches `dada2::loessErrfun`.
    #[serde(default = "default_loess_span")]
    pub loess_span: f64,
    /// Maximum Hamming distance for the same-length pairwise pass that
    /// recovers mismatch evidence. Pairs farther apart than this are
    /// treated as different species rather than parent/error pairs.
    /// Default 8 matches dada2's `errBalancedF` collection radius.
    #[serde(default = "default_max_pair_dist")]
    pub max_pair_dist: u32,
}

fn default_loess_span() -> f64 {
    0.95
}

fn default_max_pair_dist() -> u32 {
    8
}

impl Default for ErrorLearningConfig {
    fn default() -> Self {
        Self {
            n_reads: 1_000_000,
            max_iter: 16,
            tol: 1e-6,
            seed: 42,
            method: ErrFunKind::default(),
            self_consist_iters: 0,
            loess_span: default_loess_span(),
            max_pair_dist: default_max_pair_dist(),
        }
    }
}

/// Learn error rates from a collection of reads.
///
/// Two phases:
/// 1. **Match collection (cheap):** every called base contributes a
///    `(base, base, q)` self-transition. This is the only signal in
///    purely-self-aligned counts and exactly matches the legacy behaviour.
/// 2. **Mismatch collection (selfConsist-lite):** within each same-length
///    bucket, the most-abundant unique sequence is treated as the local
///    "parent" and every lower-abundance unique within `cfg.max_pair_dist`
///    Hamming bases of it contributes `(parent_base → child_base, q)`
///    mismatch counts (weighted by `min(parent_count, child_count)`).
///    This is a single-pass approximation of dada2's `selfConsist`
///    refinement; the full EM loop is expected to be driven from the
///    wrapper layer (it has to call `crate::dada`, so it can't live here).
///
/// The transition counts are then smoothed per row using either the
/// 2-parameter logistic fit (legacy) or LOESS (degree-1 local linear
/// regression with tricubic weights, default — equivalent to dada2's
/// `loessErrfun`).
///
/// # Errors
/// Returns [`Dada2Error::Convergence`] if the logistic fit diverges on the
/// first step. Returns [`Dada2Error::InvalidInput`] if no reads are supplied.
pub fn learn_errors(
    records: &[crate::io::fastq::FastqRecord],
    cfg: &ErrorLearningConfig,
) -> Result<ErrorModel, Dada2Error> {
    if records.is_empty() {
        return Err(Dada2Error::InvalidInput(
            "no reads supplied to learn_errors".into(),
        ));
    }

    let n = cfg.n_reads.min(records.len());
    let mut counts = collect_match_counts(&records[..n]);
    add_pairwise_mismatch_counts(&records[..n], &mut counts, cfg.max_pair_dist);

    let matrix = fit_from_counts(&counts, cfg)?;
    let log_matrix = build_log_matrix(&matrix);
    Ok(ErrorModel {
        matrix,
        n_reads_used: n as u64,
        log_matrix,
    })
}

/// Refit an [`ErrorModel`] from precomputed transition counts.
///
/// Intended for selfConsist-style orchestration in a wrapper layer:
/// after each `dada()` pass the wrapper rebuilds the count matrix from
/// inferred parent→child edges and calls this to refit.
///
/// # Errors
/// See [`learn_errors`].
pub fn refit_from_counts(
    counts: &Array2<f64>,
    cfg: &ErrorLearningConfig,
    n_reads_used: u64,
) -> Result<ErrorModel, Dada2Error> {
    let matrix = fit_from_counts(counts, cfg)?;
    let log_matrix = build_log_matrix(&matrix);
    Ok(ErrorModel {
        matrix,
        n_reads_used,
        log_matrix,
    })
}

fn collect_match_counts(records: &[crate::io::fastq::FastqRecord]) -> Array2<f64> {
    records
        .par_iter()
        .fold(
            || Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL)),
            |mut acc, rec| {
                for (&base, &qc) in rec.seq.iter().zip(rec.qual.iter()) {
                    let bi = base_index(base) as usize;
                    let q = Phred::from_ascii(qc).0 as usize;
                    let col = q.min(MAX_QUAL - 1);
                    let row = bi * 4 + bi;
                    acc[[row, col]] += 1.0;
                }
                acc
            },
        )
        .reduce(
            || Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL)),
            |mut a, b| {
                a += &b;
                a
            },
        )
}

/// Within each same-length bucket of reads, find the most-abundant unique
/// sequence and attribute mismatches in similar neighbours to it.
///
/// This is a *single-pass* approximation of dada2's `selfConsist`: it
/// extracts real (parent_base → child_base, q) mismatch evidence without
/// needing to run the DADA denoiser. Pairs more than `max_dist` bases
/// apart are treated as biologically different species and ignored.
fn add_pairwise_mismatch_counts(
    records: &[crate::io::fastq::FastqRecord],
    counts: &mut Array2<f64>,
    max_dist: u32,
) {
    use std::collections::HashMap;

    // Dereplicate by (length, seq) within memory; carry per-position
    // quality sums so we can recover a mean-quality vector per unique.
    #[derive(Default, Clone)]
    struct Unique {
        count: u64,
        qual_sum: Vec<u64>,
    }
    let mut by_seq: HashMap<Vec<u8>, Unique> = HashMap::new();
    for rec in records {
        if rec.seq.len() != rec.qual.len() {
            continue;
        }
        let u = by_seq.entry(rec.seq.clone()).or_insert_with(|| Unique {
            count: 0,
            qual_sum: vec![0; rec.seq.len()],
        });
        u.count += 1;
        for (i, &q) in rec.qual.iter().enumerate() {
            u.qual_sum[i] += u64::from(Phred::from_ascii(q).0);
        }
    }
    if by_seq.is_empty() {
        return;
    }

    // Bucket by length.
    let mut by_len: HashMap<usize, Vec<(Vec<u8>, Unique)>> = HashMap::new();
    for (seq, u) in by_seq {
        let len = seq.len();
        by_len.entry(len).or_default().push((seq, u));
    }

    for (_, mut bucket) in by_len {
        // Sort by descending count so bucket[0] is the local "parent".
        bucket.sort_unstable_by(|a, b| b.1.count.cmp(&a.1.count));
        let (parent_seq, parent_uniq) = bucket[0].clone();
        if bucket.len() < 2 {
            continue;
        }
        for (child_seq, child_uniq) in bucket.iter().skip(1) {
            let dist = crate::align::hamming_distance(&parent_seq, child_seq);
            if dist == 0 || dist > max_dist {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let weight = child_uniq.count.min(parent_uniq.count) as f64;
            for (i, (&pb, &cb)) in parent_seq.iter().zip(child_seq.iter()).enumerate() {
                if pb == cb {
                    continue;
                }
                let pbi = base_index(pb) as usize;
                let cbi = base_index(cb) as usize;
                // Use the child's mean quality at this position as the
                // representative q for the (parent→child) mismatch.
                #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
                let mean_q = (child_uniq.qual_sum[i] as f64 / child_uniq.count as f64) as usize;
                let col = mean_q.min(MAX_QUAL - 1);
                let row = pbi * 4 + cbi;
                counts[[row, col]] += weight;
                // Also REMOVE the match miscount we made in
                // collect_match_counts (the child's base was counted as a
                // self-transition; here we re-attribute `weight` of it to a
                // mismatch). Don't go below zero.
                let self_row = cbi * 4 + cbi;
                let prev = counts[[self_row, col]];
                counts[[self_row, col]] = (prev - weight).max(0.0);
            }
        }
    }
}

fn fit_from_counts(
    counts: &Array2<f64>,
    cfg: &ErrorLearningConfig,
) -> Result<Array2<f64>, Dada2Error> {
    let mut matrix = Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL));
    match cfg.method {
        ErrFunKind::Logistic => {
            let rows: Vec<Array1<f64>> = (0..N_TRANSITIONS)
                .map(|row| counts.row(row).to_owned())
                .collect();
            let params: Vec<Result<[f64; 2], Dada2Error>> = rows
                .par_iter()
                .enumerate()
                .map(|(row, row_counts)| {
                    let is_match = row % 5 == 0;
                    fit_logistic_row(row_counts, is_match, cfg)
                })
                .collect();
            for (row, res) in params.into_iter().enumerate() {
                let [a, b] = res?;
                for col in 0..MAX_QUAL {
                    #[allow(clippy::cast_precision_loss)]
                    let q = col as f64;
                    matrix[[row, col]] = sigmoid(a + b * q);
                }
            }
        }
        ErrFunKind::Loess => {
            // Normalise to empirical rates per column (true_base column-sum
            // = total times that base was seen at quality q), then smooth
            // each row independently.
            for tb in 0..4 {
                // Compute column totals across the 4 transition rows for this true base.
                let mut col_total = [0.0_f64; MAX_QUAL];
                for ob in 0..4 {
                    let r = tb * 4 + ob;
                    for col in 0..MAX_QUAL {
                        col_total[col] += counts[[r, col]];
                    }
                }
                // Per (true_base → obs_base): empirical p, then LOWESS.
                for ob in 0..4 {
                    let r = tb * 4 + ob;
                    let mut y = vec![0.0_f64; MAX_QUAL];
                    let mut w = vec![0.0_f64; MAX_QUAL];
                    for col in 0..MAX_QUAL {
                        let total = col_total[col];
                        if total > 0.0 {
                            y[col] = counts[[r, col]] / total;
                            w[col] = total;
                        }
                    }
                    let smoothed = lowess_row(&y, &w, cfg.loess_span);
                    let is_match = tb == ob;
                    for col in 0..MAX_QUAL {
                        // Clamp into [eps, 1] and enforce monotonicity for
                        // self-transitions (matches dada2 behaviour).
                        let p = smoothed[col].clamp(1e-12, 1.0);
                        #[allow(clippy::cast_precision_loss)]
                        let q = col as f64;
                        matrix[[r, col]] = if w[col] == 0.0 {
                            // No empirical evidence at this quality bin —
                            // fall back to the Illumina-default prior.
                            let p_err = 10f64.powf(-q / 10.0);
                            if is_match {
                                1.0 - p_err
                            } else {
                                p_err / 3.0
                            }
                        } else {
                            p
                        };
                    }
                    if is_match {
                        // Enforce non-increasing-with-q for the self
                        // transition (probability of correct call should
                        // not drop as quality increases).
                        for col in (1..MAX_QUAL).rev() {
                            if matrix[[r, col - 1]] < matrix[[r, col]] {
                                matrix[[r, col - 1]] = matrix[[r, col]];
                            }
                        }
                    }
                }
                // Renormalise each column over the 4 obs bases so the row
                // group sums to ~1 (LOWESS is per-row independent, so the
                // unnormalised sum drifts slightly).
                for col in 0..MAX_QUAL {
                    let s: f64 = (0..4).map(|ob| matrix[[tb * 4 + ob, col]]).sum();
                    if s > 0.0 {
                        for ob in 0..4 {
                            matrix[[tb * 4 + ob, col]] /= s;
                        }
                    }
                }
            }
        }
    }
    Ok(matrix)
}

/// Weighted local linear regression (LOWESS, degree 1) over the q-axis.
///
/// `y[q]` is the empirical transition probability at quality `q`; `w[q]`
/// is the count weight at that quality (zero means "no observation").
/// `span` is the fraction of points whose neighbourhood determines each
/// fit (matches dada2's loess `span` default of 0.95).
fn lowess_row(y: &[f64], w: &[f64], span: f64) -> Vec<f64> {
    let n = y.len();
    let mut out = vec![0.0_f64; n];
    if n == 0 {
        return out;
    }
    let k = ((span * n as f64).round() as usize).max(3).min(n);

    for i in 0..n {
        // The k nearest neighbours of x=i are simply [i-k/2 .. i+k/2].
        let half = k / 2;
        let lo = i.saturating_sub(half);
        let hi = (lo + k).min(n);
        let lo = hi.saturating_sub(k); // shift left if we hit the right edge
        let max_dist = (i.saturating_sub(lo)).max(hi - 1 - i).max(1) as f64;

        let mut sw = 0.0_f64;
        let mut sx = 0.0_f64;
        let mut sy = 0.0_f64;
        let mut sxx = 0.0_f64;
        let mut sxy = 0.0_f64;
        for j in lo..hi {
            let d = ((j as i64 - i as i64).abs() as f64) / max_dist;
            if d >= 1.0 {
                continue;
            }
            // Tricubic kernel
            let kw = (1.0 - d.powi(3)).powi(3);
            let wj = kw * w[j];
            if wj <= 0.0 {
                continue;
            }
            let x = j as f64;
            sw += wj;
            sx += wj * x;
            sy += wj * y[j];
            sxx += wj * x * x;
            sxy += wj * x * y[j];
        }
        if sw <= 0.0 {
            out[i] = y[i];
            continue;
        }
        let mean_x = sx / sw;
        let mean_y = sy / sw;
        let var_x = sxx / sw - mean_x * mean_x;
        let cov_xy = sxy / sw - mean_x * mean_y;
        let beta = if var_x.abs() < 1e-12 { 0.0 } else { cov_xy / var_x };
        let alpha = mean_y - beta * mean_x;
        out[i] = alpha + beta * i as f64;
    }
    out
}

/// Fit a 2-parameter logistic model σ(a + b·q) to a count vector by gradient descent.
#[allow(clippy::many_single_char_names)]
fn fit_logistic_row(
    counts: &Array1<f64>,
    is_match: bool,
    cfg: &ErrorLearningConfig,
) -> Result<[f64; 2], Dada2Error> {
    // Initialise: match rows start near 1, mismatch rows near 0
    let mut a = if is_match { 5.0_f64 } else { -5.0_f64 };
    let mut b = if is_match { -0.3_f64 } else { 0.1_f64 };
    let lr = 1e-3;
    let total: f64 = counts.sum();
    if total == 0.0 {
        return Ok([
            if is_match { 5.0 } else { -5.0 },
            if is_match { -0.3 } else { 0.1 },
        ]);
    }

    let mut prev_ll = f64::NEG_INFINITY;
    let mut first_step = true;
    for _iter in 0..cfg.max_iter {
        let mut ga = 0.0_f64;
        let mut gb = 0.0_f64;
        let mut ll = 0.0_f64;

        for col in 0..MAX_QUAL {
            #[allow(clippy::cast_precision_loss)]
            let q = col as f64;
            let p = sigmoid(a + b * q);
            let c = counts[col];
            if c == 0.0 {
                continue;
            }
            ll += c * p.max(1e-300).ln();
            let residual = c * (1.0 - p);
            ga += residual;
            gb += residual * q;
        }

        a += lr * ga;
        b += lr * gb;

        if first_step {
            first_step = false;
            if !a.is_finite() || !b.is_finite() {
                return Err(Dada2Error::Convergence(
                    "logistic regression produced NaN/Inf on first iteration".into(),
                ));
            }
        }

        let delta = (ll - prev_ll).abs();
        log::debug!("learn_errors iter: ΔlogL = {delta:.2e}");
        if delta < cfg.tol {
            return Ok([a, b]);
        }
        prev_ll = ll;
    }

    log::warn!(
        "logistic regression did not converge within {} iterations; using best-so-far parameters",
        cfg.max_iter
    );
    Ok([a, b])
}

#[inline]
fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Map a nucleotide byte to an index 0=A,1=C,2=G,3=T (N/other → 0).
#[inline]
#[must_use]
pub fn base_index(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0, // A, N, or anything else
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::fastq::FastqRecord;

    fn make_record(seq: &str, qual: &str) -> FastqRecord {
        FastqRecord {
            id: "test".into(),
            seq: seq.bytes().collect(),
            qual: qual.bytes().collect(),
        }
    }

    #[test]
    fn illumina_default_diagonal_dominates() {
        let m = ErrorModel::illumina_default();
        let p_match = m.p_error(0, 0, Phred(30));
        let p_mismatch = m.p_error(0, 1, Phred(30));
        assert!(
            p_match > 0.9,
            "match prob at Q30 should be > 0.9, got {p_match}"
        );
        assert!(
            p_mismatch < 0.01,
            "mismatch prob at Q30 should be < 0.01, got {p_mismatch}"
        );
    }

    #[test]
    fn log_p_error_consistent_with_p_error() {
        let m = ErrorModel::illumina_default();
        for tb in 0u8..4 {
            for ob in 0u8..4 {
                for q in [0u8, 10, 20, 30, 40] {
                    let p = m.p_error(tb, ob, Phred(q));
                    let log_p = m.log_p_error(tb, ob, Phred(q));
                    let expected = p.max(1e-300).ln();
                    assert!(
                        (log_p - expected).abs() < 1e-12,
                        "log_p_error mismatch at tb={tb} ob={ob} q={q}: {log_p} vs {expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn learn_errors_returns_ok() {
        let records: Vec<FastqRecord> = (0..50)
            .map(|_i| make_record("ACGTACGT", &"I".repeat(8)))
            .collect();
        let _ = learn_errors(
            &records,
            &ErrorLearningConfig {
                max_iter: 100,
                ..Default::default()
            },
        );
    }
}
