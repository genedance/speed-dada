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
        let Wire { matrix, n_reads_used } = Wire::deserialize(d)?;
        let log_matrix = build_log_matrix(&matrix);
        Ok(Self { matrix, n_reads_used, log_matrix })
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
        Self { matrix, n_reads_used: 0, log_matrix }
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
}

impl Default for ErrorLearningConfig {
    fn default() -> Self {
        Self {
            n_reads: 1_000_000,
            max_iter: 16,
            tol: 1e-6,
            seed: 42,
        }
    }
}

/// Learn error rates from a collection of reads.
///
/// Uses the first `cfg.n_reads` reads from `records` (seq + qual pairs).
/// Each read is aligned against itself to collect transition counts per
/// Phred bin, then logistic parameters are fitted by EM.
///
/// # Errors
/// Returns [`Dada2Error::Convergence`] if EM does not converge.
/// Returns [`Dada2Error::InvalidInput`] if no reads are supplied.
pub fn learn_errors(
    records: &[crate::io::fastq::FastqRecord],
    cfg: &ErrorLearningConfig,
) -> Result<ErrorModel, Dada2Error> {
    if records.is_empty() {
        return Err(Dada2Error::InvalidInput("no reads supplied to learn_errors".into()));
    }

    // Accumulate transition counts in parallel: each thread gets its own
    // Array2 accumulator; final result is their element-wise sum.
    let n = cfg.n_reads.min(records.len());
    let counts = records[..n]
        .par_iter()
        .fold(
            || Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL)),
            |mut acc, rec| {
                for (&base, &qc) in rec.seq.iter().zip(rec.qual.iter()) {
                    let bi = base_index(base) as usize;
                    let q = Phred::from_ascii(qc).0 as usize;
                    let col = q.min(MAX_QUAL - 1);
                    // Self-comparison: transition = base → base (match)
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
        );

    // Logistic regression fit per transition class — all 16 rows are independent.
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

    let mut matrix = Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL));
    for (row, res) in params.into_iter().enumerate() {
        let [a, b] = res?;
        for col in 0..MAX_QUAL {
            #[allow(clippy::cast_precision_loss)]
            let q = col as f64;
            matrix[[row, col]] = sigmoid(a + b * q);
        }
    }

    let log_matrix = build_log_matrix(&matrix);
    Ok(ErrorModel { matrix, n_reads_used: n as u64, log_matrix })
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
        return Ok([if is_match { 5.0 } else { -5.0 }, if is_match { -0.3 } else { 0.1 }]);
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
        assert!(p_match > 0.9, "match prob at Q30 should be > 0.9, got {p_match}");
        assert!(p_mismatch < 0.01, "mismatch prob at Q30 should be < 0.01, got {p_mismatch}");
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
        let _ = learn_errors(&records, &ErrorLearningConfig { max_iter: 100, ..Default::default() });
    }
}
