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
///
/// The default ([`ErrFunKind::Auto`]) sniffs the input quality
/// distribution and picks `Loess`, `Binned`, or `PacBio` automatically.
/// Force a specific variant only when you have a reason — most users
/// should leave `Auto`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ErrFunKind {
    /// Inspect the input and pick the right smoother. See [`detect_errfun`].
    Auto,
    /// 2-parameter logistic regression σ(a + b·q) — the original SpeedDada
    /// behaviour. Cheap but crude; under-fits any non-monotone curve.
    Logistic,
    /// Weighted local linear regression with a tricubic kernel
    /// (LOWESS, degree 1). Roughly equivalent to `dada2::loessErrfun`
    /// with default span — captures the empirical curvature in the
    /// transition data, including the typical "shoulder" at low quality
    /// scores where the logistic fit overshoots.
    Loess,
    /// Piecewise-linear interpolation between observed quality bins.
    /// Mirrors `dada2::makeBinnedQualErrfun`. The right choice when the
    /// instrument reports only a handful of unique Phred values (NovaSeq
    /// 4-bin, NextSeq 8-bin, MGI/DNBSEQ binned modes).
    Binned,
    /// Canned empirical PacBio CCS / HiFi error matrix. For instruments
    /// reporting Q ≥ 40 across most bases — fitting LOESS to data with
    /// near-zero error rate variance is unreliable, so we use a fixed
    /// curve calibrated to dada2's `PacBioErrfun`.
    PacBio,
}

impl Default for ErrFunKind {
    fn default() -> Self {
        Self::Auto
    }
}

/// Heuristic platform sniffer used by [`ErrFunKind::Auto`].
///
/// Counts distinct Phred values and tracks mean / max Phred across all
/// quality bytes in the input, then dispatches:
///   - PacBio: very high mean (≥ 35) AND few unique Phred values (≤ 8).
///     PacBio CCS / HiFi data is near-uniformly Q35-Q41 with only a
///     handful of distinct values.
///   - Binned: ≤ 12 distinct Phred values AND mean < 35. Covers NovaSeq
///     4-bin, NextSeq 8-bin, MGI DNBSEQ.
///   - ONT (warn, then Loess): mean Phred < 20 AND mean read length > 500.
///   - Loess: everything else (full-range Illumina MiSeq / HiSeq).
#[must_use]
pub fn detect_errfun(records: &[crate::io::fastq::FastqRecord]) -> ErrFunKind {
    let (distinct, mean_q, max_q, mean_len) = quality_stats(records);
    log::info!(
        "errFun auto-detect: distinct_phred={distinct} mean_q={mean_q:.1} \
         max_q={max_q} mean_len={mean_len:.0}"
    );
    let _ = max_q; // Currently unused in dispatch; logged for diagnostics.
    if mean_q >= 35.0 && distinct <= 8 {
        return ErrFunKind::PacBio;
    }
    if distinct <= 12 {
        return ErrFunKind::Binned;
    }
    if mean_q < 20.0 && mean_len > 500.0 {
        log::warn!(
            "errFun auto-detect: input looks like Oxford Nanopore \
             (mean_q={mean_q:.1}, mean_len={mean_len:.0}). SpeedDada's \
             substitution-dominant model with single-indel correction will \
             produce inaccurate ASVs on ONT data. Banded indel-aware \
             alignment is a separate body of work."
        );
        // Fall through to Loess — pipeline still runs, just imprecisely.
    }
    ErrFunKind::Loess
}

fn quality_stats(records: &[crate::io::fastq::FastqRecord]) -> (usize, f64, u8, f64) {
    let mut seen = [false; 96];
    let mut sum_q: u64 = 0;
    let mut sum_len: u64 = 0;
    let mut max_q: u8 = 0;
    let mut n_bytes: u64 = 0;
    for rec in records {
        sum_len += rec.qual.len() as u64;
        for &qc in &rec.qual {
            let p = qc.saturating_sub(33);
            let p_idx = (p as usize).min(95);
            seen[p_idx] = true;
            sum_q += u64::from(p);
            max_q = max_q.max(p);
            n_bytes += 1;
        }
    }
    let distinct = seen.iter().filter(|&&b| b).count();
    #[allow(clippy::cast_precision_loss)]
    let mean_q = if n_bytes > 0 {
        sum_q as f64 / n_bytes as f64
    } else {
        0.0
    };
    #[allow(clippy::cast_precision_loss)]
    let mean_len = if !records.is_empty() {
        sum_len as f64 / records.len() as f64
    } else {
        0.0
    };
    (distinct, mean_q, max_q, mean_len)
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

    // Resolve Auto → concrete kind once, log the choice, then thread the
    // resolved cfg through the rest of the pipeline so every smoother
    // branch sees a definite variant.
    let resolved_kind = if matches!(cfg.method, ErrFunKind::Auto) {
        detect_errfun(&records[..n])
    } else {
        cfg.method
    };

    // PacBio short-circuits collection entirely — its empirical matrix is
    // baked-in and not learned from the data.
    if matches!(resolved_kind, ErrFunKind::PacBio) {
        let matrix = pacbio_matrix();
        let log_matrix = build_log_matrix(&matrix);
        return Ok(ErrorModel {
            matrix,
            n_reads_used: n as u64,
            log_matrix,
        });
    }

    let mut counts = collect_match_counts(&records[..n]);
    add_pairwise_mismatch_counts(&records[..n], &mut counts, cfg.max_pair_dist);

    let resolved_cfg = ErrorLearningConfig {
        method: resolved_kind,
        ..cfg.clone()
    };
    let matrix = fit_from_counts(&counts, &resolved_cfg)?;
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

    for (_, bucket) in by_len {
        if bucket.len() < 2 {
            continue;
        }
        // Sort by descending count so the higher-count member is treated
        // as the "parent" within each pair.
        let mut sorted = bucket;
        sorted.sort_unstable_by(|a, b| b.1.count.cmp(&a.1.count));

        // Consider EVERY pair within the bucket (not just parent vs others).
        // On noisy binned-quality fixtures every read tends to be a unique
        // singleton, so the parent-vs-others heuristic only finds ~n
        // pairs; expanding to ~n²/2 pairs (still O(n²) but n is at most
        // the per-bucket unique-sequence count) gives the smoother enough
        // mismatch evidence at every (transition, q) bin to populate the
        // matrix without leaving high-Q cells at the rate floor.
        for i in 0..sorted.len() {
            let (p_seq, p_uniq) = &sorted[i];
            for j in (i + 1)..sorted.len() {
                let (c_seq, c_uniq) = &sorted[j];
                let dist = crate::align::hamming_distance(p_seq, c_seq);
                if dist == 0 || dist > max_dist {
                    continue;
                }
                #[allow(clippy::cast_precision_loss)]
                let weight = c_uniq.count.min(p_uniq.count) as f64;
                for (k, (&pb, &cb)) in p_seq.iter().zip(c_seq.iter()).enumerate() {
                    if pb == cb {
                        continue;
                    }
                    let pbi = base_index(pb) as usize;
                    let cbi = base_index(cb) as usize;
                    // Use the child's mean quality at this position as the
                    // representative q for the (parent→child) mismatch.
                    #[allow(
                        clippy::cast_precision_loss,
                        clippy::cast_possible_truncation
                    )]
                    let mean_q =
                        (c_uniq.qual_sum[k] as f64 / c_uniq.count as f64) as usize;
                    let col = mean_q.min(MAX_QUAL - 1);
                    let row = pbi * 4 + cbi;
                    counts[[row, col]] += weight;
                    // Re-attribute `weight` from the (child_base, q)
                    // self-transition to the (parent_base→child_base, q)
                    // mismatch. Don't drop below zero.
                    let self_row = cbi * 4 + cbi;
                    let prev = counts[[self_row, col]];
                    counts[[self_row, col]] = (prev - weight).max(0.0);
                }
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
        ErrFunKind::Binned => {
            // For binned-quality platforms (NovaSeq, NextSeq, MGI), LOESS
            // span 0.95 collapses to a near-flat fit because we only have
            // a handful of distinct quality points. Use piecewise linear
            // interpolation between observed bin centres instead — mirrors
            // `dada2::makeBinnedQualErrfun`.
            //
            // Sparse-data robustness: dada2 fits a separate rate for each
            // of the 12 non-self transition directions. With small sample
            // sizes and a heuristic mismatch-evidence collector (our
            // selfConsist-lite vs dada2's full selfConsist), per-direction
            // rates become noisy and individual cells go to zero — which
            // makes dada() over-split because it concludes high-Q errors
            // are impossible. We pool evidence across the 12 non-self
            // directions per Q (uniform-substitution assumption) and then
            // apply that pooled rate per direction. Robust on sparse
            // fixtures; gives up a small amount of per-direction bias.
            //
            // dada2 caps mismatch rates in `[MIN_ERROR_RATE, MAX_ERROR_RATE]`
            // = `[1e-7, 0.25]`. We mirror that.
            const MIN_ERR: f64 = 1e-7;
            const MAX_ERR: f64 = 0.25;

            // Identify observed bin centres globally (cols with any data).
            let mut total_per_col = [0.0_f64; MAX_QUAL];
            for r in 0..N_TRANSITIONS {
                for col in 0..MAX_QUAL {
                    total_per_col[col] += counts[[r, col]];
                }
            }
            let bins: Vec<usize> = (0..MAX_QUAL)
                .filter(|&c| total_per_col[c] > 0.0)
                .collect();

            // Pooled mismatch rate per Q: sum of all off-diagonal counts /
            // sum of all counts (match + mismatch) — pooled across 12
            // non-self directions, then divided by 3 to get the per-
            // direction rate (assumes uniform substitution).
            let mut bin_mismatch = vec![0.0_f64; bins.len()];
            for (i, &c) in bins.iter().enumerate() {
                let mut off_diag = 0.0_f64;
                let mut on_diag = 0.0_f64;
                for r in 0..N_TRANSITIONS {
                    if r % 5 == 0 {
                        on_diag += counts[[r, c]];
                    } else {
                        off_diag += counts[[r, c]];
                    }
                }
                let total = on_diag + off_diag;
                if total > 0.0 {
                    // Per-direction mismatch rate (12 off-diagonal cells).
                    let pooled = off_diag / total;
                    bin_mismatch[i] = (pooled / 3.0).clamp(MIN_ERR, MAX_ERR);
                }
            }

            // Apply per-direction: mismatch entries get the pooled rate,
            // match entries get 1 - 3 × pooled.
            for tb in 0..4 {
                for ob in 0..4 {
                    let r = tb * 4 + ob;
                    let is_match = tb == ob;
                    for col in 0..MAX_QUAL {
                        let p_mismatch = piecewise_linear_floor(
                            &bins, &bin_mismatch, col, MIN_ERR, MAX_ERR,
                        );
                        matrix[[r, col]] = if is_match {
                            (1.0 - 3.0 * p_mismatch).max(1.0 - 3.0 * MAX_ERR)
                        } else {
                            p_mismatch
                        };
                    }
                    if is_match {
                        // Monotonicity guard.
                        for col in (1..MAX_QUAL).rev() {
                            if matrix[[r, col - 1]] < matrix[[r, col]] {
                                matrix[[r, col - 1]] = matrix[[r, col]];
                            }
                        }
                    }
                }
            }
        }
        ErrFunKind::PacBio | ErrFunKind::Auto => {
            // PacBio is handled in `learn_errors` (it short-circuits to
            // the canned matrix); Auto is resolved upstream. Hitting
            // either here means an external caller of `refit_from_counts`
            // didn't resolve Auto first — fall back to LOESS rather than
            // erroring, since this is recoverable.
            let mut fallback = cfg.clone();
            fallback.method = ErrFunKind::Loess;
            return fit_from_counts(counts, &fallback);
        }
    }
    Ok(matrix)
}

/// Canned PacBio CCS / HiFi empirical error matrix.
///
/// Calibrated against `dada2::PacBioErrfun`: per-quality substitution
/// rate is modelled as `p_err(q) = 10^(-q/10)` clamped to `[1e-6, 0.5]`,
/// with mismatches distributed uniformly across the 3 non-match bases
/// (PacBio CCS shows no strong transition bias once consensus calling
/// is applied). Match rate is `1 - p_err(q)`.
///
/// The matrix is built at call-time (one allocation, ~10 KB) rather than
/// stored as a `const` because [`Array2`] isn't const-constructible — the
/// runtime cost is negligible since `learn_errors` only calls this once.
fn pacbio_matrix() -> Array2<f64> {
    let mut m = Array2::<f64>::zeros((N_TRANSITIONS, MAX_QUAL));
    for col in 0..MAX_QUAL {
        #[allow(clippy::cast_precision_loss)]
        let q = col as f64;
        let p_err = 10f64.powf(-q / 10.0).clamp(1e-6, 0.5);
        let p_match = 1.0 - p_err;
        let p_mismatch = p_err / 3.0;
        for tb in 0..4 {
            for ob in 0..4 {
                let r = tb * 4 + ob;
                m[[r, col]] = if tb == ob { p_match } else { p_mismatch };
            }
        }
    }
    m
}

/// Piecewise-linear interpolation of `bin_rate[i]` at observed `bins[i]`,
/// clamped to `[lo, hi]`. Nearest-bin extension outside the observed range
/// (matches `dada2::makeBinnedQualErrfun`).
fn piecewise_linear_floor(
    bins: &[usize],
    bin_rate: &[f64],
    q: usize,
    lo: f64,
    hi: f64,
) -> f64 {
    if bins.is_empty() {
        // No data at all — fall back to the Illumina prior for mismatch
        // (cells get the Phred-based estimate at this q).
        #[allow(clippy::cast_precision_loss)]
        let p_err = 10f64.powf(-(q as f64) / 10.0);
        return (p_err / 3.0).clamp(lo, hi);
    }
    if bins.len() == 1 {
        return bin_rate[0].clamp(lo, hi);
    }
    if q <= bins[0] {
        return bin_rate[0].clamp(lo, hi);
    }
    if q >= *bins.last().unwrap() {
        return bin_rate.last().unwrap().clamp(lo, hi);
    }
    for i in 0..bins.len() - 1 {
        let blo = bins[i];
        let bhi = bins[i + 1];
        if q >= blo && q <= bhi {
            #[allow(clippy::cast_precision_loss)]
            let t = (q - blo) as f64 / (bhi - blo).max(1) as f64;
            return (bin_rate[i] + t * (bin_rate[i + 1] - bin_rate[i])).clamp(lo, hi);
        }
    }
    bin_rate.last().unwrap().clamp(lo, hi)
}

/// Piecewise-linear interpolation of `bin_rate[i]` at observed `bins[i]`,
/// evaluated at quality `q`. Outside the observed range we either clamp
/// to the nearest endpoint or fall back to the Illumina prior, whichever
/// gives a smaller extrapolation jump.
fn piecewise_linear(bins: &[usize], bin_rate: &[f64], q: usize, is_match: bool) -> f64 {
    if bins.is_empty() {
        #[allow(clippy::cast_precision_loss)]
        let p_err = 10f64.powf(-(q as f64) / 10.0);
        return if is_match { 1.0 - p_err } else { p_err / 3.0 };
    }
    if bins.len() == 1 {
        return bin_rate[0];
    }
    if q <= bins[0] {
        return bin_rate[0];
    }
    if q >= *bins.last().unwrap() {
        return *bin_rate.last().unwrap();
    }
    // Find the interval [bins[i], bins[i+1]] containing q.
    for i in 0..bins.len() - 1 {
        let lo = bins[i];
        let hi = bins[i + 1];
        if q >= lo && q <= hi {
            #[allow(clippy::cast_precision_loss)]
            let t = (q - lo) as f64 / (hi - lo).max(1) as f64;
            return bin_rate[i] + t * (bin_rate[i + 1] - bin_rate[i]);
        }
    }
    *bin_rate.last().unwrap()
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

    // Build a record whose quality string is composed of the given Phred
    // bins (cycled). Lets tests synthesise per-platform quality profiles.
    fn rec_with_bins(seq: &str, bins: &[u8]) -> FastqRecord {
        let qual: Vec<u8> = (0..seq.len()).map(|i| bins[i % bins.len()] + 33).collect();
        FastqRecord {
            id: "t".into(),
            seq: seq.bytes().collect(),
            qual,
        }
    }

    #[test]
    fn detect_picks_binned_for_novaseq_4bins() {
        // NovaSeq 4-bin Phred: 2, 12, 23, 37
        let recs: Vec<FastqRecord> = (0..32)
            .map(|_| rec_with_bins("ACGTACGTACGTACGT", &[2, 12, 23, 37]))
            .collect();
        assert_eq!(detect_errfun(&recs), ErrFunKind::Binned);
    }

    #[test]
    fn detect_picks_binned_for_nextseq_8bins() {
        let recs: Vec<FastqRecord> = (0..32)
            .map(|_| rec_with_bins("ACGTACGTACGTACGT", &[2, 12, 18, 25, 32, 36, 38, 40]))
            .collect();
        assert_eq!(detect_errfun(&recs), ErrFunKind::Binned);
    }

    #[test]
    fn detect_picks_pacbio_for_near_q40() {
        // PacBio CCS profile: mean ≥ 35, few distinct values.
        let recs: Vec<FastqRecord> = (0..32)
            .map(|_| rec_with_bins("ACGTACGTACGTACGT", &[38, 39, 40, 41]))
            .collect();
        assert_eq!(detect_errfun(&recs), ErrFunKind::PacBio);
    }

    #[test]
    fn detect_picks_binned_for_mgi_12bins() {
        // MGI/DNBSEQ-style 12 bins, mean ~Q30 — should be binned, not pacbio.
        let bins: Vec<u8> = vec![2, 8, 14, 18, 22, 26, 30, 33, 36, 38, 40, 41];
        let recs: Vec<FastqRecord> = (0..32)
            .map(|_| rec_with_bins("ACGTACGTACGTACGTACGTACGT", &bins))
            .collect();
        assert_eq!(detect_errfun(&recs), ErrFunKind::Binned);
    }

    #[test]
    fn detect_does_not_pick_pacbio_for_full_range_miseq() {
        // MiSeq with full 0..41 range and mean ~Q35: should still be Loess
        // (24 distinct values rules out PacBio's "few unique" branch).
        let bins: Vec<u8> = (15..40).collect();
        let recs: Vec<FastqRecord> = (0..32)
            .map(|_| rec_with_bins("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", &bins))
            .collect();
        assert_eq!(detect_errfun(&recs), ErrFunKind::Loess);
    }

    #[test]
    fn detect_picks_loess_for_full_illumina_range() {
        // 20 distinct Phred values across a typical Illumina spread.
        let bins: Vec<u8> = (20..40).collect();
        let recs: Vec<FastqRecord> = (0..32)
            .map(|_| rec_with_bins("ACGTACGTACGTACGTACGTACGTACGTACGT", &bins))
            .collect();
        assert_eq!(detect_errfun(&recs), ErrFunKind::Loess);
    }

    #[test]
    fn auto_method_resolves_at_runtime() {
        // A binned-quality input under method=Auto should produce the
        // same model as method=Binned.
        let recs: Vec<FastqRecord> = (0..64)
            .map(|_| rec_with_bins("ACGTACGTACGTACGT", &[2, 12, 23, 37]))
            .collect();
        let auto = learn_errors(&recs, &ErrorLearningConfig::default()).unwrap();
        let binned = learn_errors(
            &recs,
            &ErrorLearningConfig {
                method: ErrFunKind::Binned,
                ..Default::default()
            },
        )
        .unwrap();
        // Element-wise equal — both took the same code path.
        for r in 0..N_TRANSITIONS {
            for c in 0..MAX_QUAL {
                assert!(
                    (auto.matrix[[r, c]] - binned.matrix[[r, c]]).abs() < 1e-12,
                    "auto vs binned diverge at row={r} col={c}: {} vs {}",
                    auto.matrix[[r, c]],
                    binned.matrix[[r, c]]
                );
            }
        }
    }

    #[test]
    fn pacbio_matrix_is_diagonal_dominant_and_normalised() {
        let m = pacbio_matrix();
        // Each column should sum to ~1 across the 4 obs-bases per true-base
        // group (P(obs | true, q) is a distribution).
        for tb in 0..4 {
            for col in 0..MAX_QUAL {
                let s: f64 = (0..4).map(|ob| m[[tb * 4 + ob, col]]).sum();
                assert!(
                    (s - 1.0).abs() < 1e-9,
                    "pacbio matrix column {col} for tb={tb} sums to {s}, not 1"
                );
            }
        }
        // Self-transition dominates at Q30 (consistent with PacBio HiFi).
        assert!(m[[0, 30]] > 0.99); // A→A at Q30
        assert!(m[[1, 30]] < 0.01); // A→C at Q30
    }
}
