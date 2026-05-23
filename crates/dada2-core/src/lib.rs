//! dada2-core — pure-Rust implementation of the DADA2 ASV pipeline.
//!
//! Pipeline stages:
//! 1. [`filter`] — quality filtering and adapter trimming
//! 2. [`error_model`] — EM-based error rate learning
//! 3. [`derep`] — dereplication of identical sequences
//! 4. [`dada`] — core DADA denoising algorithm
//! 5. [`merge`] — paired-end read merging
//! 6. [`chimera`] — bimera detection and removal
//! 7. [`taxonomy`] — naive-Bayes k-mer taxonomic classification
//! 8. [`io`] — streaming FASTQ/FASTA I/O

#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod align;
pub mod chimera;
pub mod dada;
pub mod derep;
pub mod error_model;
pub mod filter;
pub mod io;
pub mod merge;
pub mod pool;
pub mod primer;
pub mod quality_profile;
pub mod sequence_table;
pub mod taxonomy;

use thiserror::Error;

/// Unified error type for all dada2-core operations.
#[derive(Debug, Error)]
pub enum Dada2Error {
    /// I/O error reading or writing a file.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Parse error in FASTQ/FASTA/TSV data.
    #[error("parse error: {0}")]
    Parse(String),

    /// EM algorithm failed to converge within the allowed iterations.
    #[error("convergence failure: {0}")]
    Convergence(String),

    /// Invalid parameter value supplied by the caller.
    #[error("invalid input: {0}")]
    InvalidInput(String),
}

/// Phred quality score newtype — prevents confusion with raw u8 bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
pub struct Phred(pub u8);

impl Phred {
    /// Convert ASCII quality character (Phred+33 encoding) to [`Phred`].
    #[must_use]
    pub fn from_ascii(c: u8) -> Self {
        Self(c.saturating_sub(33))
    }

    /// Return the error probability P = 10^(-Q/10).
    #[must_use]
    pub fn error_prob(self) -> f64 {
        10f64.powf(-f64::from(self.0) / 10.0)
    }
}

/// k-mer hash newtype — prevents k-mer integers being mixed with counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash,
         serde::Serialize, serde::Deserialize)]
pub struct Kmer(pub u64);
