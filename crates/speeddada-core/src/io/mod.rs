//! Streaming I/O for FASTQ and FASTA formats.

pub mod fasta;
pub mod fastq;

pub use fasta::FastaRecord;
pub use fastq::FastqRecord;
