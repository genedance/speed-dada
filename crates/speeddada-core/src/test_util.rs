//! Shared test-only utilities.
//!
//! Compiled only under `#[cfg(test)]`; not part of the public API.

use std::io::Write;
use tempfile::NamedTempFile;

/// Write FASTQ records (id, seq, qual) to a fresh `NamedTempFile` and return it.
///
/// Replaces duplicated `write_temp_fastq` / `write_fastq` helpers that used
/// to live in `filter.rs::tests` and `primer.rs::tests` respectively.
pub(crate) fn write_test_fastq(records: &[(&str, &str, &str)]) -> NamedTempFile {
    let mut f = NamedTempFile::new().expect("create NamedTempFile");
    for (id, seq, qual) in records {
        writeln!(f, "@{id}\n{seq}\n+\n{qual}").expect("write FASTQ record");
    }
    f
}
