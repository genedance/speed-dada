//! FASTA reference database reader.

use crate::Dada2Error;
use std::path::Path;

/// A single FASTA record.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Sequence identifier (without leading `>`).
    pub id: String,
    /// Optional description after the first whitespace in the header.
    pub description: Option<String>,
    /// Raw sequence bytes.
    pub seq: Vec<u8>,
}

/// Read all records from a FASTA file.
///
/// # Errors
/// Returns [`Dada2Error::Io`] or [`Dada2Error::Parse`] on failure.
pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, Dada2Error> {
    use needletail::parse_fastx_file;

    let mut reader = parse_fastx_file(path)
        .map_err(|e| Dada2Error::Parse(format!("cannot open {}: {e}", path.display())))?;

    let mut records = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| Dada2Error::Parse(e.to_string()))?;
        let header = std::str::from_utf8(rec.id()).map_err(|e| Dada2Error::Parse(e.to_string()))?;
        let mut parts = header.splitn(2, ' ');
        let id = parts.next().unwrap_or("").to_owned();
        let description = parts.next().map(str::to_owned);
        let seq = rec.seq().to_vec();
        records.push(FastaRecord {
            id,
            description,
            seq,
        });
    }
    Ok(records)
}
