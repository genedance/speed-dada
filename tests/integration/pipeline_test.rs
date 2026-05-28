//! End-to-end integration test: synthetic reads → ASV recovery.

use speeddada_core::{
    dada::{dada, DadaConfig},
    derep::derep_fastq,
    error_model::{learn_errors, ErrorLearningConfig, ErrorModel},
    filter::{filter_and_trim, FilterConfig},
    io::fastq::{read_fastq, write_fastq, FastqRecord},
};
use tempfile::NamedTempFile;

/// Generate synthetic reads: `n` copies of `true_seq` with `error_rate` substitution errors.
fn generate_reads(true_seq: &[u8], n: usize, error_rate: f64, seed: u64) -> Vec<FastqRecord> {
    let mut rng = seed;
    let mut records = Vec::with_capacity(n);

    for i in 0..n {
        rng = rng
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let mut seq = true_seq.to_vec();
        let mut qual = vec![b'I'; seq.len()]; // Phred 40

        for j in 0..seq.len() {
            rng = rng
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            #[allow(clippy::cast_precision_loss, clippy::cast_lossless)]
            let p = (rng >> 33) as f64 / f64::from(u32::MAX);
            if p < error_rate {
                // Random substitution
                let bases = [b'A', b'C', b'G', b'T'];
                let alt_idx = ((rng >> 11) & 3) as usize;
                if bases[alt_idx] != seq[j] {
                    seq[j] = bases[alt_idx];
                    qual[j] = b'5'; // Phred 20
                }
            }
        }
        records.push(FastqRecord {
            id: format!("read_{i}"),
            seq,
            qual,
        });
    }
    records
}

#[test]
fn pipeline_recovers_true_asv() {
    let true_seq: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    assert_eq!(true_seq.len(), 160);

    // Generate 1000 reads with 2% error rate
    let records = generate_reads(&true_seq, 1000, 0.02, 42);

    // Write to temp file
    let input = NamedTempFile::new().unwrap();
    write_fastq(input.path(), &records).unwrap();

    // Filter
    let filtered = NamedTempFile::new().unwrap();
    let cfg = FilterConfig {
        trunc_len: 150,
        min_len: 100,
        max_ee: 5.0,
        trunc_q: 0,
        trim_left: 0,
        trim_right: 0,
    };
    let stats = filter_and_trim(&cfg, input.path(), filtered.path()).unwrap();
    assert!(stats.reads_out > 0, "all reads filtered out");

    // Learn errors
    let filtered_records = read_fastq(filtered.path()).unwrap();
    let em_cfg = ErrorLearningConfig {
        n_reads: 1000,
        max_iter: 100,
        tol: 1e-6,
        seed: 42,
        ..Default::default()
    };
    let error_model =
        learn_errors(&filtered_records, &em_cfg).unwrap_or_else(|_| ErrorModel::illumina_default());

    // Dereplicate
    let uniques = derep_fastq(&filtered_records).unwrap();
    assert!(!uniques.is_empty());

    // DADA
    let dada_cfg = DadaConfig {
        omega_a: 1e-10,
        ..Default::default()
    };
    let asvs = dada(&uniques, &error_model, &dada_cfg).unwrap();
    assert!(!asvs.is_empty(), "DADA returned no ASVs");

    // The top ASV should match the true sequence at the truncated length
    let top_asv = &asvs[0];
    let cmp_len = 150.min(true_seq.len()).min(top_asv.sequence.len());
    let matches = top_asv.sequence[..cmp_len]
        .iter()
        .zip(true_seq[..cmp_len].iter())
        .filter(|(a, b)| a == b)
        .count();
    #[allow(clippy::cast_precision_loss)]
    let identity = matches as f64 / cmp_len as f64;
    assert!(
        identity >= 0.95,
        "top ASV identity {identity:.3} < 0.95 vs true sequence"
    );
    assert!(
        top_asv.abundance >= 700,
        "top ASV abundance {} < 700",
        top_asv.abundance
    );
}
