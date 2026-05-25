use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use dada2_core::{
    dada::{Asv, dada, DadaConfig},
    derep::UniqueSeq,
    error_model::{ErrorModel, ErrorLearningConfig, learn_errors},
    io::fasta::FastaRecord,
    merge::{MergeConfig, merge_pairs},
    taxonomy::{TaxonomyConfig, TaxonomyDb},
};
use std::collections::HashMap;

// ── helpers ──────────────────────────────────────────────────────────────────

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *seed
}

fn random_seq(len: usize, seed: &mut u64) -> Vec<u8> {
    let bases = b"ACGT";
    #[allow(clippy::cast_possible_truncation)]
    (0..len).map(|_| bases[(lcg(seed) as usize) % 4]).collect()
}

fn make_asvs(n: usize, seq_len: usize, seed: u64) -> Vec<Asv> {
    let mut s = seed;
    (0..n)
        .map(|i| Asv {
            sequence: random_seq(seq_len, &mut s),
            #[allow(clippy::cast_possible_truncation)]
            abundance: (100u32.saturating_sub(i as u32)).max(1),
        })
        .collect()
}

fn make_fastq_records(n: usize, seq_len: usize, error_rate: f64, seed: u64) -> Vec<dada2_core::io::fastq::FastqRecord> {
    let mut s = seed;
    let true_seq = random_seq(seq_len, &mut s);
    (0..n)
        .map(|i| {
            let mut seq = true_seq.clone();
            let mut qual = vec![b'I'; seq_len];
            for j in 0..seq_len {
                #[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
                let p = (lcg(&mut s) >> 33) as f64 / f64::from(u32::MAX);
                if p < error_rate {
                    #[allow(clippy::cast_possible_truncation)]
                    { seq[j] = b"ACGT"[(lcg(&mut s) as usize) % 4]; }
                    qual[j] = b'5';
                }
            }
            dada2_core::io::fastq::FastqRecord { id: format!("r{i}"), seq, qual }
        })
        .collect()
}

// ── merge_pairs ───────────────────────────────────────────────────────────────

fn bench_merge_pairs(c: &mut Criterion) {
    let mut g = c.benchmark_group("merge_pairs");
    for n in [10usize, 30, 60] {
        let fwd = make_asvs(n, 200, 1);
        let rev = make_asvs(n, 200, 2);
        let cfg = MergeConfig { min_overlap: 20, ..Default::default() };
        g.throughput(Throughput::Elements((n * n) as u64));
        g.bench_with_input(BenchmarkId::new("n_asvs", n), &n, |b, _| {
            b.iter(|| merge_pairs(&fwd, &rev, &cfg).unwrap());
        });
    }
    g.finish();
}

// ── taxonomy_build ────────────────────────────────────────────────────────────

fn bench_taxonomy_build(c: &mut Criterion) {
    let mut g = c.benchmark_group("taxonomy_build");
    for n_refs in [50usize, 200, 500] {
        let mut s = 99u64;
        let records: Vec<FastaRecord> = (0..n_refs)
            .map(|i| FastaRecord {
                id: format!("ref{i}"),
                description: None,
                seq: random_seq(400, &mut s),
            })
            .collect();
        let lineages: HashMap<String, Vec<String>> = records
            .iter()
            .map(|r| (r.id.clone(), vec!["Bacteria".into(), "Phylum".into(),
                "Class".into(), "Order".into(), "Family".into(), "Genus".into(), "sp".into()]))
            .collect();
        let cfg = TaxonomyConfig::default();
        g.throughput(Throughput::Elements(n_refs as u64));
        g.bench_with_input(BenchmarkId::new("n_refs", n_refs), &n_refs, |b, _| {
            b.iter(|| TaxonomyDb::build(&records, &lineages, &cfg).unwrap());
        });
    }
    g.finish();
}

// ── taxonomy_classify ────────────────────────────────────────────────────────

fn bench_taxonomy_classify(c: &mut Criterion) {
    let mut g = c.benchmark_group("taxonomy_classify");
    let mut s = 77u64;
    let records: Vec<FastaRecord> = (0..100)
        .map(|i| FastaRecord { id: format!("ref{i}"), description: None, seq: random_seq(400, &mut s) })
        .collect();
    let lineages: HashMap<String, Vec<String>> = records
        .iter()
        .map(|r| (r.id.clone(), vec!["Bacteria".into(), "Phylum".into(),
            "Class".into(), "Order".into(), "Family".into(), "Genus".into(), "sp".into()]))
        .collect();
    let cfg = TaxonomyConfig::default();
    let db = TaxonomyDb::build(&records, &lineages, &cfg).unwrap();

    for n_queries in [10usize, 50, 200] {
        let seqs: Vec<Vec<u8>> = (0..n_queries).map(|_| random_seq(400, &mut s)).collect();
        g.throughput(Throughput::Elements(n_queries as u64));
        g.bench_with_input(BenchmarkId::new("n_queries", n_queries), &n_queries, |b, _| {
            b.iter(|| db.classify(&seqs, &cfg).unwrap());
        });
    }
    g.finish();
}

// ── dada_denoise ─────────────────────────────────────────────────────────────

fn bench_dada(c: &mut Criterion) {
    let mut g = c.benchmark_group("dada_denoise");
    for n_reads in [200usize, 500, 1000] {
        let records = make_fastq_records(n_reads, 150, 0.02, 42);
        let em_cfg = ErrorLearningConfig { n_reads, max_iter: 20, tol: 1e-4, seed: 42 };
        let error_model = learn_errors(&records, &em_cfg)
            .unwrap_or_else(|_| ErrorModel::illumina_default());
        let uniques: Vec<UniqueSeq> = {
            use std::collections::HashMap;
            let mut map: HashMap<Vec<u8>, (u32, Vec<f64>)> = HashMap::new();
            for r in &records {
                let e = map.entry(r.seq.clone()).or_insert_with(|| (0, vec![0.0; r.seq.len()]));
                e.0 += 1;
                for (acc, &q) in e.1.iter_mut().zip(r.qual.iter()) {
                    *acc += f64::from(q.saturating_sub(33));
                }
            }
            map.into_iter()
                .map(|(seq, (count, qual_sum))| UniqueSeq { seq, count, qual_sum })
                .collect()
        };
        let cfg = DadaConfig { omega_a: 1e-10, ..Default::default() };
        g.throughput(Throughput::Elements(n_reads as u64));
        g.bench_with_input(BenchmarkId::new("n_reads", n_reads), &n_reads, |b, _| {
            b.iter(|| dada(&uniques, &error_model, &cfg).unwrap());
        });
    }
    g.finish();
}

criterion_group!(benches, bench_merge_pairs, bench_taxonomy_build, bench_taxonomy_classify, bench_dada);
criterion_main!(benches);
