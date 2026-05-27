//! Standalone speeddada-core benchmark binary — no Python or R bindings.
//!
//! Internal tool, not part of the public library surface. Relax a few
//! pedantic lints that the library itself respects but that would be
//! noise here (long main fn, doc backtick nitpicks on stage names).
#![allow(
    clippy::doc_markdown,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::redundant_closure_for_method_calls,
    clippy::unnecessary_sort_by,
    clippy::struct_field_names
)]
//!
//! Runs the full pipeline (filter → `learn_errors` → derep → `dada_pseudo`
//! → merge → chimera) on N paired FASTQ samples and emits per-stage timings
//! plus a JSON summary in the same format as `bench_rust.py` /
//! `bench_speeddada.R`.
//!
//! Usage:
//!   speeddada-bench --threads 16 --in-dir <dir> --out-dir <dir>
//!               --samples stem1,stem2,stem3
//!               [--fwd-suffix .1.fq.gz] [--rev-suffix .2.fq.gz]
//!               [--prefix raw.]
//!
//! Each sample's input paths are constructed as:
//!   <in-dir>/<prefix><stem><fwd-suffix>
//!   <in-dir>/<prefix><stem><rev-suffix>

use serde::Serialize;
use speeddada_core::{
    chimera::remove_bimera_denovo,
    dada::{dada_pseudo, Asv, DadaConfig},
    derep::derep_fastq,
    error_model::{learn_errors, ErrorLearningConfig},
    filter::{filter_and_trim_paired_many, FilterConfig, FilterStatsPaired},
    io::fastq::read_fastq,
    merge::{merge_pairs, MergeConfig},
};
use std::{collections::HashMap, path::PathBuf, time::Instant};

fn parse_arg(args: &[String], name: &str) -> Option<String> {
    let key = format!("--{name}");
    args.windows(2).find_map(|w| {
        if w[0] == key {
            Some(w[1].clone())
        } else {
            None
        }
    })
}

#[derive(Serialize)]
struct StageTimings {
    filter_ms: f64,
    learn_errors_ms: f64,
    derep_ms: f64,
    dada_ms: f64,
    merge_ms: f64,
    chimera_ms: f64,
}

#[derive(Serialize)]
struct SampleStats {
    sample: String,
    reads_in: u64,
    reads_out: u64,
}

#[derive(Serialize)]
struct AsvEntry {
    sequence: String,
    abundance: u32,
}

#[derive(Serialize)]
struct Output {
    tool: &'static str,
    n_threads: usize,
    total_ms: f64,
    stages: StageTimings,
    samples: Vec<SampleStats>,
    n_asvs_before_chimera: usize,
    n_asvs_after_chimera: usize,
    total_abundance: u64,
    asvs: Vec<AsvEntry>,
}

fn ms(d: std::time::Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    let n_threads: usize = parse_arg(&args, "threads")
        .unwrap_or_else(|| "16".to_string())
        .parse()?;
    let in_dir = PathBuf::from(parse_arg(&args, "in-dir").expect("--in-dir required"));
    let out_dir = PathBuf::from(parse_arg(&args, "out-dir").expect("--out-dir required"));
    let samples: Vec<String> = parse_arg(&args, "samples")
        .expect("--samples required (comma-separated)")
        .split(',')
        .map(str::to_owned)
        .collect();
    let prefix = parse_arg(&args, "prefix").unwrap_or_else(|| "raw.".to_string());
    let fwd_suffix = parse_arg(&args, "fwd-suffix").unwrap_or_else(|| ".1.fq.gz".to_string());
    let rev_suffix = parse_arg(&args, "rev-suffix").unwrap_or_else(|| ".2.fq.gz".to_string());

    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build_global()
        .ok();

    std::fs::create_dir_all(&out_dir)?;
    println!(
        "[rust-only speeddada-bench] n_threads={n_threads}  samples={}",
        samples.len()
    );

    let fwd_in: Vec<PathBuf> = samples
        .iter()
        .map(|s| in_dir.join(format!("{prefix}{s}{fwd_suffix}")))
        .collect();
    let rev_in: Vec<PathBuf> = samples
        .iter()
        .map(|s| in_dir.join(format!("{prefix}{s}{rev_suffix}")))
        .collect();
    let fwd_filt: Vec<PathBuf> = samples
        .iter()
        .map(|s| out_dir.join(format!("{s}_R1_filt.fastq.gz")))
        .collect();
    let rev_filt: Vec<PathBuf> = samples
        .iter()
        .map(|s| out_dir.join(format!("{s}_R2_filt.fastq.gz")))
        .collect();

    let t_total = Instant::now();

    // 1. Filter — same params as bench_rust.py / bench_speeddada.R on field data.
    println!("[filter_and_trim_paired]");
    let t = Instant::now();
    let cfg_fwd = FilterConfig {
        trunc_len: 240,
        min_len: 50,
        max_ee: 2.0,
        trunc_q: 2,
        trim_left: 0,
        trim_right: 0,
    };
    let cfg_rev = FilterConfig {
        trunc_len: 180,
        min_len: 50,
        max_ee: 4.0,
        trunc_q: 2,
        trim_left: 0,
        trim_right: 0,
    };
    let pairs: Vec<(PathBuf, PathBuf, PathBuf, PathBuf)> = (0..samples.len())
        .map(|i| {
            (
                fwd_in[i].clone(),
                rev_in[i].clone(),
                fwd_filt[i].clone(),
                rev_filt[i].clone(),
            )
        })
        .collect();
    let stats: Vec<FilterStatsPaired> = filter_and_trim_paired_many(&cfg_fwd, &cfg_rev, &pairs)?;
    let t_filter = ms(t.elapsed());
    let total_in: u64 = stats.iter().map(|s| s.reads_in).sum();
    let total_out: u64 = stats.iter().map(|s| s.pairs_out).sum();
    println!("  total_in={total_in}  total_out={total_out}  ({t_filter:.1} ms)");

    // 2. Learn errors — one model fwd, one rev.
    println!("[learn_errors]");
    let t = Instant::now();
    let mut fwd_records = Vec::new();
    for p in &fwd_filt {
        fwd_records.extend(read_fastq(p)?);
    }
    let err_fwd = learn_errors(&fwd_records, &ErrorLearningConfig::default())
        .unwrap_or_else(|_| speeddada_core::error_model::ErrorModel::illumina_default());
    drop(fwd_records);
    let mut rev_records = Vec::new();
    for p in &rev_filt {
        rev_records.extend(read_fastq(p)?);
    }
    let err_rev = learn_errors(&rev_records, &ErrorLearningConfig::default())
        .unwrap_or_else(|_| speeddada_core::error_model::ErrorModel::illumina_default());
    drop(rev_records);
    let t_errors = ms(t.elapsed());
    println!("  done  ({t_errors:.1} ms)");

    // 3. Derep — per sample.
    println!("[derep_fastq]");
    let t = Instant::now();
    let derep_fwd: Vec<Vec<speeddada_core::derep::UniqueSeq>> = fwd_filt
        .iter()
        .map(|p| {
            let recs = read_fastq(p).expect("fwd derep read_fastq");
            derep_fastq(&recs).expect("fwd derep")
        })
        .collect();
    let derep_rev: Vec<Vec<speeddada_core::derep::UniqueSeq>> = rev_filt
        .iter()
        .map(|p| {
            let recs = read_fastq(p).expect("rev derep read_fastq");
            derep_fastq(&recs).expect("rev derep")
        })
        .collect();
    let t_derep = ms(t.elapsed());
    println!("  done  ({t_derep:.1} ms)");

    // 4. DADA pseudo-pool — fwd, then rev.
    println!("[dada pseudo-pool]");
    let t = Instant::now();
    let cfg = DadaConfig::default();
    let refs_fwd: Vec<&[speeddada_core::derep::UniqueSeq]> =
        derep_fwd.iter().map(|v| v.as_slice()).collect();
    let dada_fwd: Vec<Vec<Asv>> = dada_pseudo(&refs_fwd, &err_fwd, &cfg)?;
    let refs_rev: Vec<&[speeddada_core::derep::UniqueSeq]> =
        derep_rev.iter().map(|v| v.as_slice()).collect();
    let dada_rev: Vec<Vec<Asv>> = dada_pseudo(&refs_rev, &err_rev, &cfg)?;
    let t_dada = ms(t.elapsed());
    let n_asv_fwd: usize = dada_fwd.iter().map(Vec::len).sum();
    let n_asv_rev: usize = dada_rev.iter().map(Vec::len).sum();
    println!("  fwd_asvs(total)={n_asv_fwd}  rev_asvs(total)={n_asv_rev}  ({t_dada:.1} ms)");

    // 5. Merge — per sample.
    println!("[merge_pairs]");
    let t = Instant::now();
    let merge_cfg = MergeConfig {
        min_overlap: 12,
        max_mismatches: 0,
        just_concatenate: false,
    };
    let merged_per_sample: Vec<Vec<speeddada_core::merge::MergedRead>> = (0..samples.len())
        .map(|i| merge_pairs(&dada_fwd[i], &dada_rev[i], &merge_cfg).unwrap_or_default())
        .collect();
    let t_merge = ms(t.elapsed());
    let total_merged: usize = merged_per_sample.iter().map(Vec::len).sum();
    println!("  total_merged_asvs={total_merged}  ({t_merge:.1} ms)");

    // 6. Sequence table + chimera removal.
    println!("[chimera + sequence table]");
    let t = Instant::now();
    let mut counts: HashMap<Vec<u8>, u32> = HashMap::new();
    for merged_sample in &merged_per_sample {
        for m in merged_sample {
            *counts.entry(m.sequence.clone()).or_insert(0) += m.abundance;
        }
    }
    let agg: Vec<(Vec<u8>, u32)> = counts.into_iter().collect();
    let clean = remove_bimera_denovo(&agg)?;
    let n_asvs_before = agg.len();
    let n_asvs_after = clean.len();
    let t_chimera = ms(t.elapsed());
    println!("  asvs_in={n_asvs_before}  asvs_out={n_asvs_after}  ({t_chimera:.1} ms)");

    let t_total_ms = ms(t_total.elapsed());
    println!("\nTotal rust-only dada2 time: {t_total_ms:.1} ms");

    let mut asvs: Vec<AsvEntry> = clean
        .into_iter()
        .map(|(s, a)| AsvEntry {
            sequence: String::from_utf8_lossy(&s).into_owned(),
            abundance: a,
        })
        .collect();
    asvs.sort_by(|a, b| b.abundance.cmp(&a.abundance));

    let total_abundance: u64 = asvs.iter().map(|a| u64::from(a.abundance)).sum();

    let sample_stats: Vec<SampleStats> = samples
        .iter()
        .zip(stats.iter())
        .map(|(s, st)| SampleStats {
            sample: s.clone(),
            reads_in: st.reads_in,
            reads_out: st.pairs_out,
        })
        .collect();

    let out = Output {
        tool: "rust-only dada2",
        n_threads,
        total_ms: t_total_ms,
        stages: StageTimings {
            filter_ms: t_filter,
            learn_errors_ms: t_errors,
            derep_ms: t_derep,
            dada_ms: t_dada,
            merge_ms: t_merge,
            chimera_ms: t_chimera,
        },
        samples: sample_stats,
        n_asvs_before_chimera: n_asvs_before,
        n_asvs_after_chimera: n_asvs_after,
        total_abundance,
        asvs,
    };
    let out_path = out_dir.join("rust_only_output.json");
    std::fs::write(&out_path, serde_json::to_string_pretty(&out)?)?;
    println!("\nWrote {}", out_path.display());
    Ok(())
}
