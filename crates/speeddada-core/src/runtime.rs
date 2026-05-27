//! Hardware-aware runtime configuration.
//!
//! Detects available CPU cores and RAM, then configures rayon's global
//! thread pool for optimal throughput on both x86-64 and `AArch64`.

use std::thread::available_parallelism;

/// Auto-detected or manually overridden parallelism settings.
#[derive(Debug, Clone)]
pub struct RuntimeConfig {
    /// Number of rayon worker threads.
    pub n_threads: usize,
    /// Available RAM at detection time (MiB), if detectable.
    pub mem_available_mb: Option<u64>,
}

impl RuntimeConfig {
    /// Detect optimal settings from the current hardware.
    ///
    /// Thread count is capped to `available_ram_mb / 512` by default.
    /// Use [`Self::detect_with`] to override the per-thread RAM assumption.
    #[must_use]
    pub fn detect() -> Self {
        Self::detect_with(512)
    }

    /// Like [`Self::detect`] but with an explicit per-thread RAM budget in MiB.
    ///
    /// Typical values:
    /// - `512` — conservative default, suits most pipeline stages
    /// - `800` — DADA-heavy workloads where each thread holds a full logp table
    /// - `64`  — taxonomy-only or filter-only runs
    #[must_use]
    pub fn detect_with(mb_per_thread: u64) -> Self {
        let n_cpu = available_parallelism().map_or(1, std::num::NonZeroUsize::get);
        let mem_available_mb = read_available_memory_mb();

        let n_threads = match mem_available_mb {
            Some(mb) => {
                #[allow(clippy::cast_possible_truncation)]
                let by_ram = ((mb / mb_per_thread) as usize).max(1);
                n_cpu.min(by_ram)
            }
            None => n_cpu,
        };

        Self {
            n_threads,
            mem_available_mb,
        }
    }

    /// Override the thread count manually (e.g. for testing or containers).
    #[must_use]
    pub fn with_threads(mut self, n: usize) -> Self {
        self.n_threads = n.max(1);
        self
    }

    /// Apply this config to rayon's global thread pool.
    ///
    /// Must be called before any rayon work begins; subsequent calls are
    /// silently ignored by rayon (the global pool is initialised at most once).
    ///
    /// # Errors
    /// Returns an error if rayon fails to spawn the requested thread count.
    pub fn apply(&self) -> Result<(), rayon::ThreadPoolBuildError> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.n_threads)
            .build_global()
    }

    /// Split a total read budget evenly across `n_files` for error learning.
    ///
    /// Returns how many reads to sample from each file so that the total stays
    /// at most `n_reads`, with a floor of `min_per_file` to ensure statistical
    /// validity of the error model.
    #[must_use]
    pub fn reads_per_file(n_reads: usize, n_files: usize, min_per_file: usize) -> usize {
        if n_files == 0 {
            return n_reads;
        }
        (n_reads / n_files).max(min_per_file)
    }
}

impl Default for RuntimeConfig {
    fn default() -> Self {
        Self::detect()
    }
}

/// Read `MemAvailable` from `/proc/meminfo` (Linux / Raspberry Pi).
/// Returns `None` on non-Linux targets or if the file cannot be parsed.
fn read_available_memory_mb() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let content = std::fs::read_to_string("/proc/meminfo").ok()?;
        for line in content.lines() {
            if let Some(rest) = line.strip_prefix("MemAvailable:") {
                let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
                return Some(kb / 1024);
            }
        }
        None
    }
    #[cfg(not(target_os = "linux"))]
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_returns_at_least_one_thread() {
        let cfg = RuntimeConfig::detect();
        assert!(cfg.n_threads >= 1);
    }

    #[test]
    fn with_threads_clamps_to_one() {
        let cfg = RuntimeConfig::detect().with_threads(0);
        assert_eq!(cfg.n_threads, 1);
    }
}
