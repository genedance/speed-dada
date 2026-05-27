#!/usr/bin/env bash
# vendor-r-rust.sh — populate r-package/SpeedDada/src/rust/ with the Rust
# crates needed for a standalone R-package source tarball.
#
# Used before `R CMD build r-package/SpeedDada` to produce a distributable
# tarball (Bioconductor submission, CRAN-style installs). Not needed for
# in-tree development; `configure` falls back to the parent workspace then.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PKG_DIR="${REPO_ROOT}/r-package/SpeedDada"
DEST="${PKG_DIR}/src/rust"

echo "vendor-r-rust: removing any previous vendored copy at ${DEST}"
rm -rf "${DEST}"
mkdir -p "${DEST}"

echo "vendor-r-rust: copying speeddada-core and speeddada-r into ${DEST}"
cp -R "${REPO_ROOT}/crates/speeddada-core" "${DEST}/speeddada-core"
cp -R "${REPO_ROOT}/crates/speeddada-r"    "${DEST}/speeddada-r"

# Drop build artefacts that might have been copied.
find "${DEST}" -type d -name target -prune -exec rm -rf {} +

# Write a minimal workspace manifest that mirrors the root one but excludes
# speeddada-py (not needed for the R build).
cat > "${DEST}/Cargo.toml" <<'EOF'
[workspace]
resolver = "2"
members  = ["speeddada-core", "speeddada-r"]

[workspace.package]
version       = "0.99.0"
edition       = "2021"
rust-version  = "1.78"
authors       = ["Alexandre Jousset <info@genedance.com>"]
license       = "MIT"
description   = "High-performance DADA2 amplicon sequence variant pipeline — Rust core for the SpeedDada R package."
homepage      = "https://github.com/Genedance/speed-dada"
repository    = "https://github.com/Genedance/speed-dada"
documentation = "https://genedance.github.io/speed-dada/"
keywords      = ["bioinformatics", "amplicon", "dada2", "microbiome", "asv"]
categories    = ["science"]
readme        = "../../../README.md"

[workspace.dependencies]
rayon          = "1"
needletail     = "0.5"
ndarray        = { version = "0.15", features = ["serde"] }
statrs         = "0.17"
thiserror      = "1"
serde          = { version = "1", features = ["derive"] }
serde_json     = "1"
bincode        = "1"
log            = "0.4"
env_logger     = "0.11"
tempfile       = "3"
criterion      = { version = "0.5", features = ["html_reports"] }
extendr-api    = { version = "0.9", default-features = false }
speeddada-core = { path = "speeddada-core", version = "0.99.0" }

[profile.release]
opt-level     = 3
lto           = "thin"
codegen-units = 1
EOF

echo "vendor-r-rust: done."
echo "Next: R CMD build ${PKG_DIR}"
