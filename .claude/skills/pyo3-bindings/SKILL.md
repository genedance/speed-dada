# Skill: pyo3-bindings

## Invoke with
"use the pyo3-bindings skill"

## pyproject.toml template
```toml
[build-system]
requires = ["maturin>=1.5,<2"]
build-backend = "maturin"

[tool.maturin]
features = ["pyo3/extension-module"]
module-name = "dada2"
```

## GIL release pattern
```rust
let result = py.allow_threads(|| {
    heavy_rust_computation()
});
```

## Error mapping pattern
```rust
heavy_fn().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?
```

## Build & install
```bash
cd crates/dada2-py
maturin develop          # editable install into current venv
maturin build --release  # build wheel
```
