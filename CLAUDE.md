# Project Instructions

## Python Environment

**Always activate the venv before running any Python command:**
```bash
source /Users/seanharkins/projects/espa-surface-reflectance/.venv/bin/activate
```

Use `uv` for all package management. Never use `pip`, `conda`, or bare `pip install`.
```bash
uv pip install <package>
```

The venv is at `.venv` and was created with `uv venv .venv --python 3.12`.

## Build

- Rust core: `cargo test -p lasrc_core`
- Full build with Python bindings: `maturin develop --release -m lasrc_py/Cargo.toml`
