# Project Instructions

## Python Environment

Use the `.venv` in the repository root, managed with `uv`. Create it if it doesn't exist:
```bash
uv venv .venv --python 3.12
```

**Always activate the venv before running any Python command:**
```bash
source .venv/bin/activate
```

Use `uv` for all package management. Never use `pip`, `conda`, or bare `pip install`.
```bash
uv pip install <package>
```

## Build

- Rust core: `cargo test -p lasrc_core`
- Full build with Python bindings: `maturin develop --release -m lasrc_py/Cargo.toml`
