# Mappy-rs
[![build](https://github.com/Adoni5/mappy-rs/actions/workflows/CI.yml/badge.svg)](https://github.com/Adoni5/mappy-rs/actions/workflows/CI.yml)
[![CI](https://github.com/Adoni5/mappy-rs/actions/workflows/check.yml/badge.svg)](https://github.com/Adoni5/mappy-rs/actions/workflows/check.yml)


![A map with a crab on it](img/crab_map.webp)

A multi-threaded minimap2 aligner for python.

## Developers
Start with some Docs on Py03 - https://pyo3.rs/latest/

In order to build an importable module:

```bash
python -m venv .env
source -m .env/bin/activate
pip install ".[tests]"
```

To run the tests:

```bash
# Python
pytest

# Rust
cargo t --no-default-features
```

Then in your python shell of choice:

```python
import mappy_rs
aligner = mappy_rs.Aligner("/path/to/index.mmi")
```
