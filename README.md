# Mappy-rs
![A map with a crab on it](img/crab_map.webp)

A multi-threaded minimap2 aligner for python.

## Developers
Start with some Docs on Py03 - https://pyo3.rs/v0.16.4/

In order to build an importable module - 

```
python -m venv .env
source -m .env/bin/activate
pip install maturin
maturin develop
```

_NB Any conda environments cannot be activate to run maturin develop, so make sure you `conda deactivate` beforehand_

Then in your python shell of choice - 

```python
import mappy_rs
aligner = mappy_rs.Aligner(1)
```