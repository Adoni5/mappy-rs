# Mappy-rs
[![build](https://github.com/Adoni5/mappy-rs/actions/workflows/CI.yml/badge.svg)](https://github.com/Adoni5/mappy-rs/actions/workflows/CI.yml)
[![CI](https://github.com/Adoni5/mappy-rs/actions/workflows/check.yml/badge.svg)](https://github.com/Adoni5/mappy-rs/actions/workflows/check.yml)


![A map with a crab on it](img/crab_map.webp)

A multi-threaded minimap2 aligner for python. Built for [readfish](https://github.com/LooseLab/readfish/) compatibility.

Heavily leaning on and inspired by Joeseph Guhlin's minimap2-rs [repository](https://github.com/jguhlin/minimap2-rs). They also have a more heavily featured python client, however this one simply multi threads and maps.

`pip install mappy-rs`

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
aligner = mappy_rs.Aligner("resources/test/test.mmi")
```

The current iteration of `mappy-rs` serves as a drop in for `mappy`, implementing all the same methods. However if this is the use case, you may well be better off using `mappy`, as the extra level of Rust betwene your python and C++ may well add slighly slower performance.

### Multithreading
In order to use multi threading, one must first enable it.

```python
import mappy_rs
aligner = mappy_rs.Aligner("resources/test/test.mmi")
# Use 10 threads
aligner.enable_threading(10)
```

Enabling threading makes the `map_batch` method available.
This method requires a list or iterable of dictionaries, which can have any number of keys and depth, but **must** contain the key `seq` with a string value in the top-level dictionary.
Currently, the maximum batch size tobe iterated in one call is 20000.

For example:

```python
import mappy_rs
aligner = mappy_rs.Aligner("resources/test/test.mmi")
aligner.enable_threading(10)
seqs = [
    {"seq": "ACGTAGCATCGAGACTACGA", "Other_random_key": "banter"},
    {"seq": "ACGTAGCATCGAGACTACGA", "Other_random_key": "banter"},
]
for (mapping, data) in aligner.map_batch(seqs):
    print(list(mapping))
    print(data)
```
