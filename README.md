# Mappy-rs
[![build](https://github.com/Adoni5/mappy-rs/actions/workflows/CI.yml/badge.svg)](https://github.com/Adoni5/mappy-rs/actions/workflows/CI.yml)
[![CI](https://github.com/Adoni5/mappy-rs/actions/workflows/check.yml/badge.svg)](https://github.com/Adoni5/mappy-rs/actions/workflows/check.yml)


![A map with a crab on it](https://github.com/Adoni5/mappy-rs/blob/main/img/crab_map.webp)

A multi-threaded minimap2 aligner for python. Built for [readfish](https://github.com/LooseLab/readfish/) compatibility.

Heavily leaning on and inspired by Joeseph Guhlin's minimap2-rs [repository](https://github.com/jguhlin/minimap2-rs). They also have a more heavily featured python client, which also provides multithreaded alignment. This client provides a more simple streaming interface for use in pipelines.

`pip install mappy-rs`

## Developers
Start with some Docs on Py03 - https://pyo3.rs/latest/

If you wish to contribute, have a look at [CONTRIBUTING.md](.github/CONTRIBUTING.md)

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

The current iteration of `mappy-rs` serves as a drop in for `mappy`, implementing all the same methods. However if this is the use case, you may well be better off using `mappy`, as the extra level of Rust between your python and C++ may well add slightly slower performance.

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
Currently, the maximum batch size to be iterated in one call is 20000.

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

### Benchmarks

A simple benchmark against classic mappy, and mappy_rs with incrementing numbers of threads, run on a 2018 Macbook.
#### __Device__
| Property                     | Value                      |
|------------------------------|----------------------------|
| Model Name                   | MacBook Pro                |
| Model Identifier             | MacBookPro15,2             |
| Processor Name               | Quad-Core Intel Core i7    |
| Processor Speed              | 2.7 GHz                    |
| Number of Processors         | 1                          |
| Total Number of Cores        | 4                          |
| L2 Cache (per Core)          | 256 KB                     |
| L3 Cache                     | 8 MB                       |
| Hyper-Threading Technology   | Enabled                    |
| Memory                       | 16 GB                      |

#### __Results__
Name (time in s)              |  Min           | Max            | Mean          | StdDev        | Median        | IQR          | Outliers     | OPS           | Rounds        | Iterations
------------------------------|---------------|----------------|---------------|--------------|---------------|--------------|--------------|---------------|---------------|------------
test_benchmark_multi[5]       | 26.8900 (1.0) | 30.0969 (1.0)  | 28.0622 (1.0) | 1.2614 (1.0) | 27.9017 (1.0) | 1.6081 (1.35)| 1;0          | 0.0356 (1.0)  | 5             | 1
test_benchmark_multi[4]       | 28.5573 (1.06)| 43.4543 (1.44) | 32.3371 (1.15)| 6.2815 (4.98)| 29.7480 (1.07)| 5.2148 (4.37)| 1;1          | 0.0309 (0.87) | 5             | 1
test_benchmark_multi[3]       | 31.6497 (1.18)| 36.9986 (1.23) | 33.5103 (1.19)| 2.0542 (1.63)| 32.8415 (1.18)| 1.9576 (1.64)| 1;0          | 0.0298 (0.84) | 5             | 1
test_benchmark_multi[2]       | 43.2616 (1.61)| 86.3859 (2.87) | 53.8572 (1.92)| 18.3339 (14.53)| 45.9328 (1.65)| 14.6382 (12.26)| 1;1          | 0.0186 (0.52) | 5             | 1
test_classic_mappy[mappy_al]  | 78.5566 (2.92)| 82.8876 (2.75) | 79.6177 (2.84)| 1.8343 (1.45)| 78.8350 (2.83)| 1.1938 (1.0) | 1;1          | 0.0126 (0.35) | 5             | 1
test_classic_mappy[mappy_al_rs]| 83.7239 (3.11)| 87.9675 (2.92) | 85.4424 (3.04)| 1.6806 (1.33)| 85.6335 (3.07)| 2.3310 (1.95)| 2;0          | 0.0117 (0.33) | 5             | 1
test_benchmark_multi[1]       | 84.8418 (3.16)| 94.0907 (3.13) | 86.7404 (3.09)| 4.1096 (3.26)| 84.8749 (3.04)| 2.4310 (2.04)| 1;1          | 0.0115 (0.32) | 5             | 1


# Changelog

## 0.0.6
- Lowered backoff time for `map_batch` to 50 milliseconds, with 6 attempts. Each attempt will double the previous back off time.
- Improved error handling for `map_batch`, now will raise a `RuntimeError` if the backoff time is exceeded. Also prevented logging `Internal error returning data, the receiver iterator has finished. sending on a disconnected channel 2870` to stderr excessively.
- Added tests for mapping more than 50000 reads, using `back_off=True` and `back_off=False`
