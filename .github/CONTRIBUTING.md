# How To Contribute

## Local Development Environment

You should ensure that you have the rust compiler and toolchain installed (>=v1.58.1), see https://rustup.rs.
We also require Python3 at least at version 3.7.

Clone the repository:

```console
git clone git@github.com:Adoni5/mappy-rs.git
```

Change into the newly created directory, create a virtual environment, and install `mappy-rs`:

```console
cd mappy-rs
python3 -m venv .env
source ./.env/bin/activate
pip install --upgrade pip
pip install '.[dev]'
```

After modifying `src/lib.rs` you will need to recompile/install `mappy-rs` by running:

```console
pip install .
```

## Runnning tests

```console
python -m pytest
cargo t --no-default-features
```

## Running Benchmarks

The benchmarking script is found in `tests/benchmarking.py`. The paths to the mapping index and the fastq will have to be set, these variables (`FASTQ_PATH`, `INDEX_PATH`)
are located at the top of the script. `FASTQ_PATH` can be set to a directory containing fastq or gzipped fastq, or a path to a fastq file.
`INDEX_PATH` must be set to a path to an index file.

```console
cd mappy-rs
python3 -m venv .env
source ./.env/bin/activate
pip install --upgrade pip
pip install '.[benchmark]'
pytest tests/benchmark.py
```


## Linting code

To avoid committing code that will fail automated checks, you should install `pre-commit` and its hooks:
The following will setup pre-commit in the repository, and automatically run when you checkout a branch, try to commit or merge.

```console
pip install pre-commit
pre-commit install -t pre-commit -t post-checkout -t post-merge
pre-commit run --all-files
```
