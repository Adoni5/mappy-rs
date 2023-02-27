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

## Linting code

To avoid committing code that will fail automated checks, you should install `pre-commit` and its hooks:

```console
pre-commit install
pre-commit run --all-files
```
