# How To Contribute

## Local Development Environment

You should ensure that you have the rust compiler and toolchain installed (>=v1.58.1), see https://rustup.rs.
We also require Python3 at least at version 3.8.

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

The benchmarking script is found in `tests/benchmarking.py`.
The benchmarking script requires a minimap2 index or reference fasta to map against, and a folder of fastq. Data
The paths to the mapping index and the fastq will have to be set, these variables (`FASTQ_PATH`, `INDEX_PATH`)
are located at the top of the script. `FASTQ_PATH` can be set to a directory containing fastq or gzipped fastq, or a path to a fastq file.

A recommended source of Fastq is
http://s3.amazonaws.com/nanopore-human-wgs/rel6/FASTQTars/FAB42316-216722908_Multi.tar
And a good human index can be made using [seqkit](https://bioinf.shenwei.me/seqkit/) and the [HG38 reference](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)

```console
# Keep only complete contigs, filtering out alts
seqkit -j 8 grep -o "only_primary_hg38.fna" -r -p "^NC" GCF_000001405.40_GRCh38.p14_genomic.fna.gz
minimap2 -d hg38.mmi only_primary_hg38.fna
```

`INDEX_PATH` must be set to a path to an index file.

```console
cd mappy-rs
python3 -m venv .env
source ./.env/bin/activate
pip install --upgrade pip
pip install '.[benchmark]'
pytest tests/benchmark.py
```

Benchmarks can be filtered by name, for either running just the comparison against minimap2s mappy without threading, or just benchmrking, parameterising the number of threads.

For classic/single:
```console
pytest --benchmark-warmup-iterations=1 --maxfail=1 --capture=tee-sys -vk "classic"  benchmark.py
```

For multi:
```console
pytest --benchmark-warmup-iterations=1 --maxfail=1 --capture=tee-sys -vk "multi"  benchmark.py
```

In order to save previous benchmarks for comparison just add `--benchmark-autosave` or `--benchmark-save=some-name`. Full documentation can be found in the pytest-benchmark [docs](https://pytest-benchmark.readthedocs.io/en/latest/comparing.html).

#### Comparisons
Simply run

```console
pytest-benchmark compare 0001 0002
```

Where 0001 0002 are the names of previously saved benhcmarks.

You can also get a nice plot with `--benchmark-histogram`. The result is a modified Tukey box and whisker plot where the outliers (the small bullets) are Min and Max. Note that if you do not supply a name for the plot it is recommended that -`-benchmark-histogram` is the last option passed.

## Linting code

To avoid committing code that will fail automated checks, you should install `pre-commit` and its hooks:
The following will setup pre-commit in the repository, and automatically run when you checkout a branch, try to commit or merge.

```console
pip install pre-commit
pre-commit install -t pre-commit -t post-checkout -t post-merge
pre-commit run --all-files
```
