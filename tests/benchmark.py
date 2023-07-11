import mappy_rs
import mappy as mp
import pytest
from pathlib import Path

RESOURCES = (
    Path(__file__).parent.resolve().parent.resolve() / "resources/benchmarking"
)
FASTQ_PATH = RESOURCES / "fastq"
INDEX_PATH = RESOURCES / "index/hg38_simple.mmi"
_FILE_SUFFIXES = set([".fq", ".fastq", ".fastq.gz", ".fq.qz"])


def _gen_fastq(path: Path):
    """
    Generator returning fastq records from either a directory fo Fastq or a
    single Fastq file
    """
    if not path.exists():
        raise RuntimeError("Path doesn't exist")
    if path.is_dir():
        for f in path.iterdir():
            if set(f.suffixes).intersection(_FILE_SUFFIXES):
                yield from mp.fastx_read(str(f))
    else:
        if set(path.suffixes).intersection(_FILE_SUFFIXES):
            yield from mp.fastx_read(str(f))


N_READS = sum(1 for _ in _gen_fastq(FASTQ_PATH))


@pytest.fixture(scope="module")
def fasta():
    return [s for _, s, _ in _gen_fastq(FASTQ_PATH)]


@pytest.fixture
def mappy_al_rs():
    return mappy_rs.Aligner(str(INDEX_PATH))


@pytest.fixture
def mappy_al():
    yield mp.Aligner(str(INDEX_PATH))


@pytest.fixture
def aligner(request):
    yield request.getfixturevalue(request.param)


def align_multi(al, fasta):
    """
    Parameters
    ----------
    al: mappy_rs.Aligner
        Multithreaded aligner client
    fasta : list[str]
        List of sequences to align
    """
    res = al.map_batch({"seq": seq} for seq in fasta)
    return sum(1 for _ in res)


def _align(al, seqs):
    n = 0
    for s in seqs:
        _ = list(al.map(s))
        n += 1
    return n


@pytest.mark.parametrize("aligner", ["mappy_al", "mappy_al_rs"], indirect=True)
def test_classic_mappy(benchmark, aligner, fasta):
    n = benchmark.pedantic(
        _align, args=(aligner, fasta), iterations=1, rounds=1
    )
    assert N_READS == n
    print("Finished classic mappy round")


@pytest.mark.parametrize("threads", list(range(1, 6)))
def test_benchmark_multi(threads, benchmark, mappy_al_rs, fasta):
    mappy_al_rs.enable_threading(threads)
    n = benchmark.pedantic(
        align_multi, args=(mappy_al_rs, fasta), iterations=1, rounds=1
    )
    assert N_READS == n
    print("Finished threaded mappy round")
