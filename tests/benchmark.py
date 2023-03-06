from mappy_rs import Aligner
import mappy as mp
import pytest
from pathlib import Path


def _check_path_present(path: Path):
    """
    Check that the test path is present

    Parameters
    ----------
    path: Path
        Path to the Fastq directory or file that we are using for the benchmark

    Returns
    -------
    bool
        True if Present, False otherwise
    """
    return path.exists()


_FILE_SUFFIXES = set([".fq", ".fastq", ".fastq.gz", ".fq.qz"])


def _gen_fastq(path: Path):
    """
    Generator returning fastq records from either a directory fo Fastq or a
    single Fastq file
    """
    assert _check_path_present(path)
    if path.is_dir():
        for f in path.iterdir():
            if set(f.suffixes).intersection(_FILE_SUFFIXES):
                yield from mp.fastx_read(str(f))
    else:
        if set(path.suffixes).intersection(_FILE_SUFFIXES):
            yield from mp.fastx_read(str(f))


def align_multi(al):
    """
    Parameters
    ----------
    al: mappy_rs.Aligner
        Multithreaded aligner client
    """
    _ = al.map_batch(
        {"read_id": r_id, "seq": seq}
        for r_id, seq, _ in _gen_fastq(Path("../resources/benchmarking/fastq"))
    )
    for _ in res:
        continue


def align_single(al):
    """
    Parameters
    ----------
    al: mappy.Aligner
        Single threaded aligner client
    """
    for _, seq, _ in _gen_fastq(Path("../resources/benchmarking/fastq")):
        for _ in al.map(seq):
            continue


@pytest.mark.benchmark
@pytest.mark.parametrize("i", [*list(range(1, 6))])
def test_benchmark_multi(i, benchmark):
    al = Aligner("../resources/benchmarking/index/hg38.mmi")
    al.enable_threading(i)
    benchmark.pedantic(align_multi, args=(al,), iterations=1, rounds=1)


@pytest.mark.benchmark
def test_benchmark_single(benchmark):
    al = mp.Aligner("../resources/benchmarking/index/hg38.mmi")
    benchmark.pedantic(align_single, args=(al,), iterations=1, rounds=1)


if __name__ == "__main__":
    al = Aligner("../resources/benchmarking/index/hg38.mmi")
    al.enable_threading(8)
    res = al.map_batch(
        {"read_id": r_id, "seq": seq}
        for r_id, seq, _ in _gen_fastq(Path("../resources/benchmarking/fastq"))
    )
    for x in res:
        pass
