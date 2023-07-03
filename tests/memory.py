"""Run a test mapping infinitely the same reads infinitely"""
from mappy_rs import Aligner
import mappy as mp
from itertools import cycle, islice
from pathlib import Path
from pyfastx import Fastx
import time
import minimappers2


def yield_reads(
    fastq_path: Path = Path(
        "resources/benchmarking/fastq/fastq_runid_"
        "daf898996f80d549e26d8d49742e033decf51fdb_0.fastq.gz"
    ),
) -> str:
    """
    Return an infinite generator cycling through the same reads
    provided by the given fastq
    Parameters
    ----------
    fastq_path: Path
        The path to a fastq file
    Yields
    ------
    str
        Endless fasta formatted records from the fastq path
    """
    fx = Fastx(str(fastq_path.resolve()))
    fastas = [
        {"seq": seq, "Notes": "testing", "read_id": name}
        for name, seq, qual in fx
    ]
    yield from cycle(fastas)


def batched(iterable, n):
    "Batch data into tuples of length n. The last batch may be shorter."
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def single_threaded():
    """Single threaded batches of 200 reads in mappy rs"""
    al = Aligner("resources/chr20_hg38_22.mmi")
    # al.enable_threading(8)
    for i, fa_records in enumerate(
        batched(
            yield_reads(),
            200,
        )
    ):
        for record in fa_records:
            alignments = al.map(record["seq"])
            print(len(list(alignments)))


def mappy_align():
    """Align with mappy and test for memory"""
    al = mp.Aligner("resources/chr20_hg38_22.mmi")
    # al.enable_threading(8)
    for i, fa_records in enumerate(
        batched(
            yield_reads(Path()),
            200,
        )
    ):
        for record in fa_records:
            alignments = al.map(record["seq"])
            print(len(list(alignments)))


def mappy_rs_no_op_align():
    """Align returning dummy data rather calling out to minimap2-rs"""
    al = Aligner("resources/chr20_hg38_22.mmi")
    # al.enable_threading(8)
    for i, fa_records in enumerate(
        batched(
            yield_reads(Path()),
            200,
        )
    ):
        for record in fa_records:
            alignments = al.map_no_op(record["seq"])
            print(len(list(alignments)))


def mappy_rs_threaded_no_op_align():
    """Align multithreaded returning dummy
    data rather calling out to minimap2-rs"""
    al = Aligner("resources/chr20_hg38_22.mmi")
    al.enable_threading(8)
    for i, fa_records in enumerate(
        batched(
            yield_reads(Path()),
            200,
        )
    ):
        alignments = al.map_batch(fa_records)
        for alignment in alignments:
            del alignment


def map_with_minimappers2_single():
    """Including minimapper2 for completeness"""
    aligner = minimappers2.map_ont()
    aligner.index("resources/chr20_hg38_22.mmi")
    sequences = []
    for i, fa_records in enumerate(
        batched(
            yield_reads(Path()),
            200,
        )
    ):
        sequences = [
            minimappers2.Sequence(str(j), record["seq"])
            for (j, record) in enumerate(fa_records)
        ]
        alignments = aligner.map(sequences)
        print(alignments.shape)


def map_with_minimappers2_multi():
    """AHHHHHHH"""
    aligner = minimappers2.map_ont()
    aligner.threads(4)
    aligner.index("resources/chr20_hg38_22.mmi")
    sequences = []
    for i, fa_records in enumerate(
        batched(
            yield_reads(Path()),
            200,
        )
    ):
        sequences = [
            minimappers2.Sequence(str(j), record["seq"])
            for (j, record) in enumerate(fa_records)
        ]
        alignments = aligner.map(sequences)
        print(alignments.shape)


# @profile
def mappy_rs_threaded():
    al = Aligner("resources/chr20_hg38_22.mmi")
    al.enable_threading(8)
    for i, fa_records in enumerate(
        batched(
            yield_reads(),
            250,
        )
    ):
        alignments = al.map_batch(fa_records)
        for _alignment in alignments:
            continue
        time.sleep(0.1)
        # if i == 10:
        #     break

    # qprint(len(list(alignments)))


if __name__ == "__main__":
    # main()
    mappy_rs_threaded()
