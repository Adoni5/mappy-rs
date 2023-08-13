"""python_test.py

This set of tests roughly recreate those that are in lib.rs but use the
compiled package.
"""
from pathlib import Path
import copy
from itertools import repeat

import pytest

import mappy_rs

RESOURCES = Path(__file__).parent.resolve().parent.resolve() / "resources/test"
MMI_FILE = RESOURCES / "test.mmi"
FA_FILE = RESOURCES / "test.fa"


def read_fasta(fh):
    for line in fh:
        if line.startswith(">"):
            name = line[1:].strip()
            break
    fa_lines = []
    for line in fh:
        if line.startswith(">"):
            yield name, "".join(fa_lines)
            fa_lines = []
            name = line[1:].strip()
            continue
        fa_lines.append(line.strip())
    yield name, "".join(fa_lines)


@pytest.fixture
def mmi_file():
    return str(MMI_FILE)


@pytest.fixture
def al(mmi_file):
    return mappy_rs.Aligner(mmi_file)


@pytest.fixture
def fasta(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def fasta_file():
    return str(FA_FILE)


@pytest.fixture
def fasta_list(fasta_file):
    with open(fasta_file, "rt") as fh:
        seqs = [s for _, s in read_fasta(fh)]

    seqs = [
        {"id": i, "seq": seq}
        for i, seq in enumerate(copy.copy(s) for _ in range(10) for s in seqs)
    ]
    return seqs


@pytest.fixture
def fasta_iter(fasta_list):
    return iter(fasta_list)


@pytest.fixture
def fasta_tuple(fasta_list):
    return tuple(fasta_list)


@pytest.fixture
def fasta_generator(fasta_list):
    return (item for item in fasta_list)


def test_test(al):
    assert al


def test_property_k(al):
    assert al.k == 15


def test_property_n_seq(al):
    assert al.n_seq == 4


def test_property_w(al):
    assert al.w == 10


def test_property_seq_names(al):
    expected = [
        "Bacillus_subtilis",
        "Enterococcus_faecalis",
        "Escherichia_coli_1",
        "Escherichia_coli_2",
    ]
    seq_names = al.seq_names
    seq_names.sort()
    assert seq_names == expected


def test_get_seq(al):
    contig = "Bacillus_subtilis"
    expected = (
        "AGAGTGAAGCCAATATTCCGATAACGATTGCTTTCATGATATCCCTCATTCTGGCATTATTTTTTTATA"
        "CTATACTATTCGATATCGCACAGATCAATGGAGTCGTGAGAAAATAAACATGTTTTGCGAACCGCTATG"
        "TGTGGAAGACAAAAAATGGAGGTGAAATTGATGGAAGCAAAGACACAGGCGTACTTTTTTCAGGATGAT"
        "GGCAGGATTCCGAATCACCCTGATTTTCCGCTCGTTGTGTATCAAAACGCACTCAAGGACACCGGTCAG"
        "GCAGAGCGGATCGTCAACCGGCATGGCTGGTCAAACAGCTGGTCGGGGAGTGTTTTTCCATACCATCAT"
        "TATCACAGCAATACGCATGAAGTCCTGATTGCAGTTCGGGGAGAGGCTGTGATTC"
    )
    seq = al.seq(contig)
    assert seq == expected


def test_map_one(al):
    mappings = al.map(
        "AGAGCAGGTAGGATCGTTGAAAAAAGAGTACTCAGGATTCCATTCAACTTTTACTGATTTGAAGCGTAC"
        "TGTTTATGGCCAAGAATATTTACGTCTTTACAACCAATACGCAAAAAAAGGTTCATTGAGTTTGGTTGT"
        "GATTTGATGAAAATTACTGAGAATAACAGGATTATTAAGCTGATTGATGAACTAAATCAGCTTAATAAA"
        "TATTCTTTGCAGATAGGAATATTTGGGGAAAATGATTCTTTTATGGCGATGTTGGCCCAAGTTCATGAA"
        "TTTGGGGTGACTATTCGTCCCAAAGGTCGTTTTCTTGTTATACCACTTATGAAAAAGTATAGAGGTAAA"
        "AGTCCACGTCAATTTGATTTGTTTTTTATGCAAACTAAAGAAAATCACAAGTTTT",
        cs=True,
    )
    assert len(mappings) == 1
    mapping = mappings[0]
    assert mapping.target_start == 0
    assert mapping.target_end == 400


def test_map_batch_100000(al, fasta_iter):
    al.enable_threading(4)
    iter_ = repeat(next(fasta_iter), 100000)
    mappings = al.map_batch(iter_, back_off=True)
    n = 0
    for res in mappings:
        n += 1
    assert n == 100000


def test_map_batch_100000_no_backoff(al, fasta_iter):
    al.enable_threading(4)
    iter_ = repeat(next(fasta_iter), 100000)
    with pytest.raises(RuntimeError) as excinfo:
        mappings = al.map_batch(iter_, back_off=False)
        n = 0
        for res in mappings:
            n += 1
    assert "Internal error adding data to work queue, without backoff" in str(
        excinfo
    )
    assert (
        "Is your fastq batch larger than 50000? Perhaps try"
        " `map_batch` with back_off=True?" in str(excinfo)
    )


@pytest.mark.parametrize(
    "fasta",
    ["fasta_iter", "fasta_list", "fasta_tuple", "fasta_generator"],
    indirect=True,
)
def test_map_batch(al, fasta):
    al.enable_threading(2)
    mappings = al.map_batch(fasta)
    n = 0
    for res in mappings:
        n += 1
    assert n == 40


def test_map_batch_fail_dict_single(al, fasta_iter):
    fasta = next(fasta_iter)
    al.enable_threading(2)
    with pytest.raises(TypeError) as excinfo:
        _ = al.map_batch(fasta)
    e = str(excinfo)
    assert "Unsupported batch type, pass a list, iter, generator or tuple" in e


def test_map_batch_fail_dict_many(al, fasta_iter):
    fasta = {i: dct for i, dct in enumerate(fasta_iter)}
    al.enable_threading(2)
    with pytest.raises(TypeError) as excinfo:
        _ = al.map_batch(fasta)
    e = str(excinfo)
    assert "Unsupported batch type, pass a list, iter, generator or tuple" in e


def test_map_batch_fail_list_str(al, fasta_iter):
    fasta = [dct["seq"] for dct in fasta_iter]
    al.enable_threading(2)
    with pytest.raises(TypeError) as excinfo:
        _ = al.map_batch(fasta)
    assert "Element in iterable is not a dictionary" in str(excinfo.value)


def test_map_batch_fail_no_seq_key(al, fasta_iter):
    fasta = [{"SEQ": dct["seq"]} for dct in fasta_iter]
    al.enable_threading(2)
    with pytest.raises(KeyError) as excinfo:
        _ = al.map_batch(fasta)
    assert "AHHH Key 🗝️  not found in iterated dictionary" in str(excinfo)


def test_map_batch_fail_seq_not_str(al, fasta_iter):
    fasta = [{"seq": dct["seq"].encode()} for dct in fasta_iter]
    al.enable_threading(2)
    with pytest.raises(ValueError) as excinfo:
        _ = al.map_batch(fasta)
    assert "`seq` must be a string" in str(excinfo)


def test_map_batch_fail_exhausted_iter(al, fasta_iter):
    _ = list(fasta_iter)
    al.enable_threading(2)
    mappings = al.map_batch(fasta_iter)
    assert len(list(mappings)) == 0
