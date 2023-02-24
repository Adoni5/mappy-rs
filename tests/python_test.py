"""python_test.py

This set of tests roughly recreate those that are in lib.rs but use the compiled 
package.
"""
from pathlib import Path
import copy

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
def fasta_file():
    return str(FA_FILE)


@pytest.fixture
def fasta_iter(fasta_file):
    with open(fasta_file, "rt") as fh:
        seqs = [s for _, s in read_fasta(fh)]

    seqs = [
        {"id": i, "seq": seq}
        for i, seq in enumerate(copy.copy(s) for _ in range(10) for s in seqs)
    ]
    return iter(seqs)


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
        "AGAGTGAAGCCAATATTCCGATAACGATTGCTTTCATGATATCCCTCATTCTGGCATTATTTTTTTATACTATACTATTC"
        "GATATCGCACAGATCAATGGAGTCGTGAGAAAATAAACATGTTTTGCGAACCGCTATGTGTGGAAGACAAAAAATGGAGG"
        "TGAAATTGATGGAAGCAAAGACACAGGCGTACTTTTTTCAGGATGATGGCAGGATTCCGAATCACCCTGATTTTCCGCTC"
        "GTTGTGTATCAAAACGCACTCAAGGACACCGGTCAGGCAGAGCGGATCGTCAACCGGCATGGCTGGTCAAACAGCTGGTC"
        "GGGGAGTGTTTTTCCATACCATCATTATCACAGCAATACGCATGAAGTCCTGATTGCAGTTCGGGGAGAGGCTGTGATTC"
    )
    seq = al.seq(contig)
    assert seq == expected


def test_map_one(al):
    mappings = al.map(
        "AGAGCAGGTAGGATCGTTGAAAAAAGAGTACTCAGGATTCCATTCAACTTTTACTGATTTGAAGCGTACTGTTTATGGCC"
        "AAGAATATTTACGTCTTTACAACCAATACGCAAAAAAAGGTTCATTGAGTTTGGTTGTGATTTGATGAAAATTACTGAGA"
        "ATAACAGGATTATTAAGCTGATTGATGAACTAAATCAGCTTAATAAATATTCTTTGCAGATAGGAATATTTGGGGAAAAT"
        "GATTCTTTTATGGCGATGTTGGCCCAAGTTCATGAATTTGGGGTGACTATTCGTCCCAAAGGTCGTTTTCTTGTTATACC"
        "ACTTATGAAAAAGTATAGAGGTAAAAGTCCACGTCAATTTGATTTGTTTTTTATGCAAACTAAAGAAAATCACAAGTTTT",
        cs=True,
    )
    assert len(mappings) == 1
    mapping = mappings[0]
    assert mapping.target_start == 0
    assert mapping.target_end == 400


def test_map_batch(al, fasta_iter):
    al.enable_threading(2)
    mappings = al.map_batch(fasta_iter)
    n = 0
    for res in mappings:
        n += 1
    assert n == 40
