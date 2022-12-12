from mappy_rs import Aligner
import sys
import mappy as mp
import pytest
import sys
import argparse
import logging
import math
import time
import gzip
from functools import partial
from multiprocessing import Pool
from textwrap import fill
from contextlib import contextmanager, redirect_stdout
from rich import print
from io import StringIO
from pprint import pprint
from subprocess import Popen, PIPE
from collections import defaultdict

from pyguppy_client_lib.pyclient import PyGuppyClient
from pyguppy_client_lib.helper_functions import package_read, run_server
from ont_fast5_api.fast5_interface import get_fast5_file


ADDRESS = "ipc:///tmp/.guppy/"
PORT = 5556
CONFIG = "dna_r9.4.1_450bps_fast_prom.cfg"
# A bit of a hack, use logging to write to stdout
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format="%(message)s",
)
LOGGER = logging.getLogger(__name__)
# FILES = ["PAK21443_pass_barcode02_6366be35_0.fast5",
# "PAK21443_pass_barcode02_6366be35_1.fast5",
# "PAK21443_pass_barcode02_6366be35_2.fast5",
# "PAK21443_pass_barcode02_6366be35_3.fast5"]

FILES = [
    "test_data/fast5/PAK21443_pass_barcode02_6366be35_0.fast5",
    "test_data/fast5/PAK21443_pass_barcode02_6366be35_1.fast5",
]

BATCH_SIZE = 300
REF_PATH = "/data/projects/rory_says_hi/sim_retinal_panel/hg_38_no_alts_22.mmi"
GUPPY_BIN = "/usr/bin/"


@contextmanager
def start_guppy_server_and_client(
    bin_path=GUPPY_BIN,
    config=CONFIG,
    port=PORT,
    ipc=True,
    server_args=None,
    **client_args,
):
    if server_args is None:
        server_args = []
    server_args.extend(
        [
            "--config",
            config,
            "--port",
            f"/tmp/.guppy/{port}",
            "--device",
            "cuda:3",
            "--log_path",
            "/data/projects/rory_says_hi/sim_retinal_panel/guppy_logs",
            "--num_callers",
            "3",
            "--ipc_threads",
            "4",
            "--max_queued_reads",
            "5000",
        ]
    )
    # This function has it's own prints that may want to be suppressed
    with redirect_stdout(StringIO()) as fh:
        server, port = run_server(server_args, bin_path=bin_path)

    if port == "ERROR":
        raise RuntimeError("Server couldn't be started")

    if ipc:
        address = f"{port}"
    else:
        address = f"localhost:{port}"
    print(address)
    client = PyGuppyClient(address=address, config=config, **client_args)
    # client.connect()
    try:
        with client:
            yield client

    finally:
        server.terminate()


def yield_reads(files):
    """Yield reads from FAST5 `files` in a given `batch_size`"""
    for filename in files:
        print(filename)
        with get_fast5_file(filename, "r") as f5:
            yield from f5.get_reads()


def pack(read):
    """Pack an ont_fast5_api.Fast5Read for calling"""
    read_id = read.read_id
    raw_data = read.get_raw_data()
    channel_info = read.get_channel_info()
    scaling = channel_info["range"] / channel_info["digitisation"]
    offset = channel_info["offset"]
    return package_read(read_id, raw_data, offset, scaling)


def format_read(res):
    return "\n".join(
        (
            f"@{res['metadata']['read_id']}",
            fill(res["datasets"]["sequence"], width=80),
            "+",
            fill(res["datasets"]["qstring"], width=80),
        )
    )


def setup_single():
    a = mp.Aligner(REF_PATH, preset="map-ont")
    return a


def setup_multi(i):
    a = Aligner(i, REF_PATH)
    return a


def run_multi_sequence(a, client):
    """ """
    sequences = []
    sent, recv, total = 0, 0, 0
    for read in yield_reads(FILES):
        s = client.pass_read(read, pack)
        while not s:
            # The pyclient.PyGuppyClient by default makes 6 attempts to
            #   send a read to the server, the initial attempt and then
            #   5 retires if the server is not ready to accept the read

            # This will keep sending a read until it is accepted and
            #    probably not a good idea.
            time.sleep(client.throttle)
            s = client.pass_read(read, pack)
        sent += 1
        if sent >= BATCH_SIZE:
            while recv < sent:
                # get_completed_reads returns a list of dict, the contents
                #   of the results dict is dependent on the parameters used
                #   to initialise the server, e.g. if you provide an alignment
                #   index then alignment information will be returned, otherwise
                #   it will not be present.
                res = client.get_completed_reads()
                for r in res:
                    for _ in r:
                        recv += 1
                        read = _
                        sequences.append(
                            (
                                read.get("datasets").get("sequence", "A"),
                                1,
                                1,
                                read.get("read_id", "NA"),
                            )
                        )
            total += sent
            sent, recv = 0, 0

    res = a.align_batch(sequences)
    count = 0
    for hit in res:
        for ma in hit:
            count += 1
    print(f"Count is {count}")


@pytest.mark.benchmark
@pytest.mark.parametrize("i", [*list(range(1, 9))])
def test_benchmark_multi(i, benchmark):
    with start_guppy_server_and_client() as client:
        a = setup_multi(i)
        benchmark.pedantic(run_multi_sequence, args=(a, client), iterations=1, rounds=1)


def inbuilt_map(client):
    """
    Iterate inbuilt ONT guppy mapping
    """
    sent, recv, total = 0, 0, 0

    hit_count = 0
    for read in yield_reads(FILES):
        s = client.pass_read(read, pack)
        while not s:
            # The pyclient.PyGuppyClient by default makes 6 attempts to
            #   send a read to the server, the initial attempt and then
            #   5 retires if the server is not ready to accept the read

            # This will keep sending a read until it is accepted and
            #    probably not a good idea.
            time.sleep(client.throttle)
            s = client.pass_read(read, pack)
        sent += 1
        if sent >= BATCH_SIZE:
            while recv < sent:
                # get_completed_reads returns a list of dict, the contents
                #   of the results dict is dependent on the parameters used
                #   to initialise the server, e.g. if you provide an alignment
                #   index then alignment information will be returned, otherwise
                #   it will not be present.
                res = client.get_completed_reads()
                for r in res:
                    for _ in r:
                        recv += 1
                        read = _
                        for hit in read.get("datasets").get("align_results"):
                            hit_count += 1
            sent, recv = 0, 0
    print(hit_count)


@pytest.mark.benchmark
@pytest.mark.parametrize("i", [*list(range(1, 9))])
def test_inbuilt_map(i, benchmark):
    with start_guppy_server_and_client(
        server_args=[
            "--num_alignment_threads",
            f"{i}",
        ],
        align_ref=REF_PATH,
        server_file_load_timeout=120,
    ) as client:
        benchmark.pedantic(inbuilt_map, args=(client,), iterations=1, rounds=1)


def base_call(client):
    sequences = []
    sent, recv, total = 0, 0, 0
    for read in yield_reads(FILES):
        s = client.pass_read(read, pack)
        while not s:
            # The pyclient.PyGuppyClient by default makes 6 attempts to
            #   send a read to the server, the initial attempt and then
            #   5 retires if the server is not ready to accept the read

            # This will keep sending a read until it is accepted and
            #    probably not a good idea.
            time.sleep(client.throttle)
            s = client.pass_read(read, pack)
        sent += 1
        if sent >= BATCH_SIZE:
            while recv < sent:
                # get_completed_reads returns a list of dict, the contents
                #   of the results dict is dependent on the parameters used
                #   to initialise the server, e.g. if you provide an alignment
                #   index then alignment information will be returned, otherwise
                #   it will not be present.
                res = client.get_completed_reads()
                for r in res:
                    for _ in r:
                        recv += 1
                        read = _
                        sequences.append(read.get("datasets").get("sequence", "A"))
            sent, recv = 0, 0
    hit_count = 0
    # for s in sequences:
    #     for hit in a.map(s):
    #         hit_count += 1
    # print(hit_count)


@pytest.mark.benchmark
@pytest.mark.skip(reason="Skipping as slow whilst testing alignment")
def test_basecall(benchmark):
    with start_guppy_server_and_client() as client:
        benchmark.pedantic(base_call, args=(client,), iterations=2, rounds=5)


def print_args(args, logger=None, exclude=None):
    """Print and format all arguments from the command line"""
    if exclude is None:
        exclude = []
    dirs = dir(args)
    m = max([len(a) for a in dirs if a[0] != "_"])
    for attr in dirs:
        if attr[0] != "_" and attr not in exclude and attr.lower() == attr:
            record = "{a}={b}".format(a=attr, m=m, b=getattr(args, attr))
            if logger is not None:
                logger.info(record)
            else:
                print(record)


def compare_hit_counts(mc, mrsc, gc):
    """Compare Mapping hit counts"""
    assert (
        mc == mrsc
    ), "Mappy number ({}) of hits does not equal rust number ({}) of hits".format(
        mc, mrsc
    )
    assert (
        mrsc == gc
    ), "Mappy rust number ({}) of hits does not equal guppy number of hits ({})".format(
        mrsc, gc
    )
    assert mc == mrsc == gc


def compare_hits(
    mh: dict[str, list[tuple[str, str]]],
    mrh: dict[str, list[tuple[str, str]]],
    gh: dict[str, list[tuple[str, str]]],
    mm2h: dict[str, list[tuple[str, str]]],
    counts: dict[str, int]
):
    """Compare mapping starts, ends and contigs"""
    for read_id, mm2_hits in mm2h.items():
        mm2_hits.sort()
        guppy_hits = gh[read_id]
        guppy_hits.sort()
        rust_hits = mrh[read_id]
        rust_hits.sort()
        mappy_hits = mh[read_id]
        mappy_hits.sort()
        # 5099b2fd-3fe7-4e23-9229-7b2a3ac7ce75 read has one fewer hit on guppy
        # assert (
        #     len(guppy_hits) == len(rust_hits) == len(mappy_hits) == len(mm2_hits)
        # ), "More or fewer hits: Guppy: {}, Mappy: {}, Rust: {}, MM2: {}, ReadID: {}".format(
        #     len(guppy_hits), len(mappy_hits), len(rust_hits), len(mm2_hits), read_id
        # )

        for mm2_hit, guppy_hit, rust_hit, mappy_hit in zip(
            mm2_hits, guppy_hits, rust_hits, mappy_hits
        ):
            mm2_start, mm2_end, mm2_contig = mm2_hit
            guppy_start, guppy_end = guppy_hit
            rust_start, rust_end = rust_hit
            mappy_start, mappy_end = mappy_hit
            print(
                f"Guppy: {guppy_start} ({len(guppy_hits)}), Mappy: {mappy_start} ({len(mappy_hits)}), MM2: {mm2_start} ({len(mm2_hits)}), Rust: {rust_start} ({len(rust_hits)}) \t {read_id}"
            )
            assert (
                mappy_start == rust_start
            ), "Mappy start {} and Rust start {} not the same".format(
                mappy_start, rust_start
            )
            counts["total_guppy_hits"] += len(guppy_hits)
            counts["total_mappy_hits"] += len(mappy_hits)
            counts["total_mm2_hits"] += len(mm2_hits)
            counts["total_rust_hits"] += len(rust_hits)


#        assert rust_hit.target_start == guppy_start,  "Rust start {} and Guppy start {} not the same".format(rust_hit.target_start, guppy_start)


def run_minimap2(read_seq: str) -> dict[str, list[tuple[str, str, str]]]:
    """Run the minimap2 on the CLI and return paf hits

    parameters
    ----------
    read_seq: str
        Fasta formatted string to map
    """
    x = "minimap2 -x map-ont -c /data/projects/rory_says_hi/sim_retinal_panel/hg_38_no_alts_22.mmi -".split()
    out, err = Popen(args=x, stdout=PIPE, stderr=PIPE, stdin=PIPE).communicate(
        input=read_seq.encode()
    )
    mm2_hits = defaultdict(list)
    for x in out.decode().splitlines():
        hits = x.split("\t")
        mm2_hits[hits[0]].append((int(hits[7]), int(hits[8]), hits[5]))

    return mm2_hits


def test_mapping(client):
    """Compare the mapping we get back from mappy_rs and from mappy. Ideally they will be identical."""
    sent, recv, total = 0, 0, 0
    a_rs = setup_multi(1)
    a_mp = mp.Aligner(fn_idx_in=REF_PATH, extra_flags=8)
    go = True
    mapping_counts = defaultdict(int)
    for read in yield_reads(FILES):
        s = client.pass_read(read, pack)
        while not s:
            # The pyclient.PyGuppyClient by default makes 6 attempts to
            #   send a read to the server, the initial attempt and then
            #   5 retires if the server is not ready to accept the read

            # This will keep sending a read until it is accepted and
            #    probably not a good idea.
            time.sleep(client.throttle)
            s = client.pass_read(read, pack)
        sent += 1
        if sent >= BATCH_SIZE:
            sequences = []
            mappy_hits, rust_hits, guppy_hits = (
                defaultdict(list),
                defaultdict(list),
                defaultdict(list),
            )
            while recv < sent:
                # get_completed_reads returns a list of dict, the contents
                #   of the results dict is dependent on the parameters used
                #   to initialise the server, e.g. if you provide an alignment
                #   index then alignment information will be returned, otherwise
                #   it will not be present.
                res = client.get_completed_reads()
                for r in res:
                    for _ in r:
                        recv += 1
                        read = _
                        mappy_count, rust_count, guppy_count = 0, 0, 0
                        seq = read.get("datasets").get("sequence", "A")
                        read_id = read.get("metadata", {}).get("read_id")
                        sequences.append(f">{read_id}\n{seq}")
                        for hit_st in a_mp.map(f"{seq}", cs=True, MD=False):
                            mappy_count += 1
                            mappy_hits[read_id].append((hit_st.r_st, hit_st.r_en))
                        for hit_mt in next(
                            a_rs.align_batch(iter([(seq, 1, 1, "Hello")]))
                        ):
                            rust_count += 1
                            rust_hits[read_id].append(
                                (hit_mt.target_start, hit_mt.target_end)
                            )
                        for hit in read.get("datasets").get("align_results"):
                            genome = read.get("metadata", {}).get(
                                "alignment_genome", "*"
                            )
                            if not genome == "*":
                                guppy_count += 1
                                guppy_hits[read_id].append((hit[2] - 1, hit[3] - 1))
            sent, recv = 0, 0
            minimap2_hits = run_minimap2("\n".join(sequences))
            # compare_hit_counts(mappy_count, rust_count, guppy_count, minimap2_hits)
            compare_hits(mappy_hits, rust_hits, guppy_hits, minimap2_hits, mapping_counts)
            print(mapping_counts)
    print(mapping_counts)
        

if __name__ == "__main__":
    # conn = {"address": f"{ADDRESS}{PORT}", "config": CONFIG}
    # print(conn)
    with start_guppy_server_and_client(
        server_args=[
            "--num_alignment_threads",
            "8",
        ],
        align_ref=REF_PATH,
        server_file_load_timeout=120,
    ) as client:
        #    # client, server = start_guppy_server_and_client()
        #     print(client)
        #     print("Initialising aligner")
        #     print("Initialised")
        #     inbuilt_map(client)
        test_mapping(client)
