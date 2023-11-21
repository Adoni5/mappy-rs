# Test the memory issue, not dropping Aligner on rust
# side after del without manual gc.collect()
from mappy_rs import Aligner
import gc

al = Aligner("/home/adoni5/Documents/Bioinformatics/refs/hg38_no_alts_22.mmi")
al.enable_threading(8)
del al
gc.collect()
