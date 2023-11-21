# Test the memory issue
from mappy_rs import Aligner
import gc

al = Aligner("/home/adoni5/Documents/Bioinformatics/refs/hg38_no_alts_22.mmi")
al.enable_threading(8)
del al
gc.collect()
