import argparse
import gzip
import os
import subprocess
import sys
from collections import Counter

MIN_MEAN_Q = 20.0
MAX_N_FRAC = 0.05
MIN_READ_LEN = 50
MIN_GC = 0.20
MAX_GC = 0.80
MIN_MAPPING_RATE = 0.90
MIN_MEAN_MAPQ = 30.0
MAX_DUP_RATE = 0.20
MIN_DEPTH = 8
MIN_ALT_COUNT = 2
MIN_VAF = 0.05
MIN_BQ = 20
MIN_MAPQ = 20
HOM_VAF_MIN = 0.85
HET_VAF_MIN = 0.25
HET_VAF_MAX = 0.75


def read_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            if not fh.readline():
                break
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            yield seq, qual


def layer1(fastq):
    n = q_sum = q_total = len_sum = short = 0
    gc_sum = n_sum = 0.0

    for seq, qual in read_fastq(fastq):
        n += 1
        slen = len(seq)
        q_sum += sum(ord(c) - 33 for c in qual)
        q_total += len(qual)
        len_sum += slen
        if slen < MIN_READ_LEN:
            short += 1
        gc_sum += (seq.count("G") + seq.count("C")) / max(1, slen)
        n_sum += seq.count("N") / max(1, slen)

    if n == 0:
        return [("FASTQ_EMPTY", "ERROR", "no reads found")]

    results = []
    mean_q = q_sum / q_total
    results.append(("MEAN_BASE_QUALITY",
                    "PASS" if mean_q >= MIN_MEAN_Q else "ERROR",
                    f"{mean_q:.1f} (min {MIN_MEAN_Q})"))

    short_pct = short / n * 100
    results.append(("READ_LENGTH",
                    "WARN" if short_pct > 5 else "PASS",
                    f"mean {len_sum/n:.0f} bp, {short_pct:.1f}% < {MIN_READ_LEN} bp"))
