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
