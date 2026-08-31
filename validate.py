import argparse
import gzip
import os
import subprocess
import sys
from collections import Counter

MIN_MEAN_Q       = 20.0
MAX_N_FRAC       = 0.05
MIN_READ_LEN     = 50
MIN_GC           = 0.20
MAX_GC           = 0.80
MIN_MAPPING_RATE = 0.90
MIN_MEAN_MAPQ    = 30.0
MAX_DUP_RATE     = 0.20
MIN_DEPTH        = 8
MIN_ALT_COUNT    = 2
MIN_VAF          = 0.05
MIN_BQ           = 20
MIN_MAPQ         = 20
HOM_VAF_MIN      = 0.85
HET_VAF_MIN      = 0.25
HET_VAF_MAX      = 0.75


def read_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            if not fh.readline():
                break
            seq  = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            yield seq, qual


def layer1(fastq):
    n = q_sum = q_total = 0

    for seq, qual in read_fastq(fastq):
        n      += 1
        q_sum  += sum(ord(c) - 33 for c in qual)
        q_total += len(qual)

    if n == 0:
        return [("FASTQ_EMPTY", "ERROR", "no reads found")]

    results = []
    mean_q = q_sum / q_total
    results.append(("MEAN_BASE_QUALITY",
                    "PASS" if mean_q >= MIN_MEAN_Q else "ERROR",
                    f"{mean_q:.1f} (min {MIN_MEAN_Q})"))
    return results


def run(cmd):
    try:
        return subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def layer2(bam):
    results = []

    rc = subprocess.call(["samtools", "quickcheck", bam],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    results.append(("BAM_INTEGRITY",
                    "PASS" if rc == 0 else "ERROR",
                    "ok" if rc == 0 else "failed quickcheck"))

    flagstat = run(["samtools", "flagstat", bam])
    total = mapped = dups = 0
    for line in flagstat.splitlines():
        if "in total"   in line: total  = int(line.split()[0])
        if "mapped ("   in line: mapped = int(line.split()[0])
        if "duplicates" in line: dups   = int(line.split()[0])

    if total == 0:
        results.append(("BAM_READS", "ERROR", "no reads"))
        return results

    mr = mapped / total
    results.append(("MAPPING_RATE",
                    "PASS" if mr >= MIN_MAPPING_RATE else "ERROR",
                    f"{mr*100:.2f}% (min {MIN_MAPPING_RATE*100:.0f}%)"))

    dr = dups / total
    results.append(("DUPLICATE_RATE",
                    "PASS" if dr <= MAX_DUP_RATE else "WARN",
                    f"{dr*100:.2f}% (max {MAX_DUP_RATE*100:.0f}%)"))

    mapq_out = run(["samtools", "view", "-s", "0.01", bam])
    mapqs = []
    for line in mapq_out.splitlines():
        parts = line.split("\t")
        if len(parts) > 4:
            try: mapqs.append(int(parts[4]))
            except ValueError: pass
    if mapqs:
        mean_mq = sum(mapqs) / len(mapqs)
        results.append(("MEAN_MAPQ",
                        "PASS" if mean_mq >= MIN_MEAN_MAPQ else "WARN",
                        f"{mean_mq:.1f} (min {MIN_MEAN_MAPQ:.0f})"))

    return results


def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4].split(",")[0]
            fmt    = f[8].split(":")
            sample = f[9].split(":")
            gt     = dict(zip(fmt, sample)).get("GT", "./.").replace("|", "/")
            yield chrom, pos, ref, alt, gt


def get_pileup(bam, chrom, pos):
    try:
        import pysam
        reads = []
        with pysam.AlignmentFile(bam, "rb") as bf:
            for col in bf.pileup(chrom, pos - 1, pos, truncate=True,
                                  min_base_quality=0, min_mapping_quality=0):
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip:
                        continue
                    r    = pr.alignment
                    bq   = r.query_qualities[pr.query_position] if r.query_qualities else 0
                    mq   = r.mapping_quality or 0
                    base = r.query_sequence[pr.query_position].upper()
                    reads.append((base, bq, mq, not r.is_reverse))
        return reads
    except ImportError:
        pass

    out = run(["samtools", "mpileup", "-r", f"{chrom}:{pos}-{pos}",
               "-Q", str(MIN_BQ), "-q", str(MIN_MAPQ), bam])
    reads = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        for b, q in zip(parts[4].upper(), parts[5]):
            if b in ".,":     base = "REF"
            elif b in "ACGT": base = b
            else:             continue
            reads.append((base, ord(q) - 33, 30, b == b.upper()))
    return reads


def main():
    p = argparse.ArgumentParser(description="Sequencing Evidence Validator")
    p.add_argument("--no-color", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("fastq")
    s1.add_argument("fastq", metavar="FASTQ")

    args = p.parse_args()

    if args.cmd == "fastq":
        if not os.path.exists(args.fastq):
            print(f"file not found: {args.fastq}", file=sys.stderr)
            sys.exit(1)
        for name, verdict, msg in layer1(args.fastq):
            print(f"  [{verdict}]  {name:<24}  {msg}")


if __name__ == "__main__":
    main()


# FASTQ read length validation
# FASTQ GC and N fraction validation
# variant depth and allele count checks
# genotype and strand validation
# VCF layer validation
# terminal verdict formatting
# validation result reporting
# variant results graph
# command handlers and argument parsing
# validator CLI execution flow
