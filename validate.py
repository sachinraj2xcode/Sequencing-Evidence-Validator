#!/usr/bin/env python3
"""Sequencing Evidence Validator: checks FASTQ, BAM and VCF evidence."""

import argparse
import gzip
import os
import subprocess
import sys
from collections import Counter, namedtuple

# ── thresholds ────────────────────────────────────────────────────────────

MIN_MEAN_Q       = 20.0
MAX_N_FRAC       = 0.05
MIN_READ_LEN     = 50
MAX_SHORT_PCT    = 5.0
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

HET_GTS = ("0/1", "1/0")
HOM_GTS = ("1/1",)

COLOURS = {"PASS": "\033[32m", "WARN": "\033[33m", "ERROR": "\033[31m"}
SYMBOLS = {"PASS": "✓", "WARN": "⚠", "ERROR": "✗"}
RESET   = "\033[0m"
RULE    = "─" * 60

Variant = namedtuple("Variant", "chrom pos ref alt gt verdict checks")


# ── helpers ───────────────────────────────────────────────────────────────

def run(cmd):
    """Run a command and return stdout, or '' if it fails."""
    try:
        return subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def worst(verdicts):
    verdicts = set(verdicts)
    if "ERROR" in verdicts:
        return "ERROR"
    return "WARN" if "WARN" in verdicts else "PASS"


def verdict_of(checks):
    return worst(c[1] for c in checks)


def paint(verdict, text, colour):
    return f"{COLOURS.get(verdict, '')}{text}{RESET}" if colour else text


def abort_if_missing(*paths):
    for path in paths:
        if not os.path.exists(path):
            sys.exit(f"file not found: {path}")


# ── readers ───────────────────────────────────────────────────────────────

def read_fastq(path):
    """Yield (sequence, quality) for each FASTQ record."""
    with open_text(path) as fh:
        while fh.readline():
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            yield seq, qual


def read_vcf(path):
    """Yield (chrom, pos, ref, alt, gt) for each variant line."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 10:
                continue
            fmt = f[8].split(":")
            sample = f[9].split(":")
            gt = dict(zip(fmt, sample)).get("GT", "./.").replace("|", "/")
            yield f[0], int(f[1]), f[3], f[4].split(",")[0], gt


def pileup_pysam(bam, chrom, pos):
    import pysam
    reads = []
    with pysam.AlignmentFile(bam, "rb") as bf:
        for col in bf.pileup(chrom, pos - 1, pos, truncate=True,
                             min_base_quality=0, min_mapping_quality=0):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                r = pr.alignment
                bq = r.query_qualities[pr.query_position] if r.query_qualities else 0
                base = r.query_sequence[pr.query_position].upper()
                reads.append((base, bq, r.mapping_quality or 0, not r.is_reverse))
    return reads


def pileup_samtools(bam, chrom, pos):
    out = run(["samtools", "mpileup", "-r", f"{chrom}:{pos}-{pos}",
               "-Q", str(MIN_BQ), "-q", str(MIN_MAPQ), bam])
    reads = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        for b, q in zip(parts[4].upper(), parts[5]):
            if b in ".,":
                base = "REF"
            elif b in "ACGT":
                base = b
            else:
                continue
            reads.append((base, ord(q) - 33, 30, b == b.upper()))
    return reads


def get_pileup(bam, chrom, pos):
    """Return (base, base_quality, mapq, is_forward) for reads at a position."""
    try:
        return pileup_pysam(bam, chrom, pos)
    except ImportError:
        return pileup_samtools(bam, chrom, pos)


# ── validation layers ─────────────────────────────────────────────────────

def layer1(fastq):
    n = q_sum = q_total = len_sum = short = 0
    gc_sum = n_sum = 0.0

    for seq, qual in read_fastq(fastq):
        slen = max(1, len(seq))
        n += 1
        q_sum += sum(ord(c) - 33 for c in qual)
        q_total += len(qual)
        len_sum += len(seq)
        short += len(seq) < MIN_READ_LEN
        gc_sum += (seq.count("G") + seq.count("C")) / slen
        n_sum += seq.count("N") / slen

    if n == 0 or q_total == 0:
        return [("FASTQ_EMPTY", "ERROR", "no reads found")]

    mean_q = q_sum / q_total
    short_pct = short / n * 100
    mean_gc = gc_sum / n
    mean_n = n_sum / n

    return [
        ("MEAN_BASE_QUALITY", "PASS" if mean_q >= MIN_MEAN_Q else "ERROR",
         f"{mean_q:.1f} (min {MIN_MEAN_Q})"),
        ("READ_LENGTH", "WARN" if short_pct > MAX_SHORT_PCT else "PASS",
         f"mean {len_sum / n:.0f} bp, {short_pct:.1f}% < {MIN_READ_LEN} bp"),
        ("GC_CONTENT", "PASS" if MIN_GC <= mean_gc <= MAX_GC else "WARN",
         f"{mean_gc * 100:.1f}% (expected {MIN_GC * 100:.0f}–{MAX_GC * 100:.0f}%)"),
        ("N_FRACTION", "WARN" if mean_n > MAX_N_FRAC else "PASS",
         f"{mean_n * 100:.2f}% (max {MAX_N_FRAC * 100:.0f}%)"),
    ]


def parse_flagstat(text):
    total = mapped = dups = 0
    for line in text.splitlines():
        count = int(line.split()[0])
        if "in total" in line:
            total = count
        elif "mapped (" in line and "primary" not in line:
            mapped = count
        elif "duplicates" in line and "primary" not in line:
            dups = count
    return total, mapped, dups


def mean_mapq(bam):
    values = []
    for line in run(["samtools", "view", "-s", "0.01", bam]).splitlines():
        parts = line.split("\t")
        if len(parts) > 4 and parts[4].isdigit():
            values.append(int(parts[4]))
    return sum(values) / len(values) if values else None


def layer2(bam):
    rc = subprocess.call(["samtools", "quickcheck", bam],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    results = [("BAM_INTEGRITY", "PASS" if rc == 0 else "ERROR",
                "ok" if rc == 0 else "failed quickcheck")]

    total, mapped, dups = parse_flagstat(run(["samtools", "flagstat", bam]))
    if total == 0:
        return results + [("BAM_READS", "ERROR", "no reads")]

    mapping_rate = mapped / total
    dup_rate = dups / total
    results += [
        ("MAPPING_RATE", "PASS" if mapping_rate >= MIN_MAPPING_RATE else "ERROR",
         f"{mapping_rate * 100:.2f}% (min {MIN_MAPPING_RATE * 100:.0f}%)"),
        ("DUPLICATE_RATE", "PASS" if dup_rate <= MAX_DUP_RATE else "WARN",
         f"{dup_rate * 100:.2f}% (max {MAX_DUP_RATE * 100:.0f}%)"),
    ]

    mq = mean_mapq(bam)
    if mq is not None:
        results.append(("MEAN_MAPQ", "PASS" if mq >= MIN_MEAN_MAPQ else "WARN",
                        f"{mq:.1f} (min {MIN_MEAN_MAPQ:.0f})"))
    return results


def check_variant(ref, alt, gt, reads):
    hq = [r for r in reads if r[1] >= MIN_BQ and r[2] >= MIN_MAPQ]
    alts = [r for r in hq if r[0] == alt[0]]
    refs = [r for r in hq if r[0] == ref[0]]
    depth, n_alt = len(hq), len(alts)
    vaf = n_alt / depth if depth else 0.0

    if depth < MIN_DEPTH:
        return [("DEPTH", "ERROR", f"{depth} reads (min {MIN_DEPTH})")]
    checks = [("DEPTH", "PASS", str(depth))]

    if n_alt < MIN_ALT_COUNT:
        return checks + [("ALT_COUNT", "ERROR",
                          f"{n_alt} alt reads (min {MIN_ALT_COUNT})")]
    checks.append(("ALT_COUNT", "PASS", str(n_alt)))

    if vaf < MIN_VAF:
        return checks + [("VAF", "ERROR", f"{vaf:.3f} (min {MIN_VAF})")]
    checks.append(("VAF", "PASS", f"{vaf:.3f}"))

    if gt in HET_GTS and not refs:
        checks.append(("REF_SUPPORT", "WARN", "no ref reads for het call"))

    fwd = sum(1 for r in alts if r[3])
    rev = n_alt - fwd
    checks.append(("STRAND_BALANCE", "PASS" if fwd and rev else "WARN",
                   f"{fwd} fwd / {rev} rev"))

    if gt in HOM_GTS:
        checks.append(("GENOTYPE", "PASS" if vaf >= HOM_VAF_MIN else "WARN",
                       f"hom-alt vaf {vaf:.3f} (min {HOM_VAF_MIN})"))
    elif gt in HET_GTS:
        ok = HET_VAF_MIN <= vaf <= HET_VAF_MAX
        checks.append(("GENOTYPE", "PASS" if ok else "WARN",
                       f"het vaf {vaf:.3f} (expected {HET_VAF_MIN}–{HET_VAF_MAX})"))
    return checks


def in_region(chrom, pos, region):
    """region is 'CHR' or 'CHR:START-END'; None means no filter."""
    if not region:
        return True
    region_chrom, _, coords = region.partition(":")
    if chrom != region_chrom:
        return False
    if not coords:
        return True
    start, end = map(int, coords.split("-"))
    return start <= pos <= end


def layer3(vcf, bam, region=None):
    variants = []
    for chrom, pos, ref, alt, gt in read_vcf(vcf):
        if not in_region(chrom, pos, region):
            continue
        checks = check_variant(ref, alt, gt, get_pileup(bam, chrom, pos))
        variants.append(Variant(chrom, pos, ref, alt, gt, verdict_of(checks), checks))
    return variants


# ── terminal output ───────────────────────────────────────────────────────

def print_layer(title, results, colour):
    print(f"\n  {title}\n  {RULE}")
    for name, verdict, msg in results:
        tag = paint(verdict, f"{SYMBOLS[verdict]} {verdict:<7}", colour)
        print(f"  [{tag}]  {name:<24}  {msg}")


def print_variants(variants, verbose, colour):
    if verbose:
        print(f"\n  {'CHROM':<8} {'POS':>10}  {'REF':<5} {'ALT':<5} {'GT':<6}  VERDICT")
        print(f"  {RULE}")
        for v in variants:
            tag = paint(v.verdict, f"{SYMBOLS[v.verdict]} {v.verdict}", colour)
            print(f"  {v.chrom:<8} {v.pos:>10}  {v.ref:<5} {v.alt:<5} {v.gt:<6}  {tag}")
        return

    flagged = [v for v in variants if v.verdict != "PASS"]
    if flagged:
        print("\n  variants with warnings / errors:")
    for v in flagged:
        tag = paint(v.verdict, f"{SYMBOLS[v.verdict]} {v.verdict}", colour)
        note = v.checks[-1][0] if v.checks else ""
        print(f"  {tag:<20}  {v.chrom}:{v.pos}  {v.ref}>{v.alt}  {v.gt}  [{note}]")


def print_summary(l1, l2, variants, colour):
    layer_verdicts = []
    if l1:
        layer_verdicts.append(("layer 1", verdict_of(l1), ""))
    if l2:
        layer_verdicts.append(("layer 2", verdict_of(l2), ""))
    if variants:
        counts = Counter(v.verdict for v in variants)
        detail = (f"  ({counts['PASS']}/{len(variants)} pass, "
                  f"{counts['WARN']} warn, {counts['ERROR']} error)")
        layer_verdicts.append(("layer 3", worst(counts), detail))

    print(f"\n  {RULE}")
    for name, verdict, detail in layer_verdicts:
        print(f"  {name}  {paint(verdict, verdict, colour)}{detail}")
    overall = worst(v for _, v, _ in layer_verdicts)
    print(f"\n  overall  {paint(overall, overall, colour)}\n")


# ── graph ─────────────────────────────────────────────────────────────────

def make_graph(variants, out_path):
    if not variants:
        print("no variants — skipping graph")
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib/numpy not installed — skipping graph")
        return

    labels = ["PASS", "WARN", "ERROR"]
    colors = ["#2ecc71", "#f39c12", "#e74c3c"]
    edge, bg, grid = "#495057", "#f8f9fa", "#dee2e6"
    total = len(variants)

    counts = Counter(v.verdict for v in variants)
    snv = Counter(v.verdict for v in variants if len(v.ref) == len(v.alt))
    indel = Counter(v.verdict for v in variants if len(v.ref) != len(v.alt))

    fail_reasons = Counter(v.checks[-1][0] for v in variants
                           if v.verdict != "PASS" and v.checks)
    vafs = {label: [] for label in labels}
    for v in variants:
        for name, _, msg in v.checks:
            if name == "VAF":
                vafs[v.verdict].append(float(msg.split()[0]))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor=bg)
    fig.suptitle("Sequencing Evidence Validator — Layer 3 Results",
                 fontsize=12, fontweight="bold")

    def style(ax, title, ygrid=True):
        ax.set_facecolor(bg)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
        ax.set_axisbelow(True)
        ax.yaxis.grid(ygrid, color=grid, linewidth=0.6, zorder=0)

    # overall verdict
    ax = axes[0, 0]
    style(ax, "Overall verdict")
    vals = [counts[l] for l in labels]
    for bar, cnt in zip(ax.bar(labels, vals, color=colors, edgecolor=edge,
                               linewidth=0.7, zorder=3), vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{cnt}\n{cnt / total * 100:.1f}%", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, total * 1.18)

    # snv vs indel
    ax = axes[0, 1]
    style(ax, "SNV vs Indel")
    x, w = np.arange(3), 0.33
    ax.bar(x - w / 2, [snv[l] for l in labels], w, color=colors,
           edgecolor=edge, linewidth=0.7, zorder=3)
    ax.bar(x + w / 2, [indel[l] for l in labels], w, color=colors,
           edgecolor=edge, linewidth=0.7, hatch="///", alpha=0.8, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(handles=[
        mpatches.Patch(facecolor="white", edgecolor=edge, label="SNV"),
        mpatches.Patch(facecolor="white", edgecolor=edge, hatch="///", label="Indel"),
    ], fontsize=8, frameon=False)

    # fail reasons
    ax = axes[1, 0]
    style(ax, "Fail reasons", ygrid=False)
    ax.xaxis.grid(True, color=grid, linewidth=0.6, zorder=0)
    if fail_reasons:
        items = sorted(fail_reasons.items(), key=lambda kv: kv[1])
        names, nums = [k for k, _ in items], [n for _, n in items]
        for bar, cnt in zip(ax.barh(names, nums, color="#f39c12", edgecolor=edge,
                                    linewidth=0.7, zorder=3), nums):
            ax.text(cnt + 0.1, bar.get_y() + bar.get_height() / 2,
                    str(cnt), va="center", fontsize=9)
        ax.set_xlim(0, max(nums) * 1.3)
    else:
        ax.text(0.5, 0.5, "no failures", ha="center", va="center",
                transform=ax.transAxes)

    # vaf distribution
    ax = axes[1, 1]
    style(ax, "VAF distribution", ygrid=False)
    bins = np.linspace(0, 1, 21)
    for label, color in zip(labels, colors):
        if vafs[label]:
            ax.hist(vafs[label], bins=bins, color=color, alpha=0.7,
                    edgecolor=edge, linewidth=0.4,
                    label=f"{label} (n={len(vafs[label])})", zorder=3)
    ax.axvline(MIN_VAF, ls="--", lw=1, color="#6c757d", label=f"min VAF {MIN_VAF}")
    ax.axvline(HOM_VAF_MIN, ls=":", lw=1, color="#343a40", label=f"hom-alt {HOM_VAF_MIN}")
    ax.set_xlabel("VAF")
    ax.legend(fontsize=7.5, frameon=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=bg)
    plt.close(fig)
    print(f"graph saved → {out_path}")


# ── commands ──────────────────────────────────────────────────────────────

def report_fastq(path, colour):
    results = layer1(path)
    print_layer("fastq", results, colour)
    return results


def report_bam(path, colour):
    results = layer2(path)
    print_layer("bam", results, colour)
    return results


def report_vcf(args, colour):
    variants = layer3(args.vcf, args.bam, region=args.region)
    passed = sum(v.verdict == "PASS" for v in variants)
    summary = [("VARIANTS", worst(v.verdict for v in variants),
                f"{passed}/{len(variants)} pass")]
    print_layer("vcf", summary, colour)
    print_variants(variants, args.verbose, colour)
    return variants


def cmd_fastq(args):
    colour = not args.no_color
    abort_if_missing(args.fastq)
    print_summary(report_fastq(args.fastq, colour), [], [], colour)


def cmd_bam(args):
    colour = not args.no_color
    abort_if_missing(args.bam)
    print_summary([], report_bam(args.bam, colour), [], colour)


def cmd_vcf(args):
    colour = not args.no_color
    abort_if_missing(args.vcf, args.bam)
    print_summary([], [], report_vcf(args, colour), colour)


def cmd_all(args):
    colour = not args.no_color
    abort_if_missing(args.fastq, args.bam, args.vcf)
    l1 = report_fastq(args.fastq, colour)
    l2 = report_bam(args.bam, colour)
    variants = report_vcf(args, colour)
    print_summary(l1, l2, variants, colour)
    make_graph(variants, args.plot_out)


# ── entry point ───────────────────────────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        description="Sequencing Evidence Validator",
        epilog=("examples:\n"
                "  validate.py fastq reads.fastq.gz\n"
                "  validate.py bam   aln.bam\n"
                "  validate.py vcf   calls.vcf.gz aln.bam\n"
                "  validate.py all   reads.fastq.gz aln.bam calls.vcf.gz"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--no-color", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_vcf_options(sp):
        sp.add_argument("--region", metavar="CHR:START-END")
        sp.add_argument("--verbose", action="store_true")

    sp = sub.add_parser("fastq", help="validate a FASTQ file")
    sp.add_argument("fastq", metavar="FASTQ")
    sp.set_defaults(func=cmd_fastq)

    sp = sub.add_parser("bam", help="validate a BAM file")
    sp.add_argument("bam", metavar="BAM")
    sp.set_defaults(func=cmd_bam)

    sp = sub.add_parser("vcf", help="validate VCF against BAM read evidence")
    sp.add_argument("vcf", metavar="VCF")
    sp.add_argument("bam", metavar="BAM")
    add_vcf_options(sp)
    sp.set_defaults(func=cmd_vcf)

    sp = sub.add_parser("all", help="validate all three files and save a plot")
    sp.add_argument("fastq", metavar="FASTQ")
    sp.add_argument("bam", metavar="BAM")
    sp.add_argument("vcf", metavar="VCF")
    add_vcf_options(sp)
    sp.add_argument("--plot-out", default="results_graph.png", metavar="FILE")
    sp.set_defaults(func=cmd_all)

    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
