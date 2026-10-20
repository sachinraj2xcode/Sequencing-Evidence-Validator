import argparse
import gzip
import os
import subprocess
import sys
from collections import Counter

# thresholds
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

# colours
GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
RESET  = "\033[0m"


# ── helpers ───────────────────────────────────────────────────────────────

def run(cmd):
    try:
        return subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def colourise(verdict, text, use_colour):
    if not use_colour:
        return text
    return {"PASS": GREEN, "WARN": YELLOW, "ERROR": RED}.get(verdict, "") + text + RESET


def sym(verdict):
    return {"PASS": "✓", "WARN": "⚠", "ERROR": "✗"}.get(verdict, "?")


def verdict_of(checks):
    vs = [c[1] for c in checks]
    return "ERROR" if "ERROR" in vs else "WARN" if "WARN" in vs else "PASS"


def abort_if_missing(*paths):
    for f in paths:
        if not os.path.exists(f):
            print(f"file not found: {f}", file=sys.stderr)
            sys.exit(1)


# ── readers ───────────────────────────────────────────────────────────────

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


def read_vcf(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 10:
                continue
            chrom = f[0]
            pos   = int(f[1])
            ref   = f[3]
            alt   = f[4].split(",")[0]
            fmt   = f[8].split(":")
            samp  = f[9].split(":")
            gt    = dict(zip(fmt, samp)).get("GT", "./.").replace("|", "/")
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
            if b in ".,":
                base = "REF"
            elif b in "ACGT":
                base = b
            else:
                continue
            reads.append((base, ord(q) - 33, 30, b == b.upper()))
    return reads


# ── layers ────────────────────────────────────────────────────────────────

def layer1(fastq):
    n = q_sum = q_total = len_sum = short = 0
    gc_sum = n_sum = 0.0

    for seq, qual in read_fastq(fastq):
        n       += 1
        slen     = len(seq)
        q_sum   += sum(ord(c) - 33 for c in qual)
        q_total += len(qual)
        len_sum += slen
        if slen < MIN_READ_LEN:
            short += 1
        gc_sum += (seq.count("G") + seq.count("C")) / max(1, slen)
        n_sum  += seq.count("N") / max(1, slen)

    if n == 0:
        return [("FASTQ_EMPTY", "ERROR", "no reads found")]

    mean_q    = q_sum / q_total
    short_pct = short / n * 100
    mean_gc   = gc_sum / n
    mean_n    = n_sum / n

    return [
        ("MEAN_BASE_QUALITY", "PASS" if mean_q >= MIN_MEAN_Q else "ERROR",
         f"{mean_q:.1f} (min {MIN_MEAN_Q})"),
        ("READ_LENGTH", "WARN" if short_pct > 5 else "PASS",
         f"mean {len_sum/n:.0f} bp, {short_pct:.1f}% < {MIN_READ_LEN} bp"),
        ("GC_CONTENT", "WARN" if mean_gc < MIN_GC or mean_gc > MAX_GC else "PASS",
         f"{mean_gc*100:.1f}% (expected {MIN_GC*100:.0f}–{MAX_GC*100:.0f}%)"),
        ("N_FRACTION", "WARN" if mean_n > MAX_N_FRAC else "PASS",
         f"{mean_n*100:.2f}% (max {MAX_N_FRAC*100:.0f}%)"),
    ]


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
    dr = dups / total
    results.append(("MAPPING_RATE",
                    "PASS" if mr >= MIN_MAPPING_RATE else "ERROR",
                    f"{mr*100:.2f}% (min {MIN_MAPPING_RATE*100:.0f}%)"))
    results.append(("DUPLICATE_RATE",
                    "PASS" if dr <= MAX_DUP_RATE else "WARN",
                    f"{dr*100:.2f}% (max {MAX_DUP_RATE*100:.0f}%)"))

    mapqs = []
    for line in run(["samtools", "view", "-s", "0.01", bam]).splitlines():
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


def check_variant(chrom, pos, ref, alt, gt, reads):
    hq    = [r for r in reads if r[1] >= MIN_BQ and r[2] >= MIN_MAPQ]
    depth = len(hq)
    alts  = [r for r in hq if r[0] == alt[0]]
    refs  = [r for r in hq if r[0] == ref[0]]
    n_alt = len(alts)
    vaf   = n_alt / depth if depth > 0 else 0.0

    if depth < MIN_DEPTH:
        return [("DEPTH", "ERROR", f"{depth} reads (min {MIN_DEPTH})")]
    if n_alt < MIN_ALT_COUNT:
        return [("DEPTH", "PASS", str(depth)),
                ("ALT_COUNT", "ERROR", f"{n_alt} alt reads (min {MIN_ALT_COUNT})")]
    if vaf < MIN_VAF:
        return [("DEPTH", "PASS", str(depth)), ("ALT_COUNT", "PASS", str(n_alt)),
                ("VAF", "ERROR", f"{vaf:.3f} (min {MIN_VAF})")]

    checks = [
        ("DEPTH",     "PASS", str(depth)),
        ("ALT_COUNT", "PASS", str(n_alt)),
        ("VAF",       "PASS", f"{vaf:.3f}"),
    ]

    if gt in ("0/1", "1/0") and not refs:
        checks.append(("REF_SUPPORT", "WARN", "no ref reads for het call"))

    fwd = sum(1 for r in alts if r[3])
    rev = n_alt - fwd
    checks.append(("STRAND_BALANCE",
                   "PASS" if fwd > 0 and rev > 0 else "WARN",
                   f"{fwd} fwd / {rev} rev"))

    if gt in ("1/1", "1|1"):
        checks.append(("GENOTYPE",
                       "PASS" if vaf >= HOM_VAF_MIN else "WARN",
                       f"hom-alt vaf {vaf:.3f} (min {HOM_VAF_MIN})"))
    elif gt in ("0/1", "0|1", "1/0", "1|0"):
        checks.append(("GENOTYPE",
                       "PASS" if HET_VAF_MIN <= vaf <= HET_VAF_MAX else "WARN",
                       f"het vaf {vaf:.3f} (expected {HET_VAF_MIN}–{HET_VAF_MAX})"))

    return checks


def layer3(vcf, bam, region=None):
    variants = []
    for chrom, pos, ref, alt, gt in read_vcf(vcf):
        if region:
            rc, coords = (region.split(":") + [None])[:2]
            if chrom != rc:
                continue
            if coords:
                s, e = map(int, coords.split("-"))
                if not (s <= pos <= e):
                    continue
        checks = check_variant(chrom, pos, ref, alt, gt, get_pileup(bam, chrom, pos))
        variants.append((chrom, pos, ref, alt, gt, verdict_of(checks), checks))
    return variants


# ── display ───────────────────────────────────────────────────────────────

def print_layer(title, results, use_colour):
    print(f"\n  {title}")
    print(f"  {'─'*60}")
    for name, verdict, msg in results:
        tag = colourise(verdict, f"{sym(verdict)} {verdict:<7}", use_colour)
        print(f"  [{tag}]  {name:<24}  {msg}")


def print_variants(variants, verbose, use_colour):
    non_pass = [v for v in variants if v[5] != "PASS"]
    if verbose:
        print(f"\n  {'CHROM':<8} {'POS':>10}  {'REF':<5} {'ALT':<5} {'GT':<6}  VERDICT")
        print(f"  {'─'*60}")
        for chrom, pos, ref, alt, gt, verdict, _ in variants:
            tag = colourise(verdict, f"{sym(verdict)} {verdict}", use_colour)
            print(f"  {chrom:<8} {pos:>10}  {ref:<5} {alt:<5} {gt:<6}  {tag}")
    elif non_pass:
        print("\n  variants with warnings / errors:")
        for chrom, pos, ref, alt, gt, verdict, checks in non_pass:
            tag  = colourise(verdict, f"{sym(verdict)} {verdict}", use_colour)
            note = checks[-1][0] if checks else ""
            print(f"  {tag:<20}  {chrom}:{pos}  {ref}>{alt}  {gt}  [{note}]")


def print_summary(l1, l2, variants, use_colour):
    total   = len(variants)
    passed  = sum(1 for v in variants if v[5] == "PASS")
    warned  = sum(1 for v in variants if v[5] == "WARN")
    errored = sum(1 for v in variants if v[5] == "ERROR")
    lv3     = "ERROR" if errored else "WARN" if warned else "PASS"
    all_lv  = ([verdict_of(l1)] if l1 else []) + \
              ([verdict_of(l2)] if l2 else []) + \
              ([lv3] if variants else [])
    overall = "ERROR" if "ERROR" in all_lv else "WARN" if "WARN" in all_lv else "PASS"

    print(f"\n  {'─'*60}")
    if l1:       print(f"  layer 1  {colourise(verdict_of(l1), verdict_of(l1), use_colour)}")
    if l2:       print(f"  layer 2  {colourise(verdict_of(l2), verdict_of(l2), use_colour)}")
    if variants: print(f"  layer 3  {colourise(lv3, lv3, use_colour)}  "
                       f"({passed}/{total} pass, {warned} warn, {errored} error)")
    print(f"\n  overall  {colourise(overall, overall, use_colour)}\n")


# ── graph ─────────────────────────────────────────────────────────────────

def make_graph(variants, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
    except ImportError:
        print("matplotlib/numpy not installed — skipping graph")
        return

    verdicts = [v[5] for v in variants]
    counts   = Counter(verdicts)
    total    = len(verdicts)
    colors   = ["#2ecc71", "#f39c12", "#e74c3c"]
    labels   = ["PASS", "WARN", "ERROR"]

    is_indel = [len(v[2]) != len(v[3]) for v in variants]
    snv_c    = Counter(v[5] for v, ind in zip(variants, is_indel) if not ind)
    ind_c    = Counter(v[5] for v, ind in zip(variants, is_indel) if ind)

    fail_reasons = Counter()
    vafs = {"PASS": [], "WARN": [], "ERROR": []}
    for _, _, _, _, _, verdict, checks in variants:
        if verdict != "PASS" and checks:
            fail_reasons[checks[-1][0]] += 1
        for name, _, msg in checks:
            if name == "VAF":
                try: vafs[verdict].append(float(msg.split()[0]))
                except (ValueError, IndexError): pass

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor="#f8f9fa")
    fig.suptitle("Sequencing Evidence Validator — Layer 3 Results",
                 fontsize=12, fontweight="bold")

    def style(ax, title):
        ax.set_facecolor("#f8f9fa")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
        ax.yaxis.grid(True, color="#dee2e6", linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)

    style(axes[0, 0], "Overall verdict")
    vals = [counts.get(l, 0) for l in labels]
    bars = axes[0, 0].bar(labels, vals, color=colors, edgecolor="#495057", linewidth=0.7, zorder=3)
    for bar, cnt in zip(bars, vals):
        axes[0, 0].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                        f"{cnt}\n{cnt/total*100:.1f}%", ha="center", va="bottom", fontsize=9)
    axes[0, 0].set_ylim(0, total * 1.18)

    style(axes[0, 1], "SNV vs Indel")
    x, w = np.arange(3), 0.33
    axes[0, 1].bar(x - w/2, [snv_c.get(l, 0) for l in labels], w,
                   color=colors, edgecolor="#495057", linewidth=0.7, zorder=3)
    axes[0, 1].bar(x + w/2, [ind_c.get(l, 0) for l in labels], w,
                   color=colors, edgecolor="#495057", linewidth=0.7, hatch="///", alpha=0.8, zorder=3)
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels(labels)
    axes[0, 1].legend(handles=[
        mpatches.Patch(facecolor="white", edgecolor="#495057", label="SNV"),
        mpatches.Patch(facecolor="white", edgecolor="#495057", hatch="///", label="Indel"),
    ], fontsize=8, frameon=False)

    style(axes[1, 0], "Fail reasons")
    axes[1, 0].yaxis.grid(False)
    axes[1, 0].xaxis.grid(True, color="#dee2e6", linewidth=0.6, zorder=0)
    if fail_reasons:
        items = sorted(fail_reasons.items(), key=lambda x: x[1])
        hbars = axes[1, 0].barh([i[0] for i in items], [i[1] for i in items],
                                 color="#f39c12", edgecolor="#495057", linewidth=0.7, zorder=3)
        for bar, cnt in zip(hbars, [i[1] for i in items]):
            axes[1, 0].text(cnt + 0.1, bar.get_y() + bar.get_height() / 2,
                            str(cnt), va="center", fontsize=9)
        axes[1, 0].set_xlim(0, max(i[1] for i in items) * 1.3)
    else:
        axes[1, 0].text(0.5, 0.5, "no failures", ha="center", va="center",
                        transform=axes[1, 0].transAxes)

    style(axes[1, 1], "VAF distribution")
    axes[1, 1].yaxis.grid(False)
    bins = np.linspace(0, 1, 21)
    for key, color in zip(labels, colors):
        if vafs.get(key):
            axes[1, 1].hist(vafs[key], bins=bins, color=color, alpha=0.7,
                            edgecolor="#495057", linewidth=0.4,
                            label=f"{key} (n={len(vafs[key])})", zorder=3)
    axes[1, 1].axvline(MIN_VAF,     ls="--", lw=1, color="#6c757d", label=f"min VAF {MIN_VAF}")
    axes[1, 1].axvline(HOM_VAF_MIN, ls=":",  lw=1, color="#343a40", label=f"hom-alt {HOM_VAF_MIN}")
    axes[1, 1].set_xlabel("VAF")
    axes[1, 1].legend(fontsize=7.5, frameon=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="#f8f9fa")
    print(f"graph saved → {out_path}")


# ── command handlers ──────────────────────────────────────────────────────

def cmd_fastq(args):
    abort_if_missing(args.fastq)
    l1 = layer1(args.fastq)
    print_layer("fastq", l1, not args.no_color)
    print_summary(l1, [], [], not args.no_color)


def cmd_bam(args):
    abort_if_missing(args.bam)
    l2 = layer2(args.bam)
    print_layer("bam", l2, not args.no_color)
    print_summary([], l2, [], not args.no_color)


def cmd_vcf(args):
    abort_if_missing(args.vcf, args.bam)
    variants = layer3(args.vcf, args.bam, region=args.region)
    passed   = sum(1 for v in variants if v[5] == "PASS")
    verdict  = verdict_of([(None, v[5], None) for v in variants]) if variants else "PASS"
    print_layer("vcf", [("VARIANTS", verdict, f"{passed}/{len(variants)} pass")], not args.no_color)
    print_variants(variants, args.verbose, not args.no_color)
    print_summary([], [], variants, not args.no_color)


def cmd_all(args):
    abort_if_missing(args.fastq, args.bam, args.vcf)
    l1 = layer1(args.fastq)
    print_layer("fastq", l1, not args.no_color)
    l2 = layer2(args.bam)
    print_layer("bam", l2, not args.no_color)
    variants = layer3(args.vcf, args.bam, region=args.region)
    passed   = sum(1 for v in variants if v[5] == "PASS")
    verdict  = verdict_of([(None, v[5], None) for v in variants]) if variants else "PASS"
    print_layer("vcf", [("VARIANTS", verdict, f"{passed}/{len(variants)} pass")], not args.no_color)
    print_variants(variants, args.verbose, not args.no_color)
    print_summary(l1, l2, variants, not args.no_color)
    make_graph(variants, args.plot_out)


# ── entry point ───────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Sequencing Evidence Validator")
    p.add_argument("--no-color", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    s1 = sub.add_parser("fastq")
    s1.add_argument("fastq", metavar="FASTQ")

    s2 = sub.add_parser("bam")
    s2.add_argument("bam", metavar="BAM")

    s3 = sub.add_parser("vcf")
    s3.add_argument("vcf",  metavar="VCF")
    s3.add_argument("bam",  metavar="BAM")
    s3.add_argument("--region",  metavar="CHR:START-END")
    s3.add_argument("--verbose", action="store_true")

    sa = sub.add_parser("all")
    sa.add_argument("fastq", metavar="FASTQ")
    sa.add_argument("bam",   metavar="BAM")
    sa.add_argument("vcf",   metavar="VCF")
    sa.add_argument("--region",   metavar="CHR:START-END")
    sa.add_argument("--verbose",  action="store_true")
    sa.add_argument("--plot-out", default="results_graph.png", metavar="FILE")

    args = p.parse_args()
    {"fastq": cmd_fastq, "bam": cmd_bam, "vcf": cmd_vcf, "all": cmd_all}[args.cmd](args)


if __name__ == "__main__":
    main()
