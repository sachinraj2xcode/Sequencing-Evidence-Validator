# Sequencing Evidence Validator

Three-layer framework for read-level trust in variant calls.
> Raj (2026) — *End-to-End Sequencing Evidence Validation*

## Install

```bash
pip install -r requirements.txt
```

`samtools` must be on your `PATH`.

## Usage

```bash
# Demo (no files needed)
python validate.py --demo --plot

# Real files
python validate.py --fastq reads.fastq.gz --bam aln.bam --vcf calls.vcf.gz \
                   --region chr20:1000000-1100000 --plot
```

## Layers

| Layer | Input | Checks |
|-------|-------|--------|
| 1 | FASTQ | Base quality, read length, GC content, N-fraction |
| 2 | BAM | Integrity, mapping rate, duplicate rate, MAPQ |
| 3 | VCF + BAM | Depth, alt count, VAF, strand balance, genotype consistency |

## Licence

MIT
