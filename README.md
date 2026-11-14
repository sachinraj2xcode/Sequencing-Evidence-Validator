# Sequencing Evidence Validator

A small command-line tool that checks sequencing data in three layers and
reports a PASS / WARN / ERROR verdict for each.

| Layer | Input | What it checks |
|-------|-------|----------------|
| 1 | FASTQ | mean base quality, read length, GC content, N fraction |
| 2 | BAM   | file integrity, mapping rate, duplicate rate, mean MAPQ |
| 3 | VCF + BAM | per-variant depth, alt count, VAF, strand balance, genotype fit |

## Requirements

- Python 3.8+
- [samtools](https://www.htslib.org/) on your `PATH` (needed for BAM checks)
- Optional Python packages (see `requirements.txt`):
  - `pysam` for pileups (falls back to `samtools mpileup` if missing)
  - `matplotlib` and `numpy` for the results graph

```bash
pip install -r requirements.txt
```

## Usage

```bash
python validate.py fastq reads.fastq.gz
python validate.py bam   aln.bam
python validate.py vcf   calls.vcf.gz aln.bam
python validate.py all   reads.fastq.gz aln.bam calls.vcf.gz
```

Gzipped FASTQ and VCF files are detected by the `.gz` extension.

### Options

| Option | Applies to | Description |
|--------|-----------|-------------|
| `--no-color` | all (put before the command) | disable coloured output |
| `--region CHR:START-END` | `vcf`, `all` | only check variants in a region (`CHR` alone also works) |
| `--verbose` | `vcf`, `all` | list every variant, not just warnings and errors |
| `--plot-out FILE` | `all` | where to save the graph (default `results_graph.png`) |

Example:

```bash
python validate.py --no-color vcf calls.vcf.gz aln.bam --region chr1:100000-200000 --verbose
```

## Output

Each layer prints its checks, followed by a summary and an overall verdict:

- **PASS**: all checks met
- **WARN**: something looks off but isn't disqualifying
- **ERROR**: a check failed

The `all` command also saves a four-panel graph of the variant results
(overall verdicts, SNV vs indel, fail reasons, VAF distribution).

## Thresholds

All cut-offs are constants at the top of `validate.py` (for example
`MIN_MEAN_Q`, `MIN_MAPPING_RATE`, `MIN_DEPTH`, `MIN_VAF`). Edit them there to
suit your data.
