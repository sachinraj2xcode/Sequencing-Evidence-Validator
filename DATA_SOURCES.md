# Finding Data

Needed FASTQ, BAM, and VCF from the same sample. Used HG002 since it's the standard reference sample for this kind of testing.

## FASTQ

Don't download the whole 7 GB file. Just stream enough reads for testing.

```bash
wget -O - "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_1.fastq.gz" \
  2>/dev/null | zcat | head -40000 > test_R1.fastq
```

That's about 10k reads, which is enough.

HG002 accession: `SRR622461`

`fasterq-dump` also works, but it's more steps than I need here.

## BAM

BAM is huge, 100+ GB, so just grab the region I need.

First make sure `samtools` has `libcurl`:

```bash
samtools --version
```

Then:

```bash
samtools view -b \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam \
  20:1000000-1100000 > test.bam

samtools index test.bam
```

The old BAM URL changed. If this one breaks, check the GIAB data index repo:

```text
https://github.com/genome-in-a-bottle/giab_data_indexes
```

## VCF

Same idea. Just grab the region with `tabix`.

```bash
tabix -h \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz \
  20:1000000-1100000 > test.vcf
```

Important: use `-h`, otherwise the VCF headers are missing and the validator won't like it.

## Things I almost missed

**`20` vs `chr20`**

BAM and VCF need to use the same chromosome naming.

```bash
samtools view -H test.bam | grep "^@SQ" | head -3
grep -v "^#" test.vcf | head -1 | cut -f1
```

Both should agree.

**EOFError**

Usually means I ended up with an incomplete gzip file. Delete it and run the streaming command again.

**SRA Toolkit**

If Debian asks about Debian Med during installation, just accept the default.

## Test region

Using:

```text
20:1000000-1100000
```

This is the region used in the paper. Around 156 variants and only takes a few seconds.

If I want a bigger test:

```text
20:1000000-10000000
```

That's all I need for now.