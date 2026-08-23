# read-evidence-validator
Layer 4: Read-evidence validator — checks whether variant calls in a VCF file are supported by read-level evidence in a BAM file.


Goal Build a Python tool that takes a VCF and a BAM file, and for each variant call in the VCF checks whether the aligned reads in the BAM actually support that call at that position. Output the results in the same format as vcf-validator.

figure out pysam
open a BAM and get reads at a position
connect that to the VCF records
check depth
check ALT is in the reads
check base quality
check mapping quality
outputs same as vcf-validator
test on GIAB
paper
