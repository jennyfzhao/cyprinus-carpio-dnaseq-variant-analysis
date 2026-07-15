#!/usr/bin/env bash
set -euo pipefail

echo "SRA Toolkit:"
prefetch --version
fasterq-dump --version

echo
echo "FastQC:"
fastqc --version

echo
echo "AdapterRemoval:"
AdapterRemoval --version

echo
echo "Bowtie2:"
bowtie2 --version | head -n 1

echo
echo "SAMtools:"
samtools --version | head -n 1

echo
echo "BCFtools:"
bcftools --version | head -n 1

echo
echo "HTSlib/tabix:"
tabix --version | head -n 1

echo
echo "Qualimap:"
qualimap --version || true

echo
echo "Picard:"
picard MarkDuplicates --version || true

echo
echo "GATK:"
gatk --version

echo
echo "MultiQC:"
multiqc --version

echo
echo "R:"
R --version | head -n 1
