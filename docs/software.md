# DNA-seq Software Stack

This project uses free, open-source, or public-domain software for every
analysis step. The preferred installation route is Bioconda because it gives a
reproducible environment without manually compiling each tool.

For a step-by-step virtual environment guide, see
`docs/virtual_environment.md`.

## Install

Install Miniforge or Mambaforge, then create the environment:

```bash
mamba env create -f envs/dnaseq.yml
mamba activate koi-dnaseq
bash scripts/check_software.sh
```

On an Apple Silicon Mac, if the GATK or Java-heavy packages fail to solve, use
an x86_64 Conda/Miniconda install under Rosetta or run the environment on a
Linux server/HPC. Broad's GATK documentation specifically notes ARM Macs as a
special case for Conda-based installs.

If `mamba` is not available, `conda` also works, but dependency solving is
usually slower:

```bash
conda env create -f envs/dnaseq.yml
conda activate koi-dnaseq
bash scripts/check_software.sh
```

## Tools

| Step | Tool | Role | License/source status |
|---|---|---|---|
| SRA download | SRA Toolkit | Download SRA runs and convert to FASTQ | NCBI public-domain notice |
| Raw QC | FastQC | Per-sample FASTQ quality reports | GPL-3.0, with bundled Apache-2.0 dependency notice |
| Adapter trimming | AdapterRemoval | Adapter/quality trimming and paired-end handling | GPL-3.0 |
| Post-trim QC | FastQC / MultiQC | Confirm trimming improved read quality | FastQC GPL-3.0; MultiQC GPL-3.0-or-later |
| Alignment | Bowtie2 | Align reads to the carp reference genome | GPL-3.0 |
| BAM processing | SAMtools | Sort, index, flagstat, stats | MIT |
| Alignment QC | Qualimap | Mapping/coverage quality reports | GPL-2.0-or-later; available through Bioconda |
| Duplicate marking | Picard | MarkDuplicates and sequence dictionary creation | MIT |
| Variant calling | GATK4 | HaplotypeCaller, GenotypeGVCFs, VariantFiltration | Apache-2.0 |
| Variant summaries | BCFtools / R | VCF summaries and project figures | BCFtools MIT; R GPL |

## Notes for This Project

- Use GATK hard filtering, not VQSR. A koi/common carp project will not have a
  trusted training set of variants like human projects do.
- Keep the dataset small by using two to four samples and/or downsampling reads
  during the first complete run.
- For full-genome carp alignment, expect large intermediate BAM files. A
  targeted region or downsampled workflow is better for a class-size project.
- Picard and GATK need Java. The environment pins `openjdk=17`, which is
  compatible with current Picard and GATK4 packages.

## Smoke Test

After activating the environment, this should finish without errors:

```bash
bash scripts/check_software.sh
```

The script only checks that the software is installed and callable. It does not
download data.
