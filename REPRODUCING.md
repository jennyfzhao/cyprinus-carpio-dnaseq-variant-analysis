# Reproducing The Red vs Black Carp DNA-seq Project

This repository keeps the lightweight project files needed to understand and reproduce the analysis: metadata, scripts, environment definition, documentation, final chromosome A20 tables, final figures, and the final filtered A20 SNP VCF.

Large files such as FASTQ reads, BAM files, reference FASTA files, Bowtie2 indexes, local conda environments, and large GATK intermediates are intentionally regenerated from public sources.

## 1. Create The Software Environment

The analysis was designed around a conda/micromamba environment named `koi-dnaseq`.

```bash
conda env create -f envs/dnaseq.yml
conda activate koi-dnaseq
```

If using micromamba:

```bash
micromamba create -f envs/dnaseq.yml
micromamba activate koi-dnaseq
```

Check software availability:

```bash
bash scripts/check_software.sh
```

## 2. Download The SRA Reads

The selected BioProject is `PRJNA824207`, with 5 red and 5 black *Cyprinus carpio* samples listed in `config/sra_runs_red_black_10samples.txt` and `config/samples_red_black_10samples.tsv`.

```bash
RUN_LIST=config/sra_runs_red_black_10samples.txt \
THREADS=8 \
MAX_SIZE=20G \
bash scripts/download_sra.sh
```

This writes reads under:

```text
data/sra/fastq/
data/sra/cache/
data/sra/logs/
```

## 3. Download The Reference Genome

The project used RefSeq assembly `GCF_018340385.1_ASM1834038v1` for *Cyprinus carpio*.

```bash
mkdir -p references/carp_refseq references/bowtie2

curl -L -o references/carp_refseq/GCF_018340385.1_ASM1834038v1_genomic.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/340/385/GCF_018340385.1_ASM1834038v1/GCF_018340385.1_ASM1834038v1_genomic.fna.gz

curl -L -o references/carp_refseq/GCF_018340385.1_ASM1834038v1_genomic.gff.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/340/385/GCF_018340385.1_ASM1834038v1/GCF_018340385.1_ASM1834038v1_genomic.gff.gz

gunzip -c references/carp_refseq/GCF_018340385.1_ASM1834038v1_genomic.fna.gz > references/carp.fa
```

Build indexes:

```bash
bowtie2-build references/carp.fa references/bowtie2/carp
samtools faidx references/carp.fa
gatk CreateSequenceDictionary -R references/carp.fa -O references/carp.dict
```

## 4. Run QC, Alignment, And BAM Processing

Typical workflow:

```text
FastQC -> AdapterRemoval -> FastQC/MultiQC -> Bowtie2 -> SAMtools -> Qualimap -> Picard
```

The detailed software and command patterns are documented in:

```text
docs/software_resource_table.md
docs/output_files_table.md
docs/picard_dup_metrics_explained.md
```

## 5. Run GATK On Chromosome A20

The final mini-project used chromosome A20 only:

```text
NC_056591.1
```

The GATK workflow was:

```text
HaplotypeCaller per sample in GVCF mode
CombineGVCFs
GenotypeGVCFs
VariantFiltration
BCFtools filtering to PASS biallelic SNPs
```

The lightweight final VCF kept in this repository is:

```text
results/gatk_chrA20/filtered/red_black_A20.PASS.biallelic_snps.vcf.gz
```

## 6. Recreate Final Figures And Tables

Run:

```bash
Rscript scripts/plot_a20_variants.R
```

Primary outputs:

```text
results/gatk_chrA20/plots/
results/gatk_chrA20/tables/
results/gatk_chrA20/fst/
```

The main interpretation is in:

```text
docs/a20_plot_findings.md
```

## Main Reproducibility Notes

- The dataset is red-vs-black *Cyprinus carpio*, not confirmed ornamental koi.
- The full genome was too large for this mini-project, so downstream variant analysis was restricted to chromosome A20.
- BQSR was not performed because there was no trusted known-sites resource for this nonmodel-style project.
- The strongest A20 FST peak is a candidate region, not a validated causal color mutation.
