# Red vs Black Carp DNA-seq Variant Analysis

This repository documents a small DNA-seq project comparing red and black color
groups in *Cyprinus carpio* using whole-genome sequencing data from NCBI SRA.
The project was designed to mirror the structure of a previous RNA-seq
reproduction project, but with a DNA-seq workflow focused on read QC, alignment,
BAM quality control, duplicate marking, variant calling, and population-genetic
summary plots.

The biological question is:

> Can chromosome-level DNA-seq variant analysis identify candidate genomic
> regions that differ between red and black carp color groups?

The final analysis focused on chromosome A20 (`NC_056591.1`) because full
whole-genome GATK variant calling was too large and slow for a mini-project.
This made the project more manageable while still preserving a complete,
realistic DNA-seq analysis path.

## What This Project Contains

- A reproducible DNA-seq workflow using open-source tools.
- A project-local virtual environment specification in `envs/dnaseq.yml`.
- Metadata for 5 red and 5 black *Cyprinus carpio* WGS samples.
- Scripts for SRA download, environment execution, software checking, and A20 plotting.
- Documentation tables explaining software, resources, output files, and Picard metrics.
- Chromosome A20 variant summary tables and figures.
- A written interpretation of the A20 plots and candidate red-vs-black differentiation signal.

## Project Thought Process

The original idea was to do a project similar to an RNA-seq workflow, but with
DNA-seq and a nonmodel animal. Koi carp color inheritance was the preferred
biological theme, but confirmed ornamental koi datasets with clean red, white,
and black sample labels were difficult to find in SRA.

The project therefore shifted to a better documented public dataset:
BioProject `PRJNA824207`, titled `Genetic mechanism of the reddish color of
Hebao red carp`. This dataset is not confirmed ornamental koi, but it is still
*Cyprinus carpio* and has explicit red and black sample labels. That made it a
stronger choice for a reproducible red-vs-black DNA-seq comparison.

The initial goal was to compare 5 red and 5 black samples across the genome.
However, whole-genome GATK `HaplotypeCaller` took too long for the available
time and storage. To keep the project realistic and finishable, the workflow was
narrowed to one chromosome: A20 (`NC_056591.1`). This allowed the full downstream
logic to be completed: GVCF calling, cohort genotyping, SNP filtering, allele
frequency comparison, PCA, SNP density, and windowed FST.

This is a useful compromise for a mini-project: the pipeline remains authentic,
but the data size stays manageable.

## Dataset

The selected dataset is from NCBI BioProject `PRJNA824207`.

| Color group | Sample labels | SRA runs |
|---|---|---|
| Red | `red_19`, `red_18`, `red_17`, `red_16`, `red_15` | `SRR18681855`, `SRR18681856`, `SRR18681857`, `SRR18681858`, `SRR18681859` |
| Black | `black_9`, `black_8`, `black_7`, `black_6`, `black_5` | `SRR18681860`, `SRR18681861`, `SRR18681862`, `SRR18681863`, `SRR18681864` |

Project metadata files:

```text
config/samples_red_black_10samples.tsv
config/sra_runs_red_black_10samples.txt
```

Important wording: this project should be described as a red-vs-black
*Cyprinus carpio* color dataset, not as confirmed ornamental koi. It is relevant
to koi/color genetics because koi and common carp are the same species complex,
but the SRA dataset labels are Hebao red carp rather than ornamental koi.

## Reference Genome

The reference genome used was the RefSeq *Cyprinus carpio* assembly:

```text
GCF_018340385.1_ASM1834038v1
```

Important local reference files:

```text
references/carp.fa
references/carp.fa.fai
references/carp.dict
references/bowtie2/carp.*.bt2
references/carp_refseq/GCF_018340385.1_ASM1834038v1_genomic.gff.gz
```

The chromosome analyzed downstream was:

```text
NC_056591.1
```

This corresponds to chromosome A20 in the assembly report.

## Software And Environment

All major analysis tools are free and open source. The project used a local
micromamba environment named `koi-dnaseq`.

Environment file:

```text
envs/dnaseq.yml
```

Main tools:

| Step | Software | Purpose |
|---|---|---|
| Download reads | SRA Toolkit | Download SRA runs and convert to FASTQ |
| Raw read QC | FastQC | Inspect sequencing quality before processing |
| Read trimming | AdapterRemoval | Remove adapters and trim low-quality sequence |
| QC aggregation | MultiQC | Combine QC outputs into one report |
| Alignment | Bowtie2 | Map reads to the carp reference genome |
| BAM processing | SAMtools | Convert, sort, index, and summarize alignments |
| Alignment QC | Qualimap | Assess BAM quality and coverage |
| Read groups/duplicates | Picard | Add read groups and mark duplicate reads |
| Variant calling | GATK | Call per-sample GVCFs and joint cohort variants |
| VCF filtering/summaries | BCFtools | Filter, index, and summarize VCFs |
| Plotting | R / tidyverse / ggplot2 | Make A20 summary plots |

Run commands inside the environment with:

```bash
bash scripts/run_in_env.sh <command>
```

Example:

```bash
bash scripts/run_in_env.sh fastqc --version
bash scripts/run_in_env.sh bash scripts/check_software.sh
```

## Workflow Summary

The intended full DNA-seq workflow was:

```text
SRA Toolkit
  -> FastQC
  -> AdapterRemoval
  -> FastQC / MultiQC
  -> Bowtie2
  -> SAMtools
  -> Qualimap
  -> Picard AddOrReplaceReadGroups
  -> Picard MarkDuplicates
  -> GATK HaplotypeCaller
  -> GATK CombineGVCFs / GenotypeGVCFs
  -> GATK / BCFtools filtering
  -> R plots and interpretation
```

The completed downstream analysis used chromosome A20 only.

## Key Commands

Download SRA reads:

```bash
RUN_LIST=config/sra_runs_red_black_10samples.txt \
THREADS=8 \
MAX_SIZE=20G \
bash scripts/run_in_env.sh bash scripts/download_sra.sh
```

Run FastQC on FASTQ files:

```bash
mkdir -p results/fastqc/raw
bash scripts/run_in_env.sh fastqc -t 8 -o results/fastqc/raw data/sra/fastq/*.fastq.gz
```

Run MultiQC:

```bash
mkdir -p results/multiqc/raw
bash scripts/run_in_env.sh multiqc results/fastqc/raw -o results/multiqc/raw
```

Build the Bowtie2 index:

```bash
mkdir -p references/bowtie2
bash scripts/run_in_env.sh bowtie2-build references/carp.fa references/bowtie2/carp
```

Create reference indexes for SAMtools and GATK/Picard:

```bash
bash scripts/run_in_env.sh samtools faidx references/carp.fa
bash scripts/run_in_env.sh gatk CreateSequenceDictionary \
  -R references/carp.fa \
  -O references/carp.dict
```

Make A20 plots from the filtered SNP VCF:

```bash
bash scripts/run_in_env.sh Rscript scripts/plot_a20_variants.R
```

More detailed command examples and output explanations are in:

```text
docs/software_resource_table.md
docs/output_files_table.md
docs/picard_dup_metrics_explained.md
```

## Main Results

The final filtered A20 SNP dataset contained:

| Category | Value |
|---|---:|
| PASS biallelic SNPs | 109,849 |
| Samples | 10 |
| Red samples | 5 |
| Black samples | 5 |

Main result files:

```text
results/gatk_chrA20/filtered/red_black_A20.PASS.biallelic_snps.vcf.gz
results/gatk_chrA20/fst/a20_windowed_fst_100kb.tsv
results/gatk_chrA20/tables/a20_site_allele_frequencies.tsv
results/gatk_chrA20/tables/a20_pca_samples.tsv
results/gatk_chrA20/tables/a20_snp_density_100kb.tsv
results/gatk_chrA20/tables/a20_variant_summary.tsv
results/gatk_chrA20/tables/a20_fst_threshold_summary.tsv
results/gatk_chrA20/tables/a20_top_fst_window_informative_snps.tsv
results/gatk_chrA20/tables/a20_candidate_regions_fst_ge_0_45.tsv
results/gatk_chrA20/tables/a20_candidate_snps_fst_ge_0_45_abs_af_diff_ge_0_75.tsv
```

Main figures:

```text
results/gatk_chrA20/plots/a20_windowed_fst_100kb.png
results/gatk_chrA20/plots/a20_windowed_af_difference_100kb.png
results/gatk_chrA20/plots/a20_top_af_differences.png
results/gatk_chrA20/plots/a20_pca_red_black.png
results/gatk_chrA20/plots/a20_snp_density_100kb.png
results/gatk_chrA20/plots/a20_variant_summary.png
results/gatk_chrA20/plots/a20_fst_threshold_comparison.png
```

Detailed plot interpretation:

```text
docs/a20_plot_findings.md
results/gatk_chrA20/plots/a20_plot_findings.md
```

## Main Findings

The A20 analysis found one standout 100 kb window with FST above 0.6:

| Chromosome | Window start | Window end | Informative SNPs | Mean absolute allele-frequency difference | FST |
|---|---:|---:|---:|---:|---:|
| `NC_056591.1` / A20 | 14,500,001 | 14,600,000 | 3 | 0.583 | 0.750 |

This was the strongest candidate red-vs-black differentiated region on
chromosome A20. However, the window only had 3 informative SNPs contributing to
the FST calculation, so the signal is interesting but fragile.

The updated FST plot labels this peak directly and includes two visual guide
thresholds. The strict `FST >= 0.75` threshold captures the single top window,
while a broader `FST >= 0.45` threshold captures four candidate windows, 95
informative SNPs, and 16 strong SNPs with `|red allele frequency - black allele
frequency| >= 0.75`.

The most notable SNP in that region was:

```text
NC_056591.1:14557124 T>A
```

At this SNP, the called genotypes showed:

| Group | Alternate allele frequency |
|---|---:|
| Red | 0.00 |
| Black | 1.00 |

This looks like a strong red-vs-black allele-frequency difference, but it was
based on very few called alleles: 4 red alleles and 2 black alleles. For that
reason, it should be treated as a candidate marker, not a confirmed causal
mutation.

Based on the RefSeq GFF annotation, this SNP is not inside a confirmed
protein-coding exon. It is approximately 3.8 kb upstream of pseudogene
`LOC122149054`, and nearby SNPs fall in a pseudogene-rich region. The same
100 kb window also contains a protein-coding gene, `LOC109045109`, annotated as
forkhead box protein F2-like. Forkhead-box genes can be biologically interesting
because they encode transcription factors, but the differentiated SNPs do not
directly hit the coding sequence of this gene.

The most defensible interpretation is:

> Chromosome A20 contains several candidate red-vs-black differentiated regions.
> The strongest peak is at 14.5-14.6 Mb, where the strongest candidate SNP is
> `NC_056591.1:14557124 T>A`, but the evidence is limited by missing genotype
> calls and a small number of informative SNPs. A broader threshold also points
> to candidate regions near 4.4-4.6 Mb and 13.9-14.0 Mb. These signals may tag
> linked variation near color-associated regions, but they should not be claimed
> as proven pigmentation mutations.

## Plot-Level Findings

The windowed FST plot identified the strongest local differentiation on A20,
with the top peak at 14.5-14.6 Mb. Additional high-FST windows occurred near
13.9-14.0 Mb and 4.4-4.6 Mb.

The allele-frequency difference plot supported the FST result because the same
14.5-14.6 Mb region had the highest mean absolute allele-frequency difference.

The top allele-frequency difference plot highlighted individual SNPs with large
red-vs-black differences, but several of these sites had limited genotype calls.
This is why window-level and annotation-level interpretation are important.

The updated variant summary plot now separates SNPs into red-higher,
black-higher, similar, and missing-one-color-group categories. This makes the
summary figure reflect the red-vs-black research question instead of only
showing total dataset counts.

The PCA plot did not separate all red samples from all black samples. Some red
and black fish were close together in PCA space. This suggests that color group
is not the only major source of genetic variation on A20. The red-vs-black
signal is more local than chromosome-wide.

The SNP density plot showed that SNPs are not evenly distributed across A20. It
also made the top FST peak easier to interpret: the strongest peak is based on a
small number of informative SNPs rather than a broad SNP-rich region.

## Limitations

This project is intentionally small, so several limitations matter:

- Only chromosome A20 was analyzed downstream.
- The sample size was 5 red and 5 black fish.
- The dataset is red-vs-black carp, not confirmed ornamental koi.
- BQSR was not performed because this nonmodel-style dataset does not have a
  trusted known-sites database.
- The high-FST peak is based on only 3 informative SNPs.
- The strongest candidate SNP has few called alleles.
- The analysis identifies candidate regions, not experimentally validated color
  mutations.

These limitations do not make the project invalid. They simply shape the
conclusion: this is a reproducible mini DNA-seq workflow that finds candidate
red-vs-black differentiation on A20, not a definitive genetic proof of color
inheritance.

## Repository Structure

```text
config/                 sample metadata and SRA run lists
docs/                   detailed project notes and interpretation
envs/                   conda/micromamba environment YAML
scripts/                reusable workflow and plotting scripts
results/gatk_chrA20/    final A20 variant tables, plots, and selected outputs
```

## Reproducibility Notes

This project can be reproduced by:

1. Creating the `koi-dnaseq` environment from `envs/dnaseq.yml`.
2. Downloading the selected SRA runs listed in `config/sra_runs_red_black_10samples.txt`.
3. Downloading the RefSeq carp reference genome `GCF_018340385.1_ASM1834038v1`.
4. Building the Bowtie2, SAMtools, and Picard/GATK reference indexes.
5. Running QC, alignment, BAM processing, duplicate marking, and GATK variant calling.
6. Restricting GATK to chromosome A20 if keeping the project small.
7. Running `scripts/plot_a20_variants.R` to regenerate the figures and tables.

The most important final interpretation file is:

```text
docs/a20_plot_findings.md
```

That file explains each figure and gives the cautious biological interpretation
of the FST peak and candidate SNP.
