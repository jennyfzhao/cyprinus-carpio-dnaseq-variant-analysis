# Virtual Environment Setup

Use a Conda/Mamba environment for this project. A Python `venv` is not enough
because the pipeline depends on compiled command-line tools such as Bowtie2,
SAMtools, Picard, GATK, FastQC, and AdapterRemoval.

## Recommended: Miniforge + Mamba

Install Miniforge:

https://github.com/conda-forge/miniforge

Then open a new terminal and run:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq

mamba env create -f envs/dnaseq.yml
mamba activate koi-dnaseq

bash scripts/check_software.sh
```

If `mamba` is not available but `conda` is available:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq

conda env create -f envs/dnaseq.yml
conda activate koi-dnaseq

bash scripts/check_software.sh
```

## Confirm You Are Inside the Environment

After activation, your terminal prompt should usually include:

```text
(koi-dnaseq)
```

You can also check:

```bash
which fastqc
which bowtie2
which samtools
which gatk
which picard
```

The paths should point inside a Conda environment named `koi-dnaseq`.

## Run the SRA Download Inside the Environment

For the 9-sample dataset:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq
conda activate koi-dnaseq

RUN_LIST=config/sra_runs_9samples.txt THREADS=8 MAX_SIZE=20G bash scripts/download_sra.sh
```

Use `mamba activate koi-dnaseq` only if your shell is configured for mamba.
Most installations still use `conda activate koi-dnaseq` after the environment
has been created.

## If You Are on Apple Silicon

If solving or running GATK/Picard fails on an Apple Silicon Mac, use one of
these approaches:

- Install the x86_64 version of Miniforge/Miniconda under Rosetta.
- Run the same Conda environment on a Linux server/HPC.
- Use a Linux container such as Apptainer/Singularity if your class or cluster
  supports it.

## Remove the Environment

If you need to start over:

```bash
conda deactivate
conda env remove -n koi-dnaseq
```

Then recreate it from `envs/dnaseq.yml`.

## Project-local Environment Already Installed

This workspace also has a project-local micromamba installation:

```text
.local/bin/micromamba
.micromamba/envs/koi-dnaseq
```

You can run any command inside that environment without activating Conda:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq

bash scripts/run_in_env.sh bash scripts/check_software.sh
bash scripts/run_in_env.sh fastqc --version
bash scripts/run_in_env.sh samtools --version
```

To download the red/black SRA dataset using the project-local environment:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq

bash scripts/run_in_env.sh bash scripts/download_sra.sh
```

To pass the 10-sample red/black run list:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq

RUN_LIST=config/sra_runs_red_black_10samples.txt THREADS=8 MAX_SIZE=20G \
  bash scripts/run_in_env.sh bash scripts/download_sra.sh
```
