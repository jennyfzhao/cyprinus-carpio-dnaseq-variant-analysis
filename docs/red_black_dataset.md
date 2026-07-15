# Red vs Black Carp WGS Dataset

This project uses 10 paired-end WGS samples from BioProject `PRJNA824207`:
`Genetic mechanism of the reddish color of Hebao red carp`.

This is a confirmed red-vs-black color dataset in *Cyprinus carpio*. It is not
ornamental koi specifically, but it is a better fit for the DNA-seq variant
calling pipeline because the SRA labels explicitly include red and black sample
groups.

## Selected Samples

| Color | Samples | Runs |
|---|---|---|
| Red | red_19, red_18, red_17, red_16, red_15 | SRR18681855, SRR18681856, SRR18681857, SRR18681858, SRR18681859 |
| Black | black_9, black_8, black_7, black_6, black_5 | SRR18681860, SRR18681861, SRR18681862, SRR18681863, SRR18681864 |

Metadata table:

```text
config/samples_red_black_10samples.tsv
```

Run list:

```text
config/sra_runs_red_black_10samples.txt
```

## Download

Activate the virtual environment first:

```bash
cd /Users/jennyfzhao/Work/task/2026-06-15_dnaseq
conda activate koi-dnaseq
```

Then download and convert the SRA files:

```bash
RUN_LIST=config/sra_runs_red_black_10samples.txt THREADS=8 MAX_SIZE=20G bash scripts/download_sra.sh
```

The FASTQ files will be written to:

```text
data/sra/fastq/
```

The downloaded SRA cache files will be written to:

```text
data/sra/cache/
```

Logs will be written to:

```text
data/sra/logs/
```

## Storage Note

Each SRA file is about 4.9-5.1 GB compressed. Ten samples will require about
50 GB for the SRA cache alone, plus substantially more space for uncompressed
FASTQ files, trimmed reads, BAMs, duplicate-marked BAMs, and GATK outputs.

For a full 10-sample run, plan for at least 300-500 GB free disk space. If that
is too large, run the pipeline first on one red and one black sample, or
downsample reads before alignment.
