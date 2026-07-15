# Picard MarkDuplicates Metrics Explained

Example file:

```text
results/picard/SRR18681855.dup_metrics.txt
```

This file was created by Picard `MarkDuplicates`. It summarizes how many reads were examined and how many were marked as duplicate reads in `SRR18681855`.

## Run Information

| Part of file | What it says | Meaning |
|---|---|---|
| `MarkDuplicates INPUT=[results/picard/SRR18681855.rg.bam]` | Picard used `SRR18681855.rg.bam` as input | This is the BAM after read groups were added |
| `OUTPUT=results/picard/SRR18681855.dedup.bam` | Picard wrote `SRR18681855.dedup.bam` | This is the BAM file to use for GATK variant calling |
| `METRICS_FILE=results/picard/SRR18681855.dup_metrics.txt` | Picard wrote this metrics file | This is the text file being interpreted here |
| `CREATE_INDEX=true` | Picard created a BAM index | This should produce `SRR18681855.dedup.bai` |
| `REMOVE_DUPLICATES=false` | Duplicates were not removed | Duplicate reads stay in the BAM but are marked with the duplicate flag |
| `DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES` | Picard chooses the best representative read using base qualities | When duplicate reads exist, the read with stronger base quality evidence is kept as the representative |
| `VALIDATION_STRINGENCY=STRICT` | Picard used strict validation | The BAM had to follow expected formatting rules |
| `Started on: Fri Jun 26 11:38:23 EDT 2026` | Time the Picard job started | Useful for record keeping |

## Main Metrics

| Metric | Value for `SRR18681855` | What it means | Interpretation |
|---|---:|---|---|
| `LIBRARY` | `SRR18681855` | The sequencing library/sample label | This sample was treated as one library |
| `UNPAIRED_READS_EXAMINED` | `2,165,481` | Mapped reads that Picard examined as single/unpaired reads | These may come from mates where only one read mapped or from reads not usable as a proper pair |
| `READ_PAIRS_EXAMINED` | `24,808,704` | Proper read pairs Picard examined for duplication | This is the main paired-end duplicate-detection pool |
| `SECONDARY_OR_SUPPLEMENTARY_RDS` | `0` | Reads marked as secondary or supplementary alignments | None were present in this BAM, which is fine |
| `UNMAPPED_READS` | `3,887,107` | Reads that did not map to the reference genome | These are not useful for duplicate marking |
| `UNPAIRED_READ_DUPLICATES` | `59,641` | Single/unpaired reads marked as duplicates | A small number of single-read duplicates were found |
| `READ_PAIR_DUPLICATES` | `0` | Paired-end read pairs marked as duplicates | Picard did not identify duplicate read pairs in this sample |
| `READ_PAIR_OPTICAL_DUPLICATES` | `0` | Duplicate pairs likely caused by optical/imaging artifacts | None were detected; this may also depend on whether read names encode optical-position information in a recognizable way |
| `PERCENT_DUPLICATION` | `0.001152` | Fraction of examined reads/pairs marked duplicate | This is about `0.1152%`, which is very low |
| `ESTIMATED_LIBRARY_SIZE` | blank | Picard's estimate of the number of unique molecules in the library | Blank because Picard did not have enough paired duplicate information to estimate it |

## Plain-English Summary

Picard found **very low duplication** in `SRR18681855`.

The most important value is:

```text
PERCENT_DUPLICATION = 0.001152
```

That means approximately:

```text
0.1152% duplicate reads
```

This is good. It suggests that most reads represent independent DNA fragments rather than repeated sequencing of the same original molecule.

Because `REMOVE_DUPLICATES=false`, Picard did **not** delete duplicates. It only marked them in:

```text
results/picard/SRR18681855.dedup.bam
```

That duplicate-marked BAM is the correct file to use for GATK.
