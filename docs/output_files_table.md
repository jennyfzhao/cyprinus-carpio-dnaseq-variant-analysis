# Output Files Table

This table explains the main files created so far in the red vs black *Cyprinus carpio* DNA-seq project.

| File or pattern | Created by | What it means | What it looks like inside / how to peek | What it is used for | Important for final analysis? |
|---|---|---|---|---|---|
| `data/sra/cache/<RUN>/<RUN>.sra` | SRA Toolkit `prefetch` | NCBI's downloaded SRA archive file for one sequencing run | Binary/container-like SRA format; view indirectly with `fasterq-dump` or `vdb-dump`. Example: `bash scripts/run_in_env.sh vdb-dump --info data/sra/cache/SRR18681855/SRR18681855.sra` | Stored raw source data; used by `fasterq-dump` to make FASTQ files | Keep if you want to regenerate FASTQs without re-downloading |
| `data/sra/logs/<RUN>.prefetch.log` | SRA Toolkit `prefetch` | Log of the SRA download step | Plain text. Looks like download/progress/status messages from `prefetch` | Troubleshooting download failures or confirming the download ran | Useful for methods/troubleshooting |
| `data/sra/logs/<RUN>.fasterq-dump.log` | SRA Toolkit `fasterq-dump` | Log of SRA-to-FASTQ conversion | Plain text. Looks like conversion status and number of spots/reads written | Troubleshooting FASTQ conversion | Useful for methods/troubleshooting |
| `data/sra/fastq/<RUN>_1.fastq.gz` | SRA Toolkit + `pigz` | Compressed FASTQ file for read 1 of a paired-end sample | Gzipped text. Uncompressed FASTQ has 4 lines/read: `@read_id`, DNA sequence, `+`, quality string. Peek: `zcat data/sra/fastq/SRR18681855_1.fastq.gz \| head -n 8`. FASTQ format reference: [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) | Input for FastQC and Bowtie2; contains raw sequencing reads and quality scores | Yes, raw input data |
| `data/sra/fastq/<RUN>_2.fastq.gz` | SRA Toolkit + `pigz` | Compressed FASTQ file for read 2 of a paired-end sample | Same FASTQ structure as read 1. The reads correspond to the opposite end of each DNA fragment | Input for FastQC and Bowtie2; paired with `<RUN>_1.fastq.gz` | Yes, raw input data |
| `results/fastqc/raw/<RUN>_1_fastqc.html` | FastQC | Human-readable QC report for read 1 | HTML file. Open in a browser; inside it contains sections like `Per base sequence quality`, `Per sequence GC content`, and `Adapter Content` | Lets you inspect read quality, GC content, adapter contamination, duplication, and overrepresented sequences | Yes, report/figure source |
| `results/fastqc/raw/<RUN>_2_fastqc.html` | FastQC | Human-readable QC report for read 2 | Same HTML report structure as read 1, but for mate/read 2 | Same as read 1, but for the second read in each pair | Yes, report/figure source |
| `results/fastqc/raw/<RUN>_1_fastqc.zip` | FastQC | Machine-readable FastQC result archive for read 1 | Zip archive. Unzip to see files like `fastqc_data.txt`, `summary.txt`, and image files. Peek: `unzip -l results/fastqc/raw/SRR18681855_1_fastqc.zip` | Used by MultiQC to build a combined report | Yes, MultiQC input |
| `results/fastqc/raw/<RUN>_2_fastqc.zip` | FastQC | Machine-readable FastQC result archive for read 2 | Same zip structure as read 1. `fastqc_data.txt` contains the module-level values MultiQC parses | Used by MultiQC to build a combined report | Yes, MultiQC input |
| `results/multiqc/raw/multiqc_report.html` | MultiQC | Combined QC report across all FastQC outputs | HTML file. Open in a browser; contains combined plots/tables across all samples | Main summary report for raw read quality across all red and black samples | Yes, final QC report |
| `references/carp_refseq/GCF_018340385.1_ASM1834038v1_genomic.fna.gz` | NCBI RefSeq download | Compressed reference genome FASTA | Gzipped FASTA. Uncompressed FASTA has headers beginning with `>` followed by genome sequence lines. Peek: `zcat references/carp_refseq/GCF_018340385.1_ASM1834038v1_genomic.fna.gz \| head` | Original downloaded reference genome archive | Keep as original reference source |
| `references/carp_refseq/GCF_018340385.1_ASM1834038v1_assembly_report.txt` | NCBI RefSeq download | NCBI assembly metadata and chromosome/scaffold naming report | Plain text table with assembly metadata, sequence names, GenBank/RefSeq accessions, and chromosome/scaffold labels | Documents the reference genome version and sequence names | Yes, useful for methods |
| `references/carp.fa` | `gunzip -c` from downloaded RefSeq FASTA | Uncompressed working copy of the carp reference genome | FASTA text. Example shape: `>NC_... Cyprinus carpio chromosome ...` followed by many lines of `A/C/G/T/N` bases | Input for Bowtie2 indexing, SAMtools indexing, Picard dictionary, and GATK | Yes, core reference |
| `references/carp.fa.fai` | SAMtools `faidx` | FASTA index showing sequence names, lengths, and byte offsets | Plain tab-delimited text. Columns are sequence name, length, offset, bases per line, and bytes per line | Required by many tools for fast random access to the reference genome | Yes, required downstream |
| `references/carp.dict` | Picard `CreateSequenceDictionary` | Sequence dictionary for the reference genome | SAM-style header text. Contains lines like `@HD` and `@SQ SN:<sequence> LN:<length> M5:<checksum>` | Required by Picard/GATK so BAM and VCF files match the reference contigs | Yes, required for GATK |
| `references/bowtie2/carp.*.bt2` | Bowtie2 `bowtie2-build` | Bowtie2 index files for the reference genome | Binary Bowtie2 index files; not meant to be read manually. They are consumed by `bowtie2 -x references/bowtie2/carp` | Used by Bowtie2 to align reads quickly to the carp genome | Yes, required for Bowtie2 |
| `results/alignment/<RUN>.bowtie2.log` | Bowtie2 | Alignment summary for one sample | Plain text. Looks like counts for paired reads, concordant alignments, discordant alignments, and final alignment rate | Reports total paired reads, concordant/discordant alignment, and overall alignment rate | Yes, QC/report source |
| `results/alignment/<RUN>.unsorted.bam` | Bowtie2 piped to SAMtools `view` | BAM alignment file in the order reads were produced, not coordinate-sorted | Binary BAM. View as SAM text with `bash scripts/run_in_env.sh samtools view results/alignment/SRR18681855.unsorted.bam \| head`. SAM/BAM reference: [SAM/BAM spec](https://samtools.github.io/hts-specs/SAMv1.pdf) | Temporary/intermediate file used to make the sorted BAM | Usually no after sorted BAM is verified |
| `results/alignment/<RUN>.sorted.bam` | SAMtools `sort` | Coordinate-sorted BAM file | Binary BAM, but sorted by reference coordinate. View header with `samtools view -H`; view alignments with `samtools view FILE \| head` | Main alignment file for Qualimap and Picard; reads are sorted by genomic position | Yes, key downstream input |
| `results/alignment/<RUN>.sorted.bam.bai` | SAMtools `index` | BAM index for the sorted BAM | Binary index file; not human-readable. Test it with `samtools idxstats results/alignment/SRR18681855.sorted.bam` | Allows tools to jump to regions of the BAM quickly | Yes, needed with sorted BAM |
| `results/alignment/<RUN>.flagstat.txt` | SAMtools `flagstat` | Simple alignment statistics | Plain text. Contains lines like `55669996 + 0 in total` and `51782889 + 0 mapped (93.02% : N/A)` | Reports total reads, mapped reads, primary mapped reads, and other mapping categories | Yes, QC/report source |
| `results/qualimap/<RUN>/qualimapReport.html` | Qualimap `bamqc` | Human-readable BAM QC report | HTML file. Open in a browser; includes plots/tables for coverage, GC content, insert size, and mapping quality | Lets you inspect coverage, mapping quality, insert size, GC bias, and mapping distribution | Yes, final QC report |
| `results/qualimap/<RUN>/genome_results.txt` | Qualimap `bamqc` | Text summary of BAM QC metrics | Plain text. Contains sections with number of reads, coverage, mapping quality, insert size, and genome fraction metrics | Useful for extracting numeric coverage/mapping metrics into tables or figures | Yes, report/figure source |
| `results/picard/<RUN>.rg.bam` | Picard `AddOrReplaceReadGroups` | BAM file with read-group metadata added | Binary BAM. Header includes read group line like `@RG ID:<RUN> SM:<RUN> PL:DNBSEQ`. View with `samtools view -H results/picard/SRR18681855.rg.bam` | Intermediate input for duplicate marking and GATK compatibility | Intermediate, but useful until dedup BAM is verified |
| `results/picard/<RUN>.dedup.bam` | Picard `MarkDuplicates` | BAM file where duplicate reads have been marked | Binary BAM. Duplicate reads have SAM flag `1024`; inspect with `samtools view` or count with `samtools flagstat` | Main GATK input for variant calling | Yes, key downstream input |
| `results/picard/<RUN>.dedup.bai` | Picard `MarkDuplicates` with `CREATE_INDEX=true` | Index for duplicate-marked BAM | Binary index file; not human-readable. Test with `samtools idxstats results/picard/SRR18681855.dedup.bam` | Required for GATK to access the BAM efficiently | Yes, needed with dedup BAM |
| `results/picard/<RUN>.dup_metrics.txt` | Picard `MarkDuplicates` | Duplicate metrics summary | Plain text. Contains a metrics table with columns such as `LIBRARY`, `UNPAIRED_READS_EXAMINED`, `READ_PAIRS_EXAMINED`, and `PERCENT_DUPLICATION` | Reports duplicate rate and library complexity information | Yes, QC/report source |
| `results/gatk_chrA20/plots/a20_windowed_fst_100kb.png` | `scripts/plot_a20_variants.R` | Windowed FST plot across chromosome A20 with `0.45` and `0.75` guide thresholds | PNG image. Open in a browser or image viewer | Shows candidate red-vs-black differentiation regions and labels the top peak with `NC_056591.1:14557124 T>A` | Yes, main result figure |
| `results/gatk_chrA20/plots/a20_fst_with_snp_density_100kb.png` | `scripts/plot_a20_variants.R` | Combined chromosome A20 plot with windowed FST and mirrored top/bottom SNP-density strips | PNG image. Open in a browser or image viewer | Lets you compare whether high-FST windows are supported by many SNPs or only a few informative SNPs | Yes, main result figure |
| `results/gatk_chrA20/plots/a20_variant_summary.png` | `scripts/plot_a20_variants.R` | Red-vs-black SNP category summary | PNG image. Open in a browser or image viewer | Shows whether A20 SNPs have higher alternate allele frequency in red, higher in black, similar frequencies, or missing one color group | Yes, main summary figure |
| `results/gatk_chrA20/plots/a20_fst_threshold_comparison.png` | `scripts/plot_a20_variants.R` | Compares strict `FST >= 0.75` and broader `FST >= 0.45` candidate thresholds | PNG image. Open in a browser or image viewer | Shows how many candidate windows, informative SNPs, and strong SNPs are recovered at each threshold | Yes, main result figure |
| `results/gatk_chrA20/tables/a20_variant_summary.tsv` | `scripts/plot_a20_variants.R` | Table behind the updated red-vs-black SNP category plot | Tab-delimited text with category, SNP count, and percent columns | Documents the data behind `a20_variant_summary.png` | Yes, figure source |
| `results/gatk_chrA20/tables/a20_dataset_summary.tsv` | `scripts/plot_a20_variants.R` | Basic dataset counts for A20 | Tab-delimited text with total PASS biallelic SNPs and sample counts | Preserves the old dataset-count summary separately from the red-vs-black biological summary | Useful supporting table |
| `results/gatk_chrA20/tables/a20_top_fst_window_informative_snps.tsv` | `scripts/plot_a20_variants.R` | The 3 informative SNPs that produced the top `FST = 0.75` window | Tab-delimited text with position, alleles, red/black called alleles, and red/black allele frequencies | Identifies the specific SNPs behind the high FST peak | Yes, main result table |
| `results/gatk_chrA20/tables/a20_fst_threshold_summary.tsv` | `scripts/plot_a20_variants.R` | Summary of candidate yield at `FST >= 0.75` vs `FST >= 0.45` | Tab-delimited text with threshold, number of windows, informative SNPs, and strong SNPs | Supports the threshold comparison plot | Yes, figure source |
| `results/gatk_chrA20/tables/a20_candidate_regions_fst_ge_0_45.tsv` | `scripts/plot_a20_variants.R` | Candidate A20 regions found using the broader `FST >= 0.45` threshold | Tab-delimited text with region, FST, SNP count, nearby annotation, and interpretation | Connects candidate windows back to possible biological meaning for red vs black coloration | Yes, interpretation table |
| `results/gatk_chrA20/tables/a20_candidate_snps_fst_ge_0_45_abs_af_diff_ge_0_75.tsv` | `scripts/plot_a20_variants.R` | Strong candidate SNPs from broader high-FST windows | Tab-delimited text with SNP position, ref/alt alleles, red/black allele frequencies, region annotation, and interpretation | Lists candidate SNP markers with large red-vs-black allele-frequency differences | Yes, main result table |
| `.DS_Store` | macOS Finder | macOS metadata file | Binary macOS folder-view metadata; not part of the analysis | Not part of the analysis | No, ignore |

## Which Files Matter Most Going Forward

For GATK variant calling, the important files are:

```text
references/carp.fa
references/carp.fa.fai
references/carp.dict
results/picard/<RUN>.dedup.bam
results/picard/<RUN>.dedup.bai
```

For reporting your methods and quality control, the important files are:

```text
results/fastqc/raw/*_fastqc.html
results/multiqc/raw/multiqc_report.html
results/alignment/*.bowtie2.log
results/alignment/*.flagstat.txt
results/qualimap/<RUN>/qualimapReport.html
results/picard/*.dup_metrics.txt
```

For saving disk space after verifying downstream files, the most disposable files are:

```text
results/alignment/*.unsorted.bam
```

Only remove them after confirming that the matching `.sorted.bam`, `.sorted.bam.bai`, and later `.dedup.bam` files are complete.
