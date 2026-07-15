#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/plot_a20_variants.R <genotypes.tsv> <samples.tsv> <outdir>")
}

geno_path <- args[[1]]
sample_path <- args[[2]]
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "..", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "..", "fst"), recursive = TRUE, showWarnings = FALSE)

samples <- fread(sample_path)
setnames(samples, names(samples), tolower(names(samples)))
samples <- samples[, .(run_accession, sample_name, color)]

geno <- fread(geno_path)
fixed_cols <- c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER")
sample_cols <- setdiff(names(geno), fixed_cols)

if (!all(sample_cols %in% samples$run_accession)) {
  missing <- setdiff(sample_cols, samples$run_accession)
  stop("Sample IDs missing from metadata: ", paste(missing, collapse = ", "))
}

red_samples <- samples[color == "red", run_accession]
black_samples <- samples[color == "black", run_accession]

gt_to_dosage <- function(x) {
  x <- gsub("\\|", "/", x)
  x[x %in% c(".", "./.", ".|.")] <- NA_character_
  parts <- tstrsplit(x, "/", fixed = TRUE)
  a <- suppressWarnings(as.integer(parts[[1]]))
  b <- suppressWarnings(as.integer(parts[[2]]))
  a + b
}

dosage <- as.data.table(lapply(geno[, ..sample_cols], gt_to_dosage))
setnames(dosage, sample_cols)

variant_summary <- data.table(
  category = c("PASS biallelic SNPs", "samples", "red samples", "black samples"),
  value = c(nrow(geno), length(sample_cols), length(red_samples), length(black_samples))
)
fwrite(variant_summary, file.path(outdir, "..", "tables", "a20_variant_summary.tsv"), sep = "\t")

window_size <- 100000L
variant_windows <- geno[, .(n_snps = .N), by = .(CHROM, window_start = ((POS - 1L) %/% window_size) * window_size + 1L)]
variant_windows[, window_end := window_start + window_size - 1L]
fwrite(variant_windows, file.path(outdir, "..", "tables", "a20_snp_density_100kb.tsv"), sep = "\t")

group_stats <- function(mat, cols) {
  sub <- as.matrix(mat[, ..cols])
  called_alleles <- rowSums(!is.na(sub)) * 2
  alt_alleles <- rowSums(sub, na.rm = TRUE)
  p <- alt_alleles / called_alleles
  p[called_alleles == 0] <- NA_real_
  data.table(called_alleles = called_alleles, alt_alleles = alt_alleles, af = p)
}

red <- group_stats(dosage, red_samples)
black <- group_stats(dosage, black_samples)

site <- data.table(
  CHROM = geno$CHROM,
  POS = geno$POS,
  REF = geno$REF,
  ALT = geno$ALT,
  red_called_alleles = red$called_alleles,
  black_called_alleles = black$called_alleles,
  red_af = red$af,
  black_af = black$af
)
site[, af_diff := red_af - black_af]
site[, abs_af_diff := abs(af_diff)]

red_pi <- 2 * site$red_af * (1 - site$red_af)
black_pi <- 2 * site$black_af * (1 - site$black_af)
between_pi <- site$red_af * (1 - site$black_af) + site$black_af * (1 - site$red_af)
site[, fst_num := between_pi - ((red_pi + black_pi) / 2)]
site[, fst_den := between_pi]
site[is.na(fst_den) | fst_den <= 0, `:=`(fst_num = NA_real_, fst_den = NA_real_)]
site[, window_start := ((POS - 1L) %/% window_size) * window_size + 1L]
site[, window_end := window_start + window_size - 1L]

fwrite(site, file.path(outdir, "..", "tables", "a20_site_allele_frequencies.tsv"), sep = "\t")

fst_windows <- site[!is.na(fst_den), .(
  n_snps = .N,
  mean_abs_af_diff = mean(abs_af_diff, na.rm = TRUE),
  fst = sum(fst_num, na.rm = TRUE) / sum(fst_den, na.rm = TRUE)
), by = .(CHROM, window_start, window_end)]
fst_windows[, fst := pmax(fst, 0)]
fwrite(fst_windows, file.path(outdir, "..", "fst", "a20_windowed_fst_100kb.tsv"), sep = "\t")

png(file.path(outdir, "a20_snp_density_100kb.png"), width = 1400, height = 700, res = 150)
print(
  ggplot(variant_windows, aes(window_start / 1e6, n_snps)) +
    geom_col(fill = "#4C78A8", width = 0.09) +
    labs(x = "Position on A20 (Mb)", y = "PASS biallelic SNPs per 100 kb", title = "SNP Density Across A20 in Red and Black Carp Samples") +
    theme_minimal(base_size = 12)
)
dev.off()

png(file.path(outdir, "a20_windowed_fst_100kb.png"), width = 1400, height = 700, res = 150)
print(
  ggplot(fst_windows, aes(window_start / 1e6, fst)) +
    geom_line(color = "#D55E00", linewidth = 0.4) +
    geom_point(color = "#D55E00", size = 0.8, alpha = 0.75) +
    labs(x = "Position on A20 (Mb)", y = "Windowed FST", title = "Red vs Black Carp Differentiation Across A20") +
    theme_minimal(base_size = 12)
)
dev.off()

png(file.path(outdir, "a20_windowed_af_difference_100kb.png"), width = 1400, height = 700, res = 150)
print(
  ggplot(fst_windows, aes(window_start / 1e6, mean_abs_af_diff)) +
    geom_line(color = "#009E73", linewidth = 0.4) +
    geom_point(color = "#009E73", size = 0.8, alpha = 0.75) +
    labs(x = "Position on A20 (Mb)", y = "Mean |Red AF - Black AF|", title = "Windowed Red vs Black Allele Frequency Difference Across A20") +
    theme_minimal(base_size = 12)
)
dev.off()

top_diff <- site[order(-abs_af_diff)][1:min(.N, 5000)]
png(file.path(outdir, "a20_top_af_differences.png"), width = 1400, height = 700, res = 150)
print(
  ggplot(top_diff, aes(POS / 1e6, abs_af_diff)) +
    geom_point(color = "#009E73", alpha = 0.7, size = 0.8) +
    labs(x = "Position on A20 (Mb)", y = "|Red AF - Black AF|", title = "Largest Red vs Black Allele Frequency Differences on A20") +
    theme_minimal(base_size = 12)
)
dev.off()

png(file.path(outdir, "a20_variant_summary.png"), width = 900, height = 650, res = 150)
print(
  ggplot(variant_summary[category %in% c("PASS biallelic SNPs", "samples", "red samples", "black samples")],
         aes(reorder(category, value), value)) +
    geom_col(fill = "#7A5195") +
    coord_flip() +
    labs(x = NULL, y = "Count", title = "A20 Red vs Black Variant Dataset Summary") +
    theme_minimal(base_size = 12)
)
dev.off()

pca_mat <- as.matrix(dosage[, ..sample_cols])
pca_mat <- t(pca_mat)
col_means <- colMeans(pca_mat, na.rm = TRUE)
for (j in seq_len(ncol(pca_mat))) {
  pca_mat[is.na(pca_mat[, j]), j] <- col_means[j]
}
keep <- apply(pca_mat, 2, sd) > 0
pca_mat <- pca_mat[, keep, drop = FALSE]
if (ncol(pca_mat) > 50000) {
  set.seed(1)
  pca_mat <- pca_mat[, sample(seq_len(ncol(pca_mat)), 50000), drop = FALSE]
}
pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
pca_df <- data.table(run_accession = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_df <- merge(pca_df, samples, by = "run_accession", all.x = TRUE)
fwrite(pca_df, file.path(outdir, "..", "tables", "a20_pca_samples.tsv"), sep = "\t")

var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
png(file.path(outdir, "a20_pca_red_black.png"), width = 900, height = 700, res = 150)
print(
  ggplot(pca_df, aes(PC1, PC2, color = color, label = sample_name)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    scale_color_manual(values = c(black = "#222222", red = "#D55E00")) +
    labs(
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      title = "PCA of A20 Genotypes Colored by Red vs Black Carp Phenotype"
    ) +
    theme_minimal(base_size = 12)
)
dev.off()

message("Wrote plots to: ", normalizePath(outdir))
