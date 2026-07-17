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
fwrite(variant_summary, file.path(outdir, "..", "tables", "a20_dataset_summary.tsv"), sep = "\t")

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

site[, red_black_category := fifelse(
  is.na(red_af) | is.na(black_af),
  "Missing one color group",
  fifelse(
    af_diff >= 0.25,
    "Higher alternate allele frequency in red",
    fifelse(
      af_diff <= -0.25,
      "Higher alternate allele frequency in black",
      "Similar red and black allele frequency"
    )
  )
)]

fwrite(site, file.path(outdir, "..", "tables", "a20_site_allele_frequencies.tsv"), sep = "\t")

fst_windows <- site[!is.na(fst_den), .(
  n_snps = .N,
  mean_abs_af_diff = mean(abs_af_diff, na.rm = TRUE),
  fst = sum(fst_num, na.rm = TRUE) / sum(fst_den, na.rm = TRUE)
), by = .(CHROM, window_start, window_end)]
fst_windows[, fst := pmax(fst, 0)]
fwrite(fst_windows, file.path(outdir, "..", "fst", "a20_windowed_fst_100kb.tsv"), sep = "\t")

variant_categories <- site[, .N, by = red_black_category]
setnames(variant_categories, c("red_black_category", "N"), c("category", "n_snps"))
variant_categories[, percent := round(100 * n_snps / nrow(site), 2)]
variant_categories[, category := factor(category, levels = c(
  "Higher alternate allele frequency in red",
  "Higher alternate allele frequency in black",
  "Similar red and black allele frequency",
  "Missing one color group"
))]
setorder(variant_categories, category)
fwrite(variant_categories, file.path(outdir, "..", "tables", "a20_red_black_variant_categories.tsv"), sep = "\t")
fwrite(variant_categories, file.path(outdir, "..", "tables", "a20_variant_summary.tsv"), sep = "\t")

high_window <- fst_windows[which.max(fst)]
high_window_snps <- site[
  window_start == high_window$window_start &
    !is.na(fst_den),
  .(
    CHROM, POS, REF, ALT,
    red_called_alleles, black_called_alleles,
    red_af, black_af, af_diff, abs_af_diff,
    window_start, window_end
  )
][order(-abs_af_diff)]
fwrite(high_window_snps, file.path(outdir, "..", "tables", "a20_top_fst_window_informative_snps.tsv"), sep = "\t")

threshold_summary <- rbindlist(lapply(c(0.75, 0.45), function(threshold) {
  windows <- fst_windows[fst >= threshold]
  snps <- site[window_start %in% windows$window_start & !is.na(fst_den)]
  data.table(
    fst_threshold = threshold,
    candidate_windows = nrow(windows),
    informative_snps_in_candidate_windows = nrow(snps),
    strong_snps_abs_af_diff_ge_0_75 = snps[abs_af_diff >= 0.75, .N]
  )
}))
fwrite(threshold_summary, file.path(outdir, "..", "tables", "a20_fst_threshold_summary.tsv"), sep = "\t")

candidate_windows <- fst_windows[fst >= 0.45]
candidate_snps <- site[
  window_start %in% candidate_windows$window_start &
    !is.na(fst_den) &
    abs_af_diff >= 0.75,
  .(
    CHROM, POS, REF, ALT,
    red_called_alleles, black_called_alleles,
    red_af, black_af, af_diff, abs_af_diff,
    window_start, window_end
  )
][order(window_start, POS)]

region_annotations <- data.table(
  window_start = c(4400001L, 4500001L, 13900001L, 14500001L),
  candidate_region = c("A20:4.4-4.5 Mb", "A20:4.5-4.6 Mb", "A20:13.9-14.0 Mb", "A20:14.5-14.6 Mb"),
  nearby_annotation = c(
    "Within/near LOC109112681, annotated as AT-rich interactive domain-containing protein 1B-like",
    "Overlaps LOC109112681 and nearby pseudogene-rich sequence",
    "Near LOC109112719, uncharacterized protein-coding gene, and LOC122148973, cytochrome P450 2J3-like",
    "Near forkhead box protein Q1-like and F2-like genes; strongest SNP is near pseudogene LOC122149054"
  ),
  project_interpretation = c(
    "Candidate differentiation region with many informative SNPs; possible regulatory or linked variation, not a proven color mutation.",
    "Adjacent candidate region that may represent the same broader differentiated block as 4.4-4.5 Mb.",
    "Candidate region with more SNP support than the top 14.5-14.6 Mb peak; gene products are not known from this analysis to be pigmentation genes.",
    "Highest FST peak, but based on only 3 informative SNPs; strongest marker is NC_056591.1:14557124 T>A."
  )
)

candidate_regions <- merge(candidate_windows, region_annotations, by = "window_start", all.x = TRUE)
setcolorder(candidate_regions, c(
  "CHROM", "window_start", "window_end", "candidate_region",
  "n_snps", "mean_abs_af_diff", "fst",
  "nearby_annotation", "project_interpretation"
))
fwrite(candidate_regions, file.path(outdir, "..", "tables", "a20_candidate_regions_fst_ge_0_45.tsv"), sep = "\t")

candidate_snps <- merge(
  candidate_snps,
  candidate_regions[, .(window_start, candidate_region, nearby_annotation, project_interpretation)],
  by = "window_start",
  all.x = TRUE
)
setcolorder(candidate_snps, c(
  "candidate_region", "CHROM", "POS", "REF", "ALT",
  "red_called_alleles", "black_called_alleles",
  "red_af", "black_af", "af_diff", "abs_af_diff",
  "window_start", "window_end",
  "nearby_annotation", "project_interpretation"
))
fwrite(candidate_snps, file.path(outdir, "..", "tables", "a20_candidate_snps_fst_ge_0_45_abs_af_diff_ge_0_75.tsv"), sep = "\t")

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
    geom_hline(yintercept = 0.45, linetype = "dashed", color = "#4C78A8", linewidth = 0.4) +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "#C44E52", linewidth = 0.4) +
    geom_line(color = "#D55E00", linewidth = 0.4) +
    geom_point(aes(color = fst >= 0.45), size = 0.9, alpha = 0.8) +
    geom_point(data = high_window, aes(window_start / 1e6, fst), color = "#C44E52", size = 2.4) +
    annotate("text", x = high_window$window_start / 1e6 + 0.55, y = 0.69,
             label = "Top peak: FST=0.75\n3 informative SNPs\nstrongest: 14,557,124 T>A",
             hjust = 0, vjust = 0.5, size = 3.2, color = "#333333") +
    annotate("text", x = 0.25, y = 0.45, label = "FST >= 0.45 candidate windows", hjust = 0, vjust = -0.6, size = 3, color = "#4C78A8") +
    annotate("text", x = 0.25, y = 0.75, label = "Top-window threshold", hjust = 0, vjust = -0.6, size = 3, color = "#C44E52") +
    scale_color_manual(values = c(`FALSE` = "#D55E00", `TRUE` = "#C44E52"), guide = "none") +
    labs(x = "Position on A20 (Mb)", y = "Windowed FST", title = "Red vs Black Carp Differentiation Across A20") +
    coord_cartesian(ylim = c(0, 0.82), clip = "off") +
    theme_minimal(base_size = 12)
)
dev.off()

combined_windows <- merge(
  fst_windows,
  variant_windows[, .(CHROM, window_start, density_snps = n_snps)],
  by = c("CHROM", "window_start"),
  all.x = TRUE
)

png(file.path(outdir, "a20_fst_with_snp_density_100kb.png"), width = 1500, height = 820, res = 150)
print(
  ggplot(combined_windows, aes(window_start / 1e6)) +
    geom_tile(aes(y = -0.035, fill = density_snps), height = 0.055, width = 0.095) +
    geom_tile(aes(y = 0.895, fill = density_snps), height = 0.055, width = 0.095) +
    geom_hline(yintercept = 0, color = "#777777", linewidth = 0.25) +
    geom_hline(yintercept = 0.45, linetype = "dashed", color = "#4C78A8", linewidth = 0.4) +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "#C44E52", linewidth = 0.4) +
    geom_line(aes(y = fst), color = "#D55E00", linewidth = 0.45) +
    geom_point(aes(y = fst, color = fst >= 0.45), size = 0.9, alpha = 0.8) +
    geom_point(data = high_window, aes(window_start / 1e6, fst), color = "#C44E52", size = 2.4) +
    annotate("text", x = high_window$window_start / 1e6 + 0.55, y = 0.69,
             label = "Top peak: FST=0.75\n3 informative SNPs\nstrongest: 14,557,124 T>A",
             hjust = 0, vjust = 0.5, size = 3.1, color = "#333333") +
    annotate("text", x = 0, y = 0.945, label = "SNP density track", hjust = 0, size = 3.1, color = "#4C78A8") +
    annotate("text", x = 0, y = -0.075, label = "SNP density track", hjust = 0, size = 3.1, color = "#4C78A8") +
    scale_color_manual(values = c(`FALSE` = "#D55E00", `TRUE` = "#C44E52"), guide = "none") +
    scale_fill_gradient(
      low = "#DCE6F2",
      high = "#2F6FA3",
      name = "PASS biallelic\nSNPs per 100 kb"
    ) +
    scale_y_continuous(
      name = "Windowed FST",
      limits = c(-0.09, 0.96),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8),
      labels = c("0", "0.2", "0.4", "0.6", "0.8")
    ) +
    labs(
      x = "Position on A20 (Mb)",
      title = "A20 Red vs Black FST With SNP Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y.left = element_text(color = "#D55E00"),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.3, "cm")
    )
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

png(file.path(outdir, "a20_variant_summary.png"), width = 1400, height = 760, res = 150)
print(
  ggplot(variant_categories[!is.na(category)],
         aes(category, n_snps, fill = category)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = paste0(n_snps, " (", percent, "%)")), hjust = -0.05, size = 3.1) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 115000), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c(
      "Higher alternate allele frequency in red" = "#D55E00",
      "Higher alternate allele frequency in black" = "#222222",
      "Similar red and black allele frequency" = "#4C78A8",
      "Missing one color group" = "#9A9A9A"
    ), guide = "none") +
    labs(x = NULL, y = "Number of A20 SNPs", title = "Red vs Black A20 SNP Categories") +
    theme_minimal(base_size = 12)
)
dev.off()

threshold_long <- melt(
  threshold_summary,
  id.vars = "fst_threshold",
  measure.vars = c("candidate_windows", "informative_snps_in_candidate_windows", "strong_snps_abs_af_diff_ge_0_75"),
  variable.name = "metric",
  value.name = "count"
)
threshold_long[, metric := factor(metric, levels = c(
  "candidate_windows",
  "informative_snps_in_candidate_windows",
  "strong_snps_abs_af_diff_ge_0_75"
), labels = c(
  "Candidate windows",
  "Informative SNPs",
  "Strong SNPs (AF diff >= 0.75)"
))]

png(file.path(outdir, "a20_fst_threshold_comparison.png"), width = 1300, height = 760, res = 150)
print(
  ggplot(threshold_long, aes(factor(fst_threshold), count, fill = metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    geom_text(aes(label = count), position = position_dodge(width = 0.75), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("#C44E52", "#4C78A8", "#009E73")) +
    labs(
      x = "Window FST threshold",
      y = "Count",
      fill = NULL,
      title = "Candidate SNP Yield at Strict vs Broader FST Thresholds"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
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
