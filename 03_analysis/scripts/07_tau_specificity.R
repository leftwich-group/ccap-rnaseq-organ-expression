#!/usr/bin/env Rscript
# Tissue-specificity scores (tau index)
# Output: 04_results/07_Tau_specificity/
# Run from the repository root:  Rscript 03_analysis/scripts/07_tau_specificity.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

# `lcpm_adult` was a shared in-memory object in the monolithic document;
# rebuild it here so this script stands alone.
lcpm_adult <- build_adult_lcpm(counts_short_raw, meta_raw)

outdir_tau <- file.path(outdir_raw, "07_Tau_specificity")
dir.create(outdir_tau, recursive = TRUE, showWarnings = FALSE)

# Compute mean CPM per tissue (adult only)
mean_cpm_per_tissue <- lapply(names(tissues), function(tc) {
  ids <- meta_raw$id[meta_raw$tissue == tc]
  ids <- intersect(ids, colnames(lcpm_adult))
  if (length(ids) == 0) return(NULL)
  rowMeans(2^lcpm_adult[, ids, drop = FALSE] - 1)  # back to CPM
}) %>%
  setNames(names(tissues)) %>%
  Filter(Negate(is.null), .)

# Keep only tissues with data and genes present in all tissues
common_genes <- Reduce(intersect, lapply(mean_cpm_per_tissue, names))
tau_mat      <- do.call(cbind,
                        lapply(mean_cpm_per_tissue, function(x) x[common_genes]))
colnames(tau_mat) <- names(mean_cpm_per_tissue)

tau_scores <- calc_tau(tau_mat)

tau_df <- tibble(gene = names(tau_scores), tau = tau_scores) %>%
  filter(!is.na(tau)) %>%
  arrange(desc(tau))

write_tsv(tau_df, file.path(outdir_tau, "tau_scores.tsv"))

# Distribution plot
p_tau_dist <- ggplot(tau_df, aes(tau)) +
  geom_histogram(bins = 50, fill = "#0072B2", colour = "white") +
  geom_vline(xintercept = 0.85, linetype = "dashed",
             colour = "#D55E00") +
  annotate("text", x = 0.87, y = Inf, hjust = 0, vjust = 1.5,
           label = "tau > 0.85\n(tissue-specific)", size = 3.5,
           colour = "#D55E00") +
  labs(title = "Tissue-specificity distribution - tau index",
       x = "Tau (0 = ubiquitous, 1 = tissue-specific)",
       y = "Gene count") +
  plot_theme
p_tau_dist
ggsave(file.path(outdir_tau, "tau_distribution.pdf"), p_tau_dist,
       width = 7, height = 5)

# Top 20 tissue-specific genes per tissue
for (tc in colnames(tau_mat)) {
  tissue_specific <- tau_df %>%
    filter(tau > 0.85) %>%
    arrange(desc(tau_mat[gene, tc])) %>%
    slice_head(n = 20)
  message(tc, ": ", nrow(tissue_specific),
          " high-tau genes (top 20 shown)")
}

# Heatmap of top 200 tissue-specific genes
top200 <- tau_df %>% slice_head(n = 200) %>% pull(gene)
mat_tau <- log2(tau_mat[top200, ] + 1)

ann_col_tau <- data.frame(
  tissue = colnames(mat_tau), row.names = colnames(mat_tau))

pheatmap::pheatmap(t(mat_tau),
         annotation_row = ann_col_tau,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         color    = colorRampPalette(c("white", "#D55E00"))(100),
         filename = file.path(outdir_tau, "tau_top200_heatmap.pdf"),
         width = 12, height = 5)
