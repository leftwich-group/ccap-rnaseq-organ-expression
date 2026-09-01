#!/usr/bin/env Rscript
# Cross-tissue germline-soma correlation
# Output: 04_results/10_CrossTissue_correlation/
# Run from the repository root:  Rscript 03_analysis/scripts/10_crosstissue_correlation.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

# `lcpm_adult` was a shared in-memory object in the monolithic document;
# rebuild it here so this script stands alone.
lcpm_adult <- build_adult_lcpm(counts_short_raw, meta_raw)

outdir_cor <- file.path(outdir_raw, "10_CrossTissue_correlation")
dir.create(outdir_cor, recursive = TRUE, showWarnings = FALSE)

# Mean logCPM per tissue (all samples, short genome)
get_tissue_mean <- function(tc) {
  ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
  ids <- intersect(ids, colnames(lcpm_adult))
  if (length(ids) == 0) return(NULL)
  rowMeans(lcpm_adult[, ids, drop = FALSE])
}

tissue_means <- lapply(names(tissues), get_tissue_mean)
names(tissue_means) <- names(tissues)
tissue_means <- Filter(Negate(is.null), tissue_means)

# Restrict to common genes
cg <- Reduce(intersect, lapply(tissue_means, names))
mean_mat <- do.call(cbind, lapply(tissue_means, function(x) x[cg]))
colnames(mean_mat) <- tissues[colnames(mean_mat)]

# Pairwise Spearman correlation
cor_mat <- cor(mean_mat, method = "spearman")

# Heatmap
pdf(file.path(outdir_cor, "tissue_spearman_correlation.pdf"),
    width = 7, height = 6)
pheatmap(cor_mat,
         color    = colorRampPalette(c("#0072B2", "white", "#D55E00"))(100),
         breaks   = seq(-1, 1, length.out = 101),
         display_numbers = TRUE,
         number_format   = "%.2f",
         fontsize_number = 10,
         clustering_method = "ward.D2",
         main = "Tissue–tissue Spearman correlation (mean logCPM)")
dev.off()

# Scatter: testes vs distal ovary (germline comparison)
if ("TE" %in% names(tissue_means) && "YO" %in% names(tissue_means)) {
  scatter_df <- tibble(
    gene         = cg,
    Testes       = mean_mat[, "Testes"],
    Distal_Ovary = mean_mat[, "Distal_Ovary"]
  ) %>%
    mutate(diff = abs(Testes - Distal_Ovary),
           label_gene = if_else(rank(-diff) <= 20, gene, ""))

  p_germ_scatter <- ggplot(scatter_df, aes(Testes, Distal_Ovary)) +
    geom_point(alpha = 0.3, size = 1, colour = "grey60") +
    geom_smooth(method = "lm", se = FALSE,
                colour = "#D55E00") +
    geom_abline(linetype = "dashed") +
    geom_text_repel(aes(label = label_gene), size = 2.5, max.overlaps = 20) +
    labs(title = "Germline comparison: testes vs distal ovary",
         subtitle = paste0("Spearman ρ = ",
                           round(cor_mat["Testes","Distal_Ovary"], 3)),
         x = "Testes (mean log₂CPM)",
         y = "Distal ovary (mean log₂CPM)") +
    plot_theme
  p_germ_scatter
  ggsave(file.path(outdir_cor, "scatter_testes_vs_distal_ovary.pdf"),
         p_germ_scatter, width = 7, height = 6)
}
