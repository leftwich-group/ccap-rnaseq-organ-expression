#!/usr/bin/env Rscript
# Developmental timepoint confirmation (embryo QC)
# Output: 04_results/11_developmental_confirmation/
# Run from the repository root:  Rscript 03_analysis/scripts/11_developmental_confirmation.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir <- here("04_results/11_developmental_confirmation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load raw counts & metadata
expr <- read_tsv(
  here("02_preprocessing/05_results_summary/short_read_genome/salmon.merged.gene_counts.tsv"),
  show_col_types = FALSE) %>%
  column_to_rownames(var = names(.)[1]) %>%
  as.matrix()
mode(expr) <- "numeric"

meta_dev <- read_csv(here("01_data/metadata.csv"), show_col_types = FALSE) %>%
  filter(tissue == "EG", experiment == "2") %>%
  filter(id %in% colnames(expr))

expr_dev <- expr[, meta_dev$id, drop = FALSE]
expr_dev <- expr_dev[rowSums(expr_dev) > 0, ]
expr_log  <- log2(expr_dev + 1)

# Top 5,000 variable genes for PCA
top5k <- names(sort(rowVars(expr_log), decreasing = TRUE))[1:min(5000, nrow(expr_log))]
pca   <- prcomp(t(expr_log[top5k, ]), scale. = TRUE)
var_e <- pca$sdev^2 / sum(pca$sdev^2)

meta_dev$age <- factor(meta_dev$age,
                       levels = sort(unique(meta_dev$age)))
age_cols <- setNames(okabe_ito_palette(nlevels(meta_dev$age)),
                     levels(meta_dev$age))

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("id") %>%
  left_join(meta_dev, by = "id")

p_pca_dev <- ggplot(pca_df, aes(PC1, PC2, colour = factor(age), label = id)) +
  geom_point(size = 4) +
  geom_text_repel(size = 2.5, max.overlaps = 15) +
  scale_colour_manual(values = age_cols) +
  labs(x = paste0("PC1 (", round(var_e[1]*100,1), "%)"),
       y = paste0("PC2 (", round(var_e[2]*100,1), "%)"),
       colour = "Timepoint (h)",
       title = "Embryo developmental timepoints - PCA") +
  plot_theme

p_pca_dev
ggsave(file.path(outdir, "PCA_embryo_timepoints.pdf"), p_pca_dev,
       width = 8, height = 6)

# Sample correlation heatmap
cor_mat <- cor(expr_log[top5k, ], method = "spearman")
ann_row  <- data.frame(age = meta_dev$age, row.names = meta_dev$id)

         
pdf(
  file.path(outdir, "sample_correlation_heatmap.pdf"),
  width = 9,
  height = 8
)
pheatmap(cor_mat,
         annotation_row = ann_row,
         annotation_col = ann_row,
         annotation_colors = list(age = age_cols),
         clustering_method = "ward.D2",
         color    = colorRampPalette(c("white", "#0072B2"))(100),
         fontsize = 8,
)

dev.off()
# Mean expression trajectory for top 20 variable genes
gene_var_dev <- rowVars(expr_log)
top20_genes  <- names(sort(gene_var_dev, decreasing = TRUE))[1:50]

traj_df <- expr_log[top20_genes, ] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "id", values_to = "log2cpm") %>%
  left_join(meta_dev, by = "id")

p_traj <- ggplot(traj_df, aes(as.numeric(as.character(age)), log2cpm)) +
  ggdist::stat_lineribbon(.width = c(0.5, 0.8, 0.95),
                          fill = "#0072B2", colour = "#0072B2", alpha = 0.25) +
  scale_x_continuous(breaks = sort(unique(as.numeric(as.character(traj_df$age))))) +
  labs(x = "Timepoint (h)", y = "log2(counts + 1)") +
  plot_theme +
  theme(legend.position = "none")

p_traj
ggsave(file.path(outdir, "top20_trajectory.pdf"), p_traj,
       width = 8, height = 5)

message("Developmental confirmation complete.")
