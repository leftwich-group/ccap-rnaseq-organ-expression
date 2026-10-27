#!/usr/bin/env Rscript
# Global tissue atlas - PCA / UMAP
# Output: 04_results/06_Global_atlas/
# Run from the repository root:  Rscript 03_analysis/scripts/06_global_atlas.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir_atlas <- file.path(outdir_raw, "06_Global_atlas")
dir.create(outdir_atlas, recursive = TRUE, showWarnings = FALSE)

# Use adult tissues only (exclude embryos for this comparison)
adult_tissues <- c("AG", "AB", "TE", "YO", "OO")
meta_adult    <- meta_raw %>%
  filter(tissue %in% adult_tissues) %>%
  mutate(
    tissue_label = recode(tissue, !!!setNames(as.character(tissues[adult_tissues]),
                                              adult_tissues)),
    group_label  = paste0(tissue_label,
                          if_else(tissue == "AB",
                                  paste0("_", sex,
                                         if_else(!is.na(mating_status),
                                                 paste0("_", mating_status), "")),
                                  ""))
  )

counts_adult <- counts_short_raw[, meta_adult$id, drop = FALSE]

dge_adult  <- DGEList(counts = counts_adult)
keep_adult <- rowSums(cpm(dge_adult) >= 1) >= 2
dge_adult  <- dge_adult[keep_adult, , keep.lib.sizes = FALSE]
dge_adult  <- calcNormFactors(dge_adult)
lcpm_adult <- cpm(dge_adult, log = TRUE, prior.count = 2)

# ---- PCA (top 5,000 variable genes) ---------------------------------------
top5k_a <- names(sort(rowVars(lcpm_adult), decreasing = TRUE))[1:5000]
pca_a   <- prcomp(t(lcpm_adult[top5k_a, ]), scale. = TRUE)
vexp_a  <- pca_a$sdev^2 / sum(pca_a$sdev^2)

pca_df_a <- as.data.frame(pca_a$x[, 1:3]) %>%
  rownames_to_column("id") %>%
  left_join(meta_adult, by = "id")

# Scree plot
scree_df <- tibble(PC = 1:10, var_exp = vexp_a[1:10] * 100)
p_scree  <- ggplot(scree_df, aes(PC, var_exp)) +
  geom_col(fill = "#0072B2") +
  geom_line(colour = "#D55E00") +
  geom_point(colour = "#D55E00") +
  labs(title = "Global atlas - PCA scree", x = "PC", y = "% variance explained") +
  plot_theme
p_scree 
ggsave(file.path(outdir_atlas, "PCA_scree.pdf"), p_scree, width = 6, height = 4)

# PC1 vs PC2 coloured by tissue
tissue_pal <- setNames(okabe_ito_palette(length(unique(meta_adult$tissue_label))),
                       unique(meta_adult$tissue_label))

p_pca_tissue <- ggplot(pca_df_a, aes(PC1, PC2,
                                      colour = tissue_label,
                                      shape  = sex)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_colour_manual(values = tissue_pal) +
  guides(shape = guide_legend(override.aes = list(colour = "grey40"))) +
  labs(x = paste0("PC1 (", round(vexp_a[1]*100,1), "%)"),
       y = paste0("PC2 (", round(vexp_a[2]*100,1), "%)"),
       colour = "Tissue", shape = "Sex",
       title = "Global tissue atlas - PCA") +
  plot_theme
p_pca_tissue
ggsave(file.path(outdir_atlas, "PCA_tissue_sex.pdf"), p_pca_tissue,
       width = 8, height = 6)

# PC1 vs PC2 coloured by sex (to compare)
p_pca_sex <- ggplot(pca_df_a, aes(PC1, PC2, colour = sex, shape = tissue_label)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_colour_manual(values = annotation_colours$sex) +
  labs(title = "Global atlas - sex overlay") +
  plot_theme
p_pca_sex
ggsave(file.path(outdir_atlas, "PCA_sex.pdf"), p_pca_sex, width = 8, height = 6)

# ---- UMAP -----------------------------------------------------------------
set.seed(42)
umap_out <- umap(t(lcpm_adult[top5k_a, ]),
                 n_neighbors = min(15, nrow(meta_adult) - 1),
                 min_dist = 0.3)

umap_df <- as.data.frame(umap_out$layout) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  rownames_to_column("id") %>%
  left_join(meta_adult, by = "id")

p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2,
                               colour = tissue_label, shape = sex)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_colour_manual(values = tissue_pal) +
  labs(title = "Global tissue atlas - UMAP",
       colour = "Tissue", shape = "Sex") +
  plot_theme
p_umap
ggsave(file.path(outdir_atlas, "UMAP_tissue_sex.pdf"), p_umap, width = 8, height = 6)
