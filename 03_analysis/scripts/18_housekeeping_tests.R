#!/usr/bin/env Rscript
# Housekeeping orthologue tests against constitutive expression
# Tests whether Drosophila housekeeping orthologues are enriched among genes
# constitutive in each tissue and in all five adult tissues, and whether they
# vary less across tissues than other genes.
# Input:  04_results/05_Tissue_landscapes/Tissues/<tissue>/short_gene_landscape_summary.tsv
#         04_results/05_Tissue_landscapes/housekeeping/Ccap_FlyAtlas2_HK_mapping.tsv
# Output: 04_results/18_Housekeeping_tests/
# Run from the repository root:  Rscript 03_analysis/scripts/18_housekeeping_tests.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

tissue_dirs <- c("Accessory_Gland", "Abdomen", "Testes", "Distal_Ovary", "Proximal_Ovary")
land_files  <- res_dir(RES$landscapes, "Tissues", tissue_dirs, "short_gene_landscape_summary.tsv")
hk_file     <- res_dir(RES$landscapes, "housekeeping", "Ccap_FlyAtlas2_HK_mapping.tsv")
needed <- c(land_files, hk_file)
if (!all(file.exists(needed))) stop(
  "Missing: ", paste(needed[!file.exists(needed)], collapse = ", "),
  "\n  Run 03_analysis/scripts/05_tissue_landscapes.R first.", call. = FALSE)

outdir <- file.path(outdir_raw, RES$hk_tests)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hk_genes <- clean_ids(read_tsv(hk_file, show_col_types = FALSE)$gene_id)
land <- bind_rows(lapply(seq_along(tissue_dirs), function(i)
  read_tsv(land_files[i], show_col_types = FALSE) %>%
    mutate(tissue = tissue_dirs[i], gene = clean_ids(gene)))) %>%
  mutate(is_hk = gene %in% hk_genes)

fisher_row <- function(test, universe, df, flag) {
  tab <- table(factor(df$is_hk, c(FALSE, TRUE)), factor(df[[flag]], c(FALSE, TRUE)))
  ft  <- fisher.test(tab)
  tibble(test = test, universe = universe, n_universe = nrow(df),
         n_hk = sum(df$is_hk), n_hk_constitutive = sum(df$is_hk & df[[flag]]),
         odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
         odds_ratio_haldane = prod(diag(tab) + 0.5) / prod(c(tab[1, 2], tab[2, 1]) + 0.5),
         p_value = ft$p.value)
}

# ---- Enrichment among constitutive genes ----------------------------------
per_tissue <- bind_rows(lapply(tissue_dirs, function(tc)
  fisher_row("HK orthologue x constitutive", tc, land %>% filter(tissue == tc), "is_constitutive")))

all_tissue <- land %>%
  group_by(gene) %>%
  summarise(n_tissues = n(), n_const = sum(is_constitutive, na.rm = TRUE),
            mean_cv = mean(cv, na.rm = TRUE), is_hk = any(is_hk), .groups = "drop") %>%
  filter(n_tissues == length(tissue_dirs)) %>%
  mutate(constitutive_all = n_const == length(tissue_dirs))
all_label <- "all five tissues"
all_row <- fisher_row("HK orthologue x constitutive in all five tissues",
                      all_label, all_tissue, "constitutive_all")

# ---- Variation across tissues ---------------------------------------------
wilcox_row <- function(universe, df) {
  w <- wilcox.test(mean_cv ~ is_hk, data = df)
  tibble(test = "Mean cross-tissue CV, HK orthologue vs other", universe = universe,
         n_universe = nrow(df), n_hk = sum(df$is_hk),
         median_cv_hk = median(df$mean_cv[df$is_hk]), median_cv_other = median(df$mean_cv[!df$is_hk]),
         p_value = w$p.value)
}
cv_any <- land %>% group_by(gene) %>%
  summarise(mean_cv = mean(cv, na.rm = TRUE), is_hk = any(is_hk), .groups = "drop") %>%
  filter(is.finite(mean_cv))
cv_rows <- bind_rows(wilcox_row("all expressed genes", cv_any),
                     wilcox_row("genes expressed in all five tissues", all_tissue %>% filter(is.finite(mean_cv))))

tests <- bind_rows(per_tissue, all_row, cv_rows)
write_tsv(tests, file.path(outdir, "housekeeping_tests.tsv"))
print(as.data.frame(tests), digits = 3)

row_levels <- rev(gsub("_", " ", c(tissue_dirs, all_label)))
p <- bind_rows(per_tissue, all_row) %>%
  mutate(universe = factor(gsub("_", " ", universe), levels = row_levels)) %>%
  ggplot(aes(odds_ratio_haldane, universe)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(colour = "#0072B2", size = 3) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100)) +
  labs(x = "Odds ratio, housekeeping orthologues in constitutive set (log scale)", y = NULL) +
  plot_theme
ggsave(file.path(outdir, "housekeeping_enrichment_odds_ratios.pdf"), p, width = 8, height = 4)
saveRDS(list(odds_ratios = p), file.path(outdir, "plots.rds"))
message("Done: ", outdir)
