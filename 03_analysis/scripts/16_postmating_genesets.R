#!/usr/bin/env Rscript
# Q3b - Post-mating immune response as gene-set tests
# Tests whether immune-annotated genes shift as a set between mated and virgin
# females in abdomen, proximal ovary and distal ovary, using all tested genes
# (no significance threshold).
# Input:  04_results/03_PostMating/contrast_{1,2,3}/short_*_edgeR.tsv
# Output: 04_results/16_PostMating_genesets/
# Run from the repository root:  Rscript 03_analysis/scripts/16_postmating_genesets.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir <- file.path(outdir_raw, RES$mating_sets)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tissue_files <- c(
  Abdomen        = res_dir(RES$postmating, "contrast_1", "short_mated_vs_virgin_abdomen_edgeR.tsv"),
  Proximal_Ovary = res_dir(RES$postmating, "contrast_2", "short_mated_vs_virgin_OO_edgeR.tsv"),
  Distal_Ovary   = res_dir(RES$postmating, "contrast_3", "short_mated_vs_virgin_YO_edgeR.tsv"))
if (!all(file.exists(tissue_files))) stop(
  "Missing: ", paste(basename(tissue_files[!file.exists(tissue_files)]), collapse = ", "),
  "\n  Run 03_analysis/scripts/03_postmating.R first.", call. = FALSE)

min_set_size <- 5

# ---- Immune gene sets -----------------------------------------------------
immune_term_names <- c(
  "immune response", "innate immune response", "defense response to other organism",
  "Toll signaling pathway", "antimicrobial humoral response",
  "defense response to Gram-positive bacterium", "defense response to Gram-negative bacterium",
  "peptidoglycan recognition protein signaling pathway")

go_bp <- AnnotationDbi::select(GO.db, keys = keys(GO.db, "GOID"),
                               columns = c("TERM", "ONTOLOGY"), keytype = "GOID") %>%
  filter(ONTOLOGY == "BP", TERM %in% immune_term_names)
if (nrow(go_bp) != length(immune_term_names)) stop(
  "Immune GO term resolution incomplete: expected ", length(immune_term_names),
  ", matched ", nrow(go_bp), ".", call. = FALSE)

offspring <- as.list(GOBPOFFSPRING)
genes_with <- function(go_ids) names(geneID2GO_short)[
  vapply(geneID2GO_short, function(t) any(t %in% go_ids), logical(1))]

term_sets <- setNames(
  lapply(go_bp$GOID, function(id) genes_with(c(id, offspring[[id]]))),
  go_bp$TERM)
gene_sets <- c(list("All immune terms" = unique(unlist(term_sets))), term_sets)
write_tsv(tibble(set = names(gene_sets), n_genes = lengths(gene_sets)),
          file.path(outdir, "gene_set_sizes.tsv"))

# ---- Gene-set tests -------------------------------------------------------
# Signed square root of the quasi-likelihood F statistic is the per-gene
# statistic (positive = higher in mated); limma::geneSetTest uses ranks.
tests <- list(); de_genes <- list(); stats <- list()
for (tissue in names(tissue_files)) {
  res  <- read_tsv(tissue_files[[tissue]], show_col_types = FALSE) %>% mutate(id = clean_ids(gene))
  stat <- sign(res$logFC) * sqrt(res$F)
  stats[[tissue]] <- list(stat = stat, res = res)
  for (set in names(gene_sets)) {
    idx <- which(res$id %in% gene_sets[[set]])
    if (length(idx) < min_set_size) next
    tests[[paste(tissue, set)]] <- tibble(
      tissue = tissue, set = set, n_genes_tested = length(idx),
      mean_stat_set = mean(stat[idx]), mean_stat_other = mean(stat[-idx]),
      p_higher_in_mated = geneSetTest(idx, stat, alternative = "up"),
      p_higher_in_virgin = geneSetTest(idx, stat, alternative = "down"),
      p_either_direction = geneSetTest(idx, stat, alternative = "mixed"))
  }
  de_genes[[tissue]] <- res %>%
    filter(FDR < fdr_thresh, abs(logFC) >= lfc_thresh) %>%
    transmute(tissue = tissue, gene, logFC, FDR,
              higher_in = if_else(logFC > 0, "mated", "virgin"),
              in_immune_set = id %in% gene_sets[["All immune terms"]])
}
tests <- bind_rows(tests) %>%
  mutate(across(starts_with("p_"), ~ p.adjust(.x, "BH"), .names = "FDR_{.col}"))
write_tsv(tests, file.path(outdir, "immune_geneset_tests.tsv"))
write_tsv(bind_rows(de_genes), file.path(outdir, "mating_de_genes.tsv"))
print(as.data.frame(tests %>% filter(set == "All immune terms")), digits = 3)

# ---- Distribution of the statistic, immune genes against all other genes --
dist_df <- bind_rows(lapply(names(tissue_files), function(tissue) {
  d <- stats[[tissue]]
  tibble(tissue = tissue, stat = d$stat,
         group = if_else(d$res$id %in% gene_sets[["All immune terms"]], "Immune genes", "Other genes"))
})) %>% mutate(tissue = factor(gsub("_", " ", tissue), levels = gsub("_", " ", names(tissue_files))))
p_dist <- ggplot(dist_df, aes(stat, colour = group, fill = group)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_density(alpha = 0.25, adjust = 1.2) +
  coord_cartesian(xlim = c(-5, 5)) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(`Immune genes` = "#D55E00", `Other genes` = "grey40")) +
  scale_fill_manual(values = c(`Immune genes` = "#D55E00", `Other genes` = "grey40")) +
  labs(x = "Signed statistic (positive = higher in mated)", y = "Density", colour = NULL, fill = NULL) +
  plot_theme + theme(legend.position = "bottom")
ggsave(file.path(outdir, "immune_statistic_distribution.pdf"), p_dist, width = 11, height = 4)

# ---- Mean shift per set ---------------------------------------------------
p_shift <- tests %>%
  mutate(set = factor(set, levels = rev(names(gene_sets))),
         tissue = factor(gsub("_", " ", tissue), levels = gsub("_", " ", names(tissue_files))),
         sig = FDR_p_either_direction < fdr_thresh) %>%
  ggplot(aes(mean_stat_set - mean_stat_other, set, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_point(size = 2.6) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#D55E00"),
                      labels = c("FDR >= 0.05", "FDR < 0.05")) +
  labs(x = "Mean shift of the set relative to other genes (positive = higher in mated)", y = NULL, colour = NULL) +
  plot_theme + theme(legend.position = "bottom")
ggsave(file.path(outdir, "immune_set_shift.pdf"), p_shift, width = 11, height = 4.2)

saveRDS(list(shift = p_shift, distribution = p_dist), file.path(outdir, "plots.rds"))
message("Done: ", outdir)
