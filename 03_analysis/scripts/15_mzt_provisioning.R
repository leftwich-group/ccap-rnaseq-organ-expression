#!/usr/bin/env Rscript
# Q1c - Unfertilised-egg transcripts: library composition, sensitivity of the
# unfertilised-versus-fertilised contrast to normalisation, and overlap with
# genes high in ovary.
# Output: 04_results/15_MZT_provisioning/
# Run from the repository root:  Rscript 03_analysis/scripts/15_mzt_provisioning.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir <- file.path(outdir_raw, RES$mzt_provision)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

norm_methods <- c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
top_n_genes  <- 100

# Quasi-likelihood test of one contrast between named groups.
de_table <- function(counts, group, contrast, norm_method = "TMM") {
  group  <- factor(group)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  dge  <- DGEList(counts = counts)
  keep <- rowSums(dge$counts >= min_count) >= min_samples
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- if (norm_method == "none") dge else calcNormFactors(dge, method = norm_method)
  dge  <- estimateDisp(dge, design, robust = TRUE)
  fit  <- glmQLFit(dge, design, robust = TRUE)
  qlf  <- glmQLFTest(fit, contrast = makeContrasts(contrasts = contrast, levels = design))
  topTags(qlf, n = Inf)$table %>% rownames_to_column("gene")
}
n_dir <- function(res) tibble(n_higher_first = sum(res$FDR < fdr_thresh & res$logFC >=  lfc_thresh),
                              n_higher_second = sum(res$FDR < fdr_thresh & res$logFC <= -lfc_thresh),
                              median_logFC = median(res$logFC),
                              n_tested = nrow(res))

meta_eg <- meta_raw %>%
  filter(tissue == "EG") %>%
  mutate(age = as.numeric(age),
         stage = case_when(age == 0 ~ "Unfertilised", age <= 4 ~ "Early", TRUE ~ "Late"))
counts_eg <- counts_short_raw[, meta_eg$id]

# ---- Library composition --------------------------------------------------
dge_eg <- DGEList(counts_eg)
dge_eg <- dge_eg[rowSums(dge_eg$counts >= min_count) >= min_samples, , keep.lib.sizes = FALSE]
dge_eg <- calcNormFactors(dge_eg)
composition <- meta_eg %>%
  mutate(top_gene_share = apply(counts_eg, 2, function(x) sum(sort(x, decreasing = TRUE)[seq_len(top_n_genes)]) / sum(x)),
         tmm_factor     = setNames(dge_eg$samples$norm.factors, rownames(dge_eg$samples))[id],
         library_size   = colSums(counts_eg))
write_tsv(composition %>% dplyr::select(id, age, stage, library_size, top_gene_share, tmm_factor),
          file.path(outdir, "library_composition_by_sample.tsv"))

p_comp <- composition %>%
  pivot_longer(c(top_gene_share, tmm_factor), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, top_gene_share = paste0("Share of counts in top ", top_n_genes, " genes"),
                         tmm_factor = "TMM normalisation factor")) %>%
  ggplot(aes(factor(age), value, colour = stage)) +
  geom_point(size = 2.4) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = c(Unfertilised = "#E69F00", Early = "#56B4E9", Late = "#D55E00")) +
  labs(x = "Hours post-oviposition", y = NULL, colour = NULL) +
  plot_theme
ggsave(file.path(outdir, "library_composition_by_timepoint.pdf"), p_comp, width = 8.5, height = 4)

# ---- Sensitivity of unfertilised versus fertilised to normalisation -------
contrast_sets <- list(
  unfert_vs_1h   = list(meta = meta_eg %>% filter(experiment == "1"),
                        group = function(m) ifelse(m$age == 0, "Unfertilised", "Fertilised"),
                        contrast = "Unfertilised - Fertilised"),
  unfert_vs_early = list(meta = meta_eg %>% filter(stage != "Late"),
                         group = function(m) m$stage,
                         contrast = "Unfertilised - Early"))

sens <- list(); tmm_res <- list()
for (cs in names(contrast_sets)) {
  m <- contrast_sets[[cs]]$meta
  for (nm in norm_methods) {
    res <- de_table(counts_short_raw[, m$id], contrast_sets[[cs]]$group(m),
                    contrast_sets[[cs]]$contrast, nm)
    sens[[paste(cs, nm)]] <- n_dir(res) %>% mutate(contrast = cs, normalisation = nm, .before = 1)
    if (nm == "TMM") tmm_res[[cs]] <- res
    write_tsv(res, file.path(outdir, paste0(cs, "_", nm, "_edgeR.tsv")))
  }
}
sens <- bind_rows(sens) %>%
  dplyr::rename(n_higher_in_unfertilised = n_higher_first, n_higher_in_fertilised = n_higher_second)
write_tsv(sens, file.path(outdir, "normalisation_sensitivity.tsv"))
print(as.data.frame(sens), digits = 3)

p_norm <- sens %>%
  filter(contrast == "unfert_vs_1h") %>%
  pivot_longer(c(n_higher_in_unfertilised, n_higher_in_fertilised), names_to = "direction", values_to = "n_genes") %>%
  mutate(normalisation = factor(normalisation, levels = norm_methods),
         direction = factor(direction, levels = c("n_higher_in_unfertilised", "n_higher_in_fertilised"),
                            labels = c("Higher in unfertilised", "Higher in fertilised"))) %>%
  ggplot(aes(normalisation, n_genes, fill = direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("Higher in unfertilised" = "#E69F00", "Higher in fertilised" = "#56B4E9")) +
  labs(x = "Normalisation", y = "Genes (FDR < 0.05, at least 2-fold)", fill = NULL) +
  plot_theme + theme(legend.position = "bottom")
ggsave(file.path(outdir, "normalisation_sensitivity.pdf"), p_norm, width = 6.5, height = 4)

# ---- Similarity of unfertilised eggs to adult tissues ---------------------
# Group-mean log-CPM of every adult tissue and the 0 h eggs, normalised together.
meta_sim <- meta_raw %>%
  filter(tissue %in% c("TE", "AG", "AB", "YO", "OO") | (tissue == "EG" & age == 0)) %>%
  mutate(group = case_when(tissue == "EG" ~ "Unfertilised egg",
                           tissue == "TE" ~ "Testes",
                           tissue == "AG" ~ "Accessory gland",
                           tissue == "YO" ~ "Distal ovary",
                           tissue == "OO" ~ "Proximal ovary",
                           tissue == "AB" & sex == "M" ~ "Abdomen (male)",
                           TRUE ~ "Abdomen (female)"))
dge_sim <- DGEList(counts_short_raw[, meta_sim$id])
dge_sim <- dge_sim[rowSums(cpm(dge_sim) >= 1) >= min_samples, , keep.lib.sizes = FALSE]
dge_sim <- calcNormFactors(dge_sim)
lcpm_sim  <- cpm(dge_sim, log = TRUE, prior.count = 2)
group_mean <- sapply(split(meta_sim$id, meta_sim$group), function(ids) rowMeans(lcpm_sim[, ids, drop = FALSE]))
similarity <- tibble(group = setdiff(colnames(group_mean), "Unfertilised egg")) %>%
  mutate(spearman = sapply(group, function(g) cor(group_mean[, "Unfertilised egg"], group_mean[, g], method = "spearman")),
         pearson  = sapply(group, function(g) cor(group_mean[, "Unfertilised egg"], group_mean[, g])),
         n_genes  = nrow(group_mean)) %>%
  arrange(desc(spearman))
write_tsv(similarity, file.path(outdir, "unfertilised_egg_similarity_to_tissues.tsv"))
print(as.data.frame(similarity), digits = 3)

p_sim <- ggplot(similarity, aes(spearman, reorder(group, spearman))) +
  geom_col(fill = "#E69F00", width = 0.65) +
  labs(x = "Spearman correlation of mean log-CPM with unfertilised eggs", y = NULL) +
  plot_theme
ggsave(file.path(outdir, "unfertilised_egg_similarity_to_tissues.pdf"), p_sim, width = 6.5, height = 3.8)

# ---- Overlap with genes higher in mature than in younger oocytes ----------
meta_oo <- meta_raw %>% filter(mating_status == "V", tissue %in% c("OO", "YO"))
oo   <- de_table(counts_short_raw[, meta_oo$id], meta_oo$tissue, "OO - YO")
unfert <- tmm_res$unfert_vs_1h %>% dplyr::select(gene, logFC_unfert = logFC, FDR_unfert = FDR)
mrg  <- inner_join(unfert, oo %>% dplyr::select(gene, logFC_oo = logFC, FDR_oo = FDR), by = "gene") %>%
  mutate(unfert_high = FDR_unfert < fdr_thresh & logFC_unfert >= lfc_thresh,
         oo_high     = FDR_oo < fdr_thresh & logFC_oo >= lfc_thresh)
ft <- fisher.test(table(mrg$unfert_high, mrg$oo_high))
overlap <- tibble(n_universe = nrow(mrg), n_unfert_high = sum(mrg$unfert_high),
                  n_oo_high = sum(mrg$oo_high), n_overlap = sum(mrg$unfert_high & mrg$oo_high),
                  n_expected = sum(mrg$unfert_high) * sum(mrg$oo_high) / nrow(mrg),
                  odds_ratio = unname(ft$estimate), p_value = ft$p.value,
                  spearman_logFC = cor(mrg$logFC_unfert, mrg$logFC_oo, method = "spearman"))
write_tsv(overlap, file.path(outdir, "unfertilised_high_vs_mature_oocyte_overlap.tsv"))
print(as.data.frame(overlap), digits = 3)

saveRDS(list(composition = p_comp, normalisation = p_norm, similarity = p_sim),
        file.path(outdir, "plots.rds"))
message("Done: ", outdir)
