#!/usr/bin/env Rscript
# Q1 - Maternal-to-zygotic transition (MZT)
# Output: 04_results/01_MZT/
# Run from the repository root:  Rscript 03_analysis/scripts/01_mzt.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir <- file.path(outdir_raw, "01_MZT/contrast_1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>% filter(tissue == "EG", experiment == "1")
counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long  <- counts_long_raw[,  meta$id, drop = FALSE]

design <- model.matrix(~ 0 + egg_status, data = meta)
colnames(design) <- sub("egg_status", "", colnames(design))

contrast_list <- list(unfert_vs_fert = "unfertilised - fertilised")

  #positive_label <- "unfertilised"
  #negative_label <- "fertilised"
  
results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "egg_status",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
  long  = run_reference_analysis(counts_long,  meta, design, contrast_list,
                                 "long",  outdir, "egg_status",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)

# Concordance scatter
comparison_all <- lapply(names(contrast_list), function(cn)
  compare_references_per_contrast(results_all$short, results_all$long,
                                  cn, outdir, gene_map_file, fdr_thresh))

saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))

for (g in c("short","long"))
  append_de_stats(results_all, "Q1", "Fertilised vs Unfertilised Egg",
                  contrast_list[[1]], genome = g)

for (genome in c("short", "long")) {
  geneID2GO  <- get(paste0("geneID2GO_",  genome))
  gene_univ  <- get(paste0("gene_universe_", genome))
  res        <- results_all[[genome]]$contrasts$unfert_vs_fert$table

  up_genes   <- res %>% filter(FDR < fdr_thresh, logFC  >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_genes <- res %>% filter(FDR < fdr_thresh, logFC  < -lfc_thresh) %>% pull(gene) %>% clean_ids()

  go_up   <- run_topgo_all_ont(up_genes,   gene_univ, geneID2GO)
  go_down <- run_topgo_all_ont(down_genes, gene_univ, geneID2GO)

  save_go_panel(go_up, go_down,
                outdir, paste0(genome, "_unfert_vs_fert"))
  
  # outputGO 
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <- go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)
  
  go_merged$sub_contrast <- "unfert_vs_fert"
  write_tsv(go_merged, paste0(outdir,"/",genome,"GO_terms_merged.tsv"))
}

outdir_c2 <- file.path(outdir_raw, "01_MZT/contrast_2")
dir.create(outdir_c2, recursive = TRUE, showWarnings = FALSE)

# All egg samples (unfertilised = age 0; fertilised = ages 1–24)
meta_all_egg <- meta_raw %>% filter(tissue == "EG")
counts_short_egg <- counts_short_raw[, meta_all_egg$id, drop = FALSE]

dge_egg <- DGEList(counts = counts_short_egg)
keep_egg <- rowSums(dge_egg$counts >= min_count) >= min_samples
dge_egg  <- dge_egg[keep_egg, , keep.lib.sizes = FALSE]
dge_egg  <- calcNormFactors(dge_egg)
logcpm_egg <- cpm(dge_egg, log = TRUE, prior.count = 2)

# ---- Hierarchical clustering of fertilised samples to define early/late ---
meta_fert <- meta_all_egg %>% filter(age > 0)
mat_fert  <- logcpm_egg[, meta_fert$id, drop = FALSE]
top_var   <- names(sort(rowVars(mat_fert), decreasing = TRUE))[1:2000]

hc <- hclust(dist(t(mat_fert[top_var, ])), method = "ward.D2")

cluster_assign      <- cutree(hc, k = 2)
meta_fert$cluster   <- paste0("C", cluster_assign[meta_fert$id])

# Label clusters by mean age
cluster_mean_age <- meta_fert %>%
  group_by(cluster) %>%
  summarise(mean_age = mean(as.numeric(age)), .groups = "drop")
early_cluster <- cluster_mean_age$cluster[which.min(cluster_mean_age$mean_age)]

meta_fert <- meta_fert %>%
  mutate(time_group = if_else(cluster == early_cluster, "early", "late"))

meta_grouped <- meta_all_egg %>%
  left_join(meta_fert %>% dplyr::select(id, time_group), by = "id") %>%
  mutate(time_group = coalesce(time_group,
                               if_else(age == 0, "unfertilised", NA_character_)))

# ggplot dendrogram
dend_data  <- dendro_data(hc, type = "rectangle")
labels_df  <- dend_data$labels %>%
  left_join(meta_fert %>% dplyr::select(id, age), by = c("label" = "id"))

age_cols_fert <- setNames(
  okabe_ito_palette(length(unique(labels_df$age))),
  sort(unique(labels_df$age)))

p_dend <- ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x, y, xend = xend, yend = yend)) +
  geom_text(data = labels_df,
            aes(x, y - 0.5, label = label, colour = factor(age)),
            angle = 90, hjust = 1, size = 2.8) +
  scale_colour_manual(values = age_cols_fert, name = "Timepoint") +
  theme_minimal(base_size = 11) +
  labs(title = "Hierarchical clustering - fertilised egg samples") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_dend
ggsave(file.path(outdir_c2, "fertilised_hclust.pdf"), p_dend,
       width = 9, height = 7)

# ---- Three-group contrasts (unfert | early | late) -------------------------
counts_grp <- counts_short_raw[, meta_grouped$id[!is.na(meta_grouped$time_group)],
                                 drop = FALSE]
meta_grp   <- meta_grouped %>% filter(!is.na(time_group))

design_grp <- model.matrix(~ 0 + time_group, data = meta_grp)
colnames(design_grp) <- sub("time_group", "", colnames(design_grp))

contrast_list_grp <- list(
  early_vs_unfert = "early - unfertilised",
  late_vs_unfert  = "late  - unfertilised",
  late_vs_early   = "late  - early"
)

results_grp <- list(
  short = run_reference_analysis(counts_grp, meta_grp, design_grp,
                                 contrast_list_grp, "short", outdir_c2,
                                 "time_group", min_count, min_samples,
                                 fdr_thresh, lfc_thresh)
)

saveRDS(results_grp, file.path(outdir_c2, "results_all_nested.rds"))

# Record stats for each of the three temporal contrasts (short genome only)
contrast_list_grp1 <- list(
  "Unfertilised vs Early" = "early_vs_unfert",
  "Unfertilised vs Late"  = "late_vs_unfert",
  "Early vs Late"   = "late_vs_early"  
)

for (cname in names(contrast_list_grp1)) {
  append_de_stats(results_grp, "Q1",
                  paste("MZT :", cname),
                  contrast_list_grp1[[cname]],
                  genome = "short")
}


go_list <- list()
for (cname in names(contrast_list_grp)) {
  res      <- results_grp$short$contrasts[[cname]]$table
  up_g     <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_g   <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
  go_up    <- run_topgo_all_ont(up_g,   gene_universe_short, geneID2GO_short)
  go_down  <- run_topgo_all_ont(down_g, gene_universe_short, geneID2GO_short)
  save_go_panel(go_up, go_down, outdir_c2, paste0("short_", cname))
  
   # outputGO 
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <- go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)

  go_merged$sub_contrast <- cname
  
  go_list[[cname]] <- go_merged
}
go_list_merged <- go_list %>% bind_rows()
  write_tsv(go_list_merged, paste0(outdir_c2,"/short_GO_terms_merged.tsv"))

outdir_c3 <- file.path(outdir_raw, "01_MZT/contrast_3")
dir.create(outdir_c3, recursive = TRUE, showWarnings = FALSE)

# Use fertilised egg samples only
meta_fert2   <- meta_raw %>% filter(tissue == "EG", age > 0)
logcpm_fert2 <- logcpm_egg[, meta_fert2$id, drop = FALSE]

# Select top 3,000 variable genes
top3k <- names(sort(rowVars(logcpm_fert2), decreasing = TRUE))[1:3000]
mat   <- logcpm_fert2[top3k, ]

# k-means with k = 3 (M1, M2, M3)
set.seed(42)
km <- kmeans(mat, centers = 3, nstart = 25, iter.max = 100)

module_df <- tibble(gene = rownames(mat), module = paste0("M", km$cluster))

# Mean expression per module per timepoint
mean_expr <- mat %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  left_join(module_df, by = "gene") %>%
  pivot_longer(cols = -c(gene, module), names_to = "id",
               values_to = "logcpm") %>%
  left_join(meta_fert2 %>% dplyr::select(id, age), by = "id")# %>%
  #group_by(module, age) %>%
  #summarise(mean_logcpm = mean(logcpm), .groups = "drop")

# Relabel modules by expression trajectory direction
# M ordered by mean age at max expression
module_order <- mean_expr %>%
  group_by(module) %>%
  filter(logcpm == max(logcpm)) %>%
  summarise(peak_age = mean(as.numeric(age)), .groups = "drop") %>%
  arrange(peak_age)

module_df <- module_df %>%
  mutate(module = factor(module, levels = module_order$module))

p_modules <- ggplot(mean_expr,
                    aes(factor(age), logcpm,
                        group = module, colour = module)) +
   geom_jitter(size = 1.5, alpha = 0.01, position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "line",  size = 0.8)+
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_colour_manual(values = okabe_ito_palette(3)) +
  labs(x = "Timepoint (h)", y = "log2CPM",
       colour = "Module") +
  plot_theme

p_modules
ggsave(file.path(outdir_c3, "MZT_temporal_modules.pdf"), p_modules,
       width = 8, height = 5)

write_tsv(module_df, file.path(outdir_c3, "gene_modules.tsv"))

# GO per module
for (mod in unique(module_df$module)) {
  mod_genes <- module_df %>% filter(module == mod) %>% pull(gene) %>% clean_ids()
  go_mod    <- run_topgo_all_ont(mod_genes, gene_universe_short, geneID2GO_short)
  p_mod     <- plot_topgo_panel(go_mod, paste("Module", mod),
                                letters = c("A", "B", "C"))
  if (!is.null(p_mod))
    ggsave(file.path(outdir_c3, paste0("GO_module_", mod, ".pdf")),
           p_mod, width = 8, height = 10)
  
       summary_cls <- plot_topgo_bar(go_mod)
       if (!is.null(summary_cls))
        ggsave(file.path(outdir_c3, paste0(genome, "_GO_Summary_", mod, ".pdf")),
               summary_cls, width = 7, height = 6)
      
      # save summary file 
       topgo_df <- go_mod$BP
  if (is.null(topgo_df) || nrow(topgo_df) == 0) return(NULL)
 
  df <- topgo_df %>%
    mutate(p_value     = parse_topgo_p(weight01_fisher),
           Significant = as.numeric(Significant)) %>%
    filter(!is.na(0.05), !is.na(Significant)) %>%
    arrange(p_value) %>%
    slice_head(n = 12) 
 
  if (nrow(df) == 0) return(NULL)
 
  # save summary file 
      write_tsv(df,file.path(outdir_c3, paste0(genome, "_GO_Summary_", mod, ".tsv")))
  
}

# ---------------------------------------------------------------------------
# Q1 SUMMARY: Maternal-to-Zygotic Transition
# Three-panel figure:
#   A. PCA of all egg timepoints (unfertilised + fertilised, coloured by age)
#   B. Temporal module trajectories (M1 / M2 / M3)
#   C. DE gene count bar: unfertilised vs fertilised (short genome)
# ---------------------------------------------------------------------------
outdir_summary <- file.path(outdir_raw, "01_MZT/summary_panels")
dir.create(outdir_summary, recursive = TRUE, showWarnings = FALSE)

# ---- A: PCA ----------------------------------------------------------------
meta_egg_all <- meta_raw %>% filter(tissue == "EG")
counts_egg   <- counts_short_raw[, meta_egg_all$id, drop = FALSE]
dge_egg_sum  <- DGEList(counts = counts_egg)
keep_egg_sum <- rowSums(cpm(dge_egg_sum) >= 1) >= 2
dge_egg_sum  <- dge_egg_sum[keep_egg_sum, , keep.lib.sizes = FALSE]
dge_egg_sum  <- calcNormFactors(dge_egg_sum)
lcpm_egg_sum <- cpm(dge_egg_sum, log = TRUE, prior.count = 2)

top2k_egg  <- names(sort(rowVars(lcpm_egg_sum), decreasing = TRUE))[1:2000]
pca_egg    <- prcomp(t(lcpm_egg_sum[top2k_egg, ]), scale. = TRUE)
vexp_egg   <- pca_egg$sdev^2 / sum(pca_egg$sdev^2)

pca_egg_df <- as.data.frame(pca_egg$x[, 1:2]) %>%
  rownames_to_column("id") %>%
  left_join(meta_egg_all, by = "id") %>%
  mutate(
    stage = case_when(
      age == 0  ~ "Unfertilised",
      age %in% c(1,2,3,4)  ~ "Early (1-4 h)",
      TRUE      ~ "Late (8-24 h)"
    ),
    stage = factor(stage, levels = c("Unfertilised","Early (1-4 h)","Late (8-24 h)"))
  )

stage_cols <- c("Unfertilised" = "#E69F00",
                "Early (1-4 h)" = "#56B4E9",
                "Late (8-24 h)"  = "#009E73")

pA <- ggplot(pca_egg_df, aes(PC1, PC2)) +

  # ---- 1. Convex hull shading per stage (works with n >= 1) ----------------
  geom_polygon(
    data = pca_egg_df %>%
      group_by(stage) %>%
      slice(chull(PC1, PC2)) %>%        # chull silently returns 1 point when n=1
      ungroup(),
    aes(fill = stage, group = stage),
    alpha = 0.12, colour = NA
  ) +

  # ---- 2. Trajectory: connect mean position of each timepoint in order -----
  # Compute per-timepoint centroid so the arrow tracks the average trajectory
  {
    traj <- pca_egg_df %>%
      group_by(age) %>%
      summarise(PC1 = mean(PC1), PC2 = mean(PC2),
                stage = dplyr::first(stage), .groups = "drop") %>%
      arrange(as.numeric(as.character(age)))

    list(
      geom_path(data = traj,
                aes(PC1, PC2, colour = stage, group = 1),
                linewidth = 0.7,  alpha = 0.6),
      geom_segment(
        data = traj %>% slice(n()-1, n()) %>%
          summarise(x    = PC1[1], y    = PC2[1],
                    xend = PC1[2], yend = PC2[2]),
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
        colour = "grey40", linewidth = 0.6
      )
    )
  } +

  # ---- 3. Points coloured by stage, shaped by age class --------------------
  geom_point(aes(colour = stage, shape = stage), size = 3.5, alpha = 0.9) +

  # ---- 4. Timepoint labels (age in hours) ----------------------------------
  geom_text_repel(
    aes(label = paste0(age, "h"), colour = stage),
    size = 2.8, show.legend = FALSE,
    box.padding = 0.3, max.overlaps = 20
  ) +

  scale_colour_manual(values = stage_cols) +
  scale_fill_manual(  values = stage_cols) +
  scale_shape_manual(values  = c(16, 17, 15)) +
  guides(fill = "none") +
  labs(x      = paste0("PC1 (", round(vexp_egg[1]*100, 1), "%)"),
       y      = paste0("PC2 (", round(vexp_egg[2]*100, 1), "%)"),
       colour = "Stage", shape = "Stage",
       title  = "A - Developmental stage separation") +
  plot_theme +
  theme(legend.position = "bottom")
ggsave(file.path(outdir_summary, "PCA.png"),
       pA, width = 5, height = 5)
# ---- B: Temporal modules ---------------------------------------------------
# Load from Q1 contrast_3 output (module_df written there)
mod_file <- file.path(outdir_raw, "01_MZT", "contrast_3", "gene_modules.tsv")

if (file.exists(mod_file)) {
  module_df_sum <- read_tsv(mod_file, show_col_types = FALSE)

  meta_fert_sum <- meta_raw %>%
    filter(tissue == "EG", egg_status == "fertilised") %>%
    arrange(as.numeric(age))

  mat_fert_sum <- lcpm_egg_sum[
    intersect(module_df_sum$gene, rownames(lcpm_egg_sum)),
    meta_fert_sum$id, drop = FALSE]

  mod_traj <- mat_fert_sum %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "id", values_to = "logcpm") %>%
    left_join(module_df_sum, by = "gene") %>%
    left_join(meta_fert_sum %>% dplyr::select(id, age), by = "id") %>%
    group_by(module, age) %>%
    summarise(mean_logcpm = mean(logcpm, na.rm = TRUE),
              se_logcpm   = sd(logcpm, na.rm = TRUE) /
                              sqrt(n()),
              .groups = "drop")

mod_traj1 <- mat_fert_sum %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "id", values_to = "logcpm") %>%
  left_join(module_df_sum, by = "gene") %>%
  left_join(meta_fert_sum %>% dplyr::select(id, age), by = "id")


# Module label lookup — annotate with dominant GO category from your results
mod_labels <- c(
  "M1" = "M1: Maternal decay\n(metabolism / cell cycle)",
  "M2" = "M2: Transient activation\n(RNA processing / ribosomes)",
  "M3" = "M3: Late activation\n(development / nervous system)"
)

pB <- ggplot( mod_traj1 ,
             aes(factor(age),logcpm,
                 colour = module, group = module)) +
  stat_summary(fun      = mean,geom     = "point", position = position_dodge(0.9), size     = 3, show.legend = FALSE)+
  stat_summary(fun.data = mean_cl_normal, geom     = "errorbar", width    = 0.3, linewidth = 0.8, position = position_dodge(0.9)) +
  stat_summary(fun = mean, geom     = "line", linewidth = 0.5, position = position_dodge(0.9)) +
  scale_colour_manual(values = okabe_ito_palette(3),
                      labels = mod_labels) +
  scale_fill_manual(  values = okabe_ito_palette(3)) +
  guides(fill = "none",
         colour = guide_legend(nrow = 3)) +
  labs(x = "Timepoint (h post-oviposition)",
       y = "Mean log2CPM",
       colour = NULL,
       title  = "B - Temporal expression modules") +
  plot_theme +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 12),
        legend.key.height = unit(0.4, "cm"))
pB
ggsave(file.path(outdir_summary, "Temporal_modules.png"),
       pB, width = 5, height = 5)
  } else {
  pB <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Run Q1 C3 first to generate\ngene_modules.tsv",
             size = 4, colour = "grey50") +
    theme_void() +
    labs(title = "B - Temporal expression modules")
}

# ---- C: DE gene count bar --------------------------------------------------
# Load all Q1 contrast summaries (short genome only)
q1_contrast_dirs <- c(
  "unfert_vs_fert"     = file.path(outdir_raw, "01_MZT/contrast_1"),
  "early_vs_unfert"    = file.path(outdir_raw, "01_MZT/contrast_2"),
  "late_vs_unfert"     = file.path(outdir_raw, "01_MZT/contrast_2"),
  "late_vs_early"      = file.path(outdir_raw, "01_MZT/contrast_2")
)

de_counts_q1 <- lapply(names(q1_contrast_dirs), function(label) {
  f <- file.path(q1_contrast_dirs[[label]], "results_all_nested.rds")
  if (!file.exists(f)) return(NULL)
  list <- readRDS(f) 
  df <- list$short$contrasts[[label]]$summary
  df %>%
    filter(reference == "short") %>%
    mutate(contrast = case_when(
      label == "unfert_vs_fert" ~ "Fertilised vs Unfertilised",
      label == "early_vs_unfert"  ~"Unfertilised vs Early ",
      label == "late_vs_unfert" ~"Unfertilised vs Late",
      label == "late_vs_early" ~"Early vs Late"
    ))
}) %>%
  bind_rows() %>%
  filter(!is.na(contrast))

if (nrow(de_counts_q1) > 0) {
  de_long_q1 <- de_counts_q1 %>%
    pivot_longer(c(n_up, n_down),
                 names_to  = "direction",
                 values_to = "count") %>%
    mutate(
      direction     = recode(direction,
                             n_up   = "Up regulated",
                             n_down = "Down regulated"),
      count_plot    = if_else(direction == "Down regulated",
                              -count, count)
    )

  pC <- ggplot(de_long_q1,
               aes(x = contrast, y = count_plot, fill = direction)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_richtext(
  aes(
    label = abs(count),
    y = count_plot + sign(count_plot) * 50
  ),
  fill = "white",
  label.color = "darkgrey",  
  size = 3
)+
    scale_fill_manual(values = c(
      "Up regulated"   = "#D55E00",
      "Down regulated" =  "#0072B2"
    )) +
       guides(fill="none")+
    scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
    coord_flip() +
    labs(x = NULL, y = "DE gene count",
         fill  = NULL,
         title = "C - Differentially expressed genes per contrast\n( FDR < 0.05, |log2FC| > 1)") +
    plot_theme +
    theme(legend.position = "bottom")
  
  ggsave(file.path(outdir_summary, "DE_counts.png"),
       pC, width = 5, height = 5)
  
} else {
  pC <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Run Q1 contrasts first to generate\nall_contrasts_summary.tsv",
             size = 4, colour = "grey50") +
    theme_void() + labs(title = "C - DE gene counts")
}

# ---- Assemble and save -----------------------------------------------------
p_q1_summary <- (pA | pB) / pC +
  plot_annotation(
    title    = "Maternal-to-Zygotic Transition",
    tag_levels = NULL,
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )


ggsave(file.path(outdir_summary, "Q1_MZT_summary.pdf"),
       p_q1_summary, width = 14, height = 14)
ggsave(file.path(outdir_summary, "Q1_MZT_summary.png"),
       p_q1_summary, width = 14, height = 14, dpi = 200)

p_q1_summary

outdir_go_summary <- file.path(outdir_raw, "01_MZT/summary_panels", "GO_heatmap")
dir.create(outdir_go_summary, recursive = TRUE, showWarnings = FALSE)

# ---- Load all short-genome GO summary files from Q1 -----------------------
q1_go_files <- list.files(
  path       = file.path(outdir_raw, "01_MZT"),
  pattern    = "GO_terms_merged\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
) %>%
  .[str_detect(., "short")]

message("Q1 GO files found: ", length(q1_go_files))
print(q1_go_files)
# should be 2 because C2 has all 3 merged 
if (length(q1_go_files) == 0) stop("No Q1 short GO files found - check paths")

# Read and label by contrast (parent folder name)
q1_go <- map_dfr(q1_go_files, function(f) {
  contrast <- basename(dirname(f))   # e.g. "contrast_1", "contrast_2"
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      contrast       = contrast,
      p_value        = parse_topgo_p(weight01_fisher),
      neg_log10p     = -log10(pmax(p_value, 1e-10)),
      contrast_group = paste0(sub_contrast," | ", group)
    )
  
  
})

message("Total GO rows loaded: ", nrow(q1_go))
message("Contrast-groups found:")
print(unique(q1_go$contrast_group))

# ============================================================
# Helper: make heatmap from long GO data
#   x-axis (columns) = contrast_group
#   y-axis (rows)    = GO Term
#   fill             = -log10(p)
# ============================================================
make_contrast_go_heatmap <- function(go_data,
                                     title_text,
                                     n_terms     = 30,
                                     min_groups  = 2,
                                     p_threshold = 0.05) {
  
  # Keep terms significant in at least min_groups contrast-groups
  sig_terms <- go_data %>%
    filter(p_value < p_threshold) %>%
    dplyr::count(Term) %>%
    filter(n >= min_groups) %>%
    pull(Term)
  
  # If nothing is shared, take top terms by mean enrichment instead
  if (length(sig_terms) == 0) {
    message("No shared terms across groups - taking top terms by mean enrichment")
    sig_terms <- go_data %>%
      filter(p_value < p_threshold) %>%
      group_by(Term) %>%
      summarise(mean_p = mean(neg_log10p), .groups = "drop") %>%
      arrange(desc(mean_p)) %>%
      slice_head(n = n_terms) %>%
      pull(Term)
  }
  
  if (length(sig_terms) == 0) {
    message("No significant GO terms found")
    return(NULL)
  }
  
  # Build matrix: best p-value per term × contrast_group
  mat_long <- go_data %>%
    filter(Term %in% sig_terms) %>%
    group_by(Term, contrast_group) %>%
    slice_min(p_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(Term, contrast_group, neg_log10p)
  
  mat_wide <- mat_long %>%
    pivot_wider(names_from  = contrast_group,
                values_from = neg_log10p,
                values_fill = 0)
  
  # Select top n_terms by mean enrichment across all groups
  top_terms <- mat_wide %>%
    mutate(mean_p = rowMeans(across(-Term))) %>%
    arrange(desc(mean_p)) %>%
    slice_head(n = n_terms) %>%
    pull(Term)
  
  mat <- mat_wide %>%
    filter(Term %in% top_terms) %>%
    column_to_rownames("Term") %>%
    as.matrix()
  
  mat[mat > 10] <- 10
  
  # Wrap long term names
  rownames(mat) <- str_trunc(rownames(mat), 50)
  
  # Row clustering
  if (nrow(mat) > 2) {
    row_order <- hclust(dist(mat), method = "ward.D2")$order
    mat       <- mat[row_order, , drop = FALSE]
  }
  
  pheatmap(
    mat,
    color             = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
    cluster_rows      = FALSE,
    cluster_cols      = TRUE,
    clustering_method = "ward.D2",
    fontsize_row      = 9,
    fontsize_col      = 9,
    angle_col         = "45",
    main              = title_text,
    legend_breaks     = c(0, 2, 4, 6, 8, 10),
    legend_labels     = c("0","2","4","6","8",">=10")
  )
}

# ============================================================
# PLOT 1: All Q1 contrasts in one heatmap
#   Columns = contrast | group (e.g. "contrast_1 | unfert")
#   Shows which GO terms are shared across MZT contrasts
#   and which are specific to one comparison
# ============================================================
p_q1_heatmap <- make_contrast_go_heatmap(
  go_data    = q1_go,
  title_text = "MZT: GO enrichment across all contrasts\n(-log10 p)",
  n_terms    = 30,
  min_groups = 2
)
p_q1_heatmap<- wrap_elements(
  full = grid.grabExpr(ComplexHeatmap::draw( p_q1_heatmap)))

if (!is.null(p_q1_heatmap)) {

  ggsave(file.path(outdir_go_summary, "Q1_GO_heatmap_all_contrasts.pdf"), p_q1_heatmap, width = 16, height = 10, dpi = 200)
  

}

# ============================================================
# PLOT 2: Faceted bar — top 8 terms per contrast-group
#   One panel per contrast-group
#   Good for seeing what each specific comparison is driven by
# ============================================================
p_q1_facet <- q1_go %>%
  filter(p_value < 0.05) %>%
  group_by(contrast_group) %>%
  slice_min(p_value, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Term_short = str_trunc(Term, 40)) %>%
  ggplot(aes(x      = reorder_within(Term_short, neg_log10p, contrast_group),
             y      = neg_log10p,
             fill   = neg_log10p)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ contrast_group, scales = "free_y", ncol = 2) +
  scale_x_reordered() +
  scale_fill_gradientn(colours = enrichment_palette(100),
                       name    = expression(-log[10](p))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "MZT: top GO terms per contrast group",
    subtitle = "Top 8 significant BP terms ( p < 0.05)",
    x = NULL, y = expression(-log[10](p))
  ) +
  plot_theme +
  theme(
    strip.text      = element_text(face = "bold", size = 9),
    axis.text.y     = element_text(size = 8),
    legend.position = "right"
  )

p_q1_facet
ggsave(file.path(outdir_go_summary, "Q1_GO_faceted_per_contrast.pdf"),
       p_q1_facet, width = 14, height = 0.6 * length(unique(q1_go$contrast_group)) + 4,
       limitsize = FALSE)
