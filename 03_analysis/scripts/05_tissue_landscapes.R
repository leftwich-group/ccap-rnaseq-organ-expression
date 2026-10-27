#!/usr/bin/env Rscript
# Tissue-specific transcriptional landscapes
# Output: 04_results/05_Tissue_landscapes/
# Run from the repository root:  Rscript 03_analysis/scripts/05_tissue_landscapes.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

# Shared build-landscape function used for every tissue section
build_tissue_landscape <- function(tissue_code, tissue_label,
                                   min_cpm = 1, min_samples = 2,
                                   n_mixtures = 3,
                                   outdir_raw ) {

  outdir <- file.path(outdir_raw, "05_Tissue_landscapes/Tissues", tissue_label)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  for (genome in c("short", "long")) {
    counts_raw <- if (genome == "short") counts_short_raw else counts_long_raw
    meta       <- meta_raw %>% filter(tissue == tissue_code)
    if (nrow(meta) == 0) { message("No samples for ", tissue_label, " (", genome, ")"); next }

    counts <- counts_raw[, meta$id, drop = FALSE]
    dge    <- DGEList(counts = counts)
    keep   <- rowSums(cpm(dge) >= min_cpm) >= min_samples
    dge    <- dge[keep, , keep.lib.sizes = FALSE]
    dge    <- calcNormFactors(dge)
    dge    <- estimateDisp(dge, model.matrix(~ 1, data = meta), robust = TRUE)
    lcpm   <- cpm(dge, log = TRUE, prior.count = 2)

    mean_cpm <- rowMeans(cpm(dge))

    # Expression class via Gaussian mixture
    fit_mix <- tryCatch(
      Mclust(log2(mean_cpm + 1), G = n_mixtures, verbose = FALSE),
      error = function(e) NULL
    )
    expr_class <- if (!is.null(fit_mix)) {
      cls <- fit_mix$classification
      mu  <- fit_mix$parameters$mean
      tibble(gene = names(cls),
             expr_class = case_when(
               cls == which.max(mu) ~ "high",
               cls == which.min(mu) ~ "low",
               TRUE                 ~ "medium"
             ))
    } else {
      tibble(gene = rownames(lcpm),
             expr_class = cut(mean_cpm,
                              breaks = quantile(mean_cpm, c(0, 1/3, 2/3, 1)),
                              labels = c("low","medium","high"),
                              include.lowest = TRUE))
    }

    # Overdispersion (log2 observed / fitted BCV)
    fitted_disp <- dge$trended.dispersion
    names(fitted_disp) <- rownames(dge)
    obs_disp    <- dge$tagwise.dispersion
    names(obs_disp) <- rownames(dge)
    overdisp <- log2((obs_disp + 0.001) / (fitted_disp + 0.001))

    # Constitutive: expressed in >= 80% of samples at >= 1 log-CPM
    is_const <- rowSums(lcpm >= 1) >= 0.8 * ncol(lcpm)

    cv <- apply(lcpm, 1, function(x) sd(x) / (mean(x) + 0.1))

    landscape <- as.data.frame(tibble(
      gene          = rownames(lcpm),
      mean_logcpm   = rowMeans(lcpm),
      cv            = cv,
      overdispersion = overdisp[rownames(lcpm)],
      is_constitutive = is_const[rownames(lcpm)],
      is_overdispersed = overdisp[rownames(lcpm)] > 0.5
    ) %>%
      left_join(expr_class, by = "gene"))
    
    landscape$expr_class <- factor(landscape$expr_class, levels = c("low", "medium", "high"))

    write_tsv(landscape,
              file.path(outdir, paste0(genome, "_gene_landscape_summary.tsv")))

    # ---- Expression class plot -----------------------------------------------
    p_class <- landscape %>%
      dplyr::count(expr_class) %>%
      ggplot(aes(expr_class, n, fill = expr_class)) +
      geom_col() +
      scale_fill_manual(values = okabe_ito_palette(3)) +
      labs(title = paste(tissue_label, genome, "- expression classes"),
           x = NULL, y = "Gene count") +
      plot_theme + theme(legend.position = "none")

    ggsave(file.path(outdir, paste0(genome, "_expression_classes.pdf")),
           p_class, width = 5, height = 4)

    # ---- GO per class --------------------------------------------------------
    for (cls_name in c("high", "medium")) {
      cls_genes <- landscape %>% filter(expr_class == cls_name) %>%
        pull(gene) %>% clean_ids()
      gene_univ <- get(paste0("gene_universe_", genome))
      geneID2GO <- get(paste0("geneID2GO_",     genome))
      go_cls    <- run_topgo_all_ont(cls_genes, gene_univ, geneID2GO)
      p_cls     <- plot_topgo_panel(go_cls, paste(tissue_label, cls_name))
      if (!is.null(p_cls))
        ggsave(file.path(outdir, paste0(genome, "_GO_", cls_name, ".pdf")),
               p_cls, width = 8, height = 10)
      
      summary_cls <- plot_topgo_bar(go_cls, paste(tissue_label, cls_name))
       if (!is.null(summary_cls))
        ggsave(file.path(outdir, paste0(genome, "_GO_Summary_", cls_name, ".png")),
               summary_cls, width = 7, height = 6)
      saveRDS(  summary_cls, file.path(outdir, paste0(genome, "_GO_Summary_", cls_name, ".RDS")))
      # save summary file 
       topgo_df <- go_cls$BP
  if (is.null(topgo_df) || nrow(topgo_df) == 0) return(NULL)
 
  df <- topgo_df %>%
    mutate(p_value     = parse_topgo_p(weight01_fisher),
           Significant = as.numeric(Significant)) %>%
    filter(!is.na(0.05), !is.na(Significant)) %>%
    arrange(p_value) %>%
    slice_head(n = 12) 
 
  if (nrow(df) == 0) return(NULL)
 
  # save summary file 
      write_tsv(df,file.path(outdir, paste0(genome, "_GO_Summary_", cls_name, ".tsv")))
  
    }

    # Short vs long comparison
    if (genome == "long") {
      short_lscp <- read_tsv(
        file.path(outdir, "short_gene_landscape_summary.tsv"),
        show_col_types = FALSE)

      if (!is.null(gene_map) && nrow(gene_map) > 0) {
        comp <- landscape %>%
          rename(long_gene_name = gene) %>%
          left_join(gene_map, by = "long_gene_name") %>%
          left_join(short_lscp %>% rename(short_gene_id = gene,
                                          mean_logcpm_short = mean_logcpm),
                    by = "short_gene_id") %>%
          filter(!is.na(mean_logcpm_short))

        if (nrow(comp) > 10) {
          p_comp <- ggplot(comp, aes(mean_logcpm_short, mean_logcpm)) +
            geom_point(alpha = 0.3, size = 1) +
            geom_smooth(method = "lm", se = FALSE,
                        colour = "#D55E00") +
            geom_abline(linetype = "dashed") +
            labs(title = paste(tissue_label, "- short vs long mean logCPM"),
                 x = "Short read", y = "Long read") +
            plot_theme
          ggsave(file.path(outdir, "short_vs_long_comparison.pdf"),
                 p_comp, width = 6, height = 5)
        }
      }
    }
  }
  message(tissue_label, " landscape complete.")
}

tissues <- c(AG = "Accessory_Gland",
             AB = "Abdomen",
             TE = "Testes",
             YO = "Distal_Ovary",
             OO = "Proximal_Ovary",
             EG = "Egg")

gene_map <- if (file.exists(gene_map_file))
  read_tsv(gene_map_file, show_col_types = FALSE) else NULL


for (code in names(tissues))
  build_tissue_landscape(code, tissues[code],min_cpm = 1, min_samples = 2,n_mixtures = 3,outdir_raw)




# MErging all sumamry go terms 
base_dir <- paste0(outdir_raw, "/05_Tissue_landscapes/Tissues")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

merged_go <- list.files(
  path = base_dir ,
  pattern = "^short_GO_Summary_(high|medium)\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  map_dfr(function(f) {

    tissue <- basename(dirname(f))

    classification <- str_extract(
      basename(f),
      "(?<=short_GO_Summary_)(high|medium)(?=\\.tsv)"
    )

    read_tsv(f, show_col_types = FALSE) |>
      mutate(
        tissue = tissue,
        classification = classification,
        weight01_fisher = as.character(weight01_fisher)
      )
  })

write_tsv(merged_go, file.path(base_dir , "merged_all_tissues_GO_summary.tsv"))

outdir_rank <- file.path(outdir_raw, "05_Tissue_landscapes", "rank_ecdf")
dir.create(outdir_rank, recursive = TRUE, showWarnings = FALSE)

adult_tissues <- c("AG", "AB", "TE", "YO", "OO")

# Build per-tissue mean logCPM table (short genome, adult tissues)
tissue_mean_lcpm <- lapply(adult_tissues, function(tc) {
  ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
  ids <- intersect(ids, colnames(counts_short_raw))
  if (length(ids) == 0) return(NULL)
  dge <- DGEList(counts = counts_short_raw[, ids, drop = FALSE])
  keep <- rowSums(cpm(dge) >= 1) >= 2
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)
  lcpm <- cpm(dge, log = TRUE, prior.count = 2)
  tibble(gene = rownames(lcpm),
         mean_logcpm = rowMeans(lcpm),
         tissue = tissues[tc])
}) %>%
  bind_rows()

# ---- Rank plot: each tissue as a curve ------------------------------------
rank_df <- tissue_mean_lcpm %>%
  group_by(tissue) %>%
  arrange(desc(mean_logcpm)) %>%
  mutate(rank = row_number(),
         n    = n(),
         rank_pct = rank / n * 100) %>%
  ungroup()

tissue_pal_adult <- setNames(
  okabe_ito_palette(length(unique(rank_df$tissue))),
  unique(rank_df$tissue))

p_rank <- ggplot(rank_df, aes(rank_pct, mean_logcpm,
                               colour = tissue, group = tissue)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = tissue_pal_adult) +
  labs(title  = "Expression rank plot ",
       x      = "Gene rank (percentile, highest to lowest)",
       y      = "Mean log2CPM",
       colour = "Tissue") +
  plot_theme
p_rank
ggsave(file.path(outdir_rank, "expression_rank_plot.pdf"), p_rank,
       width = 8, height = 5)

# ---- ECDF plot ------------------------------------------------------------
p_ecdf <- ggplot(tissue_mean_lcpm, aes(mean_logcpm,
                                        colour = tissue, group = tissue)) +
  stat_ecdf(linewidth = 0.9) +
  scale_colour_manual(values = tissue_pal_adult) +
  labs(title  = "Cumulative distribution of mean expression",
       x      = "Mean log2CPM",
       y      = "Cumulative fraction of genes",
       colour = "Tissue") +
  plot_theme
p_ecdf
ggsave(file.path(outdir_rank, "expression_ecdf.pdf"), p_ecdf,
       width = 8, height = 5)

# ---- Violin: expression distribution per tissue --------------------------

# Dunn post-hoc test
dunn_res <- tissue_mean_lcpm %>%
  dunn_test(mean_logcpm ~ tissue,
            p.adjust.method = "BH")

# Named vector of adjusted p-values
pvals <- setNames(
  dunn_res$p.adj,
  paste(dunn_res$group1, dunn_res$group2, sep = "-")
)

# Compact letter display
letters <- multcompLetters(
  pvals,
  threshold = 0.05
)$Letters

letter_df <- tibble(
  tissue = names(letters),
  letter = letters
)

# Position letters above each violin
letter_df <- letter_df %>%
  left_join(
    tissue_mean_lcpm %>%
      group_by(tissue) %>%
      summarise(
        y = max(mean_logcpm, na.rm = TRUE) + 1.2,
        .groups = "drop"
      ),
    by = "tissue"
  )

p_violin <- ggplot(tissue_mean_lcpm,
                   aes(tissue, mean_logcpm, fill = tissue, colour = tissue)) +
  ggdist::stat_halfeye(
    adjust = 0.6, width = 0.7, .width = c(0.66, 0.95),
    slab_alpha = 0.35, slab_colour = NA,
    point_size = 2, interval_size_range = c(0.5, 1.4),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = tissue_pal_adult) +
  scale_colour_manual(values = tissue_pal_adult) +
  labs(title = "Gene level expression distribution per tissue",
       x = NULL, y = "Mean log2CPM") +
  plot_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
p_violin
ggsave(file.path(outdir_rank, "expression_violin.pdf"), p_violin,
       width = 7, height = 5)

# with letters 
p_violin <- p_violin +
 geom_text(
  data = letter_df,
  aes(x = tissue, y = y, label = letter),
  inherit.aes = FALSE,
  size = 5,
  fontface = "bold"
)
ggsave(file.path(outdir_rank, "expression_violin_letters.pdf"), p_violin,
       width = 7, height = 5)

outdir_breadth <- file.path(outdir_raw, "05_Tissue_landscapes", "breadth")
dir.create(outdir_breadth, recursive = TRUE, showWarnings = FALSE)

# Expressed = mean logCPM > 0 (i.e. CPM > 1 after prior count)
breadth_mat <- tissue_mean_lcpm %>%
  mutate(expressed = mean_logcpm > 0) %>%
  group_by(gene) %>%
  summarise(n_tissues_expressed = sum(expressed),
            max_logcpm = max(mean_logcpm),
            .groups = "drop")

n_tissues_total <- length(adult_tissues)

p_breadth <- ggplot(breadth_mat, aes(factor(n_tissues_expressed))) +
  geom_bar(fill = "#0072B2", colour = "white") +
  geom_text(stat = "count",
            aes(label = scales::comma(after_stat(count))),
            vjust = -0.4, size = 3.2) +
  scale_x_discrete(labels = c(
    "0" = "None", as.character(seq_len(n_tissues_total)))) +
  labs(title = "Breadth of expression across adult tissues",
       x     = "Number of tissues expressed in",
       y     = "Gene count") +
  plot_theme
p_breadth
ggsave(file.path(outdir_breadth, "breadth_barplot.pdf"), p_breadth,
       width = 7, height = 5)

# Composition of single-tissue vs multi-tissue vs ubiquitous genes
breadth_summary <- breadth_mat %>%
  mutate(category = case_when(
    n_tissues_expressed == 0 ~ "Not expressed",
    n_tissues_expressed == 1 ~ "Single-tissue",
    n_tissues_expressed == n_tissues_total ~ "Ubiquitous",
    TRUE ~ "Multi-tissue (2-4)"
  )) %>%
  dplyr::count(category) %>%
  mutate(pct = n / sum(n) * 100,
         category = factor(category,
                           levels = c("Not expressed","Single-tissue",
                                      "Multi-tissue (2-4)","Ubiquitous")))

p_donut <- ggplot(breadth_summary, aes(x = "", y = n, fill = category)) +
  geom_col(colour = "white", linewidth = 0.5, width = 0.6) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = okabe_ito_palette(5)) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Expression breadth categories",
       x = NULL, y = "Gene count", fill = NULL) +
  plot_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())
p_donut
ggsave(file.path(outdir_breadth, "breadth_donut.pdf"), p_donut,
       width = 7, height = 3)

write_tsv(breadth_mat, file.path(outdir_breadth, "gene_breadth.tsv"))

outdir_bcv <- file.path(outdir_raw, "05_Tissue_landscapes", "BCV")
dir.create(outdir_bcv, recursive = TRUE, showWarnings = FALSE)

for (tc in adult_tissues) {
  ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
  ids <- intersect(ids, colnames(counts_short_raw))
  if (length(ids) < 2) next

  dge  <- DGEList(counts = counts_short_raw[, ids, drop = FALSE])
  keep <- rowSums(cpm(dge) >= 1) >= 2
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)

  # BCV plot requires a design; use intercept-only for single-group
  design_bc <- model.matrix(~ 1, data = data.frame(id = ids))
  dge <- tryCatch(estimateDisp(dge, design_bc, robust = TRUE),
                  error = function(e) estimateDisp(dge))

  # Extract into a tidy data frame for ggplot
  bcv_df <- tibble(
    logCPM       = aveLogCPM(dge),
    sqrt_disp_tagged = sqrt(dge$tagwise.dispersion),
    sqrt_disp_trended = sqrt(dge$trended.dispersion),
    sqrt_disp_common  = sqrt(dge$common.dispersion)
  )

  p_bcv <- ggplot(bcv_df, aes(logCPM, sqrt_disp_tagged)) +
    geom_point(alpha = 0.2, size = 0.8, colour = "grey50") +
    geom_line(aes(y = sqrt_disp_trended),
              colour = "#D55E00", linewidth = 1) +
    geom_hline(aes(yintercept = sqrt_disp_common),
               colour = "#0072B2", linetype = "dashed",
               linewidth = 0.8) +
    annotate("text", x = min(bcv_df$logCPM) + 0.5,
             y = unique(bcv_df$sqrt_disp_common) + 0.02,
             label = paste0("Common BCV = ",
                            round(unique(bcv_df$sqrt_disp_common), 3)),
             hjust = 0, size = 3, colour = "#0072B2") +
    labs(title = paste0(tissues[tc], " - BCV plot"),
         x = "Average log2CPM",
         y = "Biological coefficient of variation") +
    plot_theme
  ggsave(file.path(outdir_bcv, paste0(tc, "_BCV.pdf")), p_bcv,
         width = 6, height = 5)
}

message("BCV plots complete.")

outdir_hk <- file.path(outdir_raw, "05_Tissue_landscapes", "housekeeping")
dir.create(outdir_hk, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# STEP 1: Load all landscape summaries
# ===========================================================================
adult_tissues <- c("AG", "AB", "TE", "YO", "OO")

all_lscps <- lapply(adult_tissues, function(tc) {
  f <- file.path(outdir_raw, "05_Tissue_landscapes/Tissues",
                 tissues[tc], "short_gene_landscape_summary.tsv")
  if (!file.exists(f)) return(NULL)
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(tissue = tissues[tc],
           gene   = clean_ids(gene))
}) %>% bind_rows()

# ===========================================================================
# STEP 2: Expression-based constitutive gene set
#   Genes that are high-expression, low-CV, constitutive in ALL adult tissues
# ===========================================================================
constitutive_all <- all_lscps %>%
  group_by(gene) %>%
  summarise(
    n_tissues_constitutive = sum(is_constitutive, na.rm = TRUE),
    n_tissues_total        = n(),
    mean_expr_all          = mean(mean_logcpm, na.rm = TRUE),
    mean_cv_all            = mean(cv,          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    n_tissues_constitutive == length(adult_tissues),
    mean_expr_all > 2,
    mean_cv_all   < 0.3
  ) %>%
  arrange(desc(mean_expr_all))

message("Empirically constitutive genes: ", nrow(constitutive_all))
write_tsv(constitutive_all, file.path(outdir_hk, "constitutive_all_tissues.tsv"))

# ===========================================================================
# STEP 3: GO enrichment on the empirical constitutive set
#   Expect: ribosome, translation, mitochondria, proteasome
#   If these dominate the classification is working
# ===========================================================================
go_const <- run_topgo_all_ont(
  constitutive_all$gene,
  gene_universe_short,
  geneID2GO_short
)

# Save full panel
p_go_const <- plot_topgo_panel(go_const,
                               paste0("Constitutive genes across all adult tissues (n = ",
                                      nrow(constitutive_all), ")"))
if (!is.null(p_go_const))
  ggsave(file.path(outdir_hk, "GO_constitutive_all_tissues.pdf"),
         p_go_const, width = 9, height = 12)

# Compact bar for inline display
go_const_bp <- go_const$BP
if (!is.null(go_const_bp) && nrow(go_const_bp) > 0) {
  p_const_bar <- go_const_bp %>%
    mutate(p_value     = parse_topgo_p(weight01_fisher),
           Significant = as.numeric(Significant)) %>%
    filter(!is.na(p_value)) %>%
    arrange(p_value) %>%
    slice_head(n = 12) %>%
    mutate(Term = str_wrap(Term, 45),
           Term = fct_reorder(Term, Significant)) %>%
    ggplot(aes(Term, Significant, fill = -log10(p_value))) +
    geom_col() +
    coord_flip() +
    scale_fill_gradientn(colours = enrichment_palette(100)) +
    labs(title    = "GO: Biological Process - empirical constitutive genes",
         subtitle = paste0("n = ", nrow(constitutive_all),
                           " genes constitutive across all adult tissues"),
         x = NULL, y = "Significant genes",
         fill = expression(-log[10](p))) +
    plot_theme
  p_const_bar
  ggsave(file.path(outdir_hk, "GO_constitutive_BP_bar.pdf"),
         p_const_bar, width = 9, height = 6)
}

# ===========================================================================
# STEP 4: Orthology-based prong
#   Map eggNOG Preferred_name to FlyAtlas2 ubiquitous gene list
#   Then get gene IDs matching the landscape table via the GFF
#   Finally: category bar chart showing which HK functional classes
#   are represented in C. capitata
# ===========================================================================

# ---- Reference list -------------------------------------------------------
# FlyAtlas2 ubiquitous gene symbols (Leader et al. 2018)
# Covers Drosophila AND human naming conventions (both seen in eggNOG)
# Everything uppercased before matching so capitalisation doesn't matter

flyatlas2_upper <- toupper(c(
  # Cytoplasmic ribosome — large subunit
  paste0("RpL", c(3:15, 17:19, 21:24, 26:32, 35, 36, 38:41)),
  paste0("RPL", c(3:15, 17:19, 21:24, 26:32, 35, 36, 38:41)),
  "RPL10A","RPL13A","RPL18A","RPL23A","RPL27A","RPL36AL","RPL22L1",
  # Cytoplasmic ribosome — small subunit
  paste0("RpS", c(2:11, 13:21, 23:30)),
  paste0("RPS", c(2:11, 13:21, 23:30)),
  "RPS15A","RPS27A","RPS27L","FAU",
  # Mitochondrial ribosome
  paste0("mRpL", c(1:4, 9:17, 19, 20, 22:24, 27, 28,
                   32, 33, 35:40, 44, 46, 48:55)),
  paste0("mRpS", c(5, 7, 9:12, 14:18, 21:31, 33:36)),
  paste0("MRPL", c(1:4, 9:17, 19, 20, 22:24, 27, 28,
                   32, 33, 35:40, 44, 46, 48:55)),
  paste0("MRPS", c(5, 7, 9:12, 14:18, 21:31, 33:36)),
  # Translation factors
  "eIF-2alpha","eIF2beta","eIF2gamma",
  paste0("eIF3-S", c(4, 6, 8, 9, 10)),
  "eIF4A","eIF4E-1","eIF4G","eIF5","eIF5B",
  "EF1alpha48D","EF1alpha100E","EF1beta","EF1gamma","EF2",
  "EIF2S1","EIF2S2","EIF2S3",
  paste0("EIF3", c("A","B","C","D","E","F","G","H","I","J","K","L","M")),
  "EIF4A1","EIF4A2","EIF4E","EIF4G1","EIF4G2","EIF5B",
  "EEF1A1","EEF1A2","EEF1B2","EEF1G","EEF2",
  # Oxidative phosphorylation
  "ND2","ND4","ND5",
  paste0("NDUFS", 1:8), paste0("NDUFV", 1:3),
  paste0("NDUFA", 1:13), paste0("NDUFB", 1:11), paste0("NDUFC", 1:2),
  "SdhA","SdhB","SdhC","SdhD","SDHA","SDHB","SDHC","SDHD",
  "CG3621","CG4169","UQCRC1","UQCRC2","UQCRB","UQCRFS1","UQCR10","UQCR11",
  "CoVa","CoVb","levy","cype",
  "COX4I1","COX4I2","COX5A","COX5B","COX6A1","COX6A2","COX6B1","COX6C",
  "COX7A1","COX7A2","COX7B","COX7C","COX8A",
  "ATPsynalpha","ATPsynbeta","ATPsyngamma","ATPsynd","ATPsyne",
  "ATPsynG","ATPsynCF6","ATPsyneps",
  "ATP5A1","ATP5B","ATP5C1","ATP5D","ATP5E",
  "ATP5F1","ATP5G1","ATP5G2","ATP5G3","ATP5H","ATP5I",
  "ATP5J","ATP5J2","ATP5L","ATP5O",
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
  "ATP5MC1","ATP5MC2","ATP5MC3","ATP5MF","ATP5MG",
  "ATP5PB","ATP5PD","ATP5PF","ATP5PO",
  # Core metabolism — glycolysis
  "Gapdh1","Gapdh2","GAPDH","GAPD",
  "Pgk","PGK1","PGK2",
  "Eno","ENO1","ENO2","ENO3",
  "Pfk","PFKL","PFKM","PFKP",
  "Pgi","GPI","Pgm","PGM1","PGM2",
  "Ald","ALDOA","ALDOB","ALDOC",
  "Tpi","TPI1","Pyk","PKLR","PKM",
  # Core metabolism — TCA cycle
  "Cs","CS","Aco","ACO1","ACO2",
  "Idh","IDH1","IDH2","IDH3A","IDH3B","IDH3G",
  "sucA","sucB","SUCLG1","SUCLG2","SUCLA2",
  "Scs-alfa","Mdh1","Mdh2","MDH1","MDH2","Men","ME1","ME2",
  # Cytoskeleton
  "Act5C","Act42A","Act57B","Act79B","Act87E","Act88F",
  "ACTB","ACTG1","ACTA1","ACTA2",
  "alpha-Tub84B","alpha-Tub84D","beta-Tub56D","beta-Tub60D",
  "TUBA1A","TUBA1B","TUBA1C","TUBA4A",
  "TUBB","TUBB2A","TUBB2B","TUBB3","TUBB4A","TUBB4B","TUBB6",
  # Proteasome — 20S core
  paste0("Prosalpha", 1:7), paste0("Prosbeta", 1:7),
  paste0("PSMA", 1:7),      paste0("PSMB", 1:7),
  # Proteasome — 19S regulatory
  paste0("Rpn", c(1:3, 5:12)), paste0("Rpt", 1:6),
  paste0("PSMC", 1:6),
  paste0("PSMD", c(1:4, 6:8, 11:14)),
  "PSMD5","PSMD9","PSMD10",
  # Ubiquitin pathway
  "Ubi-p63E","Ubi-p5E","UBB","UBC","UBA52",
  "Ubc84D","Ubc87F","Ubc2","Ubc4","Ubc6",
  "UBE2A","UBE2B","UBE2C","UBE2D1","UBE2D2","UBE2D3",
  "UBE2E1","UBE2I","UBE2K","UBE2L3","UBE2N",
  # Chaperones
  "Hsc70-1","Hsc70-2","Hsc70-3","Hsc70-4","Hsp83","CG5525",
  "HSPA1A","HSPA1B","HSPA2","HSPA5","HSPA8","HSPA9",
  "HSP90AA1","HSP90AB1","HSP90B1","HSPD1","HSPE1",
  paste0("CCT", 2:8), "TCP1",
  # RNA processing / spliceosome
  "SmB","SmD1","SmD2","SmD3","SmE","SmF","SmG",
  "SNRPB","SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF","SNRPG",
  "U2af38","U2af50","U2AF1","U2AF2",
  "Prp8","Prp31","PRPF8","PRPF31",
  "sf3a1","sf3b1","SF3A1","SF3B1","SF3B2","SF3B3",
  "CG6905",paste0("LSM", 2:8)
)) %>%
  unique()

# Category lookup — for the bar chart
# Each symbol maps to one functional category
hk_categories <- c(
  setNames(rep("Cytoplasmic ribosome (LSU)", length(grep("^RPL|^RpL", flyatlas2_upper, value = TRUE))),
           grep("^RPL|^RpL", flyatlas2_upper, value = TRUE)),
  setNames(rep("Cytoplasmic ribosome (SSU)", length(grep("^RPS|^RpS", flyatlas2_upper, value = TRUE))),
           grep("^RPS|^RpS", flyatlas2_upper, value = TRUE)),
  setNames(rep("Mitochondrial ribosome", length(grep("^MRPL|^MRPS|^MRpL|^MRpS", flyatlas2_upper, value = TRUE))),
           grep("^MRPL|^MRPS|^MRpL|^MRpS", flyatlas2_upper, value = TRUE)),
  setNames(rep("Translation factors", length(grep("^EIF|^EEF|^EF1|^eIF|^EF2", flyatlas2_upper, value = TRUE))),
           grep("^EIF|^EEF|^EF1|^eIF|^EF2", flyatlas2_upper, value = TRUE)),
  setNames(rep("Oxidative phosphorylation", length(grep("^NDUF|^COX|^ATP|^SDH|^UQCR|^ND[245]|ATPSYN|COVA|COVB|LEVY|CYPE|CG3621|CG4169", flyatlas2_upper, value = TRUE))),
           grep("^NDUF|^COX|^ATP|^SDH|^UQCR|^ND[245]|ATPSYN|COVA|COVB|LEVY|CYPE|CG3621|CG4169", flyatlas2_upper, value = TRUE)),
  setNames(rep("Core metabolism", length(grep("GAPDH|GAPD|^PGK|^ENO|^PFK|^GPI|^PGM|^ALD|^TPI|^PKM|PKLR|^CS$|^ACO|^IDH|^MDH|^MDH|SUCLG|SUCLA|SCS|^ME[12]$|^MEN$|^ACO|^CS$|^Cs$|^Aco$|^Pgk$|^Eno$|^Pfk$|^Pgi$|^Pgm$|^Ald$|^Tpi$|^Pyk$|GAPDH|Gapdh|^Mdh|^Idh", flyatlas2_upper, value = TRUE))),
           grep("GAPDH|GAPD|^PGK|^ENO|^PFK|^GPI|^PGM|^ALD|^TPI|^PKM|PKLR|^CS$|^ACO|^IDH|^MDH|SUCLG|SUCLA|SCS|^ME[12]$|^MEN$|^ACO|^CS$|^Aco$|^Pgk$|^Eno$|^Pfk$|^Pgi$|^Pgm$|^Ald$|^Tpi$|^Pyk$|GAPDH|Gapdh|^Mdh|^Idh", flyatlas2_upper, value = TRUE)),
  setNames(rep("Cytoskeleton", length(grep("^ACT|^TUB|^TUBA|^TUBB|TUB84|TUB56|TUB60", flyatlas2_upper, value = TRUE))),
           grep("^ACT|^TUB|^TUBA|^TUBB|TUB84|TUB56|TUB60", flyatlas2_upper, value = TRUE)),
  setNames(rep("Proteasome / Ubiquitin", length(grep("^PSM|^PROS|^RPN|^RPT|^UBB|^UBC|^UBA|^UBE|^UBI", flyatlas2_upper, value = TRUE))),
           grep("^PSM|^PROS|^RPN|^RPT|^UBB|^UBC|^UBA|^UBE|^UBI", flyatlas2_upper, value = TRUE)),
  setNames(rep("Chaperones", length(grep("^HSP|^CCT|^TCP1|^CG5525", flyatlas2_upper, value = TRUE))),
           grep("^HSP|^CCT|^TCP1|^CG5525", flyatlas2_upper, value = TRUE)),
  setNames(rep("RNA processing", length(grep("^SNRP|^SM[BDEFG]|^U2AF|^PRPF|^SF3|^LSM|^CG6905", flyatlas2_upper, value = TRUE))),
           grep("^SNRP|^SM[BDEFG]|^U2AF|^PRPF|^SF3|^LSM|^CG6905", flyatlas2_upper, value = TRUE))
)

# ---- Read eggNOG and match ------------------------------------------------
eggnog_raw_file <- here("01_data/GO_terms/eggNOG-mapper/medfly_short.emapper.annotations")

hk_prot_map <- NULL

if (file.exists(eggnog_raw_file)) {
  egg <- tryCatch(
    read_tsv(eggnog_raw_file, comment = "##",
             show_col_types = FALSE, skip = 3),
    error = function(e) { message("Could not read eggNOG: ", e); NULL }
  )
  
  if (!is.null(egg)) {
    query_col <- names(egg)[1]
    name_col  <- names(egg)[
      str_detect(names(egg), regex("preferred.name", ignore_case = TRUE))
    ][1]
    
    hk_prot_map <- egg %>%
      dplyr::select(protein        = all_of(query_col),
                    preferred_name = all_of(name_col)) %>%
      filter(!is.na(preferred_name),
             preferred_name != "-",
             preferred_name != "") %>%
      mutate(name_upper = toupper(preferred_name),
             is_hk      = name_upper %in% flyatlas2_upper) %>%
      filter(is_hk) %>%
      mutate(category = hk_categories[name_upper]) %>%
      dplyr::select(protein, preferred_name, name_upper, category)
    
    message("Protein accessions matching FlyAtlas2 HK list: ", nrow(hk_prot_map))
  }
}

# ---- Map protein accessions to gene IDs via GFF ----------------------------
gff_file <- here("01_data/Genomic_files/GCF_000347755.3_Ccap_2.1_genomic.gff.gz")

ccap_hk_genes <- character(0)
hk_full_map   <- tibble()

if (!is.null(hk_prot_map) && nrow(hk_prot_map) > 0 && file.exists(gff_file)) {
  
  prot2gene <- read_tsv(
    gff_file,
    comment   = "#",
    col_names = c("seqid","source","type","start","end",
                  "score","strand","phase","attributes"),
    show_col_types = FALSE
  ) %>%
    filter(type == "CDS") %>%
    mutate(
      protein = str_extract(attributes, "(?<=protein_id=)[^;]+"),
      gene    = str_extract(attributes, "(?<=gene=)[^;]+"),
      locus   = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
      gene_id = coalesce(gene, locus)
    ) %>%
    filter(!is.na(protein), !is.na(gene_id)) %>%
    distinct(protein, gene_id)
  
  hk_full_map <- hk_prot_map %>%
    left_join(prot2gene, by = "protein") %>%
    filter(!is.na(gene_id),
           gene_id %in% all_lscps$gene) %>%
    distinct(gene_id, preferred_name, category)
  
  ccap_hk_genes <- hk_full_map$gene_id
  
  message("HK genes successfully mapped to landscape gene IDs: ",
          length(ccap_hk_genes))
  
  write_tsv(hk_full_map,
            file.path(outdir_hk, "Ccap_FlyAtlas2_HK_mapping.tsv"))
  
} else {
  message("GFF mapping failed - falling back to GO-term HK proxy")
  gene2go_tbl <- tibble(
    gene  = rep(names(geneID2GO_short), lengths(geneID2GO_short)),
    go_id = unlist(geneID2GO_short)
  )
  ccap_hk_genes <- gene2go_tbl %>%
    filter(go_id %in% c("GO:0005840","GO:0006412","GO:0005739",
                        "GO:0003735","GO:0006511")) %>%
    pull(gene) %>%
    unique()
  message("GO-proxy HK genes: ", length(ccap_hk_genes))
}

# ---- Category bar chart ---------------------------------------------------
if (nrow(hk_full_map) > 0) {
  
  cat_tbl <- table(hk_full_map$category)
  cat_df  <- data.frame(
    category = names(cat_tbl),
    n        = as.integer(cat_tbl),
    stringsAsFactors = FALSE
  )
  cat_df  <- cat_df[order(cat_df$n), ]           # ascending so coord_flip puts highest at top
  cat_df$category <- factor(cat_df$category, levels = cat_df$category)
  
  p_cat_summary <- ggplot(cat_df, aes(category, n, fill = category)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    
    geom_richtext(
      aes(
        label = n),
      fill = "white",
      label.color = "darkgrey",  
      size = 3
    ) + 
    scale_fill_manual(values = setNames(
      okabe_ito_palette(nrow(cat_df)),
      levels(cat_df$category)
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    coord_flip() +
    labs(
      title    = "Drosophila FlyAtlas2 HK orthologues in C. capitata",
      subtitle = paste0(length(ccap_hk_genes),
                        " Ccap genes with a ubiquitous Dmel orthologue",
                        " (FlyAtlas2; Leader et al. 2018)"),
      x = NULL, y = "Gene count"
    ) +
    plot_theme
  
  p_cat_summary
  ggsave(file.path(outdir_hk, "FlyAtlas2_HK_ortho_categories.pdf"),
         p_cat_summary, width = 8, height = 6)
}

# ===========================================================================
# STEP 5: Are the HK orthologues constitutive in C. capitata?
#   Fisher's exact test: HK orthologue × is_constitutive
# ===========================================================================
# One row per gene — use first tissue (all tissues have the same constitutive flag)
hk_landscape <- all_lscps %>%
  filter(tissue == tissues[adult_tissues[1]]) %>%
  mutate(is_hk_ortho = gene %in% ccap_hk_genes)

message("HK genes found in landscape: ", sum(hk_landscape$is_hk_ortho))

if (sum(hk_landscape$is_hk_ortho) > 5) {
  ct <- table(hk_ortho        = hk_landscape$is_hk_ortho,
              is_constitutive = hk_landscape$is_constitutive)
  ft <- fisher.test(ct)
  message("Fisher OR = ", round(ft$estimate, 2),
          "  p = ", signif(ft$p.value, 3))
  
  p_hk_enrich <- hk_landscape %>%
    mutate(gene_class = if_else(is_hk_ortho, "HK orthologue", "Background"),
           expr_class  = factor(expr_class, c("low","medium","high"))) %>%
    filter(!is.na(expr_class)) %>%
    dplyr::count(gene_class, expr_class) %>%
    group_by(gene_class) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(gene_class, frac, fill = expr_class)) +
    geom_col(colour = "lightgrey", linewidth = 0.4) +
    scale_fill_manual(values = c(low="#56B4E9" ,
                                 medium="white",
                                 high= "#E69F00")) +
    scale_y_continuous(labels = scales::percent) +
    annotate("text", x = 1.5, y = 1.08,
             label = paste0("Fisher OR = ", round(ft$estimate, 2),
                            "\np = ", signif(ft$p.value, 3)),
             size = 3.5, hjust = 0.5) +
    labs(title    = "A - Dmel HK orthologues enriched in high-expression class",
         subtitle = paste0(sum(hk_landscape$is_hk_ortho),
                           " HK orthologues vs ",
                           sum(!hk_landscape$is_hk_ortho),
                           " background genes"),
         x = NULL, y = "Fraction", fill = "Expression class") +
    plot_theme
  
  p_hk_enrich
  ggsave(file.path(outdir_hk, "HK_ortho_expression_class.pdf"),
         p_hk_enrich, width = 6, height = 5)

  # =========================================================================
  # STEP 6: Do HK orthologues cluster in low-CV space?
  #   Wilcoxon rank-sum test on mean cross-tissue CV, HK orthologue vs
  #   background.
  # =========================================================================
  cv_all_genes <- all_lscps %>%
    group_by(gene) %>%
    summarise(mean_cv_all = mean(cv, na.rm = TRUE), .groups = "drop") %>%
    mutate(is_hk_ortho = gene %in% ccap_hk_genes) %>%
    filter(is.finite(mean_cv_all))

  # Test on the full, unclipped data -- Wilcoxon is rank-based and robust to
  # the small number of extreme outliers described below.
  wt_cv <- wilcox.test(mean_cv_all ~ is_hk_ortho, data = cv_all_genes)
  message("HK ortho vs background mean CV -- Wilcoxon p = ",
          signif(wt_cv$p.value, 3),
          "  (HK median = ", signif(median(cv_all_genes$mean_cv_all[cv_all_genes$is_hk_ortho]), 3),
          ", background median = ", signif(median(cv_all_genes$mean_cv_all[!cv_all_genes$is_hk_ortho]), 3), ")")

  # CV = SD(log2CPM) / (mean(log2CPM) + 0.1) is unstable for genes whose mean
  # log2CPM sits near -0.1 (denominator near zero) -- a handful of very
  # lowly-expressed genes get an extreme or even negative mean_cv_all that
  # is a formula artifact, not a real variability estimate. Excluded from
  # the *plot* only (a log10 axis needs positive values regardless); the
  # Wilcoxon test above already ran on the complete, unfiltered data.
  plot_df_cv <- cv_all_genes %>%
    filter(mean_cv_all > 0) %>%
    mutate(gene_class = if_else(is_hk_ortho, "HK orthologue", "Background"))

  # Half-eye: density with median and interval.
  p_hk_cv <- ggplot(plot_df_cv, aes(x = gene_class, y = mean_cv_all,
                                     colour = gene_class, fill = gene_class)) +
    ggdist::stat_halfeye(
      adjust = 0.5, width = 0.6, .width = c(0.5, 0.95),
      justification = -0.15, point_colour = "black", alpha = 0.6
    ) +
    scale_colour_manual(values = c(Background = "grey50", `HK orthologue` = "#009E73"), guide = "none") +
    scale_fill_manual(values = c(Background = "grey80", `HK orthologue` = "#009E73"), guide = "none") +
    scale_y_log10() +
    annotate("text", x = 1.5, y = max(plot_df_cv$mean_cv_all) * 0.9,
             label = paste0("Wilcoxon p = ", signif(wt_cv$p.value, 3)),
             size = 3.5, hjust = 0.5) +
    labs(title    = "B - HK orthologues show lower cross-tissue expression variability",
         subtitle = paste0(sum(cv_all_genes$is_hk_ortho), " HK orthologues vs ",
                           sum(!cv_all_genes$is_hk_ortho), " background genes"),
         x = NULL, y = "Mean CV across adult tissues (log scale)") +
    plot_theme

  p_hk_cv
  ggsave(file.path(outdir_hk, "HK_ortho_CV_comparison.pdf"),
         p_hk_cv, width = 6, height = 5)

  # Combined two-panel Figure S3: expression-class enrichment + CV comparison.
  # Panel letters are baked into each panel's own title ("A - ...", "B - ...")
  # rather than added via plot_annotation(tag_levels), so no separate tag is
  # layered on top.
  p_hk_summary <- p_hk_enrich + p_hk_cv
  ggsave(file.path(outdir_hk, "HK_ortho_summary.pdf"),
         p_hk_summary, width = 11, height = 5)
}

outdir_bb <- file.path(outdir_raw, "05_Tissue_landscapes", "battenberg")
dir.create(outdir_bb, recursive = TRUE, showWarnings = FALSE)

# Build a single normalised expression matrix across all adult samples
meta_adult_bb <- meta_raw %>%
  filter(tissue %in% adult_tissues)

counts_bb <- counts_short_raw[, meta_adult_bb$id, drop = FALSE]
dge_bb    <- DGEList(counts = counts_bb)
keep_bb   <- rowSums(cpm(dge_bb) >= 1) >= 3
dge_bb    <- dge_bb[keep_bb, , keep.lib.sizes = FALSE]
dge_bb    <- calcNormFactors(dge_bb)
lcpm_bb   <- cpm(dge_bb, log = TRUE, prior.count = 2)

# Select top 1000 variable genes (by inter-sample variance)
gene_var_bb <- rowVars(lcpm_bb)
top1k_bb    <- names(sort(gene_var_bb, decreasing = TRUE))[1:min(1000, length(gene_var_bb))]
mat_bb      <- lcpm_bb[top1k_bb, ]

# Z-score rows
mat_z <- t(scale(t(mat_bb)))
mat_z[mat_z >  3] <-  3
mat_z[mat_z < -3] <- -3

# Annotation bars
ann_col_bb <- meta_adult_bb %>%
  dplyr::select(id, tissue, sex, mating_status) %>%
  column_to_rownames("id") %>%
  mutate(tissue = tissues[tissue])

# Colour lists for annotation
ann_pal_bb <- list(
  tissue        = annotation_colours$tissue,
  sex           = annotation_colours$sex,
  mating_status = annotation_colours$mating_status
)


# ComplexHeatmap Battenberg
ha_bb <- HeatmapAnnotation(
  tissue         = ann_col_bb$tissue,
  sex            = ann_col_bb$sex,
  mating_status  = ann_col_bb$mating_status,
  col            = ann_pal_bb,
  annotation_name_gp = gpar(fontsize = 9),
  na_col         = "grey70"
)

ht_bb <- ComplexHeatmap::Heatmap(
  mat_z,
  name             = "Z-score",
  top_annotation   = ha_bb,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names   = FALSE,
  show_column_names = TRUE,
  column_names_rot  = 45,
  column_names_gp   = gpar(fontsize = 8),
  col = circlize::colorRamp2(
    c(-3, -1.5, 0, 1.5, 3),
    c("#0072B2", "#56B4E9", "white", "#E69F00", "#D55E00")
  ),
  row_title = paste0("Top ", nrow(mat_z), " variable genes"),
  column_title = "All adult tissue samples",
  use_raster = TRUE, raster_quality = 2
)

pdf(file.path(outdir_bb, "battenberg_heatmap.pdf"), width = 14, height = 10)
draw(ht_bb)
dev.off()

message("Battenberg heatmap complete.")

outdir_sp <- file.path(outdir_raw, "05_Tissue_landscapes", "shared_private")
dir.create(outdir_sp, recursive = TRUE, showWarnings = FALSE)

# ---- Build expressed gene sets --------------------------------------------
# First get per-compartment sets
expressed_raw <- lapply(adult_tissues, function(tc) {
  ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
  ids <- intersect(ids, colnames(lcpm_bb))
  if (length(ids) == 0) return(character(0))
  rownames(lcpm_bb)[rowMeans(lcpm_bb[, ids, drop = FALSE]) > 0]
}) %>% setNames(tissues[adult_tissues])

# Merge distal and proximal ovary into one "Ovary" set
# A gene counts as expressed in the ovary if it is expressed in EITHER compartment
expressed_sets <- list(
  "Accessory gland" = expressed_raw[["Accessory_Gland"]],
  "Abdomen"         = expressed_raw[["Abdomen"]],
  "Testes"          = expressed_raw[["Testes"]],
  "Ovary"           = union(expressed_raw[["Distal_Ovary"]],
                            expressed_raw[["Proximal_Ovary"]])
)

message("Expressed genes per tissue/organ:")
sapply(expressed_sets, length) %>% print()

# ---- UpSet plot ------------------------------------------------------------
pdf(file.path(outdir_sp, "UpSet_expressed_genes.pdf"), width = 12, height = 7)
upset(
  fromList(expressed_sets),
  order.by        = "freq",
  sets            = names(expressed_sets),
  keep.order      = FALSE,
  sets.bar.color  = okabe_ito_palette(length(expressed_sets)),
  main.bar.color  = "grey30",
  text.scale      = 1.3,
  mainbar.y.label = "Gene count",
  sets.x.label    = "Expressed genes per organ",
  mb.ratio        = c(0.6, 0.4)
)
dev.off()

# ---- Private vs shared counts ---------------------------------------------
private_counts <- sapply(names(expressed_sets), function(tc) {
  others <- unname(unlist(expressed_sets[names(expressed_sets) != tc]))
  sum(!expressed_sets[[tc]] %in% others)
})

shared_counts <- sapply(names(expressed_sets), function(tc) {
  all_others <- Reduce(intersect,
                       expressed_sets[names(expressed_sets) != tc])
  sum(expressed_sets[[tc]] %in% all_others)
})

summary_sp <- tibble(
  organ           = names(expressed_sets),
  total_expressed = sapply(expressed_sets, length),
  private         = private_counts,
  shared_all      = shared_counts,
  other           = total_expressed - private - shared_all,
  pct_private     = private / total_expressed * 100
) %>%
  mutate(organ = factor(organ, levels = organ[order(total_expressed)]))

write_tsv(summary_sp, file.path(outdir_sp, "shared_private_summary.tsv"))

# Also save the split version for reference
tibble(
  compartment     = names(expressed_raw),
  total_expressed = sapply(expressed_raw, length)
) %>%
  write_tsv(file.path(outdir_sp, "expressed_per_compartment_split.tsv"))

# ---- Stacked bar ----------------------------------------------------------
category_cols <- c(
  "Tissue-private"   = "#D55E00",
  "Partially shared" = "#56B4E9",
  "Shared by all"    = "#0072B2"
)

p_sp_bar <- summary_sp %>%
  pivot_longer(c(private, other, shared_all),
               names_to  = "category",
               values_to = "count") %>%
  mutate(category = factor(category,
                           levels = c("private", "other", "shared_all"),
                           labels = c("Tissue-private",
                                      "Partially shared",
                                      "Shared by all"))) %>%
  ggplot(aes(organ, count, fill = category)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = category_cols) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title    = "Expressed gene sharing across adult organs",
       subtitle = "Ovary = union of distal and proximal compartments",
       x = NULL, y = "Gene count", fill = NULL) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

p_sp_bar
ggsave(file.path(outdir_sp, "shared_private_stacked.pdf"),
       p_sp_bar, width = 7, height = 5)

# ---- Uniqueness dot plot --------------------------------------------------
organ_cols <- c(
  "Accessory gland" = "#C46E7A",
  "Abdomen"         = "#D79A4A",
  "Testes"          = "#4A8D8A",
  "Ovary"           = "#0072B2"
)

p_sp_dot <- ggplot(summary_sp,
                   aes(total_expressed, pct_private,
                       colour = organ, label = organ)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3.8, show.legend = FALSE,
                  box.padding = 0.5, max.overlaps = 10) +
  scale_colour_manual(values = organ_cols) +
  scale_x_continuous(labels = scales::comma) +
  labs(title    = "Transcriptional uniqueness per organ",
       x = "Total expressed genes",
       y = "% tissue-private genes") +
  plot_theme +
  theme(legend.position = "none")

p_sp_dot
ggsave(file.path(outdir_sp, "uniqueness_dot.pdf"),
       p_sp_dot, width = 6, height = 5)

# ---- Side-by-side: split vs merged ovary ----------------------------------
# Unique and total expressed genes with the ovary regions split and merged.
split_df <- tibble(
  organ           = c("Distal ovary", "Proximal ovary"),
  total_expressed = c(length(expressed_raw$Distal_Ovary),
                      length(expressed_raw$Proximal_Ovary)),
  private         = sapply(c("Distal_Ovary","Proximal_Ovary"), function(tc) {
    others <- unname(unlist(expressed_raw[names(expressed_raw) != tc]))
    sum(!expressed_raw[[tc]] %in% others)
  }),
  approach        = "Split compartments"
) %>%
  mutate(pct_private = private / total_expressed * 100)

merged_df <- summary_sp %>%
  filter(organ == "Ovary") %>%
  mutate(approach = "Merged ovary")

compare_df <- bind_rows(
  split_df %>% dplyr::select(organ, total_expressed, pct_private, approach),
  merged_df %>% dplyr::select(organ, total_expressed, pct_private, approach)
)

p_compare <- ggplot(compare_df,
                    aes(total_expressed, pct_private,
                        colour = approach, label = organ)) +
  geom_point(size = 5) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Split compartments" = "#56B4E9",
    "Merged ovary"       =  "#0072B2"
  )) +
  scale_x_continuous(labels = scales::comma) +
  labs(title    = "Effect of merging ovary compartments",
       x = "Total expressed genes",
       y = "% tissue-private genes",
       colour = NULL) +
  plot_theme +
  theme(legend.position = "bottom")

p_compare
ggsave(file.path(outdir_sp, "ovary_split_vs_merged.pdf"),
       p_compare, width = 6, height = 5)

# Summary table
summary_sp %>%
  mutate(across(where(is.numeric), round, 1)) %>%
  knitr::kable(caption = "Expressed gene sharing - organs (ovary merged)")

outdir_gomap <- file.path(outdir_raw, "05_Tissue_landscapes", "GO_comparison")
dir.create(outdir_gomap, recursive = TRUE, showWarnings = FALSE)

# For each tissue, run topGO on the high-expressed gene set (short genome)
# Collect the top 8 BP terms per tissue, union the term set, then score all tissues

tissue_list <- c(
  "Abdomen",
 "Accessory_Gland",
 "Distal_Ovary",
 "Egg",
 "Proximal_Ovary",
 "Testes"
)

collect_go_terms <- function(tc, n_top = 30) {
  lscp_f <- file.path(outdir_raw, "05_Tissue_landscapes/Tissues",
                       tc, "short_gene_landscape_summary.tsv")
  if (!file.exists(lscp_f)) return(NULL)
  lscp <- read_tsv(lscp_f, show_col_types = FALSE)

  high_genes <- lscp %>%
    filter(expr_class == "high") %>%
    pull(gene) %>% clean_ids()

  res <- run_topgo(high_genes, gene_universe_short, geneID2GO_short, "BP")
  if (is.null(res) || nrow(res) == 0) return(NULL)

  res %>%
    mutate(
      p_value = parse_topgo_p(weight01_fisher),
      tissue  = tc
    ) %>%
    arrange(p_value) %>%
    slice_head(n = n_top) %>%
    dplyr::select(GO.ID, Term, p_value, tissue)
}

go_per_tissue <- lapply(tissue_list, collect_go_terms) %>%
  bind_rows()

if (nrow(go_per_tissue) > 0) {
  # Union of all terms that appear in >= 2 tissues
  term_freq <- go_per_tissue %>%
    dplyr::count(Term, sort = TRUE) %>%
    filter(n >= 2)

  # Re-score all tissues for those shared terms
  rescore_all <- function(tc) {
    lscp_f <- file.path(outdir_raw, "05_Tissue_landscapes/Tissues",
                       tc, "short_gene_landscape_summary.tsv")
    if (!file.exists(lscp_f)) return(NULL)
    lscp     <- read_tsv(lscp_f, show_col_types = FALSE)
    high_genes <- lscp %>% filter(expr_class == "high") %>%
      pull(gene) %>% clean_ids()
    res <- run_topgo(high_genes, gene_universe_short, geneID2GO_short, "BP")
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res %>%
      filter(Term %in% term_freq$Term) %>%
      mutate(p_value = parse_topgo_p(weight01_fisher),
             neg_log10p = -log10(pmax(p_value, 1e-10)),
             tissue = tc) %>%
      dplyr::select(Term, neg_log10p, tissue)
  }

  go_matrix_long <- lapply(tissue_list, rescore_all) %>%
    bind_rows()

  
  go_matrix_wide <- go_matrix_long %>%
    pivot_wider(names_from = tissue, values_from = neg_log10p, values_fill = 0)

  mat_go <- go_matrix_wide %>%
    column_to_rownames("Term") %>%
    as.matrix()

  # Clamp and scale for display
  mat_go[mat_go > 10] <- 10

  pdf(file.path(outdir_gomap, "GO_tissue_comparison_heatmap.pdf"),
      width = 10, height = 0.35 * nrow(mat_go))
  pheatmap(mat_go,
           color             = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
           clustering_method = "ward.D2",
           fontsize_row      = 9,
           fontsize_col      = 10,
           main              = "GO term enrichment across tissues\n(high-expressed genes, -log10 p)",
           angle_col         = "45")
  dev.off()

  write_tsv(go_matrix_wide, file.path(outdir_gomap, "GO_matrix_all_tissues.tsv"))
  saveRDS(mat_go,file.path(outdir_gomap, "GO_matrix_all_tissues.RDS"))
}

outdir_traj <- file.path(outdir_raw, "05_Tissue_landscapes", "category_trajectories")
dir.create(outdir_traj, recursive = TRUE, showWarnings = FALSE)

# Each category is a named list of GO IDs, mapped to genes via geneID2GO_short.
gene2go_tbl <- tibble(
  gene  = rep(names(geneID2GO_short), sapply(geneID2GO_short, length)),
  go_id = unlist(geneID2GO_short)
)

# Category definitions: GO terms associated with key biological processes
category_terms <- list(
  "Ribosome/Translation"   = c("GO:0005840","GO:0006412","GO:0003735"),
  "Mitochondria/Respiration" = c("GO:0005739","GO:0006119","GO:0005747"),
  "Reproduction"           = c("GO:0000003","GO:0007276","GO:0048477"),
  "Immunity/Defence"       = c("GO:0006955","GO:0045087","GO:0098542"),
  "Cytoskeleton/Motility"  = c("GO:0005856","GO:0007018","GO:0003777"),
  "Proteostasis/UPS"       = c("GO:0006508","GO:0016567","GO:0006511"),
  "Chromatin/Gene regulation" = c("GO:0006325","GO:0006355","GO:0003682")
)

# Map each category to Ccap genes
category_genes <- lapply(category_terms, function(terms) {
  gene2go_tbl %>%
    filter(go_id %in% terms) %>%
    pull(gene) %>%
    unique()
})

# For each tissue, get mean logCPM for genes in each category
traj_list <- lapply(adult_tissues, function(tc) {
  ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
  ids <- intersect(ids, colnames(lcpm_bb))
  if (length(ids) == 0) return(NULL)
  tissue_lcpm <- lcpm_bb[, ids, drop = FALSE]
rownames(tissue_lcpm) <- gsub("gene-","", rownames(tissue_lcpm))
  
lapply(names(category_genes), function(cat) {
    cat_genes <- intersect(category_genes[[cat]], rownames(tissue_lcpm))
    if (length(cat_genes) == 0) return(NULL)
    tibble(
      tissue   = tissues[tc],
      category = cat,
      mean_logcpm     = mean(rowMeans(tissue_lcpm[cat_genes, , drop = FALSE])),
      median_logcpm   = median(rowMeans(tissue_lcpm[cat_genes, , drop = FALSE])),
      n_genes         = length(cat_genes)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

if (nrow(traj_list) > 0) {
  # Heatmap: categories × tissues (mean logCPM)
  cat_mat <- traj_list %>%
    dplyr::select(tissue, category, mean_logcpm) %>%
    pivot_wider(names_from = tissue, values_from = mean_logcpm, values_fill = NA) %>%
    column_to_rownames("category") %>%
    as.matrix()

  # Scaled so each category is comparable across tissues
  cat_mat_z <- t(scale(t(cat_mat)))

  pdf(file.path(outdir_traj, "category_expression_heatmap.pdf"),
      width = 8, height = 5)
  pheatmap(cat_mat_z,
           color    = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
           border_color = NA,
           clustering_method = "ward.D2",
           fontsize_row = 10, fontsize_col = 10,
           main     = "Functional category mean expression across tissues\n(row z-score)",
           angle_col = "45")
  dev.off()

  # Line plot: each category as a line across tissues
  tissue_order <- c("Testes","Accessory_Gland","Abdomen",
                    "Distal_Ovary","Proximal_Ovary")
  tissue_order <- tissue_order[tissue_order %in% traj_list$tissue]

  p_cat_lines <- traj_list %>%
    filter(tissue %in% tissue_order) %>%
    mutate(tissue = factor(tissue, levels = tissue_order)) %>%
    ggplot(aes(tissue, mean_logcpm,
               group = category, colour = category)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_colour_manual(values = okabe_ito_palette(length(category_genes))) +
    labs(title  = "Functional category expression profiles across tissues",
         x      = NULL,
         y      = "Mean log2CPM",
         colour = "Category") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right")
  p_cat_lines
  ggsave(file.path(outdir_traj, "category_line_plot.pdf"),
         p_cat_lines, width = 9, height = 6)

  write_tsv(traj_list,
            file.path(outdir_traj, "category_expression_by_tissue.tsv"))
}

# ===========================================================================
# TISSUE LANDSCAPE SUMMARY PANEL
# ---------------------------------------------------------------------------
# Modelled on Graveley et al. 2011 (modENCODE) and Yanai et al. 2005:
# the canonical tissue atlas summary combines breadth of expression,
# expression class distribution, tissue uniqueness, and the tau distribution
# into one figure.
#
# Panel layout (2×2 patchwork):
#   A. Expression breadth histogram — how many tissues each gene is expressed in
#      (Graveley et al. 2011 Fig 2; GTEx pilot Fig 2)
#   B. Stacked bar — expression class (high / medium / low) per tissue
#      (modENCODE; Cardoso-Moreira et al. 2019 Nature)
#   C. Transcriptional uniqueness dot — total expressed genes vs % private
#      (GTEx; Brawand et al. 2011 Nature)
#   D. Tau distribution — tissue-specificity score histogram
#      (Yanai et al. 2005 Bioinformatics; Kryuchkova-Mostacci & Robinson-Rechavi 2017)
# ===========================================================================

outdir_landscape_summary <- file.path(outdir_raw, "05_Tissue_landscapes", "summary")
dir.create(outdir_landscape_summary, recursive = TRUE, showWarnings = FALSE)

adult_tissues_ls <- c("AG", "AB", "TE", "YO", "OO")
tissue_labels_ls  <- c(
  AG = "Accessory\nGland", AB = "Abdomen",
  TE = "Testes",  YO = "Distal\nOvary", OO = "Proximal\nOvary"
)

# ---- A: Breadth of expression histogram ------------------------------------
# Reload tissue_mean_lcpm if not in environment (landscape section may have run first)
if (!exists("tissue_mean_lcpm")) {
  tissue_mean_lcpm <- lapply(adult_tissues_ls, function(tc) {
    ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
    ids <- intersect(ids, colnames(counts_short_raw))
    if (length(ids) == 0) return(NULL)
    dge  <- DGEList(counts = counts_short_raw[, ids, drop = FALSE])
    keep <- rowSums(cpm(dge) >= 1) >= 2
    dge  <- dge[keep, , keep.lib.sizes = FALSE]
    dge  <- calcNormFactors(dge)
    lcpm <- cpm(dge, log = TRUE, prior.count = 2)
    tibble(gene = rownames(lcpm), mean_logcpm = rowMeans(lcpm),
           tissue = tissue_labels_ls[tc])
  }) %>% bind_rows()
}

breadth_df <- tissue_mean_lcpm %>%
  mutate(expressed = mean_logcpm > 0) %>%
  group_by(gene) %>%
  summarise(n_tissues = sum(expressed), .groups = "drop")

n_tiss <- length(adult_tissues_ls)

breadth_counts <- breadth_df %>%
  dplyr::count(n_tissues) %>%
  mutate(
    label    = as.character(n_tissues),
    category = case_when(
      n_tissues == 0 ~ "None",
      n_tissues == 1 ~ "Single-tissue",
      n_tissues == n_tiss ~ "Ubiquitous",
      TRUE ~ "Multi-tissue"
    ),
    category = factor(category,
                      c("None","Single-tissue","Multi-tissue","Ubiquitous"))
  )

cat_cols <- c(
  "Single-tissue" = "#D55E00",
  "Multi-tissue"  = "#56B4E9" ,
  "Ubiquitous"    = "grey85"
)

pA_ls <- breadth_counts |>
  filter(category != "None") |>
  ggplot(
                aes(factor(n_tissues), n, fill = category)) +
  geom_col(colour = "white", width = 0.75) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = cat_cols) +
  scale_x_discrete(
    labels = c("0" = "0\n(not expr.)", as.character(seq_len(n_tiss)))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "No. tissues expressed in",
       y = "Gene count",
       fill = NULL,
       title = "A - Expression breadth") +
  plot_theme +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))
ggsave(file.path(outdir_landscape_summary, "PA_ls.pdf"),
       pA_ls, width = 6, height = 6)

# ---- B: Expression class stacked bar per tissue ----------------------------
# Load per-tissue landscape summaries
lscp_files <- file.path(outdir_raw, "05_Tissue_landscapes/Tissues",
                         tissues[adult_tissues_ls],
                         "short_gene_landscape_summary.tsv")
names(lscp_files) <- adult_tissues_ls

class_df <- lapply(adult_tissues_ls, function(tc) {
  f <- lscp_files[tc]
  if (!file.exists(f)) return(NULL)
  read_tsv(f, show_col_types = FALSE) %>%
    dplyr::count(expr_class) %>%
    mutate(tissue = tissue_labels_ls[tc])
}) %>%
  bind_rows() %>%
  filter(!is.na(expr_class)) %>%
  mutate(expr_class = factor(expr_class, c("low","medium","high")))

class_cols <- c(
  low    = "grey85",
  medium = "#56B4E9" ,
  high   = "#D55E00"
)

pB_ls <- if (nrow(class_df) > 0) {
  class_df %>%
    group_by(tissue) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    ggplot(aes(tissue, pct, fill = expr_class)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = class_cols,
                      labels = c(low="Low",medium="Medium",high="High")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(x = NULL, y = "% of genes", fill = "Expression class",
         title = "Expression class per tissue") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
} else {
  ggplot() + annotate("text",x=0.5,y=0.5,
                      label="Run landscape section first",size=4) +
    theme_void() + labs(title="Expression class per tissue")
}
ggsave(file.path(outdir_landscape_summary, "PB_ls.pdf"),
       pB_ls, width = 6, height = 6)
# ---- C: Transcriptional uniqueness dot plot --------------------------------
sp_file <- file.path(outdir_raw, "05_Tissue_landscapes",
                     "shared_private", "shared_private_summary.tsv")

pC_ls <- if (file.exists(sp_file)) {
  sp_df <- read_tsv(sp_file, show_col_types = FALSE)

  tissue_pal_ls <- setNames(
    okabe_ito_palette(nrow(sp_df)),
    sp_df$organ
  )

  ggplot(sp_df, aes(total_expressed, pct_private,
                    colour = organ, label = organ)) +
    geom_point(size = 5, alpha = 0.9) +

        geom_label_repel(
  fill = "white",
  label.color = "darkgrey",  
  size = 3
) +
    scale_colour_manual(values = tissue_pal_ls) +
    labs(x = "Total expressed genes",
         y = "% tissue-private genes",
         colour = NULL,
         title = "B - Transcriptional uniqueness") +
    plot_theme +
    theme(legend.position = "none")
} else {
  ggplot() + annotate("text",x=0.5,y=0.5,
                      label="Run shared_private section first",size=4) +
    theme_void() + labs(title="B - Transcriptional uniqueness")
}
ggsave(file.path(outdir_landscape_summary, "PC_ls.pdf"),
       pC_ls, width = 6, height = 6)

# ---- D: Tau distribution ---------------------------------------------------
tau_file <- file.path(outdir_raw, "07_Tau_specificity", "tau_scores.tsv")

pD_ls <- if (file.exists(tau_file)) {
  tau_df_s <- read_tsv(tau_file, show_col_types = FALSE)

  # Proportion with tau > 0.85 (tissue-specific threshold from Yanai 2005)
  pct_ts <- round(mean(tau_df_s$tau > 0.85, na.rm = TRUE) * 100, 1)

  ggplot(tau_df_s, aes(tau)) +
    geom_histogram(bins = 60,
                   fill = "#0072B2",
                   colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.85,
               linetype = "dashed",
               colour = "#D55E00", linewidth = 0.8) +
    annotate("text", x = 0.87, y = Inf, vjust = 1.5, hjust = 0,
             label = paste0("tau > 0.85\n(tissue-specific)\n",
                            pct_ts, "% of genes"),
             size = 3.2, colour = "#D55E00") +
    labs(x = "Tau (0 = ubiquitous,  1 = tissue-specific)",
         y = "Gene count",
         title = "Tau tissue-specificity distribution") +
    plot_theme
} else {
  ggplot() + annotate("text",x=0.5,y=0.5,
                      label="Run tau section first",size=4) +
    theme_void() + labs(title="D  Tau distribution")
}
ggsave(file.path(outdir_landscape_summary, "PD_ls.pdf"),
       pD_ls, width = 6, height = 6)

# ---- Assemble --------------------------------------------------------------
p_landscape_summary <- (pA_ls | pB_ls) / (pC_ls | pD_ls) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, colour = "grey40")
    )
  )

ggsave(file.path(outdir_landscape_summary, "Tissue_landscape_summary.pdf"),
       p_landscape_summary, width = 13, height = 10)
ggsave(file.path(outdir_landscape_summary, "Tissue_landscape_summary.png"),
       p_landscape_summary, width = 13, height = 10, dpi = 200)

p_landscape_summary

## Another final figure for tissues 

mat_go <- readRDS(file.path(outdir_raw, "05_Tissue_landscapes/GO_comparison/GO_matrix_all_tissues.RDS"))
  

 pC2_ls <-  pheatmap(mat_go,
           color             = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
           clustering_method = "ward.D2",
           fontsize_row      = 9,
           fontsize_col      = 10,
           main              = "C - GO term enrichment across tissues\n(high-expressed genes, -log10 p)",
           angle_col         = "45")
    pC2_ls <- wrap_elements(
  full = grid.grabExpr(ComplexHeatmap::draw( pC2_ls)))
  
p_landscape_summary2 <- ((pA_ls /pC_ls) |  pC2_ls )+
    plot_annotation(
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, colour = "grey40")
    )
  ) +
  plot_layout(widths = c(0.5,1))

ggsave(file.path(outdir_landscape_summary, "Tissue_landscape_summary2.pdf"),
       p_landscape_summary2, width = 13, height = 9)
ggsave(file.path(outdir_landscape_summary, "Tissue_landscape_summary2.png"),
       p_landscape_summary2, width = 13, height = 9, dpi = 200)

p_landscape_summary2

outdir_gomap <- file.path(outdir_raw, "05_Tissue_landscapes", "GO_comparison")
dir.create(outdir_gomap, recursive = TRUE, showWarnings = FALSE)

# Read the merged GO summary (built at end of build_tissue_landscape loop)
base_dir <- paste0(outdir_raw, "/05_Tissue_landscapes/Tissues")
merged_go <- read_tsv(
  file.path(base_dir, "merged_all_tissues_GO_summary.tsv"),
  show_col_types = FALSE
)

# Clean up p-values
merged_go <- merged_go %>%
  mutate(
    p_value    = parse_topgo_p(weight01_fisher),
    neg_log10p = -log10(pmax(p_value, 1e-10))
  )

# ============================================================
# PLOT 1: Heatmap — high expressed genes across all tissues
#   Rows = GO terms, Columns = tissue
#   Fill = -log10(p)
#   Shows which functional programmes are shared vs tissue-specific
# ============================================================

make_go_heatmap <- function(go_data, cls, title_text,
                             n_terms = 20, min_tissues = 2) {

  df <- go_data %>% filter(classification == cls)
  if (nrow(df) == 0) return(NULL)

  # Keep terms that appear in at least min_tissues tissues
  term_freq <- df %>%
    filter(p_value < 0.05) %>%
    dplyr::count(Term) %>%
    filter(n >= min_tissues) %>%
    pull(Term)

  if (length(term_freq) == 0) {
    message("No shared terms found for ", cls)
    return(NULL)
  }

  # For terms that appear in multiple tissues, take the most significant
  # p-value per tissue (some tissues have duplicates from both classes)
  mat_long <- df %>%
    filter(Term %in% term_freq) %>%
    group_by(Term, tissue) %>%
    slice_min(p_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(Term, tissue, neg_log10p)

  # Pivot to matrix: rows = terms, columns = tissues
  mat_wide <- mat_long %>%
    pivot_wider(names_from  = tissue,
                values_from = neg_log10p,
                values_fill = 0)

  # Select top n_terms by mean enrichment across tissues
  top_terms <- mat_wide %>%
    mutate(mean_p = rowMeans(across(-Term))) %>%
    arrange(desc(mean_p)) %>%
    slice_head(n = n_terms) %>%
    pull(Term)

  mat <- mat_wide %>%
    filter(Term %in% top_terms) %>%
    column_to_rownames("Term") %>%
    as.matrix()

  # Cap at 10 for display
  mat[mat > 10] <- 10

  # Cluster rows by similarity so related terms group together
  row_order <- hclust(dist(mat), method = "ward.D2")$order
  mat       <- mat[row_order, , drop = FALSE]

  # Clean up tissue names for display
  colnames(mat) <- str_replace_all(colnames(mat), "_", " ")

  pheatmap::pheatmap(
    mat,
    color             = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
    cluster_rows      = FALSE,   # already ordered above
    cluster_cols      = TRUE,
    clustering_method = "ward.D2",
    fontsize_row      = 9,
    fontsize_col      = 10,
    angle_col         = "45",
    main              = title_text,
    legend_breaks     = c(0, 2, 4, 6, 8, 10),
    legend_labels     = c("0", "2", "4", "6", "8", ">=10")
  )
}

# ---- High expressed genes -------------------------------------------------
p_high <- make_go_heatmap(
  merged_go,
  cls        = "high",
  title_text = "GO enrichment - high expressed genes\n(-log10 p, shared terms across >= 2 tissues)",
  n_terms    = 25,
  min_tissues = 2
)

if (!is.null(p_high)) {
  pdf(file.path(outdir_gomap, "GO_heatmap_high_expressed.pdf"),
      width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_high$gtable)
  dev.off()
  grid::grid.newpage()
  grid::grid.draw(p_high$gtable)
}

# ---- Medium expressed genes -----------------------------------------------
p_med <- make_go_heatmap(
  merged_go,
  cls        = "medium",
  title_text = "GO enrichment - medium expressed genes\n(-log10 p, shared terms across >= 2 tissues)",
  n_terms    = 25,
  min_tissues = 2
)

if (!is.null(p_med)) {
  pdf(file.path(outdir_gomap, "GO_heatmap_medium_expressed.pdf"),
      width = 10, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_med$gtable)
  dev.off()
  grid::grid.newpage()
  grid::grid.draw(p_med$gtable)
}

# ============================================================
# PLOT 2: Faceted bar chart — top GO terms per tissue
#   One panel per tissue, bars = top 8 terms, coloured by -log10(p)
#   Good for showing the tissue-specific signature at a glance
# ============================================================

top_per_tissue <- merged_go %>%
  filter(classification == "high",
         p_value < 0.05) %>%
  group_by(tissue) %>%
  slice_min(p_value, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    tissue     = str_replace_all(tissue, "_", " "),
    Term_short = str_trunc(Term, 35),
    Term_short = factor(Term_short,
                        levels = rev(unique(Term_short[order(tissue, p_value)])))
  )

p_facet <- ggplot(top_per_tissue,
                  aes(x = reorder_within(Term_short, neg_log10p, tissue),
                      y = neg_log10p,
                      fill = neg_log10p)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ tissue, scales = "free_y", ncol = 2) +
  scale_x_reordered() +
  scale_fill_gradientn(colours = enrichment_palette(100),
                       name    = expression(-log[10](p))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Top GO terms per tissue - high expressed genes",
    subtitle = "Top 8 significant BP terms per tissue (p < 0.05)",
    x = NULL, y = expression(-log[10](p))
  ) +
  plot_theme +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(size = 8),
    legend.position  = "right"
  )

p_facet
ggsave(file.path(outdir_gomap, "GO_faceted_per_tissue_high.pdf"),
       p_facet, width = 14, height = 12)

# ---- Same for medium -------------------------------------------------------
top_per_tissue_med <- merged_go %>%
  filter(classification == "medium",
         p_value < 0.05) %>%
  group_by(tissue) %>%
  slice_min(p_value, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    tissue     = str_replace_all(tissue, "_", " "),
    Term_short = str_trunc(Term, 35)
  )

p_facet_med <- ggplot(top_per_tissue_med,
                      aes(x = reorder_within(Term_short, neg_log10p, tissue),
                          y = neg_log10p,
                          fill = neg_log10p)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~ tissue, scales = "free_y", ncol = 2) +
  scale_x_reordered() +
  scale_fill_gradientn(colours = enrichment_palette(100),
                       name    = expression(-log[10](p))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Top GO terms per tissue - medium expressed genes",
    subtitle = "Top 8 significant BP terms per tissue (p < 0.05)",
    x = NULL, y = expression(-log[10](p))
  ) +
  plot_theme +
  theme(
    strip.text      = element_text(face = "bold", size = 10),
    axis.text.y     = element_text(size = 8),
    legend.position = "right"
  )

p_facet_med
ggsave(file.path(outdir_gomap, "GO_faceted_per_tissue_medium.pdf"),
       p_facet_med, width = 14, height = 12)
