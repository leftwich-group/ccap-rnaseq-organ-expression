# ===========================================================================
# R/setup.R
# Package loads, colour palettes, ggplot themes, and every helper function
# used across the analysis scripts. Sourced (via R/init.R) at the top of each
# script in scripts/. No data is read here.
# ===========================================================================

# ---- Libraries ------------------------------------------------------------
suppressPackageStartupMessages({
  # Core data wrangling
  library(tidyverse)
  library(here)
  library(tibble)

  # RNA-seq / DE
  library(edgeR)
  library(limma)
  library(DESeq2)          # used for VST in global atlas section
  library(matrixStats)

  # Dimensionality reduction
  library(umap)            # NEW: Section 10

  # Clustering & heatmaps
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggdendro)
  library(dynamicTreeCut)  # for MZT WGCNA-style modules
  library(mclust)          # Gaussian mixture for landscape classification

  # GO enrichment
  library(topGO)
  library(clusterProfiler)
  library(GO.db)

  # Visualisation
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(ggh4x)
  library(patchwork)
  library(forcats)
  library(EnhancedVolcano)
  library(UpSetR)
  library(gghalves)
  library(ggforce)
  library(ggdist)
  library(ggtext)
  library(cowplot)
  library(tidytext)

  # Plot statistics
  library(dplyr)
  library(rstatix)
  library(multcompView)
  # Tissue specificity
  # tau() implemented manually below — no extra package needed

  # Motif enrichment (Section 13)
  # Requires HOMER installed externally; R calls via system()
})

# ---- Colour palettes ----------------------------------------------------
# ---------------------------------------------------------------------------
# Okabe-Ito colour-blind safe palette (primary palette throughout)
# ---------------------------------------------------------------------------
okabe_ito <- c(
  orange     = "#E69F00",
  sky_blue   = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  pink       = "#CC79A7",
  black      = "#000000"
)
okabe_ito_palette <- colorRampPalette(unname(okabe_ito))

# Diverging palette for enrichment bubble plots
enrichment_palette <- colorRampPalette(c(
  "#0072B2", "#56B4E9", "#FFFFFF", "#E69F00", "#D55E00"
))

# Per-variable colour maps used in annotation bars
annotation_colours <- list(
  sex = c(M = "#F0A3FF", F = "#A6761D"),
  mating_status = c(M = "#4A8D8A", V = "#C46E7A"),
  age = c(
    "0"  = "#E69F00", "1"  = "#56B4E9", "2"  = "#009E73",
    "3"  = "#F0E442", "4"  = "#0072B2", "8"  = "#D55E00",
    "12" = "#CC79A7", "24" = "#000000"
  ),
  tissue = c(
    Testes          = "#4A8D8A",
    Accessory_Gland = "#C46E7A",
    Abdomen_M       = "#D79A4A",
    Abdomen_F       = "#E69F00",
    Distal_Oocyte   = "#0072B2",
    Proximal_Oocyte = "#56B4E9",
    Egg             = "#009E73",
    Abdomen       = "#D79A4A",
    Distal_Ovary   = "#0072B2",
    Proximal_Ovary = "#56B4E9"
  )
)

# ---- Themes -------------------------------------------------------------
plot_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title       = element_text(size = 11, face = "bold"),
    axis.text        = element_text(size = 10),
    strip.text       = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    strip.background = element_rect(fill = "grey95"),
    plot.title       = element_text(size = 12, face = "bold"),
     panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor.y = element_line(colour = "grey92", linewidth = 0.2)
  )

volcano_theme <- theme_bw(base_size = 14) +
  theme(
    axis.title  = element_text(face = "bold"),
    plot.title  = element_text(face = "bold")
  )

# ---- Helper functions -------------------------------------------------------
# ---------------------------------------------------------------------------
# clean_ids  — strip genome-specific prefixes for consistent gene matching
# ---------------------------------------------------------------------------
clean_ids <- function(x, strip_prefix = TRUE) {
  x <- as.character(x)
  if (strip_prefix) {
    x <- sub("^gene-", "", x)
    x <- sub("^CCAP1982_", "", x)
  }
  trimws(x)
}

# ---------------------------------------------------------------------------
# run_reference_analysis
#   Full edgeR QL pipeline: filter to normalise to dispersion to fit to contrasts
#   Returns nested list: dge, fit, logcpm, pca, contrasts[[name]]
# ---------------------------------------------------------------------------
run_reference_analysis <- function(count_matrix,
                                   metadata,
                                   design,
                                   contrast_list,
                                   prefix,
                                   outdir,
                                   plot_grouping_value = "group",
                                   min_count   = 10,
                                   min_samples = 2,
                                   fdr_thresh  = 0.05,
                                   lfc_thresh  = 1,
                                   label_top_n = 15,
                                   positive_label = "Up",
                                   negative_label = "Down",
                                   top_n = 100) {

  message("\n============================")
  message("Reference: ", prefix)
  message("============================")

  dge  <- DGEList(counts = count_matrix)
  keep <- rowSums(dge$counts >= min_count) >= min_samples
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  message("Genes retained after filtering: ", nrow(dge))

  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)

  logcpm <- cpm(dge, log = TRUE, prior.count = 2)

  # ---- PCA ----------------------------------------------------------------
  pca        <- prcomp(t(logcpm), scale. = TRUE)
  var_exp    <- pca$sdev^2 / sum(pca$sdev^2)
  pc1v       <- round(var_exp[1] * 100, 1)
  pc2v       <- round(var_exp[2] * 100, 1)

  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("id") %>%
    left_join(metadata, by = "id")

  grp_vals <- na.omit(unique(pca_df[[plot_grouping_value]]))
  pal      <- setNames(okabe_ito_palette(length(grp_vals)), grp_vals)

  p_pca <- ggplot(pca_df, aes(PC1, PC2,
                               colour = .data[[plot_grouping_value]],
                               label  = id)) +
    geom_point(size = 3) +
    geom_text_repel(size = 2.8, max.overlaps = 15) +
    scale_colour_manual(values = pal) +
    labs(x = paste0("PC1 (", pc1v, "%)"),
         y = paste0("PC2 (", pc2v, "%)"),
         colour = "") +
    plot_theme
  ggsave(file.path(outdir, paste0(prefix, "_PCA.pdf")), p_pca,
         width = 7, height = 5)

  # ---- MDS (legacy QC) ----------------------------------------------------
  pdf(file.path(outdir, paste0(prefix, "_MDS.pdf")), width = 7, height = 6)
  grp_factor <- factor(metadata[[plot_grouping_value]])
  plotMDS(logcpm,
          labels = metadata$id,
          col    = okabe_ito_palette(nlevels(grp_factor))[grp_factor])
  dev.off()

  # ---- Contrasts ----------------------------------------------------------
  contrast_results <- list()

  for (cname in names(contrast_list)) {
    message("  Contrast: ", cname)
    cmat <- makeContrasts(contrasts = contrast_list[[cname]], levels = design)
    qlf  <- glmQLFTest(fit, contrast = cmat)

    pos_label <- trimws(strsplit(contrast_list[[cname]], "-")[[1]][1])
    neg_label <- trimws(strsplit(contrast_list[[cname]], "-")[[1]][2])

    res <- topTags(qlf, n = Inf)$table %>%
      rownames_to_column("gene") %>%
      mutate(DE_status = case_when(
        FDR < fdr_thresh & logFC >=  lfc_thresh ~ pos_label,
        FDR < fdr_thresh & logFC <= -lfc_thresh ~ neg_label,
        TRUE ~ "NS"
      ))

    write_tsv(res,
              file.path(outdir, paste0(prefix, "_", cname, "_edgeR.tsv")))

    # MA plot
    pdf(file.path(outdir, paste0(prefix, "_", cname, "_MA.pdf")),
        width = 7, height = 5)
    status <- limma::decideTests(qlf, p.value = fdr_thresh)
    plotMD(qlf, status = status, values = c(-1, 1),
           col  = c("#0072B2", "#D55E00"),
           main = paste(prefix, cname))
    abline(h = c(-lfc_thresh, lfc_thresh), lty = 2)
    dev.off()

        # Volcano
      top_up <- res %>%
    filter(FDR < fdr_thresh, logFC > lfc_thresh) %>%
    arrange(FDR) %>%
    slice_head(n = label_top_n)
  
  top_down <- res %>%
    filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>%
    arrange(FDR) %>%
    slice_head(n = label_top_n)
  
  top_labels <- bind_rows(top_up, top_down)
  

    p_vol <- ggplot(res, aes(x = logFC, y = -log10(FDR), colour = DE_status)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = 2) +
    geom_hline(yintercept = -log10(fdr_thresh), linetype = 2) +
    geom_label_repel(
      data = top_labels,
      aes(label = gene),
      size = 2,
      max.overlaps = 50,
      fill = "white",
      colour = "black",
      label.size = 0.2,
      box.padding = 0.35,
      point.padding = 0.2
    ) +
    volcano_theme +
    scale_colour_manual(
  values = c(
    "NS" = "grey70", # forces NS to be grey 
    setNames(c("#D55E00", "#0072B2"),
             setdiff(unique(res$DE_status), "NS"))
  )
)+
    labs(
      x = "log2 fold change",
      y = "-log10(FDR)",
      colour = ""
    )
    ggsave(file.path(outdir, paste0(prefix, "_", cname, "_volcano.pdf")),
           p_vol, width = 8, height = 7)


      # Top genes from each direction
  top_up <- res %>%
 filter(FDR < fdr_thresh, abs(logFC) >= lfc_thresh)%>%
    arrange(FDR) %>%
    slice_head(n = top_n)
  
  top_down <- res %>%
   filter(FDR < fdr_thresh, abs(logFC) >= lfc_thresh)%>%
    arrange(FDR) %>%
    slice_head(n = top_n)
  
  sig_genes <- bind_rows(top_up, top_down) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
    sig_genes <- res %>%
      mutate(p_n = case_when(
        logFC >=  lfc_thresh ~ neg_label,
        logFC <= -lfc_thresh ~ pos_label,
        TRUE ~ "NS"
      )) %>%
      group_by(p_n) %>%
      filter(FDR < fdr_thresh, abs(logFC) >= lfc_thresh) %>%
      arrange(FDR) %>%
      slice_head(n = 100) %>%
      pull(gene)

    if (length(sig_genes) >= 2) {
      mat_hm <- logcpm[sig_genes[sig_genes %in% rownames(logcpm)], ,
                       drop = FALSE]
      mat_hm <- t(scale(t(mat_hm)))

      
if (!plot_grouping_value %in% names(annotation_colours)) {
  grouping_levels <- unique(metadata[[plot_grouping_value]])
  grouping_levels <- grouping_levels[!is.na(grouping_levels)]
  annotation_colours[[plot_grouping_value]] <-
    setNames(okabe_ito_palette(length(grouping_levels)), grouping_levels)
}

ha <- HeatmapAnnotation(
  df = metadata[, plot_grouping_value, drop = FALSE],
  col = annotation_colours[plot_grouping_value],
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = NULL
  ),
  na_col = "grey70"
)

ht <- ComplexHeatmap::Heatmap(
  mat_hm,
  name = "Z-score",
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  column_title = NULL,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#56B4E9", "white", "#E69F00"))
)

pdf(file.path(outdir,
              paste0(prefix, "_", cname, "_heatmap.pdf")),
    width = 8, height = 10)
draw(ht)
dev.off()

  
    }

    summary_tbl <- tibble(
      contrast   = cname,
      reference  = prefix,
      n_tested   = nrow(res),
      n_up       = sum(res$DE_status == pos_label),
      n_down     = sum(res$DE_status == neg_label),
      n_sig_total = sum(res$FDR < fdr_thresh)
    )

    contrast_results[[cname]] <- list(
      contrast_string = contrast_list[[cname]],
      qlf     = qlf,
      table   = res,
      summary = summary_tbl,
      volcano = p_vol
    )
  }

  list(dge = dge, design = design, fit = fit, logcpm = logcpm,
       pca = pca, contrasts = contrast_results)
}

#---------------
# compare_references_per_contrast
#   Short vs long logFC scatter + concordance metrics
# ---------------------------------------------------------------------------
compare_references_per_contrast <- function(short_obj, long_obj,
                                            contrast_name, outdir,
                                            gene_map_file = NULL,
                                            fdr_thresh    = 0.05) {
  short_res <- short_obj$contrasts[[contrast_name]]$table
  long_res  <- long_obj$contrasts[[contrast_name]]$table

  if (!is.null(gene_map_file) && file.exists(gene_map_file)) {
    gene_map <- read_tsv(gene_map_file, show_col_types = FALSE)
    short_cmp <- rename(short_res, short_gene_id = gene,
                        logFC_short = logFC, FDR_short = FDR,
                        logCPM_short = logCPM, DE_status_short = DE_status)
    long_cmp  <- rename(long_res,  long_gene_name = gene,
                        logFC_long = logFC, FDR_long = FDR,
                        logCPM_long = logCPM, DE_status_long = DE_status)
    compare_df <- gene_map %>%
      inner_join(short_cmp, by = "short_gene_id") %>%
      inner_join(long_cmp,  by = "long_gene_name")
  } else {
    shared <- intersect(short_res$gene, long_res$gene)
    short_cmp <- short_res %>% filter(gene %in% shared) %>%
      rename(gene_shared = gene, logFC_short = logFC, FDR_short = FDR,
             logCPM_short = logCPM, DE_status_short = DE_status)
    long_cmp  <- long_res  %>% filter(gene %in% shared) %>%
      rename(gene_shared = gene, logFC_long  = logFC, FDR_long  = FDR,
             logCPM_long  = logCPM, DE_status_long  = DE_status)
    compare_df <- inner_join(short_cmp, long_cmp, by = "gene_shared")
  }

  if (nrow(compare_df) == 0) return(NULL)

  r_pearson  <- cor(compare_df$logFC_short, compare_df$logFC_long)
  r_spearman <- cor(compare_df$logFC_short, compare_df$logFC_long,
                    method = "spearman")

  p <- ggplot(compare_df, aes(logFC_short, logFC_long)) +
    geom_point(alpha = 0.4, size = 1.2, colour = "grey40") +
    geom_smooth(method = "lm", se = FALSE,
                colour = "#D55E00") +
    geom_abline(linetype = "dashed", colour = "black") +
    labs(title   = paste(contrast_name, "- short vs long"),
         subtitle = paste0("Pearson r = ", round(r_pearson, 3),
                           "  |  Spearman rho = ", round(r_spearman, 3)),
         x = "log2FC (short)", y = "log2FC (long)") +
    plot_theme
  ggsave(file.path(outdir,
                   paste0("short_vs_long_", contrast_name, "_logFC.pdf")),
         p, width = 6, height = 5)

  summary_tbl <- tibble(
    contrast           = contrast_name,
    n_compared         = nrow(compare_df),
    pearson_logFC      = r_pearson,
    spearman_logFC     = r_spearman,
    short_sig          = sum(compare_df$FDR_short < fdr_thresh),
    long_sig           = sum(compare_df$FDR_long  < fdr_thresh),
    concordant_sig     = sum(compare_df$FDR_short < fdr_thresh &
                               compare_df$FDR_long  < fdr_thresh)
  )

  list(table = compare_df, summary = summary_tbl, scatter = p)
}


# ---------------------------------------------------------------------------
# GO enrichment helpers  (topGO wrapper + clusterProfiler wrapper)
# ---------------------------------------------------------------------------

# Parse topGO p-value strings
parse_topgo_p <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("< 1e-30", "<1e-30")] <- "1e-30"
  suppressWarnings(as.numeric(x))
}

# Run topGO for one gene set / one ontology.
# The result is cached to 04_results/.go_cache/, keyed by the exact inputs
# (gene set + genome universe + ontology). Re-running a section therefore
# reuses earlier GO runs instead of recomputing them -- this is what makes a
# re-run after a downstream failure fast. Set env var GO_FORCE=1 (or FORCE=1)
# to ignore the cache, or delete 04_results/.go_cache/ to clear it (do this
# after rebuilding the GO annotations).
run_topgo <- function(target_genes, gene_universe, geneID2GO,
                      ontology = "BP") {
  target_genes <- intersect(clean_ids(target_genes), gene_universe)
  if (length(target_genes) < 2) return(NULL)

  .go_force <- nzchar(Sys.getenv("GO_FORCE")) || nzchar(Sys.getenv("FORCE"))
  .cache_dir <- here::here("04_results", ".go_cache")
  .key <- rlang::hash(list(sort(target_genes), names(geneID2GO), ontology))
  .cache_f <- file.path(.cache_dir, paste0("topgo_", .key, ".rds"))
  if (!.go_force && file.exists(.cache_f)) return(readRDS(.cache_f))

  gene_list <- factor(as.integer(gene_universe %in% target_genes),
                      levels = c(0, 1))
  names(gene_list) <- gene_universe

  GOdata <- new("topGOdata", ontology = ontology, allGenes = gene_list,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  out <- tibble::as_tibble(GenTable(GOdata, weight01_fisher = result,
                                    orderBy = "weight01_fisher", topNodes = 30))
  dir.create(.cache_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, .cache_f)
  out
}

# Run all three ontologies
run_topgo_all_ont <- function(target_genes, gene_universe, geneID2GO) {
  list(
    BP = run_topgo(target_genes, gene_universe, geneID2GO, "BP"),
    MF = run_topgo(target_genes, gene_universe, geneID2GO, "MF"),
    CC = run_topgo(target_genes, gene_universe, geneID2GO, "CC")
  )
}

# Plot one ontology result as a bubble plot
plot_one_topgo <- function(topgo_df, title_text,
                           n_show = 8, wrap_width = 40) {
  if (is.null(topgo_df) || nrow(topgo_df) == 0) return(NULL)
  df <- topgo_df %>%
    mutate(p_value = parse_topgo_p(weight01_fisher),
           Significant = as.numeric(Significant)) %>%
    filter(!is.na(p_value), !is.na(Significant)) %>%
    arrange(p_value) %>%
    slice_head(n = n_show) %>%
    mutate(Description = str_wrap(Term, wrap_width),
           Description = fct_reorder(Description, Significant),
           logp = -log10(p_value))

  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(Description, Significant, size = Significant, fill = logp)) +
    geom_point(shape = 21, alpha = 0.9, colour = "black") +
    coord_flip() +
    labs(title = title_text, x = NULL, y = "Gene count",
         fill = expression(-log[10](p))) +
    scale_fill_gradientn(colours = enrichment_palette(100)) +
    plot_theme
  
  
}

# plotting summary go term barchart
plot_topgo_bar <- function(topgo_df,  title_text = "GO term summary",
                           n_show = 12, wrap_width = 45) {
  
  topgo_df <- topgo_df$BP
  if (is.null(topgo_df) || nrow(topgo_df) == 0) return(NULL)
 
  df <- topgo_df %>%
    mutate(p_value     = parse_topgo_p(weight01_fisher),
           Significant = as.numeric(Significant)) %>%
    filter(!is.na(p_value), !is.na(Significant)) %>%
    arrange(p_value) %>%
    slice_head(n = n_show) %>%
    mutate(Term = str_wrap(Term, wrap_width),
           Term = fct_reorder(Term, Significant),
           logp = -log10(p_value))
 
  if (nrow(df) == 0) return(NULL)
 

      
  ggplot(df, aes(Term, Significant, fill = logp)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_gradientn(colours = enrichment_palette(100),
                         name = expression(-log[10](p))) +
    labs(title = title_text, x = NULL, y = "Significant genes") +
    plot_theme +
    theme(
      # slightly more compact for the summary role
      plot.title      = element_text(size = rel(0.95)),
      axis.text.y     = element_text(size = rel(0.85)),
      legend.key.height = unit(0.5, "cm")
    )
}

# Patchwork panel: BP / MF / CC stacked
# `letters` (optional, length-3, one per BP/MF/CC) bakes a panel label into
# each sub-panel's own title ("A - <title> — BP", ...) rather than a corner
# tag, so the label stays attached when combined with other letter-less panels.
plot_topgo_panel <- function(ont_list, panel_title,
                             n_show = 8, wrap_width = 40,
                             letters = c(NA, NA, NA)) {
  onts <- c("BP", "MF", "CC")
  plots <- lapply(seq_along(onts), function(i) {
    ont <- onts[i]
    title_text <- paste0(panel_title, " - ", ont)
    if (!is.na(letters[i])) title_text <- paste0(letters[i], " - ", title_text)
    plot_one_topgo(ont_list[[ont]], title_text, n_show, wrap_width)
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(NULL)
  wrap_plots(plots, ncol = 1)
}

# Up / Down side-by-side and saved to file. Letters run in on-page reading
# order (left-to-right within each ontology row): A/B = BP Up/Down,
# C/D = MF Up/Down, E/F = CC Up/Down.
save_go_panel <- function(up_list, down_list, outdir, label,
                          n_show = 8) {
  p_up   <- plot_topgo_panel(up_list,   paste(label, "Up"),   n_show,
                             letters = c("A", "C", "E"))
  p_down <- plot_topgo_panel(down_list, paste(label, "Down"), n_show,
                             letters = c("B", "D", "F"))
  plots  <- Filter(Negate(is.null), list(p_up, p_down))
  if (length(plots) == 0) return(invisible(NULL))
  p_combined <- wrap_plots(plots, ncol = length(plots))
  ggsave(file.path(outdir, paste0("GO_", label, "_combined.pdf")),
         p_combined, width = 9 * length(plots), height = 14,
         limitsize = FALSE)
}

# ---------------------------------------------------------------------------
# Tau tissue-specificity index (Yanai et al. 2005)
#   tau = 0  to  ubiquitously expressed
#   tau = 1  to  tissue-specific
# ---------------------------------------------------------------------------
calc_tau <- function(expr_matrix) {
  # expr_matrix: genes × tissues (mean CPM or TPM per tissue)
  apply(expr_matrix, 1, function(x) {
    x[x < 0] <- 0
    if (max(x) == 0) return(NA_real_)
    n <- length(x)
    (sum(1 - x / max(x))) / (n - 1)
  })
}

# ---------------------------------------------------------------------------
# Summary-panel helpers, shared by the Q2 / Q3 / Q4 summary sections
# (defined once here; in the monolithic document they lived in the Q2 chunk).
# ---------------------------------------------------------------------------

# Load one short-genome contrast table from a saved results_all_nested.rds
load_short_res <- function(rds_path, contrast_name) {
  if (!file.exists(rds_path)) return(NULL)
  readRDS(rds_path)$short$contrasts[[contrast_name]]$table
}

# Volcano plot for a summary panel; NS greyed, top 8 hits per direction labelled
make_summary_volcano <- function(res, title_text,
                                 pos_lab = "Female-biased",
                                 neg_lab = "Male-biased") {
  if (is.null(res)) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5, label="Data not found", size=4) +
             theme_void() + labs(title = title_text))
  }

  # Top labels by FDR in each direction
  labels_df <- bind_rows(
    res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>%
      arrange(FDR) %>% slice_head(n = 8),
    res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>%
      arrange(FDR) %>% slice_head(n = 8)
  )

  res <- res %>%
    mutate(group = case_when(
      FDR < fdr_thresh & logFC >  lfc_thresh ~ pos_lab,
      FDR < fdr_thresh & logFC < -lfc_thresh ~ neg_lab,
      TRUE ~ "NS"
    ))

  n_pos  <- sum(res$group == pos_lab)
  n_neg  <- sum(res$group == neg_lab)

  ggplot(res, aes(logFC, -log10(FDR), colour = group)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = -log10(fdr_thresh),
               linetype = "dashed", colour = "grey60") +
    geom_label_repel(data = labels_df,
                     aes(label = gene),
                     size = 2.2, fill = "white",
                     colour = "black", label.size = 0.15,
                     max.overlaps = 20,
                     box.padding = 0.3) +
    scale_colour_manual(values = c(
      setNames("#D55E00", pos_lab),
      setNames("#0072B2",       neg_lab),
      "NS" = "grey80"
    )) +
    annotate(
  "label",
  x = max(res$logFC, na.rm = TRUE) * 0.85,
  y = max(-log10(res$FDR), na.rm = TRUE) * 0.3,
  label = paste0(pos_lab, "\nn = ", scales::comma(n_pos)),
  size = 3,
  colour = "#D55E00",
  fill = "grey95",
  label.size = 0.15,
  hjust = 1
) +
annotate(
  "label",
  x = min(res$logFC, na.rm = TRUE) * 0.85,
  y = max(-log10(res$FDR), na.rm = TRUE) * 0.3,
  label = paste0(neg_lab, "\nn = ", scales::comma(n_neg)),
  size = 3,
  colour = "#0072B2",
  fill = "grey95",
  label.size = 0.15,
  hjust = 0
) +
    labs(x = "log2FC", y = "-log10(FDR)",
         colour = NULL, title = title_text) +
    plot_theme +
    theme(legend.position = "none")
}

# ---- Shared constants -----------------------------------------------------
# Tissue code -> label. EG (egg) is adult-excluded in most cross-tissue views.
tissues <- c(AG = "Accessory_Gland",
             AB = "Abdomen",
             TE = "Testes",
             YO = "Distal_Ovary",
             OO = "Proximal_Ovary",
             EG = "Egg")

adult_tissues <- c("AG", "AB", "TE", "YO", "OO")

# ---- build_adult_lcpm --------------------------------------------------------
# Normalised log2-CPM matrix across the adult tissue samples. Used by the
# global atlas, tau, and cross-tissue correlation scripts, which in the
# original monolithic document shared one in-memory `lcpm_adult` object.
build_adult_lcpm <- function(counts_raw, meta,
                             tissue_codes = adult_tissues,
                             min_cpm = 1, min_samples = 2) {
  ids <- meta %>% dplyr::filter(tissue %in% tissue_codes) %>% dplyr::pull(id)
  ids <- intersect(ids, colnames(counts_raw))
  dge <- edgeR::DGEList(counts = counts_raw[, ids, drop = FALSE])
  keep <- rowSums(edgeR::cpm(dge) >= min_cpm) >= min_samples
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  edgeR::cpm(dge, log = TRUE, prior.count = 2)
}
