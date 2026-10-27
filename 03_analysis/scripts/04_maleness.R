#!/usr/bin/env Rscript
# Q4 - Transcriptional determinants of phenotypic maleness
# Output: 04_results/04_Maleness/   (reads 04_results/02_Sexual_dimorphism/)
# Run from the repository root:  Rscript 03_analysis/scripts/04_maleness.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir_q4 <- file.path(outdir_raw, "04_Maleness")
dir.create(outdir_q4, recursive = TRUE, showWarnings = FALSE)

male_contrasts <- list(
  C1 = list(tissues = c("TE", "AG"), var = "organ4",
            levels = c(TE = "Testes", AG = "Accessory_Gland"),
            contrast_list = list(AG_vs_TE = "Accessory_Gland - Testes")),
  C2 = list(tissues = c("TE", "YO"), var = "organ4",
            levels = c(TE = "Testes", YO = "Distal_Oocyte"),
            contrast_list = list(YO_vs_TE = "Distal_Oocyte - Testes")),
  C3 = list(tissues = c("TE", "OO"), var = "organ4",
            levels = c(TE = "Testes", OO = "Proximal_Oocyte"),
            contrast_list = list(OO_vs_TE = "Proximal_Oocyte - Testes")),
  C4 = list(tissues = c("AG", "AB"), var = "organ4",
            levels = c(AG = "Accessory_Gland", AB = "Abdomen"),
            contrast_list = list(AB_vs_AG = "Abdomen - Accessory_Gland"))
)

results_q4 <- list()

for (cid in names(male_contrasts)) {
  cfg    <- male_contrasts[[cid]]
  outdir <- file.path(outdir_q4, cid)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  meta <- meta_raw %>%
    filter(tissue %in% cfg$tissues) %>%
    mutate(organ4 = recode(tissue, !!!cfg$levels)) %>%
    filter(!is.na(organ4))

  counts_short <- counts_short_raw[, meta$id, drop = FALSE]
  counts_long  <- counts_long_raw[,  meta$id, drop = FALSE]

  design        <- model.matrix(~ 0 + organ4, data = meta)
  colnames(design) <- sub("organ4", "", colnames(design))

  results_q4[[cid]] <- list(
    short = run_reference_analysis(counts_short, meta, design,
                                   cfg$contrast_list, "short", outdir,
                                   "organ4", min_count, min_samples,
                                   fdr_thresh, lfc_thresh),
    long  = run_reference_analysis(counts_long,  meta, design,
                                   cfg$contrast_list, "long",  outdir,
                                   "organ4", min_count, min_samples,
                                   fdr_thresh, lfc_thresh)
  )
  saveRDS(results_q4[[cid]], file.path(outdir, "results_all_nested.rds"))

  # Record DE stats for this Q4 contrast
  cname_q4  <- names(cfg$contrast_list)[1]
  cstr_q4   <- cfg$contrast_list[[cname_q4]]
  q4_label  <- paste0("Maleness: ", sub("_vs_", " vs ", cname_q4))
  for (g in c("short","long"))
    append_de_stats(results_q4[[cid]], "Q4", q4_label, cstr_q4, genome = g)
}

for (cid in names(results_q4)) {
  outdir <- file.path(outdir_q4, cid)
  for (genome in c("short", "long")) {
    geneID2GO <- get(paste0("geneID2GO_",    genome))
    gene_univ <- get(paste0("gene_universe_", genome))
    for (cname in names(male_contrasts[[cid]]$contrast_list)) {
      res    <- results_q4[[cid]][[genome]]$contrasts[[cname]]$table
      up_g   <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
      down_g <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
      save_go_panel(
        run_topgo_all_ont(up_g,   gene_univ, geneID2GO),
        run_topgo_all_ont(down_g, gene_univ, geneID2GO),
        outdir, paste0(genome, "_", genome)
      )
      
       # outputGO 
  
  go_down <- run_topgo_all_ont(down_g, gene_univ, geneID2GO)
  go_up <- run_topgo_all_ont(up_g, gene_univ, geneID2GO)
  
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <-  go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)
  write_tsv(go_merged, paste0(outdir,"/",cname,"_",genome,"GO_terms_merged.tsv"))

    }
  }
}

# ---------------------------------------------------------------------------
# Q4 SUMMARY: Transcriptional determinants of phenotypic maleness
# Four-panel figure:
#   A. Mirror bar: DE gene counts across all Q4 contrasts
#   B. Volcano: Accessory gland vs Testes (defines male-specific repertoire)
#   C. Expression dot plot: top male-specific genes across male tissues
#   D. Overlap bar: male-biased genes — germline vs soma (links to Q2)
# ---------------------------------------------------------------------------

outdir_summary <- file.path(outdir_raw, "/04_Maleness/summary_panels")
dir.create(outdir_summary, recursive = TRUE, showWarnings = FALSE)

load_short_res <- function(rds_path, contrast_name) {
  if (!file.exists(rds_path)) return(NULL)
  readRDS(rds_path)$short$contrasts[[contrast_name]]$table
}

q4_contrast_meta <- tribble(
  ~label,                              ~rds_path,                         ~contrast_name,
  "Testes vs Accessory Gland ",                      "04_Maleness/C1/results_all_nested.rds", "AG_vs_TE",
  "Testes vs Distal oocyte ",           "04_Maleness/C2/results_all_nested.rds", "YO_vs_TE",
  "Testes vs Proximal oocyte",         "04_Maleness/C3/results_all_nested.rds", "OO_vs_TE",
  " Accessory Gland vs Abdomen",                     "04_Maleness/C4/results_all_nested.rds", "AB_vs_AG"
)

q4_res_list <- pmap(q4_contrast_meta,
                    function(label, rds_path, contrast_name, ...) {
  load_short_res(file.path(outdir_raw, rds_path), contrast_name)
}) %>% setNames(q4_contrast_meta$label)

# ---- A: Mirror bar --------------------------------------------------------
de_counts_q4 <- pmap_dfr(q4_contrast_meta,
                          function(label, rds_path, contrast_name, ...) {
  res <- load_short_res(file.path(outdir_raw, rds_path), contrast_name)
  if (is.null(res)) return(NULL)
  tibble(
    contrast = label,
    n_up   = sum(res$FDR < fdr_thresh & res$logFC >  lfc_thresh),
    n_down = sum(res$FDR < fdr_thresh & res$logFC < -lfc_thresh)
  )
})

pA_q4 <- if (nrow(de_counts_q4) > 0) {
  de_counts_q4 %>%
    mutate(contrast = factor(contrast, levels = rev(contrast))) %>%
    pivot_longer(c(n_up, n_down), names_to = "dir", values_to = "n") %>%
    mutate(n_plot  = if_else(dir == "n_down", -n, n),
           dir_lab = if_else(dir == "n_up",
                             "Up regulated", "Down regulated")) %>%
    ggplot(aes(contrast, n_plot, fill = dir_lab)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, colour = "grey40") +

    geom_richtext(
  aes(
    label = scales::comma(abs(n)),
    y = n_plot + sign(n_plot) * 80),
  fill = "white",
  label.color = "darkgrey",  
  size = 3
) + 
    scale_fill_manual(values = c(
      "Up regulated"  = "#D55E00",
      "Down regulated" = "#0072B2"
    )) +
    scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
    coord_flip() +
    labs(x = NULL, y = "DE gene count", fill = NULL,
         title = "A - Differential expression gene counts \n( FDR < 0.05, |log2FC| > 1)") +
       guides(fill="none")+
    plot_theme + theme(legend.position = "bottom")

    
} else {
  ggplot() + annotate("text",x=0.5,y=0.5,label="Run Q4 first",size=4) +
    theme_void() + labs(title="A - DE counts")
}
    ggsave(file.path(outdir_summary, "Q4_PA.png"),
       pA_q4, width = 5, height = 5, dpi = 200)
# ---- B: Volcano — AG vs Testes --------------------------------------------
pB_q4 <- make_summary_volcano(
  q4_res_list[["Testes vs Accessory Gland "]],
  "B - Testes vs Accessory Gland ",
  pos_lab = "AG-biased",
  neg_lab = "Testes-biased"
)
  ggsave(file.path(outdir_summary, "Q4_PB.png"),
       pB_q4, width = 5, height = 5, dpi = 200)
  
# ---- C: Multi-tissue expression dot plot — top testes-specific genes ------
# Genes most strongly testes-enriched (vs AG, vs distal oocyte)
res_ag_te  <- q4_res_list[["Testes vs Accessory Gland "]]
res_do_te  <- q4_res_list[["Testes vs Distal oocyte "]]

if (!is.null(res_ag_te) && !is.null(res_do_te)) {
  # Testes-enriched = negative logFC in both contrasts (Testes > other)
  testes_sig <- intersect(
    res_ag_te %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene),
    res_do_te %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene)
  )

  # Expression across male tissues
  meta_male <- meta_raw %>%
    filter(tissue %in% c("TE","AG","AB"), sex == "M")

  if (nrow(meta_male) >= 2 && length(testes_sig) >= 2) {
    counts_male <- counts_short_raw[, meta_male$id, drop = FALSE]
    dge_male    <- DGEList(counts = counts_male)
    keep_m      <- rowSums(cpm(dge_male) >= 1) >= 2
    dge_male    <- dge_male[keep_m, , keep.lib.sizes = FALSE]
    dge_male    <- calcNormFactors(dge_male)
    lcpm_male   <- cpm(dge_male, log = TRUE, prior.count = 2)

    top_te <- head(intersect(testes_sig, rownames(lcpm_male)), 25)

    dot_df <- lcpm_male[top_te, , drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "id", values_to = "logcpm") %>%
      left_join(meta_male %>% dplyr::select(id, tissue), by = "id") %>%
      mutate(tissue = recode(tissue,
                             TE = "Testes", AG = "Accessory\nGland",
                             AB = "Abdomen")) %>%
      group_by(gene, tissue) %>%
      summarise(mean_logcpm = mean(logcpm), .groups = "drop")

    pC_q4 <- ggplot(dot_df, aes(tissue, fct_reorder(gene, mean_logcpm),
                                  size = mean_logcpm,
                                  colour = mean_logcpm)) +
      geom_point() +
      scale_colour_gradientn(
        colours = enrichment_palette(3),
        name = "Mean\nlog2CPM"
      ) +
      scale_size_continuous(range = c(1, 6), guide = "none") +
      labs(x = NULL, y = NULL,
           title = "C - Testes-enriched genes\nacross male tissues") +
      plot_theme +
      theme(axis.text.y = element_text(size = 7))
  } else {
    pC_q4 <- ggplot() + annotate("text",x=0.5,y=0.5,
                                   label="Insufficient data",size=4) +
      theme_void() + labs(title="C - Testes-enriched dot plot")
  }
} else {
  pC_q4 <- ggplot() + annotate("text",x=0.5,y=0.5,
                                 label="Run Q4 C1 & C2 first",size=4) +
    theme_void() + labs(title="C - Testes-enriched dot plot")
}
  ggsave(file.path(outdir_summary, "Q4_PC.png"),
       pC_q4, width = 5, height = 5, dpi = 200)
  
# ---- D: Male-biased gene overlap: germline vs soma (bridge to Y-chr) ------
# Re-use Q2 germline and Q2 soma results already loaded
# NOTE: results_all_nested.rds objects are the full nested DE result lists
# from each contrast. They are regenerated by the DE pipeline above and are
# not tracked in this repo (each is tens of MB); regenerate them by running
# the corresponding contrast section, or point these paths at your own copy.
res_germline_q4 <- readRDS(here("04_results", "02_Sexual_dimorphism", "contrast_1", "results_all_nested.rds"))
res_germline_q4 <- res_germline_q4$short$contrasts$female_gonad_vs_male_gonad$table

res_soma_mated_q4 <- readRDS(here("04_results", "02_Sexual_dimorphism", "contrast_4", "results_all_nested.rds"))
res_soma_mated_q4 <- res_soma_mated_q4$short$contrasts$AbdomenF_mated_vs_AbdomenM$table

if (!is.null(res_germline_q4) && !is.null(res_soma_mated_q4)) {
  # Male-biased = lower in female contrast = logFC < 0
  male_gl   <- res_germline_q4  %>%
    filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene)
  male_soma <- res_soma_mated_q4 %>%
    filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene)

  shared_male <- intersect(male_gl, male_soma)

  overlap_df <- tibble(
    category = c("Germline\nonly", "Soma\nonly",
                 "Germline\n& soma"),
    n = c(
      length(setdiff(male_gl,   male_soma)),
      length(setdiff(male_soma, male_gl)),
      length(shared_male)
    ),
    fill_col = c("#0072B2" , "#D55E00", "#CC79A7")
  )

  pD_q4 <- ggplot(overlap_df,
                   aes(factor(category, levels  = c(
                     "Germline\nonly",
                     "Soma\nonly",
                     "Germline\n& soma"
                   )), n, fill = category)) +
    geom_col(width = 0.55) +
    geom_richtext(
  aes(
    label = scales::comma(n)),
  fill = "white",
  label.color = "darkgrey",  
  size = 3
) + 
    scale_fill_manual(values = setNames(overlap_df$fill_col,
                                        overlap_df$category)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Gene count",
         title = "D - Male-biased gene overlap\n(germline vs soma)") +
    plot_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9))
    ggsave(file.path(outdir_summary, "Q4_PD.png"),
       pD_q4, width = 5, height = 5, dpi = 200)
  
} else {
  pD_q4 <- ggplot() +
    annotate("text",x=0.5,y=0.5,label="Run Q2 contrasts first",size=4) +
    theme_void() + labs(title="D - Male-biased overlap")
}

# ---- Assemble ---------------------------------------------------------------
p_q4_summary <- (pA_q4 | pB_q4) / (pC_q4 | pD_q4) +
  plot_annotation(
    title    = "Transcriptional Determinants of Phenotypic Maleness",
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave(file.path(outdir_summary, "Q4_Maleness_summary.pdf"),
       p_q4_summary, width = 14, height = 11)
ggsave(file.path(outdir_summary, "Q4_Maleness_summary.png"),
       p_q4_summary, width = 14, height = 11, dpi = 200)

p_q4_summary

outdir_go_summary <- file.path(outdir_raw, "04_Maleness/summary_panels", "GO_heatmap")
dir.create(outdir_go_summary, recursive = TRUE, showWarnings = FALSE)

# ---- Load all short-genome GO summary files from q4 -----------------------
q4_go_files <- list.files(
  path       = file.path(outdir_raw, "04_Maleness"),
  pattern    = "GO_terms_merged\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
) %>%
  .[str_detect(., "short")]

message("q4 GO files found: ", length(q4_go_files))
print(q4_go_files)
# should be 2 because C2 has all 3 merged 
if (length(q4_go_files) == 0) stop("No q4 short GO files found - check paths")

# Read and label by contrast (parent folder name)
q4_go <- purrr::map_dfr(q4_go_files, function(f) {
  
  contrast <- basename(dirname(f))
  
  # Skip empty files
  if (file.info(f)$size == 0) {
    return(tibble())
  }
  
  dat <- tryCatch(
    read_tsv(f, show_col_types = FALSE),
    error = function(e) NULL
  )
  
  # Skip files that failed to read
  if (is.null(dat)) {
    return(tibble())
  }
  
  # Skip files with no rows
  if (nrow(dat) == 0) {
    return(tibble())
  }
  
  dat %>%
    mutate(
      weight01_fisher = as.character(weight01_fisher),
      contrast        = contrast,
      p_value         = parse_topgo_p(weight01_fisher),
      neg_log10p      = -log10(pmax(p_value, 1e-10)),
      contrast_group  = paste0(contrast, " | ", group)
    )
  
})

message("Total GO rows loaded: ", nrow(q4_go))
message("Contrast-groups found:")
print(unique(q4_go$contrast_group))

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
# PLOT 1: All q4 contrasts in one heatmap
#   Columns = contrast | group (e.g. "contrast_1 | unfert")
#   Shows which GO terms are shared across MZT contrasts
#   and which are specific to one comparison
# ============================================================
p_q4_heatmap <- make_contrast_go_heatmap(
  go_data    = q4_go,
  title_text = "MZT: GO enrichment across all contrasts\n(-log10 p)",
  n_terms    = 30,
  min_groups = 2
)
p_q4_heatmap<- wrap_elements(
  full = grid.grabExpr(ComplexHeatmap::draw( p_q4_heatmap)))
p_q4_heatmap
if (!is.null(p_q4_heatmap)) {

  ggsave(file.path(outdir_go_summary, "q4_GO_heatmap_all_contrasts.pdf"), p_q4_heatmap, width = 16, height = 10, dpi = 200)
  

}

# ============================================================
# PLOT 2: Faceted bar — top 8 terms per contrast-group
#   One panel per contrast-group
#   Good for seeing what each specific comparison is driven by
# ============================================================
p_q4_facet <- q4_go %>%
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

p_q4_facet
ggsave(file.path(outdir_go_summary, "q4_GO_faceted_per_contrast.pdf"),
       p_q4_facet, width = 14, height = 0.6 * length(unique(q4_go$contrast_group)) + 4,
       limitsize = FALSE)
