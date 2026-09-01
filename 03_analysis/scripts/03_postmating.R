#!/usr/bin/env Rscript
# Q3 - Post-mating transcriptional remodelling in females
# Output: 04_results/03_PostMating/
# Run from the repository root:  Rscript 03_analysis/scripts/03_postmating.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir_q3 <- file.path(outdir_raw, "03_PostMating")
dir.create(outdir_q3, recursive = TRUE, showWarnings = FALSE)

outdir <- file.path(outdir_q3, "contrast_1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue == "AB", sex == "F", experiment =="1") %>%
  mutate(mating_grp = if_else(mating_status == "M",
                              "Mated_abdomen", "Virgin_abdomen"))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long <- counts_long_raw[, meta$id, drop = FALSE]
design        <- model.matrix(~ 0 + mating_grp, data = meta)
colnames(design) <- sub("mating_grp", "", colnames(design))
contrast_list <- list(mated_vs_virgin_abdomen = "Mated_abdomen - Virgin_abdomen")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "mating_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
    long = run_reference_analysis(counts_long, meta, design, contrast_list,
                                 "long", outdir, "mating_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q3", "Female somvirgin vs mated",
                  "Mated_abdomen - Virgin_abdomen", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$mated_vs_virgin_abdomen$table
  up_g      <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_g    <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
  
  if (length(up_g) == 0 && length(down_g) == 0) {
    message("No significant DE genes for ", genome, " in C3.C1; skipping GO enrichment.")
    next
  }
  save_go_panel(run_topgo_all_ont(up_g, gene_univ, geneID2GO),
                run_topgo_all_ont(down_g, gene_univ, geneID2GO),
                outdir, paste0(genome))
  
 # outputGO 
  
  go_down <- run_topgo_all_ont(down_g, gene_univ, geneID2GO)
  go_up <- run_topgo_all_ont(up_g, gene_univ, geneID2GO)
  
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <-  go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)
  write_tsv(go_merged, paste0(outdir,"/",genome,"GO_terms_merged.tsv"))
}

outdir <- file.path(outdir_q3, "contrast_2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue == "OO", experiment =="1") %>%
  mutate(mating_grp = if_else(mating_status == "M",
                              "Mated_OO", "Virgin_OO"))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long <- counts_long_raw[, meta$id, drop = FALSE]
design        <- model.matrix(~ 0 + mating_grp, data = meta)
colnames(design) <- sub("mating_grp", "", colnames(design))
contrast_list <- list(mated_vs_virgin_OO = "Mated_OO - Virgin_OO")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "mating_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
    long = run_reference_analysis(counts_long, meta, design, contrast_list,
                                 "long", outdir, "mating_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q3", "Proximal OVary : virgin vs mated",
                  "Mated_OO - Virgin_OO", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$mated_vs_virgin_OO$table
  up_g      <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_g    <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
  save_go_panel(run_topgo_all_ont(up_g, gene_univ, geneID2GO),
                run_topgo_all_ont(down_g, gene_univ, geneID2GO),
                outdir, paste0(genome))
  
   # outputGO 
  
  go_down <- run_topgo_all_ont(down_g, gene_univ, geneID2GO)
  go_up <- run_topgo_all_ont(up_g, gene_univ, geneID2GO)
  
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <-  go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)
  write_tsv(go_merged, paste0(outdir,"/",genome,"GO_terms_merged.tsv"))
}

outdir <- file.path(outdir_q3, "contrast_3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue == "YO", experiment =="1") %>%
  mutate(mating_grp = if_else(mating_status == "M",
                              "Mated_YO", "Virgin_YO"))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long <- counts_long_raw[, meta$id, drop = FALSE]
design        <- model.matrix(~ 0 + mating_grp, data = meta)
colnames(design) <- sub("mating_grp", "", colnames(design))
contrast_list <- list(mated_vs_virgin_YO = "Mated_YO - Virgin_YO")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "mating_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
    long = run_reference_analysis(counts_long, meta, design, contrast_list,
                                 "long", outdir, "mating_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q3", "Distal ovary : virgin vs mated",
                  "Mated_YO - Virgin_YO", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$mated_vs_virgin_YO$table
  up_g      <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_g    <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
  save_go_panel(run_topgo_all_ont(up_g, gene_univ, geneID2GO),
                run_topgo_all_ont(down_g, gene_univ, geneID2GO),
                outdir, paste0(genome))
  
   # outputGO 
  
  go_down <- run_topgo_all_ont(down_g, gene_univ, geneID2GO)
  go_up <- run_topgo_all_ont(up_g, gene_univ, geneID2GO)
  
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <-  go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)
  write_tsv(go_merged, paste0(outdir,"/",genome,"GO_terms_merged.tsv"))
}

for (ms in c("V", "M")) {
  label  <- if (ms == "V") "virgin" else "mated"
  outdir <- file.path(outdir_q3, paste0("contrast_4/ovary_vs_abdomen_", label))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  meta <- meta_raw %>%
    filter(sex == "F", tissue %in% c("OO", "YO", "AB"),
           mating_status == ms, experiment =="1") %>%
    mutate(tissue_grp = case_when(
      tissue %in% c("OO", "YO") ~ "Ovary",
      tissue == "AB"             ~ "Abdomen",
      TRUE ~ NA_character_
    ))

  counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long <- counts_long_raw[, meta$id, drop = FALSE]
  design        <- model.matrix(~ 0 + tissue_grp, data = meta)
  colnames(design) <- sub("tissue_grp", "", colnames(design))
  contrast_list <- setNames(
    list("Ovary - Abdomen"),
    paste0("ovary_vs_abdomen_", label)
  )

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "tissue_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
    long = run_reference_analysis(counts_long, meta, design, contrast_list,
                                 "long", outdir, "tissue_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
  saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
  for (g in c("short","long"))
    append_de_stats(results_all, "Q3",
                    paste0(label," : abdomen vs ovary"),
                    "Ovary - Abdomen", genome = g)

  contrast <- names(contrast_list)
  ## GO terms
  for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts[[contrast]]$table
  up_g      <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_g    <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
  save_go_panel(run_topgo_all_ont(up_g, gene_univ, geneID2GO),
                run_topgo_all_ont(down_g, gene_univ, geneID2GO),
                outdir, paste0(genome))
  
   # outputGO 
  
  go_down <- run_topgo_all_ont(down_g, gene_univ, geneID2GO)
  go_up <- run_topgo_all_ont(up_g, gene_univ, geneID2GO)
  
  go_up_merged <- go_up %>% bind_rows() %>% mutate(group = "up")
  go_down_merged <-  go_down %>% bind_rows()%>% mutate(group = "down")
  
  go_merged <- rbind(go_up_merged, go_down_merged)
  write_tsv(go_merged, paste0(outdir,"/",genome,"GO_terms_merged.tsv"))

}
}

# ---------------------------------------------------------------------------
# Q3 SUMMARY: Post-mating transcriptional remodelling in females
# Four-panel figure:
#   A. Mirror bar: DE gene counts across mated vs virgin contrasts per tissue
#   B. Volcano: mated vs virgin distal ovary (the strongest signal)
#   C. Dot plot: top DE genes shared across tissues post-mating
#   D. Tissue × mating_status heatmap of top DE genes (distal ovary signal)
# ---------------------------------------------------------------------------


outdir_summary <- file.path(outdir_raw, "03_PostMating/summary_panels")
dir.create(outdir_summary, recursive = TRUE, showWarnings = FALSE)

# Helper already defined above
load_short_res <- function(rds_path, contrast_name) {
  if (!file.exists(rds_path)) return(NULL)
  readRDS(rds_path)$short$contrasts[[contrast_name]]$table
}

q3_contrast_meta <- tribble(
  ~label,              ~rds_path,                              ~contrast_name,
  "Female soma - \n Virgin vs Mated",           "03_PostMating/contrast_1/results_all_nested.rds", "mated_vs_virgin_abdomen",
  
  "Proximal Ovary - \n Virgin vs Mated",    "03_PostMating/contrast_2/results_all_nested.rds", "mated_vs_virgin_OO",
  
    "Distal Ovary - \n Virgin vs Mated",      "03_PostMating/contrast_3/results_all_nested.rds", "mated_vs_virgin_YO",
  
  "Virgin:\n Abdomen vs Ovary",   "03_PostMating/contrast_4/ovary_vs_abdomen_virgin/results_all_nested.rds", "ovary_vs_abdomen_virgin",
  
  "Mated:\n Abdomen vs Ovary",    "03_PostMating/contrast_4/ovary_vs_abdomen_mated/results_all_nested.rds", "ovary_vs_abdomen_mated",

)

q3_res_list <- pmap(q3_contrast_meta,
                    function(label, rds_path, contrast_name, ...) {
  load_short_res(file.path(outdir_raw, rds_path), contrast_name)
}) %>% setNames(q3_contrast_meta$label)

# ---- A: Mirror bar of DE counts -------------------------------------------
de_counts_q3 <- pmap_dfr(q3_contrast_meta,
                          function(label, rds_path, contrast_name, ...) {
  res <- load_short_res(file.path(outdir_raw, rds_path), contrast_name)
  if (is.null(res)) return(NULL)
  tibble(
    contrast = label,
    n_up   = sum(res$FDR < fdr_thresh & res$logFC >  lfc_thresh),
    n_down = sum(res$FDR < fdr_thresh & res$logFC < -lfc_thresh)
  )
})

pA_q3 <- if (nrow(de_counts_q3) > 0) {
  de_counts_q3 %>%
    mutate(contrast = factor(contrast, levels = rev(contrast))) %>%
    pivot_longer(c(n_up, n_down), names_to = "dir", values_to = "n") %>%
    mutate(n_plot  = if_else(dir == "n_down", -n, n),
           dir_lab = if_else(dir == "n_up",
                             "Up regulated", "Down regulated")) %>%
    ggplot(aes(contrast, n_plot, fill = dir_lab)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, colour = "grey40") +
geom_label_repel(
  aes(
    label = scales::comma(abs(n)),
    y = n_plot + sign(n_plot) * 80
  ),
  fill = "white",
  color = "darkgrey",
  size = 3
)+
    scale_fill_manual(values = c(
      "Up regulated" = "#D55E00" ,
      "Down regulated" = "#0072B2"
    )) +
    scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
    coord_flip() +
    guides(fill="none")+
    labs(x = NULL, y = "DE gene count", fill = NULL,
         title = "A - Differential expression gene counts \n(FDR < 0.05, |log₂FC| > 1)") +
    plot_theme + theme(legend.position = "bottom")
  

} else {
  ggplot() + annotate("text",x=0.5,y=0.5,label="Run Q3 first",size=4) +
    theme_void() + labs(title="A - DE counts")
}
  ggsave(file.path(outdir_summary, "Q3_PA.png"),
       pA_q3, width = 5, height = 5, dpi = 200)
# ---- B: Volcano — distal ovary (largest mating effect) --------------------
pB_q3 <- make_summary_volcano(
  q3_res_list[["Distal Ovary - \n Virgin vs Mated"]],
  "B - Disytal Ovary - Virgin vs Mated  ",
  pos_lab = "Mated-biased",
  neg_lab = "Virgin-biased"
)
  ggsave(file.path(outdir_summary, "Q3_PB.png"),
       pB_q3, width = 5, height = 5, dpi = 200)
# ---- C: Shared up-regulated genes across mating-status contrasts -----------
# (abdomen, proximal ovary, distal ovary — the three tissue contrasts)
mating_tissue_contrasts <- c("Soma - \n Virgin vs Mated","Proximal Ovary - \n Virgin vs Mated","Distal Ovary - \n Virgin vs Mated")

up_sets <- lapply(mating_tissue_contrasts, function(tc) {
  res <- q3_res_list[[tc]]
  if (is.null(res)) return(character(0))
  res %>% filter(FDR < fdr_thresh, logFC > lfc_thresh) %>% pull(gene)
}) %>% setNames(mating_tissue_contrasts)

shared_up <- Reduce(intersect, Filter(function(x) length(x) > 0, up_sets))

shared_df <- tibble(
  set = c(names(up_sets), "Shared across\nall 3 tissues"),
  n   = c(sapply(up_sets, length), length(shared_up))
) 

pC_q3 <- ggplot(shared_df, aes(fct_reorder(set, -n), n, fill = set)) +
  geom_col(width = 0.6) +
  geom_richtext( aes( label = scales::comma(n)), fill = "white", label.color = "darkgrey", size = 3 ) +
  scale_fill_manual(values = c(
    "Soma - \n Virgin vs Mated" = "#0072B2",
    "Proximal Ovary - \n Virgin vs Mated" = "#D55E00",
    "Distal Ovary - \n Virgin vs Mated" = "#CC79A7",
    "Shared across\nall 3 tissues" = "#000000"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Up-regulated gene count",
       title = "C - Post-mating up-regulated genes\nper tissue (and shared)") +
  plot_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(outdir_summary, "Q3_PC.png"),
       pC_q3, width = 5, height = 5, dpi = 200)
  
# ---- D: Expression heatmap — top distal ovary DE genes across all tissues -
res_do <- q3_res_list[["Distal Ovary - \n Virgin vs Mated"]]

if (!is.null(res_do)) {
  top_do_genes <- res_do %>%
    filter(FDR < fdr_thresh, abs(logFC) >= lfc_thresh) %>%
    arrange(FDR) %>%
    slice_head(n = 40) %>%
    pull(gene)

  # Build expression matrix across ovary + abdomen samples
  meta_hm_q3 <- meta_raw %>%
    filter(tissue %in% c("YO","OO","AB"), sex == "F") %>%
    mutate(tissue_mating = paste0(
      recode(tissue, YO="Distal_OV", OO="Proximal_OV", AB="Abdomen"),
      "_", mating_status))

  if (nrow(meta_hm_q3) >= 3) {
    counts_hm_q3 <- counts_short_raw[, meta_hm_q3$id, drop = FALSE]
    dge_hm_q3 <- DGEList(counts = counts_hm_q3)
    keep_hm   <- rowSums(cpm(dge_hm_q3) >= 1) >= 2
    dge_hm_q3 <- dge_hm_q3[keep_hm, , keep.lib.sizes = FALSE]
    dge_hm_q3 <- calcNormFactors(dge_hm_q3)
    lcpm_hm_q3 <- cpm(dge_hm_q3, log = TRUE, prior.count = 2)

    genes_in_mat <- intersect(top_do_genes, rownames(lcpm_hm_q3))

    if (length(genes_in_mat) >= 2) {
      mat_q3 <- lcpm_hm_q3[genes_in_mat, , drop = FALSE]
      mat_q3_z <- t(scale(t(mat_q3)))
      mat_q3_z[mat_q3_z >  3] <-  3
      mat_q3_z[mat_q3_z < -3] <- -3

      ann_q3 <- meta_hm_q3 %>%
        dplyr::select(id, tissue_mating) %>%
        column_to_rownames("id")

    ann_colors <- list(
  tissue_mating = c(
    "Distal_OV_M"     = "#0072B2",
    "Distal_OV_V"     = "#56B4E9",
    "Proximal_OV_M"   = "#D55E00",
    "Proximal_OV_V"   = "#E69F00",
    "Abdomen_M"       = "#CC79A7",
    "Abdomen_V"       = "#E6A9C8",
    "Abdomen_NA"      = "#999999"
  )
)

      p_hm_q3 <- pheatmap(mat_q3_z,
               annotation_col    = ann_q3,
               annotation_colors = ann_colors,
               clustering_method = "ward.D2",
               show_row    = FALSE,
               show_colnames     = FALSE,
               color    = circlize::colorRamp2(
                 c(-3, 0, 3),
                 c("#0072B2","white","#D55E00"))(seq(-3, 3, length.out = 100)),
               main     = "D -  Top DE genes \nacross all female tissues",
               silent   = TRUE)
      pD_q3 <- wrap_elements(
  full = grid.grabExpr(ComplexHeatmap::draw(p_hm_q3))
)
    } else {
      pD_q3 <- NULL
    }
  } else {
    pD_q3 <- NULL
  }
} else {
  pD_q3 <- NULL
}
  ggsave(file.path(outdir_summary, "Q3_PD.png"),
       pD_q3, width = 5, height = 5, dpi = 200)
  
# ---- Assemble --------------------------------------------------------------
top_row   <- pA_q3 | pB_q3
bot_row   <- pC_q3 | pD_q3

p_q3_summary <- top_row / bot_row +
  plot_annotation(
    title    = "Post-mating Transcriptional Remodelling in Females",
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave(file.path(outdir_summary, "Q3_PostMating_summary.pdf"),
       p_q3_summary, width = 14, height = 11)
ggsave(file.path(outdir_summary, "Q3_PostMating_summary.png"),
       p_q3_summary, width = 14, height = 11, dpi = 200)

p_q3_summary

outdir_go_summary <- file.path(outdir_raw, "03_PostMating/summary_panels", "GO_heatmap")
dir.create(outdir_go_summary, recursive = TRUE, showWarnings = FALSE)

# ---- Load all short-genome GO summary files from q3 -----------------------
q3_go_files <- list.files(
  path       = file.path(outdir_raw, "03_PostMating"),
  pattern    = "GO_terms_merged\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
) %>%
  .[str_detect(., "short")]

message("q3 GO files found: ", length(q3_go_files))
print(q3_go_files)
# should be 2 because C2 has all 3 merged 
if (length(q3_go_files) == 0) stop("No q3 short GO files found — check paths")

# Read and label by contrast (parent folder name)
q3_go <- purrr::map_dfr(q3_go_files, function(f) {
  
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

message("Total GO rows loaded: ", nrow(q3_go))
message("Contrast-groups found:")
print(unique(q3_go$contrast_group))

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
    message("No shared terms across groups — taking top terms by mean enrichment")
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
    legend_labels     = c("0","2","4","6","8","≥10")
  )
}

# ============================================================
# PLOT 1: All q3 contrasts in one heatmap
#   Columns = contrast | group (e.g. "contrast_1 | unfert")
#   Shows which GO terms are shared across MZT contrasts
#   and which are specific to one comparison
# ============================================================
p_q3_heatmap <- make_contrast_go_heatmap(
  go_data    = q3_go,
  title_text = "MZT: GO enrichment across all contrasts\n(−log₁₀ p)",
  n_terms    = 30,
  min_groups = 2
)
p_q3_heatmap<- wrap_elements(
  full = grid.grabExpr(ComplexHeatmap::draw( p_q3_heatmap)))
p_q3_heatmap
if (!is.null(p_q3_heatmap)) {

  ggsave(file.path(outdir_go_summary, "q3_GO_heatmap_all_contrasts.pdf"), p_q3_heatmap, width = 16, height = 10, dpi = 200)
  

}

# ============================================================
# PLOT 2: Faceted bar — top 8 terms per contrast-group
#   One panel per contrast-group
#   Good for seeing what each specific comparison is driven by
# ============================================================
p_q3_facet <- q3_go %>%
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

p_q3_facet
ggsave(file.path(outdir_go_summary, "q3_GO_faceted_per_contrast.pdf"),
       p_q3_facet, width = 14, height = 0.6 * length(unique(q3_go$contrast_group)) + 4,
       limitsize = FALSE)
