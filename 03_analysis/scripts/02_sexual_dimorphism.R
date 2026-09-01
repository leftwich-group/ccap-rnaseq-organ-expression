#!/usr/bin/env Rscript
# Q2 - Sexual dimorphism across germline and soma
# Output: 04_results/02_Sexual_dimorphism/
# Run from the repository root:  Rscript 03_analysis/scripts/02_sexual_dimorphism.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir_q2 <- file.path(outdir_raw, "02_Sexual_dimorphism")
dir.create(outdir_q2, recursive = TRUE, showWarnings = FALSE)

outdir <- file.path(outdir_q2, "contrast_1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue %in% c("TE", "OO", "YO")) %>%
  mutate(gamete = case_when(
    tissue == "TE"              ~ "Testes",
    tissue %in% c("OO", "YO")  ~ "Ovary",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(gamete))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long  <- counts_long_raw[,  meta$id, drop = FALSE]

design <- model.matrix(~ 0 + gamete, data = meta)
colnames(design) <- sub("gamete", "", colnames(design))

contrast_list <- list(female_gonad_vs_male_gonad = "Ovary - Testes")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "gamete",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
  long  = run_reference_analysis(counts_long,  meta, design, contrast_list,
                                 "long",  outdir, "gamete",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q2", "Germline: testes vs ovary",
                  "Ovary - Testes", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$female_gonad_vs_male_gonad$table
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

outdir <- file.path(outdir_q2, "contrast_2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue %in% c("YO", "OO")) %>%
  mutate(ovary_pos = if_else(tissue == "YO", "Distal", "Proximal"))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long  <- counts_long_raw[,  meta$id, drop = FALSE]

design        <- model.matrix(~ 0 + ovary_pos, data = meta)
colnames(design) <- sub("ovary_pos", "", colnames(design))
contrast_list <- list(distal_vs_proximal = "Distal - Proximal")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "ovary_pos",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
  long  = run_reference_analysis(counts_long,  meta, design, contrast_list,
                                 "long",  outdir, "ovary_pos",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q2", "Ovary axis: proximal vs distal",
                  "Distal - Proximal", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$distal_vs_proximal$table
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

outdir <- file.path(outdir_q2, "contrast_3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue == "AB",
         !(sex == "F" & mating_status == "M")) %>%
  mutate(soma_sex = if_else(sex == "M", "Abdomen_M", "Abdomen_F_virgin"))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long  <- counts_long_raw[,  meta$id, drop = FALSE]

design        <- model.matrix(~ 0 + soma_sex, data = meta)
colnames(design) <- sub("soma_sex", "", colnames(design))
contrast_list <- list(AbdomenF_virgin_vs_AbdomenM = "Abdomen_F_virgin - Abdomen_M")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "soma_sex",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
  long  = run_reference_analysis(counts_long,  meta, design, contrast_list,
                                 "long",  outdir, "soma_sex",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q2", "Sommale vs virgin female",
                  "Abdomen_F_virgin - Abdomen_M", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$AbdomenF_virgin_vs_AbdomenM$table
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

outdir <- file.path(outdir_q2, "contrast_4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- meta_raw %>%
  filter(tissue == "AB",
         !(sex == "F" & mating_status == "V")) %>%
  mutate(soma_grp = case_when(
    sex == "M" ~ "Abdomen_M",
    sex == "F" & mating_status == "M" ~ "Abdomen_F_mated",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(soma_grp))

counts_short <- counts_short_raw[, meta$id, drop = FALSE]
counts_long  <- counts_long_raw[,  meta$id, drop = FALSE]

design        <- model.matrix(~ 0 + soma_grp, data = meta)
colnames(design) <- sub("soma_grp", "", colnames(design))
contrast_list <- list(AbdomenF_mated_vs_AbdomenM = "Abdomen_F_mated - Abdomen_M")

results_all <- list(
  short = run_reference_analysis(counts_short, meta, design, contrast_list,
                                 "short", outdir, "soma_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh),
  long  = run_reference_analysis(counts_long,  meta, design, contrast_list,
                                 "long",  outdir, "soma_grp",
                                 min_count, min_samples, fdr_thresh, lfc_thresh)
)
saveRDS(results_all, file.path(outdir, "results_all_nested.rds"))
for (g in c("short","long"))
  append_de_stats(results_all, "Q2", "Sommale vs mated female",
                  "Abdomen_F_mated - Abdomen_M", genome = g)

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$AbdomenF_mated_vs_AbdomenM$table
  up_g      <- res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene) %>% clean_ids()
  down_g    <- res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene) %>% clean_ids()
  save_go_panel(run_topgo_all_ont(up_g, gene_univ, geneID2GO),
                run_topgo_all_ont(down_g, gene_univ, geneID2GO),
                outdir, paste0(genome))
}

for (genome in c("short", "long")) {
  geneID2GO <- get(paste0("geneID2GO_",      genome))
  gene_univ <- get(paste0("gene_universe_",   genome))
  res       <- results_all[[genome]]$contrasts$AbdomenF_mated_vs_AbdomenM$table
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

# ---------------------------------------------------------------------------
# Q2 SUMMARY: Sexual dimorphism across germline and soma
# Four-panel figure:
#   A. Volcano: male vs female gonad (germline dimorphism)
#   B. Volcano: male vs female abdomen (somatic dimorphism, virgin)
#   C. Mirror bar: DE gene counts across all Q2 contrasts
#   D. Venn/UpSet: overlap of male-biased genes between germline and soma
# ---------------------------------------------------------------------------

# Helper: load short-genome results from a saved rds
load_short_res <- function(rds_path, contrast_name) {
  if (!file.exists(rds_path)) return(NULL)
  readRDS(rds_path)$short$contrasts[[contrast_name]]$table
}

outdir_summary <- file.path(outdir_raw, "02_Sexual_dimorphism/summary_panels")
dir.create(outdir_summary, recursive = TRUE, showWarnings = FALSE)

# ---- Load the four Q2 contrasts -------------------------------------------
res_germline   <- load_short_res(
  file.path(outdir_raw, "02_Sexual_dimorphism/contrast_1/results_all_nested.rds"),
  "female_gonad_vs_male_gonad"
)
res_distal_prox <- load_short_res(
  file.path(outdir_raw, "02_Sexual_dimorphism/contrast_2/results_all_nested.rds"),
  "distal_vs_proximal"
)
res_soma_virgin <- load_short_res(
  file.path(outdir_raw, "02_Sexual_dimorphism/contrast_3/results_all_nested.rds"),
  "AbdomenF_virgin_vs_AbdomenM"
)
res_soma_mated  <- load_short_res(
  file.path(outdir_raw, "02_Sexual_dimorphism/contrast_4/results_all_nested.rds"),
  "AbdomenF_mated_vs_AbdomenM"
)

# ---- A: Germline volcano ---------------------------------------------------
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
    labs(x = "log₂FC", y = "-log₁₀(FDR)",
         colour = NULL, title = title_text) +
    plot_theme +
    theme(legend.position = "none")
}

pA_q2 <- make_summary_volcano(
  res_germline,
  "Germline: Ovary vs Testes",
  pos_lab = "Ovary-biased",
  neg_lab = "Testis-biased"
)
ggsave(file.path(outdir_summary, "Q2_plotA.png"),
       pA_q2, width = 5, height = 5, dpi = 200)
pB_q2 <- make_summary_volcano(
  res_soma_virgin,
  "Soma (Abdomen): Male vs Female",
  pos_lab = "Female-biased",
  neg_lab = "Male-biased"
)
ggsave(file.path(outdir_summary, "Q2_plotB.png"),
       pB_q2, width = 5, height = 5, dpi = 200)
# ---- C: Mirror bar — DE counts across all contrasts -----------------------
q2_contrast_meta <- tribble(
  ~label,                          ~rds_path,                                            ~contrast_name,
  "Germline\n(Testes vs Ovary)",   "02_Sexual_dimorphism/contrast_1/results_all_nested.rds", "female_gonad_vs_male_gonad",
  
  "Ovary axis\n(Proximal vs Distal)",  "02_Sexual_dimorphism/contrast_2/results_all_nested.rds", "distal_vs_proximal",
  
  "Soma - \n Male vs virgin female", "02_Sexual_dimorphism/contrast_3/results_all_nested.rds", "AbdomenF_virgin_vs_AbdomenM",
  
  "Soma - \n Male vs mated female", "02_Sexual_dimorphism/contrast_4/results_all_nested.rds", "AbdomenF_mated_vs_AbdomenM",
)

de_counts_q2 <- pmap_dfr(q2_contrast_meta, function(label, rds_path, contrast_name) {
  res <- load_short_res(file.path(outdir_raw, rds_path), contrast_name)
  if (is.null(res)) return(NULL)
  tibble(
    contrast = label,
    n_up   = sum(res$FDR < fdr_thresh & res$logFC >  lfc_thresh),
    n_down = sum(res$FDR < fdr_thresh & res$logFC < -lfc_thresh)
  )
})

pC_q2 <- if (nrow(de_counts_q2) > 0) {
  de_counts_q2 %>%
    mutate(contrast = factor(contrast, levels = rev(contrast))) %>%
    pivot_longer(c(n_up, n_down), names_to = "dir", values_to = "n") %>%
    mutate(n_plot = if_else(dir == "n_down", -n, n),
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
      "Up regulated"   = "#D55E00",
      "Down regulated" = "#0072B2"
    )) +
    
    scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
    coord_flip() +
    labs(x = NULL, y = "DE gene count",
         fill  = NULL,
         title = "Differential expression gene counts \n( FDR < 0.05, |log₂FC| > 1)") +
    plot_theme +
       guides(fill="none")+
    theme(legend.position = "bottom")
} else {
  ggplot() + annotate("text", x=0.5, y=0.5,
                      label="Run Q2 contrasts first", size=4) +
    theme_void() + labs(title="C  DE gene counts")
}
ggsave(file.path(outdir_summary, "Q2_plotC.png"),
       pC_q2, width = 5, height = 5, dpi = 200)
# ---- D: Shared male-biased genes across germline and soma -----------------
if (!is.null(res_germline) && !is.null(res_soma_virgin)) {
  male_germ <- res_germline %>%
    filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene)
  male_soma <- res_soma_virgin %>%
    filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene)

  # Build a simple overlap dot chart (UpSet-style summary without the package)
  venn_df <- tibble(
    set  = c("Germline only", "Soma only", "Shared"),
    n    = c(
      length(setdiff(male_germ, male_soma)),
      length(setdiff(male_soma, male_germ)),
      length(intersect(male_germ, male_soma))
    )
  ) %>%
    mutate(set = factor(set, levels = c("Germline only","Soma only","Shared")))

  pD_q2 <- ggplot(venn_df, aes(set, n, fill = set)) +
    geom_col(width = 0.55) +
    geom_richtext(
    aes(label = scales::comma(n)),
    fill = "white",
    label.color = "darkgrey",  
    size = 3
) + 
    scale_fill_manual(values = c(
      "Germline only" = "#0072B2",
      "Shared"        = "#009E73",
      "Soma only"     = "#E69F00"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Gene count",
         title = "Overlap of male-biased\ngenes: germline vs soma") +
    plot_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
} else {
  pD_q2 <- ggplot() +
    annotate("text", x=0.5, y=0.5, label="Data not found", size=4) +
    theme_void() + labs(title="D  Overlap")
}
ggsave(file.path(outdir_summary, "Q2_plotD.png"),
       pD_q2, width = 5, height = 5, dpi = 200)
# ---- Assemble ---------------------------------------------------------------
p_q2_summary <- (pA_q2 | pB_q2) / (pC_q2 | pD_q2) +
  plot_annotation(
    title    = "Sexual Dimorphism across Germline and Soma",
    tag_levels = "A",
    theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave(file.path(outdir_summary, "Q2_SexualDimorphism_summary.pdf"),
       p_q2_summary, width = 14, height = 11)
ggsave(file.path(outdir_summary, "Q2_SexualDimorphism_summary.png"),
       p_q2_summary, width = 14, height = 11, dpi = 200)

p_q2_summary


# save just volcano 
pp <- (pA_q2 | pB_q2) + 
  plot_annotation(tag_levels = "A")
ggsave(file.path(outdir_summary, "Q2_SexualDimorphism_volcano.png"),
       pp, width = 10, height = 7, dpi = 200)

# save just bar 
pp <- (pC_q2 | pD_q2) + 
  plot_annotation(tag_levels = "A")
ggsave(file.path(outdir_summary, "Q2_SexualDimorphism_barplot.png"),
       pp, width = 10, height = 7, dpi = 200)

outdir_go_summary <- file.path(outdir_raw, "02_Sexual_dimorphism/summary_panels", "GO_heatmap")
dir.create(outdir_go_summary, recursive = TRUE, showWarnings = FALSE)

# ---- Load all short-genome GO summary files from q2 -----------------------
q2_go_files <- list.files(
  path       = file.path(outdir_raw, "02_Sexual_dimorphism"),
  pattern    = "GO_terms_merged\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
) %>%
  .[str_detect(., "short")]

message("q2 GO files found: ", length(q2_go_files))
print(q2_go_files)
# should be 2 because C2 has all 3 merged 
if (length(q2_go_files) == 0) stop("No q2 short GO files found — check paths")

# Read and label by contrast (parent folder name)
q2_go <- map_dfr(q2_go_files, function(f) {
  contrast <- basename(dirname(f))   # e.g. "contrast_1", "contrast_2"
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(
      weight01_fisher = as.character(weight01_fisher),
      contrast       = contrast,
      p_value        = parse_topgo_p(weight01_fisher),
      neg_log10p     = -log10(pmax(p_value, 1e-10)),
      contrast_group = paste0(contrast," | ", group)
    )
  
  
})

message("Total GO rows loaded: ", nrow(q2_go))
message("Contrast-groups found:")
print(unique(q2_go$contrast_group))

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
# PLOT 1: All q2 contrasts in one heatmap
#   Columns = contrast | group (e.g. "contrast_1 | unfert")
#   Shows which GO terms are shared across MZT contrasts
#   and which are specific to one comparison
# ============================================================
p_q2_heatmap <- make_contrast_go_heatmap(
  go_data    = q2_go,
  title_text = "MZT: GO enrichment across all contrasts\n(−log₁₀ p)",
  n_terms    = 30,
  min_groups = 2
)
p_q2_heatmap<- wrap_elements(
  full = grid.grabExpr(ComplexHeatmap::draw( p_q2_heatmap)))
p_q2_heatmap
if (!is.null(p_q2_heatmap)) {

  ggsave(file.path(outdir_go_summary, "q2_GO_heatmap_all_contrasts.pdf"), p_q2_heatmap, width = 16, height = 10, dpi = 200)
  

}

# ============================================================
# PLOT 2: Faceted bar — top 8 terms per contrast-group
#   One panel per contrast-group
#   Good for seeing what each specific comparison is driven by
# ============================================================
p_q2_facet <- q2_go %>%
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

p_q2_facet
ggsave(file.path(outdir_go_summary, "q2_GO_faceted_per_contrast.pdf"),
       p_q2_facet, width = 14, height = 0.6 * length(unique(q2_go$contrast_group)) + 4,
       limitsize = FALSE)
