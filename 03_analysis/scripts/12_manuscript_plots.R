#!/usr/bin/env Rscript
# Manuscript figure assembly
# Output: 04_results/12_manuscript_plots/Figure_1..6 and Figure_S1..S18 ({pdf,png}) + manifests
# Run from the repository root:  Rscript 03_analysis/scripts/12_manuscript_plots.R
#
# Collects the per-section summary panels (scripts 01-05, 13, 19) under the
# figure numbers used in the manuscript. Run those scripts first.
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

output_dir <- here("04_results/12_manuscript_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# figure number -> source panel (built by the section scripts) + title
figures <- tibble::tribble(
  ~figure, ~source_rel, ~caption,
  1L, "05_Tissue_landscapes/summary/Tissue_landscape_summary2",
      "Tissue identity structures the adult transcriptome, and the testes are the most transcriptionally specialised of the profiled tissues.",
  2L, "01_MZT/summary_panels/Q1_MZT_summary",
      "The maternal-to-zygotic transition unfolds as three sequential phases of gene expression.",
  3L, "19_New_figures/Figure_3_MZT_intron_exon",
      "Intronic reads indicate that the main wave of zygotic transcription occurs between 4 and 8 h after oviposition.",
  4L, "02_Sexual_dimorphism/summary_panels/Q2_SexualDimorphism_summary",
      "Sexual dimorphism is concentrated in the germline and is smaller in the abdomen.",
  5L, "03_PostMating/summary_panels/Q3_PostMating_summary",
      "Mating changes few genes in the ovary and almost none in the abdomen.",
  6L, "04_Maleness/summary_panels/Q4_Maleness_summary",
      "Testes are the most transcriptionally distinct male tissue, and a small set of genes is higher in males in both germline and soma."
)

missing <- character(0)
for (i in seq_len(nrow(figures))) {
  fig  <- figures$figure[i]
  base <- here("04_results", figures$source_rel[i])
  copied_any <- FALSE
  for (ext in c("pdf", "png")) {
    src <- paste0(base, ".", ext)
    if (file.exists(src)) {
      file.copy(src, file.path(output_dir, sprintf("Figure_%d.%s", fig, ext)),
                overwrite = TRUE)
      copied_any <- TRUE
    }
  }
  if (copied_any) {
    message(sprintf("  Figure %d  <-  %s", fig, figures$source_rel[i]))
  } else {
    missing <- c(missing, sprintf("Figure %d (%s.*)", fig, figures$source_rel[i]))
  }
}

readr::write_tsv(
  figures %>% dplyr::mutate(file = sprintf("Figure_%d.pdf", figure)) %>%
    dplyr::select(figure, file, source = source_rel, caption),
  file.path(output_dir, "figure_manifest.tsv")
)

# ---- Supplementary figures ---------------------------------------------------
# Numbered in order of first citation in the manuscript text.
supp_figures <- tibble::tribble(
  ~figure, ~source_rel, ~caption,
  1L, "06_Global_atlas/PCA_tissue_sex",
      "Tissue identity is the dominant axis of variation across adult samples.",
  2L, "10_CrossTissue_correlation/tissue_spearman_correlation",
      "Adult tissues form distinct clusters by pairwise correlation.",
  3L, "07_Tau_specificity/tau_distribution",
      "Tissue specificity of expressed genes.",
  4L, "05_Tissue_landscapes/shared_private/ovary_split_vs_merged",
      "Splitting the ovary into two compartments lowers the percentage of tissue-private genes.",
  5L, "13_Extended_analyses/A_per_organ_GO/per_organ_private_GO",
      "Functional enrichment of genes expressed in only one adult tissue.",
  6L, "19_New_figures/Supp_HK_odds_ratios",
      "Enrichment of housekeeping orthologues among constitutively expressed genes by tissue.",
  7L, "05_Tissue_landscapes/housekeeping/HK_ortho_summary",
      "Drosophila housekeeping orthologues are enriched among highly expressed genes and vary little across tissues.",
  8L, "19_New_figures/Supp_GO_constitutive",
      "Functional enrichment of highly expressed genes in each adult tissue.",
  9L, "01_MZT/contrast_1/GO_short_unfert_vs_fert_combined",
      "Functional enrichment of genes differentially expressed between unfertilised and fertilised eggs.",
  10L, "19_New_figures/Supp_MZT_composition",
      "Library composition, normalisation and tissue similarity of unfertilised eggs.",
  11L, "01_MZT/contrast_3/GO_module_M1",
      "Functional enrichment of co-expression module M1 (maternal decay).",
  12L, "01_MZT/contrast_3/GO_module_M2",
      "Functional enrichment of co-expression module M2 (transient activation).",
  13L, "01_MZT/contrast_3/GO_module_M3",
      "Functional enrichment of co-expression module M3 (late activation).",
  14L, "02_Sexual_dimorphism/contrast_1/GO_short_combined",
      "Functional enrichment of genes differentially expressed between testes and ovary.",
  15L, "02_Sexual_dimorphism/contrast_2/GO_short_combined",
      "Functional enrichment of genes differentially expressed between distal and proximal ovary.",
  16L, "02_Sexual_dimorphism/contrast_3/GO_short_combined",
      "Functional enrichment of genes differentially expressed between male and virgin female abdomen.",
  18L, "19_New_figures/Supp_PostMating_immune",
      "Immune gene-set tests for mated against virgin females."
)

supp_missing <- character(0)
for (i in seq_len(nrow(supp_figures))) {
  fig  <- supp_figures$figure[i]
  base <- here("04_results", supp_figures$source_rel[i])
  copied_any <- FALSE
  for (ext in c("pdf", "png")) {
    src <- paste0(base, ".", ext)
    if (file.exists(src)) {
      file.copy(src, file.path(output_dir, sprintf("Figure_S%d.%s", fig, ext)),
                overwrite = TRUE)
      copied_any <- TRUE
    }
  }
  if (copied_any) {
    message(sprintf("  Figure S%d  <-  %s", fig, supp_figures$source_rel[i]))
  } else {
    supp_missing <- c(supp_missing, sprintf("Figure S%d (%s.*)", fig, supp_figures$source_rel[i]))
  }
}

# S17: sex-biased UpSet + GO summary composite -- assembled from its RDS
# rather than copied, since no single pre-rendered file for it exists.
upset_rds <- here("04_results/08_UpSet_sex_bias/go_summary_upset.RDS")
if (file.exists(upset_rds)) {
  rl <- readRDS(upset_rds)

  # UpSetR plots have no built-in title argument (unlike the ggplot panels
  # elsewhere in this script), so the panel letter is added as a text grob
  # stacked above each captured UpSet plot instead -- g3/g4 already carry
  # their C/D letters in-title from 08_upset_sex_bias.R's plot_topgo_bar()
  # call, matching the same house style.
  titled_upset <- function(upset_obj, title_text) {
    gridExtra::arrangeGrob(
      grid::textGrob(title_text, x = 0.01, hjust = 0,
                     gp = grid::gpar(fontsize = 13, fontface = "bold")),
      grid::grid.grabExpr(print(upset_obj)),
      ncol = 1, heights = c(0.06, 0.94)
    )
  }

  g1 <- titled_upset(rl$male_biased$p,
                     "A - Male-biased gene overlap across tissues")
  g2 <- titled_upset(rl$female_biased$p,
                     "B - Female-biased gene overlap across tissues")
  g3 <- rl$male_biased$summary_cls
  g4 <- rl$female_biased$summary_cls
  png(file.path(output_dir, "Figure_S17.png"),
      width = 12, height = 9, units = "in", res = 300)
  grid::grid.draw(gridExtra::arrangeGrob(g1, g2, g3, g4, ncol = 2))
  dev.off()
  message("  Figure S17  <-  08_UpSet_sex_bias/go_summary_upset.RDS")
} else {
  supp_missing <- c(supp_missing, "Figure S17 (08_UpSet_sex_bias/go_summary_upset.RDS)")
}

supp_figures <- dplyr::bind_rows(
  supp_figures,
  tibble::tibble(figure = 17L, source_rel = "08_UpSet_sex_bias/go_summary_upset.RDS",
                  caption = "Genes with higher expression in each sex are mostly specific to the germline.")
) %>% dplyr::arrange(figure)

readr::write_tsv(
  supp_figures %>%
    dplyr::mutate(file = vapply(figure, function(fig) {
      pdf_path <- file.path(output_dir, sprintf("Figure_S%d.pdf", fig))
      if (file.exists(pdf_path)) sprintf("Figure_S%d.pdf", fig) else sprintf("Figure_S%d.png", fig)
    }, character(1))) %>%
    dplyr::select(figure, file, source = source_rel, caption),
  file.path(output_dir, "supplementary_figure_manifest.tsv")
)

if (length(supp_missing) > 0) {
  warning("Supplementary source panel(s) not found - run the corresponding ",
          "section script(s) first:\n  ", paste(supp_missing, collapse = "\n  "),
          call. = FALSE)
}

if (length(missing) > 0) {
  stop("Source panel(s) not found - run the corresponding section script(s) first:\n  ",
       paste(missing, collapse = "\n  "), call. = FALSE)
}
message("Manuscript figures written to ", output_dir)
