#!/usr/bin/env Rscript
# Manuscript figure assembly
# Output: 04_results/12_manuscript_plots/Figure_1..5.{pdf,png} + figure_manifest.tsv
# Run from the repository root:  Rscript 03_analysis/scripts/12_manuscript_plots.R
#
# The five main-text figures are the per-section summary panels built by the
# 01-05 scripts. This stage just collects and renames them to match
# de_Coriolis_et_al_Bulk_Figures. Run 01-05 first.
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

output_dir <- here("04_results/12_manuscript_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# figure number -> source panel (built by the section scripts) + draft caption
figures <- tibble::tribble(
  ~figure, ~source_rel,                                                              ~caption,
  1L, "05_Tissue_landscapes/summary/Tissue_landscape_summary2",                       "Tissue identity structures the adult transcriptome, and the testes are the most transcriptionally specialized tissue.",
  2L, "01_MZT/summary_panels/Q1_MZT_summary",                                         "The maternal-to-zygotic transition unfolds as three sequential phases of gene expression.",
  3L, "02_Sexual_dimorphism/summary_panels/Q2_SexualDimorphism_summary",              "Sexual dimorphism is concentrated in the germline, not the soma.",
  4L, "03_PostMating/summary_panels/Q3_PostMating_summary",                           "Mating remodels the ovary, not the soma, and each reproductive compartment responds independently.",
  5L, "04_Maleness/summary_panels/Q4_Maleness_summary",                               "The testes are the most male-defining tissue, and a small germline-soma-shared gene set marks candidate sex-linked loci."
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
# Matches de_Coriolis_et_al_Bulk_suppl_mat.docx's Figure S1-S6 list exactly.
# S1 is the one inference here, not an exact source-file name match like the
# rest -- it's the PCA showing tissue as the dominant axis of variation
# (matching "Splitting the ovary..." etc. for S2 onward is unambiguous).
supp_figures <- tibble::tribble(
  ~figure, ~source_rel,                                            ~caption,
  1L, "06_Global_atlas/PCA_tissue_sex",                             "Tissue identity is the dominant axis of variation across all adult samples.",
  2L, "05_Tissue_landscapes/shared_private/ovary_split_vs_merged",  "Splitting the ovary into two compartments reveals no strong proximal or distal ovary identity.",
  3L, "05_Tissue_landscapes/housekeeping/HK_ortho_expression_class","Drosophila housekeeping orthologues behave as housekeeping genes in C. capitata.",
  5L, "07_Tau_specificity/tau_distribution",                        "Most expressed genes show strong tissue specificity (tau close to 1).",
  6L, "10_CrossTissue_correlation/tissue_spearman_correlation",     "Adult tissues form distinct expression clusters by pairwise correlation."
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

# S4: sex-biased UpSet + GO summary composite -- assembled from its RDS
# rather than copied, since no single pre-rendered file for it exists.
upset_rds <- here("04_results/08_UpSet_sex_bias/go_summary_upset.RDS")
if (file.exists(upset_rds)) {
  rl <- readRDS(upset_rds)
  g1 <- grid::grid.grabExpr(print(rl$male_biased$p))
  g2 <- grid::grid.grabExpr(print(rl$female_biased$p))
  g3 <- rl$male_biased$summary_cls
  g4 <- rl$female_biased$summary_cls
  png(file.path(output_dir, "Figure_S4.png"),
      width = 12, height = 9, units = "in", res = 300)
  grid::grid.draw(gridExtra::arrangeGrob(g1, g2, g3, g4, ncol = 2))
  dev.off()
  message("  Figure S4  <-  08_UpSet_sex_bias/go_summary_upset.RDS")
} else {
  supp_missing <- c(supp_missing, "Figure S4 (08_UpSet_sex_bias/go_summary_upset.RDS)")
}

supp_figures <- dplyr::bind_rows(
  supp_figures,
  tibble::tibble(figure = 4L, source_rel = "08_UpSet_sex_bias/go_summary_upset.RDS",
                  caption = "Male- and female-biased genes are overwhelmingly germline-specific, with a small set shared with the soma.")
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
