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
# Numbered in the order each is first cited in the manuscript text (S1 = the
# first "(Figure Sx)" a reader hits, and so on) -- NOT the order they were
# added to this script. Two consequences of that renumbering worth noting:
#
#  - S11 (distal-vs-proximal ovary GO enrichment) backs a sentence that
#    previously cited the wrong figure outright ("(Figure S4)", pointing at
#    the UpSet composite below, which has nothing to do with ovary GO terms)
#    -- this is the correct figure for that claim, spot-checked against the
#    specific terms named in the text (SRP-dependent protein targeting,
#    chitin binding, etc.) before being added here.
#  - S13 (the UpSet composite) is the figure that mis-citation was pointing
#    at. Its own natural home is the "small, shared" germline-soma overlap
#    claim next to Figure 3D, which now also cites it.
supp_figures <- tibble::tribble(
  ~figure, ~source_rel,                                                     ~caption,
  1L, "06_Global_atlas/PCA_tissue_sex",                                      "Tissue identity is the dominant axis of variation across all adult samples.",
  2L, "10_CrossTissue_correlation/tissue_spearman_correlation",              "Adult tissues form distinct expression clusters by pairwise correlation.",
  3L, "07_Tau_specificity/tau_distribution",                                 "Most expressed genes show strong tissue specificity (tau close to 1).",
  4L, "05_Tissue_landscapes/shared_private/ovary_split_vs_merged",           "Splitting the ovary into two compartments reveals no strong proximal or distal ovary identity.",
  5L, "05_Tissue_landscapes/housekeeping/HK_ortho_summary",                  "Drosophila housekeeping orthologues behave as housekeeping genes in C. capitata.",
  6L, "01_MZT/contrast_1/GO_short_unfert_vs_fert_combined",                  "GO enrichment among genes differentially expressed between unfertilised and fertilised eggs.",
  7L, "01_MZT/contrast_3/GO_module_M1",                                     "GO enrichment for co-expression module M1 (maternal decay).",
  8L, "01_MZT/contrast_3/GO_module_M2",                                     "GO enrichment for co-expression module M2 (transient activation).",
  9L, "01_MZT/contrast_3/GO_module_M3",                                     "GO enrichment for co-expression module M3 (late activation).",
  10L, "02_Sexual_dimorphism/contrast_1/GO_short_combined",                 "GO enrichment among genes differentially expressed between testes and ovary (germline sex bias).",
  11L, "02_Sexual_dimorphism/contrast_2/GO_short_combined",                 "GO enrichment among genes differentially expressed between distal and proximal ovary.",
  12L, "02_Sexual_dimorphism/contrast_3/GO_short_combined",                 "GO enrichment among genes differentially expressed between male and virgin female abdomen (somatic sex bias)."
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

# S13: sex-biased UpSet + GO summary composite -- assembled from its RDS
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
  png(file.path(output_dir, "Figure_S13.png"),
      width = 12, height = 9, units = "in", res = 300)
  grid::grid.draw(gridExtra::arrangeGrob(g1, g2, g3, g4, ncol = 2))
  dev.off()
  message("  Figure S13  <-  08_UpSet_sex_bias/go_summary_upset.RDS")
} else {
  supp_missing <- c(supp_missing, "Figure S13 (08_UpSet_sex_bias/go_summary_upset.RDS)")
}

supp_figures <- dplyr::bind_rows(
  supp_figures,
  tibble::tibble(figure = 13L, source_rel = "08_UpSet_sex_bias/go_summary_upset.RDS",
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
