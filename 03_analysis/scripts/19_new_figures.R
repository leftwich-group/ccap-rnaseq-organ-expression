#!/usr/bin/env Rscript
# Assembles the manuscript figures that combine the extra analyses (14-18).
# Input:  plots.rds written by scripts 14, 15, 16 and 18; the PDF/PNG from 17
# Output: 04_results/19_New_figures/
#   Figure_3_MZT_intron_exon      main figure (intronic and exonic change)
#   Supp_MZT_composition          library composition, normalisation, tissue similarity
#   Supp_PostMating_immune        immune gene-set tests
#   Supp_GO_constitutive          GO of highly expressed genes per tissue
#   Supp_HK_odds_ratios           housekeeping orthologue enrichment
# Run from the repository root:  Rscript 03_analysis/scripts/19_new_figures.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

plot_files <- c(mzt = res_dir(RES$mzt_intron, "plots.rds"),
                prov = res_dir(RES$mzt_provision, "plots.rds"),
                sets = res_dir(RES$mating_sets, "plots.rds"),
                hk = res_dir(RES$hk_tests, "plots.rds"))
go_pdf <- res_dir(RES$supp_go, "Supp_GO_high_expressed.pdf")
go_png <- res_dir(RES$supp_go, "Supp_GO_high_expressed.png")
needed <- c(plot_files, go_pdf, go_png)
if (!all(file.exists(needed))) stop(
  "Missing: ", paste(basename(dirname(needed[!file.exists(needed)])), collapse = ", "),
  "\n  Run scripts 14 to 18 first.", call. = FALSE)

outdir <- file.path(outdir_raw, RES$new_figs)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
p <- lapply(plot_files, readRDS)

save_fig <- function(plot, name, width, height) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(outdir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
}
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 14))

# ---- Figure 3: intronic and exonic change across the MZT ------------------
fig3 <- ((p$mzt$fraction + theme(legend.position = "bottom") + ggtitle(NULL)) |
         p$mzt$onset) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(fig3, "Figure_3_MZT_intron_exon", 11, 9)

# ---- Supp: composition, normalisation and tissue similarity ---------
s_comp <- p$prov$composition /
       (p$prov$normalisation | p$prov$similarity) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(s_comp, "Supp_MZT_composition", 11, 8.5)

# ---- Supp: immune gene-set tests ------------------------------------
s_imm <- p$sets$shift / p$sets$distribution +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(s_imm, "Supp_PostMating_immune", 11, 9)

# ---- Supp: constitutive-layer GO (figure from script 17) ------------
file.copy(go_pdf, file.path(outdir, "Supp_GO_constitutive.pdf"), overwrite = TRUE)
file.copy(go_png, file.path(outdir, "Supp_GO_constitutive.png"), overwrite = TRUE)

# ---- Supp: housekeeping odds ratios ---------------------------------
save_fig(p$hk$odds_ratios, "Supp_HK_odds_ratios", 8, 4)

message("Done: ", outdir)
