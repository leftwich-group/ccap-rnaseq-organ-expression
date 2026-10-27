#!/usr/bin/env Rscript
# Supplementary GO figures for the adult tissue landscapes
# One panel per tissue: top GO biological process terms among genes in the
# medium (variable layer) and high (constitutive layer) expression classes.
# Input:  04_results/05_Tissue_landscapes/Tissues/merged_all_tissues_GO_summary.tsv
# Output: 04_results/17_Supp_GO_tissues/
# Run from the repository root:  Rscript 03_analysis/scripts/17_supp_go_tissues.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

go_file <- res_dir(RES$landscapes, "Tissues", "merged_all_tissues_GO_summary.tsv")
if (!file.exists(go_file)) stop(
  "Missing: ", go_file, "\n  Run 03_analysis/scripts/05_tissue_landscapes.R first.", call. = FALSE)

outdir <- file.path(outdir_raw, RES$supp_go)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tissue_levels <- c("Testes", "Accessory_Gland", "Abdomen", "Distal_Ovary", "Proximal_Ovary")
n_terms       <- 10
min_genes     <- 3
class_titles  <- c(medium = "Medium-expressed genes (variable layer)",
                   high   = "Highly expressed genes (constitutive layer)")

go <- read_tsv(go_file, show_col_types = FALSE) %>%
  filter(tissue %in% tissue_levels, classification %in% names(class_titles), Significant >= min_genes)

# The summary table truncates long term names; take full names from GO.db.
full_names <- AnnotationDbi::select(GO.db, keys = unique(go$GO.ID),
                                    columns = "TERM", keytype = "GOID") %>%
  dplyr::rename(GO.ID = GOID, term_full = TERM)

top <- go %>%
  left_join(full_names, by = "GO.ID") %>%
  mutate(term_full = coalesce(term_full, Term),
         tissue = factor(tissue, levels = tissue_levels)) %>%
  group_by(classification, tissue) %>%
  arrange(p_value, desc(Significant), .by_group = TRUE) %>%
  slice_head(n = n_terms) %>%
  ungroup()
write_tsv(top %>% dplyr::select(classification, tissue, GO.ID, term = term_full,
                                Annotated, Significant, p_value),
          file.path(outdir, "supp_go_terms_plotted.tsv"))

for (cls in names(class_titles)) {
  d <- top %>% filter(classification == cls) %>%
    mutate(label = str_wrap(term_full, 42),
           label_r = tidytext::reorder_within(label, -p_value, tissue))
  p <- ggplot(d, aes(Significant, label_r, fill = -log10(p_value))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(Significant, "/", Annotated)), hjust = -0.1, size = 2.8) +
    tidytext::scale_y_reordered() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
    scale_fill_gradient(low = "#56B4E9", high = "#D55E00", name = expression(-log[10](p))) +
    facet_wrap(~ tissue, ncol = 2, scales = "free",
               labeller = as_labeller(function(x) gsub("_", " ", x))) +
    labs(x = "Genes in term (labels show genes / annotated)", y = NULL,
         title = class_titles[[cls]]) +
    plot_theme + theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(outdir, paste0("Supp_GO_", cls, "_expressed.pdf")), p, width = 11, height = 13)
  ggsave(file.path(outdir, paste0("Supp_GO_", cls, "_expressed.png")), p, width = 11, height = 13, dpi = 200)
}
message("Done: ", outdir)
