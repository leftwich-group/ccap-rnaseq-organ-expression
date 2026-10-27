#!/usr/bin/env Rscript
# Master DE-statistics table - final wide table + overview figure
# Output: 04_results/DE_stats_master_wide.tsv, 04_results/DE_overview_all_contrasts.*
# Run from the repository root:  Rscript 03_analysis/scripts/99_de_stats_master.R
# Run after 01-04 (the Q scripts) have populated 04_results/DE_stats_master.tsv.
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

# ===========================================================================
# Print and save the completed DE statistics master table
# ===========================================================================

if (file.exists(de_master_file)) {
  de_master <- read_tsv(de_master_file, show_col_types = FALSE)

  # Pivot to wide: short and long side-by-side for easy comparison
  de_wide <- de_master %>%
    pivot_wider(
      id_cols     = c(section, question, contrast),
      names_from  = genome,
      values_from = c(n_tested, n_up, n_down, n_sig_total, pct_sig),
      names_sep   = "_"
    ) %>%
    arrange(section, question)

  write_tsv(de_wide,
            file.path(outdir_raw, "DE_stats_master_wide.tsv"))

  # Summary bar: total DE genes per section (short genome only)
  de_section_summary <- de_master %>%
    filter(genome == "short") %>%
    group_by(section) %>%
    summarise(
      total_contrasts = n(),
      total_up        = sum(n_up),
      total_down      = sum(n_down),
      total_DE        = sum(n_sig_total),
      .groups = "drop"
    )

  p_de_overview <- de_master %>%
    filter(genome == "short") %>%
    mutate(
      label = paste0(section, ": ", str_trunc(question, 35)),
      label = fct_inorder(label)
    ) %>%
   
    pivot_longer(c(n_up, n_down), names_to = "dir", values_to = "n") %>%
    mutate(n_plot  = if_else(dir == "n_down", -n, n),
           dir_lab = if_else(dir == "n_up", "Up", "Down")) %>%
     mutate(
  label_pos = case_when(
    abs(n_plot) < 100 ~ n_plot + sign(n_plot) * 500,
    TRUE              ~ n_plot * 0.9
  ) # if DE count is low increat the distance between labels 
) %>%
    ggplot(aes(label, n_plot, fill = dir_lab)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.5) +
geom_richtext(
  aes(
    label = scales::comma(abs(n_plot)),
    y = label_pos
  ),
  fill = "white",
  label.color = "darkgrey",
  size = 3
) +
    scale_fill_manual(values = c(Up = "#D55E00",
                                 Down = "#0072B2")) +
    scale_y_continuous(labels = function(x) scales::comma(abs(x))) +
    coord_flip() +
    facet_grid(section ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = "DE gene count",
         fill = NULL,
         title = "DE gene counts across all contrasts",
         subtitle = "FDR < 0.05  |  |log2FC| > 1") +
    plot_theme +
    theme(
      legend.position   = "bottom",
      strip.text.y      = element_text(angle = 0, size = 9),
      axis.text.y       = element_text(size = 8)
    )+
    guides(fill="none")

  ggsave(file.path(outdir_raw, "DE_overview_all_contrasts.pdf"),
         p_de_overview,
         width  = 12,
         height = 0.35 * nrow(de_master[de_master$genome=="short",]) + 3,
         limitsize = FALSE)
  ggsave(file.path(outdir_raw, "DE_overview_all_contrasts.png"),
         p_de_overview,
         width  = 12,
         height = 0.35 * nrow(de_master[de_master$genome=="short",]) + 3,
         dpi = 200, limitsize = FALSE)

  p_de_overview

  # Print wide table as knitr table
  knitr::kable(
    de_wide,
    caption = "DE statistics - all contrasts (short and long genome)",
    format  = "html"
  )

} else {
  message("DE master file not found - run contrasts first.")
}
