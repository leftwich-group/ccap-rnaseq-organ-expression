#!/usr/bin/env Rscript
# Transcription-factor motif enrichment (HOMER)
# Output: 04_results/09_HOMER_motifs/   (reads 02_Sexual_dimorphism/, 03_PostMating/)
# Run from the repository root:  Rscript 03_analysis/scripts/09_homer_motifs.R
# Promoter motif discovery needs HOMER installed and configured (see the
# archived Rmd for setup notes). Set RUN_HOMER_DISCOVERY <- TRUE once ready.
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

RUN_HOMER_DISCOVERY <- FALSE

if (RUN_HOMER_DISCOVERY) {
  # Set paths — adjust to your HOMER installation
  HOMER_PATH   <- "/path/to/your/HOMER/bin"
  GENOME_FASTA <- here("01_data/Genomic_files/GCF_000347755.3_Ccap_2.1_genomic.fna")
  GFF_FILE <- here("01_data/Genomic_files/GCF_000347755.3_Ccap_2.1_genomic.gff")
  HOMER_SIZE   <- 1000   # promoter window (bp upstream of TSS)

  outdir_homer <- file.path(outdir_raw, "09_HOMER_motifs")
  dir.create(outdir_homer, recursive = TRUE, showWarnings = FALSE)

  # Helper: extract promoter BED from gene list
  genes_to_promoter_bed <- function(gene_ids, gff_file,
                                     upstream = 1000, downstream = 100) {
    gff <- read_tsv(gff_file, comment = "#",
                    col_names = c("seqid","source","type","start","end",
                                  "score","strand","phase","attributes"),
                    show_col_types = FALSE) %>%
      filter(type == "gene") %>%
      mutate(gene = str_extract(attributes, "(?<=gene=)[^;]+"),
             locus = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
             gene_id = coalesce(gene, locus)) %>%
      filter(gene_id %in% gene_ids)

    gff %>%
      transmute(
        chrom  = seqid,
        prom_start = if_else(strand == "+",
                             pmax(0, start - upstream),
                             end - downstream),
        prom_end   = if_else(strand == "+",
                             start + downstream,
                             end + upstream),
        name   = gene_id,
        score  = 0,
        strand = strand
      )
  }

  # Helper: run HOMER findMotifs.pl
  run_homer <- function(target_genes, background_genes, label,
                        gff_file = GFF_FILE,
                        genome_fa = GENOME_FASTA) {

    outdir_h <- file.path(outdir_homer, label)
    dir.create(outdir_h, recursive = TRUE, showWarnings = FALSE)

    tgt_bed <- genes_to_promoter_bed(target_genes, gff_file)
    bg_bed  <- genes_to_promoter_bed(background_genes, gff_file)

    tgt_file <- file.path(outdir_h, "target_promoters.bed")
    bg_file  <- file.path(outdir_h, "background_promoters.bed")

    write_tsv(tgt_bed, tgt_file, col_names = FALSE)
    write_tsv(bg_bed,  bg_file,  col_names = FALSE)

    cmd <- paste0(
      file.path(HOMER_PATH, "findMotifsGenome.pl"),
      " ", tgt_file,
      " ", genome_fa,
      " ", outdir_h,
      " -bg ", bg_file,
      " -size ", HOMER_SIZE,
      " -len 6,8,10",
      " -p 4"        # threads
    )
    message("Running HOMER: ", label)
    system(cmd)
  }

  # Run for germline sex-biased genes
  germ_res    <- readRDS(file.path(outdir_raw,
                                   "02_Sexual_dimorphism/contrast_1/results_all_nested.rds"))
  germ_table  <- germ_res$short$contrasts$female_gonad_vs_male_gonad$table

  male_biased_germ   <- germ_table %>%
    filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>% pull(gene)
  female_biased_germ <- germ_table %>%
    filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene)
  background_germ    <- germ_table$gene

  run_homer(male_biased_germ,   background_germ, "germline_male_biased")
  run_homer(female_biased_germ, background_germ, "germline_female_biased")

  # Run for post-mating distal ovary (immune activation set)
  pm_res   <- readRDS(file.path(outdir_raw,
                                 "03_PostMating/contrast_3/results_all_nested.rds"))
  pm_table <- pm_res$short$contrasts$mated_vs_virgin_YO$table
  pm_up    <- pm_table %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>% pull(gene)
  run_homer(pm_up, pm_table$gene, "postmating_DO_upregulated")
}

# ---- Plot whatever HOMER results already exist on disk --------------------
outdir_homer <- file.path(outdir_raw, "09_HOMER_motifs")
dir.create(outdir_homer, recursive = TRUE, showWarnings = FALSE)

read_homer_known <- function(label) {
  f <- file.path(outdir_homer, label, "knownResults.txt")
  if (!file.exists(f)) {
    message("Not found: ", f)
    return(NULL)
  }
  read_tsv(f, show_col_types = FALSE) %>%
    setNames(c("motif_name","consensus","pvalue","log_pvalue",
               "q_value","n_target","pct_target","n_bg","pct_bg")) %>%
    mutate(
      label       = label,
      pct_target  = as.numeric(str_remove(pct_target, "%")),
      pct_bg      = as.numeric(str_remove(pct_bg,     "%")),
      log10p      = -log_pvalue / log(10),   # HOMER gives natural log
      fold_enrich = pct_target / pct_bg,
      # Clean up motif name — HOMER appends "/Motif Family/Accession"
      tf_name     = str_extract(motif_name, "^[^/]+")
    ) %>%
    arrange(log10p)
}

homer_sets <- list(
  "germline_male_biased"     = read_homer_known("germline_male_biased"),
  "germline_female_biased"   = read_homer_known("germline_female_biased"),
  "postmating_DO_upregulated" = read_homer_known("postmating_DO_upregulated")
)

# ---- Plot: top 15 known motifs per contrast --------------------------------
plot_homer <- function(df, title_text, n_top = 15, q_thresh = 0.05) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df_plot <- df %>%
    filter(q_value < q_thresh) %>%
    arrange(desc(log10p)) %>%
    slice_head(n = n_top) %>%
    mutate(tf_name = factor(tf_name, levels = rev(tf_name)))
  
  if (nrow(df_plot) == 0) {
    message("No significant motifs (q < ", q_thresh, ") for: ", title_text)
    return(NULL)
  }
  
  ggplot(df_plot, aes(tf_name, log10p, fill = fold_enrich)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    coord_flip() +
    scale_fill_gradientn(
      colours = enrichment_palette(100),
      name    = "Fold\nenrichment\n(% target /\n% background)"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = title_text,
      subtitle = paste0("Top ", nrow(df_plot),
                        " known motifs (q < ", q_thresh, ")"),
      x = NULL,
      y = expression(-log[10](p))
    ) +
    plot_theme +
    theme(legend.position = "right")
}

p_homer_male   <- plot_homer(homer_sets$germline_male_biased,
                              "Motifs enriched in testis-biased gene promoters")
p_homer_female <- plot_homer(homer_sets$germline_female_biased,
                              "Motifs enriched in ovary-biased gene promoters")
p_homer_pm     <- plot_homer(homer_sets$postmating_DO_upregulated,
                              "Motifs enriched in post-mating distal ovary gene promoters")

# Save individual plots
outdir_homer_plots <- file.path(outdir_raw, "09_HOMER_motifs", "plots")
dir.create(outdir_homer_plots, recursive = TRUE, showWarnings = FALSE)

for (nm in names(homer_sets)) {
  p <- plot_homer(homer_sets[[nm]], nm)
  if (!is.null(p))
    ggsave(file.path(outdir_homer_plots, paste0(nm, "_known_motifs.pdf")),
           p, width = 9, height = 6)
}

# ---- Combined patchwork panel ----------------------------------------------
plots_combined <- Filter(Negate(is.null),
                          list(p_homer_male, p_homer_female, p_homer_pm))

if (length(plots_combined) > 0) {
  p_homer_panel <- wrap_plots(plots_combined, ncol = 1) +
    plot_annotation(
      title    = "Transcription factor motif enrichment — HOMER",
      subtitle = "Known motifs in promoters of DE gene sets (1 kb upstream of TSS)",
      tag_levels = "A"
    )
  ggsave(file.path(outdir_homer_plots, "HOMER_motif_panel.pdf"),
         p_homer_panel,
         width  = 10,
         height = 5 * length(plots_combined),
         limitsize = FALSE)
  p_homer_panel
}
