#!/usr/bin/env Rscript
# UpSet plot - sex-biased gene overlap across tissues
# Output: 04_results/08_UpSet_sex_bias/   (reads 04_results/02_Sexual_dimorphism/)
# Run from the repository root:  Rscript 03_analysis/scripts/08_upset_sex_bias.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir_upset <- file.path(outdir_raw, "08_UpSet_sex_bias")
dir.create(outdir_upset, recursive = TRUE, showWarnings = FALSE)

# Load all Q2 contrast results (short genome)
contrast_files <- list(
  Germline   = file.path(outdir_raw, "02_Sexual_dimorphism/contrast_1/results_all_nested.rds"),
  Abdomen_V  = file.path(outdir_raw, "02_Sexual_dimorphism/contrast_3/results_all_nested.rds"),
  Abdomen_M  = file.path(outdir_raw, "02_Sexual_dimorphism/contrast_4/results_all_nested.rds")
)

de_sets <- lapply(contrast_files, function(f) {
  if (!file.exists(f)) return(list(male_biased = character(), female_biased = character()))
  res_obj <- readRDS(f)
  cname   <- names(res_obj$short$contrasts)[1]
  res     <- res_obj$short$contrasts[[cname]]$table
  list(
    male_biased   = res %>% filter(FDR < fdr_thresh, logFC < -lfc_thresh) %>%
      pull(gene),
    female_biased = res %>% filter(FDR < fdr_thresh, logFC >  lfc_thresh) %>%
      pull(gene)
  )
})

# Build UpSet input list
results_list <- list()

for (direction in c("male_biased", "female_biased")) {
  set_list <- lapply(de_sets, `[[`, direction)
  set_list <- Filter(function(x) length(x) > 0, set_list)

  if (length(set_list) < 2) next

   p<- upset(fromList(set_list),
        order.by  = "freq",
        sets.bar.color = "#D55E00",
        main.bar.color = "#0072B2",
        text.scale = 1.4,
        mainbar.y.label = paste(direction, "gene count"),
        sets.x.label    = "Total DE genes")
  
  pdf(file.path(outdir_upset, paste0("UpSet_", direction, ".pdf")),
      width = 10, height = 6)
print(p)
  dev.off()

  # Genes shared across ALL tissues (most stringent sex-biased set)
  shared_all <- Reduce(intersect, set_list)
  message(direction, " — shared across all tissues: ", length(shared_all))
  
  if (length(shared_all) > 0) {
    write_lines(shared_all,
                file.path(outdir_upset,
                          paste0(direction, "_shared_all_tissues.txt")))
    # GO for the shared core
    go_shared <- run_topgo_all_ont(clean_ids(shared_all),
                                   gene_universe_short, geneID2GO_short)
    p_shared  <- plot_topgo_panel(go_shared,
                                  paste0(direction, " — shared across tissues"))
    if (!is.null(p_shared))
      ggsave(file.path(outdir_upset,
                       paste0("GO_", direction, "_shared_all.pdf")),
             p_shared, width = 8, height = 10)
    
          summary_cls <- plot_topgo_bar(go_shared)
       if (!is.null(summary_cls))
        ggsave(file.path(outdir_upset, paste0( "_GO_Summary_",direction ,".pdf")),
               summary_cls, width = 7, height = 6)
      
      # save summary file 
       topgo_df <- go_shared$BP
  if (is.null(topgo_df) || nrow(topgo_df) == 0) return(NULL)
 
  df <- topgo_df %>%
    mutate(p_value     = parse_topgo_p(weight01_fisher),
           Significant = as.numeric(Significant)) %>%
    filter(!is.na(0.05), !is.na(Significant)) %>%
    arrange(p_value) %>%
    slice_head(n = 12) 
 
  if (nrow(df) == 0) return(NULL)
 
  # save summary file 
      write_tsv(df,file.path(outdir_upset, paste0( "_GO_Summary_", direction , ".tsv")))
  
      
      results_list[[direction]] <- list (
        p = p,
        p_shared = p_shared,
        summary_cls = summary_cls
      )
  }
}

saveRDS(results_list, file.path(outdir_upset, paste0("go_summary_upset.RDS")))
