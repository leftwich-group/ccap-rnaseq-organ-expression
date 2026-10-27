#!/usr/bin/env Rscript
# Q1b - MZT with intronic and exonic reads
# Input:  02_preprocessing/06_intron_exon_counts/results/ (exon_counts.txt,
#         intron_counts.txt, intron_lengths.tsv)
# Output: 04_results/14_MZT_intron_exon/
# Run from the repository root:  Rscript 03_analysis/scripts/14_mzt_intron_exon.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

indir <- Sys.getenv("INTRON_EXON_DIR",
                    here("02_preprocessing/06_intron_exon_counts/results/exon_intron_counts"))
inputs <- file.path(indir, c("exon_counts.txt", "intron_counts.txt", "intron_lengths.tsv"))
if (!all(file.exists(inputs))) stop(
  "Missing: ", paste(basename(inputs[!file.exists(inputs)]), collapse = ", "),
  "\n  Produced by 02_preprocessing/06_intron_exon_counts/ (see its README).",
  call. = FALSE)

outdir <- file.path(outdir_raw, RES$mzt_intron)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

min_intron_bp <- 300    # genes with less intronic sequence are not tested
qc_ratio_max  <- 0.2    # warn if the 0 h intronic fraction exceeds this share of the late-stage fraction
age_levels    <- c(1, 2, 3, 4, 8, 12, 24)

# ---- Load -----------------------------------------------------------------
read_fc <- function(path) {
  x <- read_tsv(path, comment = "#", show_col_types = FALSE)
  m <- as.matrix(x[, -(1:6)])
  rownames(m) <- clean_ids(x$Geneid)
  colnames(m) <- sub("\\.(markdup\\.)?sorted\\.bam$|\\.bam$", "", basename(colnames(m)))
  m
}
exon_all   <- read_fc(inputs[1])
intron_all <- read_fc(inputs[2])
intron_len <- read_tsv(inputs[3], show_col_types = FALSE) %>% mutate(gene_id = clean_ids(gene_id))

meta <- meta_raw %>%
  filter(tissue == "EG") %>%
  mutate(age   = as.numeric(age),
         stage = case_when(age == 0 ~ "Unfertilised",
                           age <= 4 ~ "Early",
                           TRUE     ~ "Late"),
         stage = factor(stage, levels = c("Unfertilised", "Early", "Late")),
         age_f = factor(paste0("h", age), levels = paste0("h", c(0, age_levels))))
if (!all(meta$id %in% colnames(exon_all)) || !all(meta$id %in% colnames(intron_all))) stop(
  "Embryo samples missing from the count tables: ",
  paste(setdiff(meta$id, intersect(colnames(exon_all), colnames(intron_all))), collapse = ", "),
  call. = FALSE)

common <- intersect(rownames(exon_all), rownames(intron_all))
exon   <- exon_all[common, meta$id]
intron <- intron_all[common, meta$id]

# ---- QC: intronic fraction per sample -------------------------------------
# Number of genes needed to account for half of a sample's intronic fragments.
genes_for_half <- apply(intron, 2, function(x) {
  v <- sort(x, decreasing = TRUE)
  which(cumsum(v) >= sum(v) / 2)[1]
})
frac <- meta %>%
  mutate(exon_reads   = colSums(exon),
         intron_reads = colSums(intron),
         intron_fraction = intron_reads / (intron_reads + exon_reads),
         genes_for_half  = genes_for_half[id])
write_tsv(frac %>% dplyr::select(id, age, stage, exon_reads, intron_reads, intron_fraction, genes_for_half),
          file.path(outdir, "intronic_fraction_by_sample.tsv"))

f_unfert <- median(frac$intron_fraction[frac$age == 0])
f_late   <- median(frac$intron_fraction[frac$age >= 8])
message(sprintf("Intronic fraction, median: 0 h = %.4f, 8-24 h = %.4f", f_unfert, f_late))
if (f_unfert > qc_ratio_max * f_late) warning(
  "The 0 h intronic fraction is above ", qc_ratio_max, " x the late-stage fraction. ",
  "Examine genomic DNA carry-over and unspliced maternal transcripts before using the intronic results.",
  call. = FALSE)

p_frac <- ggplot(frac, aes(factor(age), intron_fraction)) +
  ggdist::stat_halfeye(adjust = 1.5, width = 0.75, .width = 0, justification = -0.2,
                       point_colour = NA, fill = "grey80", slab_alpha = 0.8, normalize = "groups") +
  geom_point(aes(colour = genes_for_half),
             position = position_jitter(width = 0.05, height = 0, seed = 1), size = 2.8) +
  scale_colour_viridis_c(option = "plasma", trans = "log10", end = 0.9,
                         name = "Genes giving half\nof intronic fragments") +
  labs(x = "Hours post-oviposition", y = "Intronic fraction of gene-assigned fragments",
       title = "Intronic fraction by embryo timepoint") +
  plot_theme
ggsave(file.path(outdir, "intronic_fraction_by_timepoint.pdf"), p_frac, width = 6.5, height = 4.5)

# ---- DE on exon and intron counts -----------------------------------------
# Both layers use library sizes and TMM factors from exon + intron counts, so
# their fold changes share one scale.
dge_comb <- calcNormFactors(DGEList(exon + intron))
comb_lib <- setNames(dge_comb$samples$lib.size, colnames(exon))
comb_nf  <- setNames(dge_comb$samples$norm.factors, colnames(exon))

run_de <- function(counts, design, contrasts) {
  dge  <- DGEList(counts = counts)
  keep <- rowSums(dge$counts >= min_count) >= min_samples
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge$samples$lib.size     <- comb_lib[rownames(dge$samples)]
  dge$samples$norm.factors <- comb_nf[rownames(dge$samples)]
  dge  <- estimateDisp(dge, design, robust = TRUE)
  fit  <- glmQLFit(dge, design, robust = TRUE)
  lapply(contrasts, function(cn) {
    qlf <- glmQLFTest(fit, contrast = makeContrasts(contrasts = cn, levels = design))
    topTags(qlf, n = Inf)$table %>% rownames_to_column("gene")
  })
}

intron_genes <- intersect(common, intron_len$gene_id[intron_len$intron_bp >= min_intron_bp])
message("Genes with at least ", min_intron_bp, " bp of intron: ", length(intron_genes))

design_stage <- model.matrix(~ 0 + stage, data = meta)
colnames(design_stage) <- levels(meta$stage)
stage_contrasts <- c(early_vs_unfert = "Early - Unfertilised",
                     late_vs_unfert  = "Late - Unfertilised",
                     late_vs_early   = "Late - Early")

design_age <- model.matrix(~ 0 + age_f, data = meta)
colnames(design_age) <- levels(meta$age_f)
age_contrasts <- setNames(paste0("h", age_levels, " - h0"), paste0("h", age_levels))

exon_stage   <- run_de(exon[intron_genes, ],   design_stage, stage_contrasts)
intron_stage <- run_de(intron[intron_genes, ], design_stage, stage_contrasts)
exon_age     <- run_de(exon[intron_genes, ],   design_age,   age_contrasts)
intron_age   <- run_de(intron[intron_genes, ], design_age,   age_contrasts)

# ---- Intronic against exonic change per stage contrast --------------------
status <- function(fdr, lfc) case_when(fdr < fdr_thresh & lfc >=  lfc_thresh ~ "up",
                                       fdr < fdr_thresh & lfc <= -lfc_thresh ~ "down",
                                       TRUE ~ "ns")
class_levels <- c("Intron and exon up", "Intron up only", "Exon down, intron not up",
                  "Exon up, intron not up", "Other")

paired <- bind_rows(lapply(names(stage_contrasts), function(cn) {
  inner_join(exon_stage[[cn]]   %>% dplyr::select(gene, logFC_exon = logFC, FDR_exon = FDR),
             intron_stage[[cn]] %>% dplyr::select(gene, logFC_intron = logFC, FDR_intron = FDR),
             by = "gene") %>%
    mutate(contrast = factor(cn, levels = names(stage_contrasts)),
           exon_status   = status(FDR_exon,   logFC_exon),
           intron_status = status(FDR_intron, logFC_intron),
           class = case_when(
             intron_status == "up" & exon_status == "up"   ~ class_levels[1],
             intron_status == "up"                         ~ class_levels[2],
             exon_status == "down" & intron_status != "up" ~ class_levels[3],
             exon_status == "up"   & intron_status != "up" ~ class_levels[4],
             TRUE                                          ~ class_levels[5]),
           class = factor(class, levels = class_levels))
}))
write_tsv(paired, file.path(outdir, "exon_vs_intron_by_stage_contrast.tsv"))
write_tsv(paired %>% dplyr::count(contrast, class) %>% pivot_wider(names_from = class, values_from = n, values_fill = 0),
          file.path(outdir, "exon_vs_intron_class_counts.tsv"))

class_cols <- c("Intron and exon up" = "#D55E00", "Intron up only" = "#E69F00",
                "Exon down, intron not up" = "#0072B2", "Exon up, intron not up" = "#56B4E9",
                "Other" = "grey80")
p_scatter <- ggplot(paired %>% arrange(class != "Other"),
                    aes(logFC_exon, logFC_intron, colour = class)) +
  geom_hline(yintercept = c(-lfc_thresh, lfc_thresh), linetype = 2, colour = "grey60") +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = 2, colour = "grey60") +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ contrast, nrow = 1,
             labeller = as_labeller(c(early_vs_unfert = "Early vs unfertilised",
                                      late_vs_unfert  = "Late vs unfertilised",
                                      late_vs_early   = "Late vs early"))) +
  scale_colour_manual(values = class_cols) +
  labs(x = "Exonic log2 fold change", y = "Intronic log2 fold change", colour = NULL) +
  plot_theme + theme(legend.position = "bottom")
ggsave(file.path(outdir, "exon_vs_intron_logFC.pdf"), p_scatter, width = 11, height = 4.8)

# ---- Onset: genes above the 0 h level at each timepoint -------------------
count_dir <- function(res_list, layer) {
  bind_rows(lapply(names(res_list), function(cn) {
    r <- res_list[[cn]]
    tibble(layer = layer, age = as.numeric(sub("h", "", cn)),
           n_up   = sum(r$FDR < fdr_thresh & r$logFC >=  lfc_thresh),
           n_down = sum(r$FDR < fdr_thresh & r$logFC <= -lfc_thresh))
  }))
}
onset <- bind_rows(count_dir(exon_age, "Exon"), count_dir(intron_age, "Intron"))
write_tsv(onset, file.path(outdir, "genes_changed_vs_0h_by_timepoint.tsv"))

p_onset <- onset %>%
  pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n_genes") %>%
  mutate(direction = factor(direction, levels = c("n_up", "n_down"),
                            labels = c("Up vs unfertilised", "Down vs unfertilised"))) %>%
  ggplot(aes(age, n_genes, colour = layer)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  facet_wrap(~ direction) +
  scale_x_continuous(breaks = c(1, 4, 8, 12, 24)) +
  scale_colour_manual(values = c(Exon = "#0072B2", Intron = "#D55E00")) +
  labs(x = "Hours post-oviposition", y = "Genes (FDR < 0.05, at least 2-fold)", colour = NULL) +
  plot_theme
ggsave(file.path(outdir, "genes_changed_vs_0h_by_timepoint.pdf"), p_onset, width = 8.5, height = 4.2)

saveRDS(list(fraction = p_frac, onset = p_onset, scatter = p_scatter),
        file.path(outdir, "plots.rds"))
saveRDS(list(exon_stage = exon_stage, intron_stage = intron_stage,
             exon_age = exon_age, intron_age = intron_age),
        file.path(outdir, "results_nested.rds"))
message("Done: ", outdir)
