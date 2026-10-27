#!/usr/bin/env Rscript
# Extended analyses: per-organ GO enrichment, immune-gene overlap test,
# mating x tissue interaction, Drosophila tissue-programme overlap.
# Output: 04_results/13_Extended_analyses/
# Run from the repository root:  Rscript 03_analysis/scripts/13_extended_analyses.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

outdir_ext <- file.path(outdir_raw, "13_Extended_analyses")
dir.create(outdir_ext, recursive = TRUE, showWarnings = FALSE)

adult_tissues <- c(AG = "Accessory_Gland", AB = "Abdomen", TE = "Testes",
                    YO = "Distal_Ovary", OO = "Proximal_Ovary")

# ===========================================================================
# Part A - per-organ GO enrichment of tissue-private genes
# ===========================================================================
outdir_a <- file.path(outdir_ext, "A_per_organ_GO")
dir.create(outdir_a, recursive = TRUE, showWarnings = FALSE)

expressed_by_tissue <- lapply(names(adult_tissues), function(tc) {
  ids <- meta_raw %>% filter(tissue == tc) %>% pull(id)
  ids <- intersect(ids, colnames(counts_short_raw))
  if (length(ids) == 0) return(character(0))
  dge  <- DGEList(counts = counts_short_raw[, ids, drop = FALSE])
  keep <- rowSums(cpm(dge) >= 1) >= 2
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge)
  lcpm <- cpm(dge, log = TRUE, prior.count = 2)
  rownames(lcpm)[rowMeans(lcpm) > 0]
}) %>% setNames(adult_tissues)

private_by_tissue <- lapply(names(expressed_by_tissue), function(tc) {
  others <- unname(unlist(expressed_by_tissue[names(expressed_by_tissue) != tc]))
  expressed_by_tissue[[tc]][!expressed_by_tissue[[tc]] %in% others]
}) %>% setNames(names(expressed_by_tissue))

message("Tissue-private gene counts:")
print(sapply(private_by_tissue, length))

go_private_by_tissue <- lapply(names(private_by_tissue), function(tc) {
  genes <- clean_ids(private_by_tissue[[tc]])
  if (length(genes) < 2) return(NULL)
  run_topgo_all_ont(genes, gene_universe_short, geneID2GO_short)
}) %>% setNames(names(private_by_tissue))

panel_letters <- LETTERS[seq_along(go_private_by_tissue)]
private_go_plots <- Map(function(tissue_name, go_res, letter) {
  if (is.null(go_res)) return(NULL)
  plot_topgo_bar(go_res, title_text = paste0(letter, " - ", gsub("_", " ", tissue_name)),
                 n_show = 8)
}, names(go_private_by_tissue), go_private_by_tissue, panel_letters)
private_go_plots <- Filter(Negate(is.null), private_go_plots)

if (length(private_go_plots) > 0) {
  p_private_go <- wrap_plots(private_go_plots, ncol = 2) +
    plot_annotation(title = "Per-organ GO enrichment (BP) of tissue-private genes")
  ggsave(file.path(outdir_a, "per_organ_private_GO.pdf"), p_private_go,
         width = 14, height = 5 * ceiling(length(private_go_plots) / 2),
         limitsize = FALSE)
  ggsave(file.path(outdir_a, "per_organ_private_GO.png"), p_private_go,
         width = 14, height = 5 * ceiling(length(private_go_plots) / 2),
         dpi = 300, limitsize = FALSE)
} else {
  message("Part A: no tissue had >=2 private genes with GO terms.")
}

go_private_tsv <- Map(function(tissue_name, go_res) {
  if (is.null(go_res)) return(NULL)
  bind_rows(go_res, .id = "ontology") %>% mutate(tissue = tissue_name)
}, names(go_private_by_tissue), go_private_by_tissue) %>%
  Filter(Negate(is.null), .) %>%
  bind_rows()
if (nrow(go_private_tsv) > 0) {
  write_tsv(go_private_tsv, file.path(outdir_a, "per_organ_private_GO_terms.tsv"))
}

# ===========================================================================
# Part B - immune-gene overlap test (post-mating distal ovary)
# ===========================================================================
outdir_b <- file.path(outdir_ext, "B_immune_overlap")
dir.create(outdir_b, recursive = TRUE, showWarnings = FALSE)

# Terms resolved by exact name from GO.db, not hardcoded IDs.
immune_term_names <- c(
  "immune response",
  "innate immune response",
  "defense response to other organism",
  "Toll signaling pathway",
  "antimicrobial humoral response",
  "defense response to Gram-positive bacterium",
  "defense response to Gram-negative bacterium",
  "peptidoglycan recognition protein signaling pathway"
)

go_bp_terms <- AnnotationDbi::select(
  GO.db, keys = keys(GO.db, "GOID"),
  columns = c("TERM", "ONTOLOGY"), keytype = "GOID"
) %>% filter(ONTOLOGY == "BP")

immune_go <- go_bp_terms %>% filter(TERM %in% immune_term_names)
if (nrow(immune_go) != length(immune_term_names)) {
  stop("Immune GO term resolution incomplete: expected ",
       length(immune_term_names), " terms, matched ", nrow(immune_go),
       ". Check immune_term_names against GO.db.")
}
write_tsv(immune_go, file.path(outdir_b, "immune_go_terms_used.tsv"))
immune_go_ids <- immune_go$GOID

immune_genes <- names(geneID2GO_short)[
  vapply(geneID2GO_short, function(terms) any(terms %in% immune_go_ids), logical(1))
]
message("Immune gene set size: ", length(immune_genes))

distal_ov_rds <- res_dir(RES$postmating, "contrast_3", "results_all_nested.rds")
if (file.exists(distal_ov_rds)) {
  distal_ov_res <- load_short_res(distal_ov_rds, "mated_vs_virgin_YO")
  if (is.null(distal_ov_res)) {
    distal_ov_res <- readRDS(distal_ov_rds)$short$contrasts[[1]]$table
  }

  mated_up <- distal_ov_res %>%
    filter(FDR < fdr_thresh, logFC > lfc_thresh) %>%
    pull(gene) %>% clean_ids()

  universe        <- gene_universe_short
  in_de_immune    <- length(intersect(mated_up, immune_genes))
  in_de_notimmune <- length(mated_up) - in_de_immune
  notde_immune    <- length(intersect(setdiff(universe, mated_up), immune_genes))
  notde_notimmune <- length(universe) - length(mated_up) - notde_immune

  contingency <- matrix(
    c(in_de_immune, notde_immune, in_de_notimmune, notde_notimmune),
    nrow = 2,
    dimnames = list(c("Immune", "Not immune"), c("Mated-up", "Not mated-up"))
  )
  ft <- fisher.test(contingency)
  message("Immune overlap: OR = ", round(unname(ft$estimate), 2),
          ", p = ", signif(ft$p.value, 3))

  overlap_summary <- tibble(
    contrast          = "distal_ovary_mated_vs_virgin_up",
    n_mated_up        = length(mated_up),
    n_immune_universe = length(immune_genes),
    n_overlap         = in_de_immune,
    odds_ratio        = unname(ft$estimate),
    p_value           = ft$p.value,
    ci_low            = ft$conf.int[1],
    ci_high           = ft$conf.int[2]
  )
  write_tsv(overlap_summary, file.path(outdir_b, "immune_overlap_fisher_test.tsv"))
  write_tsv(tibble(gene = intersect(mated_up, immune_genes)),
            file.path(outdir_b, "immune_overlap_genes.tsv"))
} else {
  message("Part B skipped: ", distal_ov_rds, " not found. Run 03_postmating.R first.")
}

# ===========================================================================
# Part C - mating x tissue interaction (ovary vs abdomen)
# ===========================================================================
outdir_c <- file.path(outdir_ext, "C_mating_tissue_interaction")
dir.create(outdir_c, recursive = TRUE, showWarnings = FALSE)

meta_int <- meta_raw %>%
  filter(sex == "F", experiment == "1",
         tissue %in% c("AB", "YO", "OO")) %>%
  mutate(
    tissue_grp = if_else(tissue == "AB", "Abdomen", "Ovary"),
    mating_grp = if_else(mating_status == "M", "Mated", "Virgin"),
    group      = paste(tissue_grp, mating_grp, sep = "_")
  )

counts_int <- counts_short_raw[, meta_int$id, drop = FALSE]
design_int <- model.matrix(~ 0 + group, data = meta_int)
colnames(design_int) <- sub("^group", "", colnames(design_int))

interaction_contrast <- list(
  mating_x_tissue = "(Ovary_Mated - Ovary_Virgin) - (Abdomen_Mated - Abdomen_Virgin)"
)

results_interaction <- run_reference_analysis(
  counts_int, meta_int, design_int, interaction_contrast,
  "short", outdir_c, "group",
  min_count, min_samples, fdr_thresh, lfc_thresh
)
saveRDS(results_interaction, file.path(outdir_c, "results_interaction.rds"))
append_de_stats(list(short = results_interaction), "Extended",
                 "Mating x tissue interaction (ovary vs abdomen)",
                 interaction_contrast$mating_x_tissue, genome = "short")

int_res <- results_interaction$contrasts$mating_x_tissue$table
n_int_sig <- sum(int_res$FDR < fdr_thresh, na.rm = TRUE)
message("Mating x tissue interaction: ", n_int_sig, " / ", nrow(int_res),
        " genes FDR < ", fdr_thresh)

if (n_int_sig >= 2) {
  int_genes <- int_res %>% filter(FDR < fdr_thresh) %>% pull(gene) %>% clean_ids()
  go_interaction <- run_topgo_all_ont(int_genes, gene_universe_short, geneID2GO_short)
  write_tsv(bind_rows(go_interaction, .id = "ontology"),
            file.path(outdir_c, "interaction_GO_terms.tsv"))
}

# ===========================================================================
# Part D - Drosophila FlyAtlas2 tissue-programme overlap
# ===========================================================================
outdir_d <- file.path(outdir_ext, "D_flyatlas2_tissue_overlap")
dir.create(outdir_d, recursive = TRUE, showWarnings = FALSE)

# flyatlas2_tissue_enrichment.tsv (columns: fly_symbol, tissue, enriched) is
# built here from the mean-FPKM gene-data workbook, which has no enrichment
# call of its own. A gene is "enriched" in a tissue if its FPKM there is
# >= 4x the mean FPKM across the other adult tissue columns (+1 pseudocount).
flyatlas2_tissue_file <- here("01_data/flyatlas2_tissue_enrichment.tsv")
flyatlas2_xlsx        <- here("01_data/FlyAtlas2_gene_data_2025.xlsx")
flyatlas2_xlsx_url    <- "https://motif.mvls.gla.ac.uk/downloads/FlyAtlas2_gene_data_2025.xlsx"
flyatlas2_enrich_thresh <- 4

if (!file.exists(flyatlas2_tissue_file) && !file.exists(flyatlas2_xlsx)) {
  dl_ok <- tryCatch({
    download.file(flyatlas2_xlsx_url, flyatlas2_xlsx, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!dl_ok) message("Part D: could not download ", flyatlas2_xlsx_url)
}

if (!file.exists(flyatlas2_tissue_file) && file.exists(flyatlas2_xlsx)) {
  fa2 <- readxl::read_excel(flyatlas2_xlsx, sheet = "FPKMs 2025")

  fa2_adult_cols <- c(
    "Carcass M", "Carcass F", "Head M", "Head F", "Eye M", "Eye F",
    "Brain M", "Brain F", "Ganglion M", "Ganglion F", "Crop M", "Crop F",
    "Midgut M", "Midgut F", "Hindgut M", "Hindgut F", "Tubule M", "Tubule F",
    "RectalPad M", "RectalPad F", "Fat Body M", "Fat Body F", "Heart M",
    "Heart F", "Salivary M", "Salivary F", "Ovary", "V Sp.", "M Sp.",
    "Testis", "Accessory Gland"
  )
  fa2_mat <- as.matrix(fa2[, fa2_adult_cols])

  fa2_tissue_values <- list(
    testis           = fa2_mat[, "Testis"],
    `accessory gland` = fa2_mat[, "Accessory Gland"],
    midgut           = rowMeans(fa2_mat[, c("Midgut M", "Midgut F")], na.rm = TRUE)
  )

  flyatlas2_tissue <- lapply(names(fa2_tissue_values), function(ft) {
    fpkm_this  <- fa2_tissue_values[[ft]]
    other_mean <- rowMeans(fa2_mat[, setdiff(fa2_adult_cols, ft), drop = FALSE], na.rm = TRUE)
    tibble(fly_symbol = fa2$Symbol, tissue = ft,
           enriched = (fpkm_this / (other_mean + 1)) >= flyatlas2_enrich_thresh)
  }) %>% bind_rows()

  write_tsv(flyatlas2_tissue, flyatlas2_tissue_file)
}

eggnog_raw_file <- here("01_data/GO_terms/eggNOG-mapper/medfly_short.emapper.annotations")
gff_file        <- here("01_data/Genomic_files/GCF_000347755.3_Ccap_2.1_genomic.gff.gz")

# eggNOG-mapper's own gene identifier is a protein accession, mapped to the
# RNA-seq count-matrix gene ID via the GFF's protein_id/gene attributes
# (same two-step mapping as the FlyAtlas2 housekeeping comparison above).
prot2name <- tibble()
if (file.exists(eggnog_raw_file)) {
  egg <- read_tsv(eggnog_raw_file, comment = "##", show_col_types = FALSE, skip = 3)
  name_col <- names(egg)[str_detect(names(egg), regex("preferred.name", ignore_case = TRUE))][1]
  prot2name <- egg %>%
    dplyr::select(protein = 1, preferred_name = all_of(name_col)) %>%
    filter(!is.na(preferred_name), preferred_name != "-", preferred_name != "")
}

prot2gene <- tibble()
if (file.exists(gff_file)) {
  prot2gene <- read_tsv(
    gff_file, comment = "#",
    col_names = c("seqid", "source", "type", "start", "end",
                  "score", "strand", "phase", "attributes"),
    show_col_types = FALSE
  ) %>%
    filter(type == "CDS") %>%
    mutate(
      protein = str_extract(attributes, "(?<=protein_id=)[^;]+"),
      gene    = str_extract(attributes, "(?<=gene=)[^;]+"),
      locus   = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
      gene_id = coalesce(gene, locus)
    ) %>%
    filter(!is.na(protein), !is.na(gene_id)) %>%
    distinct(protein, gene_id)
}

name2gene <- prot2name %>%
  inner_join(prot2gene, by = "protein") %>%
  mutate(name_upper = toupper(preferred_name)) %>%
  distinct(name_upper, gene_id)

if (file.exists(flyatlas2_tissue_file) && nrow(name2gene) > 0) {
  flyatlas2_tissue <- read_tsv(flyatlas2_tissue_file, show_col_types = FALSE)

  tissue_fly_map <- c(Testes = "testis", Accessory_Gland = "accessory gland",
                       Abdomen = "midgut")

  overlap_rows <- lapply(names(tissue_fly_map), function(tc) {
    fly_tissue <- tissue_fly_map[[tc]]
    fly_genes_this_tissue <- flyatlas2_tissue %>%
      filter(tissue == fly_tissue, enriched) %>%
      pull(fly_symbol) %>% toupper()

    ccap_orthologues <- name2gene %>%
      filter(name_upper %in% fly_genes_this_tissue) %>%
      pull(gene_id) %>% clean_ids()

    private_this <- clean_ids(private_by_tissue[[tc]])
    overlap_n <- length(intersect(private_this, ccap_orthologues))

    universe <- gene_universe_short
    a <- overlap_n
    b <- length(private_this) - a
    c_ <- length(intersect(setdiff(universe, private_this), ccap_orthologues))
    d <- length(universe) - length(private_this) - c_
    ft <- fisher.test(matrix(c(a, c_, b, d), nrow = 2))

    tibble(tissue = tc, fly_tissue = fly_tissue,
           n_private = length(private_this),
           n_fly_tissue_orthologues = length(ccap_orthologues),
           n_overlap = a, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
  }) %>% bind_rows()

  write_tsv(overlap_rows, file.path(outdir_d, "flyatlas2_tissue_overlap.tsv"))
} else {
  message("Part D skipped: missing FlyAtlas2 enrichment data, eggNOG annotations, or GFF.")
}

# ===========================================================================
# Part E - direct gene-list overlap with published post-mating studies
# ===========================================================================
outdir_e <- file.path(outdir_ext, "E_literature_gene_overlap")
dir.create(outdir_e, recursive = TRUE, showWarnings = FALSE)

# Short et al. 2012 (J Insect Physiol) measured exactly these six
# antimicrobial peptide genes by qPCR; it is not a genome-wide gene list.
short2012_genes <- c("AttA", "AttB", "Mtk", "Def", "Drs", "DptA")

if (nrow(name2gene) > 0) {
  short2012_orthologues <- name2gene %>%
    filter(name_upper %in% toupper(short2012_genes)) %>%
    pull(gene_id) %>% clean_ids()

  postmating_rds <- c(
    Abdomen        = res_dir(RES$postmating, "contrast_1", "results_all_nested.rds"),
    Proximal_Ovary = res_dir(RES$postmating, "contrast_2", "results_all_nested.rds"),
    Distal_Ovary   = res_dir(RES$postmating, "contrast_3", "results_all_nested.rds")
  )
  postmating_contrast <- c(
    Abdomen        = "mated_vs_virgin_abdomen",
    Proximal_Ovary = "mated_vs_virgin_OO",
    Distal_Ovary   = "mated_vs_virgin_YO"
  )

  short2012_overlap <- lapply(names(postmating_rds), function(tc) {
    rds <- postmating_rds[[tc]]
    if (!file.exists(rds)) return(NULL)
    res <- load_short_res(rds, postmating_contrast[[tc]])
    if (is.null(res)) res <- readRDS(rds)$short$contrasts[[1]]$table
    mated_up <- res %>% filter(FDR < fdr_thresh, logFC > lfc_thresh) %>% pull(gene) %>% clean_ids()
    tibble(tissue = tc, n_mated_up = length(mated_up),
           n_short2012_orthologues = length(short2012_orthologues),
           n_overlap = length(intersect(mated_up, short2012_orthologues)),
           overlap_genes = paste(intersect(mated_up, short2012_orthologues), collapse = ";"))
  }) %>% bind_rows()

  write_tsv(short2012_overlap, file.path(outdir_e, "short2012_AMP_overlap.tsv"))
} else {
  message("Part E (Short et al. 2012) skipped: name2gene mapping unavailable.")
}

# McGraw et al. 2004 supplementary tables (Curr Biol), one per mating
# component, each with its own "Immune/Defense"-type category already
# assigned by the authors -- used directly rather than re-derived.
mcgraw2004_dir <- here("01_data/McGraw_2004_data")
mcgraw2004_sheets <- tibble::tribble(
  ~file,                                   ~immune_category,
  "1-s2.0-S0960982204006098-mmc2.xls",     "Immune/Defense",
  "1-s2.0-S0960982204006098-mmc3.xls",     "Immune response",
  "1-s2.0-S0960982204006098-mmc4.xls",     "Defense/Immunity"
)

extract_mcgraw_category <- function(path, category) {
  d <- readxl::read_excel(path, sheet = "Sheet1", col_names = FALSE)
  hdr_row <- which(apply(d, 1, function(r) any(grepl("^FBID$", r, ignore.case = TRUE))))[1]
  col_names <- as.character(d[hdr_row, ])
  fbid_col   <- which(col_names == "FBID")
  name_col   <- which(col_names == "name")
  updown_col <- which(col_names == "up/down regulated")

  is_header <- !is.na(d[[1]]) & is.na(d[[2]]) & is.na(d[[3]])
  header_rows <- which(is_header)
  idx <- which(d[[1]][header_rows] == category)
  if (length(idx) == 0) return(tibble())
  start <- header_rows[idx] + 1
  end   <- if (idx < length(header_rows)) header_rows[idx + 1] - 1 else nrow(d)

  tibble(fbid = d[[fbid_col]][start:end],
         name = d[[name_col]][start:end],
         direction = d[[updown_col]][start:end])
}

mcgraw2004_immune <- purrr::pmap_dfr(mcgraw2004_sheets, function(file, immune_category) {
  path <- file.path(mcgraw2004_dir, file)
  if (!file.exists(path)) return(tibble())
  extract_mcgraw_category(path, immune_category) %>% mutate(source_file = file)
})

if (nrow(mcgraw2004_immune) > 0 && nrow(name2gene) > 0) {
  write_tsv(mcgraw2004_immune, file.path(outdir_e, "mcgraw2004_immune_genes.tsv"))

  mcgraw2004_orthologues <- name2gene %>%
    filter(name_upper %in% toupper(na.omit(mcgraw2004_immune$name))) %>%
    pull(gene_id) %>% clean_ids()

  mcgraw2004_overlap <- lapply(names(postmating_rds), function(tc) {
    rds <- postmating_rds[[tc]]
    if (!file.exists(rds)) return(NULL)
    res <- load_short_res(rds, postmating_contrast[[tc]])
    if (is.null(res)) res <- readRDS(rds)$short$contrasts[[1]]$table
    mated_up <- res %>% filter(FDR < fdr_thresh, logFC > lfc_thresh) %>% pull(gene) %>% clean_ids()
    tibble(tissue = tc, n_mated_up = length(mated_up),
           n_mcgraw2004_immune_orthologues = length(mcgraw2004_orthologues),
           n_overlap = length(intersect(mated_up, mcgraw2004_orthologues)),
           overlap_genes = paste(intersect(mated_up, mcgraw2004_orthologues), collapse = ";"))
  }) %>% bind_rows()

  write_tsv(mcgraw2004_overlap, file.path(outdir_e, "mcgraw2004_immune_overlap.tsv"))
} else {
  message("Part E (McGraw et al. 2004) skipped: no immune-category rows extracted, ",
          "or name2gene mapping unavailable.")
}

message("Extended analyses complete: ", outdir_ext)
