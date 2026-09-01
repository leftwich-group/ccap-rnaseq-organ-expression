# ===========================================================================
# R/load_data.R
# Reads the count matrices, sample metadata, and GO annotations once, sets the
# global DE thresholds, and defines append_de_stats(). Sourced via R/init.R.
# ===========================================================================

# ---------------------------------------------------------------------------
# Global DE thresholds
# ---------------------------------------------------------------------------
min_count   <- 10
min_samples <- 2
fdr_thresh  <- 0.05
lfc_thresh  <- 1          # log2FC threshold

outdir_raw  <- here("04_results")
dir.create(outdir_raw, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load count matrices
# ---------------------------------------------------------------------------
counts_short_raw <- read_tsv(
  here("02_preprocessing/05_results_summary/short_read_genome/salmon.merged.gene_counts.tsv"),
  show_col_types = FALSE) %>%
  column_to_rownames(names(.)[1]) %>%
  as.matrix()
mode(counts_short_raw) <- "numeric"

counts_long_raw <- read_tsv(
  here("02_preprocessing/05_results_summary/long_read_genome/salmon.merged.gene_counts.tsv"),
  show_col_types = FALSE) %>%
  column_to_rownames(names(.)[1]) %>%
  as.matrix()
mode(counts_long_raw) <- "numeric"

# Gene name corrections for long genome (strip CCAP1982_ prefix from rownames)
rownames(counts_long_raw) <- sub("^CCAP1982_", "", rownames(counts_long_raw))

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
meta_raw <- read_csv(here("01_data/metadata.csv"), show_col_types = FALSE) %>%
  mutate(across(where(is.character), trimws)) %>%
  # experiment-2 abdomens (RN34-RN36, all female) are contaminated libraries
  # and are held out of every downstream section
  filter(!(tissue == "AB" & experiment == 2))

# ---------------------------------------------------------------------------
# Short ↔ long gene ID mapping
# ---------------------------------------------------------------------------
gene_map_file <- here("data/gene_map_short_to_long.tsv")
gene_map_file <-(here("data/Genomic_files/old/gene_mapping_rbh_raw.tsv"))
# ---------------------------------------------------------------------------
# GO annotations (pre-built in Section 2)
# ---------------------------------------------------------------------------
geneID2GO_short <- readRDS(here("02_preprocessing/04_go_annotation/geneID2GO_short.rds"))
geneID2GO_long  <- readRDS(here("02_preprocessing/04_go_annotation/geneID2GO_long.rds"))
gene_universe_short <- names(geneID2GO_short)
gene_universe_long  <- names(geneID2GO_long)

message("Genes in short count matrix: ", nrow(counts_short_raw))
message("Genes in long count matrix:  ", nrow(counts_long_raw))
message("Samples in metadata:          ", nrow(meta_raw))

# ---- DE statistics accumulator -------------------------------------------
# ===========================================================================
# RUNNING DE STATISTICS TABLE
# ---------------------------------------------------------------------------
# Every contrast in this analysis appends one row to a single master TSV.
# This chunk initialises the file (deleting any stale version from a previous
# run) and defines the helper.  The helper is called at the end of every
# summary panel chunk, so the table grows automatically as you run sections.
#
# Output: results/final_output/DE_stats_master.tsv
#
# Columns:
#   section         – e.g. "Q1", "Q2", "Tissue_landscape", "Y_chromosome"
#   question        – human-readable label
#   genome          – "short" | "long"
#   contrast        – contrast name from run_reference_analysis
#   group_A / group_B – the two groups being compared
#   n_tested        – total genes tested
#   n_up            – genes up in group_A (logFC > lfc_thresh, FDR < fdr_thresh)
#   n_down          – genes up in group_B
#   n_sig_total     – all DE genes
#   pct_sig         – % of tested genes that are DE
# ===========================================================================


de_master_file <- file.path(outdir_raw, "DE_stats_master.tsv")

# ---------------------------------------------------------------------------
# Column spec — used both for initialisation and every read-back.
# Keeping this in one place prevents the type-mismatch error that occurs when
# read_tsv() guesses columns as <chr> on an empty file.
# ---------------------------------------------------------------------------
de_master_cols <- cols(
  section     = col_character(),
  question    = col_character(),
  genome      = col_character(),
  contrast    = col_character(),
  group_A     = col_character(),
  group_B     = col_character(),
  n_tested    = col_integer(),
  n_up        = col_integer(),
  n_down      = col_integer(),
  n_sig_total = col_integer(),
  pct_sig     = col_double()
)

# Initialise with a typed empty tibble so the header row has the right types.
# Comment out the write_tsv line to accumulate across sessions instead.
de_master_empty <- tibble(
  section     = character(), question    = character(),
  genome      = character(), contrast    = character(),
  group_A     = character(), group_B     = character(),
  n_tested    = integer(),   n_up        = integer(),
  n_down      = integer(),   n_sig_total = integer(),
  pct_sig     = numeric()
)

# ---------------------------------------------------------------------------
# append_de_stats()
#   results_obj     — the list returned by run_reference_analysis
#                     (or a list with $short / $long elements)
#   section         — e.g. "Q1", "Q2"
#   question        — plain-English contrast description
#   contrast_string — the contrast formula, e.g. "Ovary - Testes"
#   genome          — "short" | "long"
# ---------------------------------------------------------------------------
append_de_stats <- function(results_obj,
                             section,
                             question,
                             contrast_string,
                             genome = "short") {

  # Defensive: accept results objects that have the genome nested directly
  # OR wrapped in $short / $long
  genome_obj <- if (!is.null(results_obj[[genome]])) {
    results_obj[[genome]]
  } else {
    results_obj   # already the genome-level object
  }

  if (is.null(genome_obj$contrasts) || length(genome_obj$contrasts) == 0) {
    message("append_de_stats: no contrasts found for ", genome, " — skipping.")
    return(invisible(NULL))
  }

  if (grepl("MZT", question)) {
    cname <- contrast_string
    
  } else {
    cname <- names(genome_obj$contrasts)[1]
  }
  
  res   <- genome_obj$contrasts[[cname]]$table

  if (is.null(res) || nrow(res) == 0) {
    message("append_de_stats: empty results table — skipping.")
    return(invisible(NULL))
  }

  # Split "GroupA - GroupB" cleanly (handles spaces around dash)
  parts <- trimws(strsplit(contrast_string, "\\s*-\\s*", perl = TRUE)[[1]])
  grp_A <- parts[1]
  grp_B <- if (length(parts) >= 2) paste(parts[-1], collapse = " - ") else ""

  n_up    <- sum(res$FDR < fdr_thresh & res$logFC >  lfc_thresh, na.rm = TRUE)
  n_down  <- sum(res$FDR < fdr_thresh & res$logFC < -lfc_thresh, na.rm = TRUE)
  n_tot   <- sum(res$FDR < fdr_thresh,                            na.rm = TRUE)
  n_test  <- nrow(res)

  new_row <- tibble(
    section     = as.character(section),
    question    = as.character(question),
    genome      = as.character(genome),
    contrast    = as.character(cname),
    group_A     = as.character(grp_A),
    group_B     = as.character(grp_B),
    n_tested    = as.integer(n_test),
    n_up        = as.integer(n_up),
    n_down      = as.integer(n_down),
    n_sig_total = as.integer(n_tot),
    pct_sig     = round(as.numeric(n_tot) / as.numeric(n_test) * 100, 2)
  )

  # Read back with explicit col_types so integers stay integers on empty files
  existing <- if (file.exists(de_master_file)) {
    read_tsv(de_master_file, col_types = de_master_cols, show_col_types = FALSE)
  } else de_master_empty

  existing <- existing %>%
    dplyr::filter(!(section  == new_row$section  &
                    genome   == new_row$genome   &
                    contrast == new_row$contrast &
                    group_A  == new_row$group_A  &
                    group_B  == new_row$group_B))

  write_tsv(bind_rows(existing, new_row), de_master_file)
  message(sprintf("  [DE stats] %s | %s | %s | up=%d down=%d",
                  section, genome, cname, n_up, n_down))
  invisible(new_row)
}

