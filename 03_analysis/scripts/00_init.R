#!/usr/bin/env Rscript
# Initialise a fresh analysis run
# Output: 04_results/DE_stats_master.tsv (truncated to an empty typed table)
# Run from the repository root:  Rscript 03_analysis/scripts/00_init.R
# ---------------------------------------------------------------------------
if (!file.exists("03_analysis/R/init.R")) stop(
  "Run from the repository root (the folder containing 01_data/, ",
  "03_analysis/, 04_results/).", call. = FALSE)
source("03_analysis/R/init.R")
# ---------------------------------------------------------------------------

# Truncate the DE-statistics master table so a full pipeline run starts clean.
# Skip this (run the section scripts directly) to append to an existing table.
write_tsv(de_master_empty, de_master_file)
message("DE stats master initialised: ", de_master_file)
