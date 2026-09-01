# ===========================================================================
# R/init.R
# Common bootstrap: locate the repo root, then load setup + paths + data.
# Every script in scripts/ sources this file first.
# ===========================================================================
suppressPackageStartupMessages(library(here))

if (!dir.exists(here::here("01_data")) ||
    !dir.exists(here::here("02_preprocessing"))) {
  stop("here() resolved to '", here::here(), "', which is not the repository ",
       "root.\n  Run these scripts with the repository root (the folder ",
       "containing 01_data/, 03_analysis/, 04_results/) as the working ",
       "directory, e.g.\n    Rscript 03_analysis/scripts/01_mzt.R", call. = FALSE)
}

source(here::here("03_analysis/R/setup.R"))
source(here::here("03_analysis/R/paths.R"))
source(here::here("03_analysis/R/load_data.R"))
