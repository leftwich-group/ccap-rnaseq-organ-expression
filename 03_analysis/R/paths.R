# ===========================================================================
# R/paths.R
# One place that names every 04_results/ sub-directory. Script code writes to
# these zero-padded names (matching the folders already on disk); the runner
# and README refer to RES for the script <-> output-folder mapping.
# ===========================================================================

res_dir <- function(...) here::here("04_results", ...)

RES <- list(
  mzt          = "01_MZT",
  dimorphism   = "02_Sexual_dimorphism",
  postmating   = "03_PostMating",
  maleness     = "04_Maleness",
  landscapes   = "05_Tissue_landscapes",
  atlas        = "06_Global_atlas",
  tau          = "07_Tau_specificity",
  upset        = "08_UpSet_sex_bias",
  homer        = "09_HOMER_motifs",
  crosstissue  = "10_CrossTissue_correlation",
  devconfirm   = "11_developmental_confirmation",
  manuscript   = "12_manuscript_plots"
)
