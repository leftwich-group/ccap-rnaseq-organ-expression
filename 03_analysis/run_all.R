#!/usr/bin/env Rscript
# =============================================================================
# run_all.R
# R wrapper around the same stages as run_pipeline.sh. Each stage runs in its
# own Rscript process; output is teed to 03_analysis/logs/<stage>.log. A stage
# that finishes cleanly writes 03_analysis/logs/<stage>.ok and is SKIPPED next
# time.
#
#   Rscript 03_analysis/run_all.R              # resume: skip stages already .ok
#   env FORCE=1 Rscript 03_analysis/run_all.R  # re-run every stage
#   (use the `env` prefix; `FORCE=1 Rscript ...` is bash-only and fails in tcsh)
#
# Within a stage that runs, GO enrichment is reused from 04_results/.go_cache/
# (FORCE=1 bypasses it; delete the folder to clear it).
#
# scripts/09_homer_motifs.R (HOMER motif discovery) is intentionally excluded —
# it needs a local HOMER install; run it by hand once configured.
# =============================================================================

# Work from the repo root so here::here() resolves correctly.
this <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
repo_root <- normalizePath(file.path(dirname(this), ".."))
setwd(repo_root)

log_dir <- file.path("03_analysis", "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

force <- nzchar(Sys.getenv("FORCE"))
if (force) cat("FORCE=1  re-running every stage, ignoring GO cache\n\n")

stages <- c(
  "00_init",
  "11_developmental_confirmation",
  "01_mzt",
  "02_sexual_dimorphism",
  "03_postmating",
  "04_maleness",
  "05_tissue_landscapes",
  "06_global_atlas",
  "07_tau_specificity",
  "08_upset_sex_bias",
  "10_crosstissue_correlation",
  "99_de_stats_master",
  "12_manuscript_plots"
)

rscript <- file.path(R.home("bin"), "Rscript")
failed <- character(0)
skipped <- character(0)

for (stage in stages) {
  script <- file.path("03_analysis", "scripts", paste0(stage, ".R"))
  log    <- file.path(log_dir, paste0(stage, ".log"))
  ok     <- file.path(log_dir, paste0(stage, ".ok"))

  if (!force && file.exists(ok)) {
    cat(sprintf("=== %-32s skip (%s exists; FORCE=1 to redo)\n", stage, basename(ok)))
    skipped <- c(skipped, stage)
    next
  }

  unlink(ok)
  cat(sprintf("=== %-32s ", stage))
  t0 <- Sys.time()
  status <- system2(rscript, script, stdout = log, stderr = log)
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (status == 0) {
    file.create(ok)
    cat(sprintf("ok   (%ds)  -> %s\n", dt, log))
  } else {
    cat(sprintf("FAIL (%ds)  -> %s\n", dt, log))
    failed <- c(failed, stage)
  }
}

cat("\n")
if (length(skipped)) cat("Skipped (already .ok):", paste(skipped, collapse = " "), "\n")
if (length(failed) == 0) {
  cat("All stages up to date. Results in 04_results/.\n")
} else {
  cat("Stages that FAILED (see their logs):\n")
  cat(paste0("  - ", failed, "\n"))
  quit(status = 1)
}
