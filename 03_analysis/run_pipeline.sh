#!/bin/bash
# =============================================================================
# run_pipeline.sh
# Runs the full downstream analysis, one stage per R script, in dependency
# order. Each stage is a separate Rscript process; console output is teed to
# 03_analysis/logs/<stage>.log. A stage that finishes cleanly writes
# 03_analysis/logs/<stage>.ok and is SKIPPED on the next run.
#
#   bash 03_analysis/run_pipeline.sh              # resume: skip stages already .ok
#   env FORCE=1 bash 03_analysis/run_pipeline.sh  # re-run every stage
#   (the `env` prefix works in bash/zsh/tcsh alike; `FORCE=1 bash ...` is
#    bash-only syntax and errors in tcsh)
#
# Within a stage that does run, GO enrichment results are reused from
# 04_results/.go_cache/ (keyed by input), so a re-run after a downstream
# failure is minutes, not hours. FORCE=1 also bypasses that cache; delete
# 04_results/.go_cache/ to clear it entirely.
#
# Independent stages keep going if an earlier one fails; the summary lists any
# failures and the script exits non-zero if there were any.
#
# Run a single stage directly instead with, e.g.:
#   Rscript 03_analysis/scripts/02_sexual_dimorphism.R
#
# HOMER motif discovery (scripts/09_homer_motifs.R) is NOT run here: it needs a
# local HOMER install. Configure it, set RUN_HOMER_DISCOVERY <- TRUE in that
# script, and run it by hand.
# =============================================================================
set -u

# Work from the repo root so here::here() in the R scripts resolves there.
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
REPO_ROOT="$(pwd)"
LOG_DIR="03_analysis/logs"
mkdir -p "$LOG_DIR"

FORCE="${FORCE:-}"

echo "Repository root: $REPO_ROOT"
echo "Logs:            $LOG_DIR/"
[ -n "$FORCE" ] && echo "FORCE=1          re-running every stage, ignoring GO cache"
echo

# Stage list, in run order. 11 (embryo QC) is independent; 01-03 are
# independent of each other; 04 needs 02; 07/10 rebuild what they need; 08
# needs 02; 99 needs 01-04; 12 needs 08.
STAGES=(
  "00_init"
  "11_developmental_confirmation"
  "01_mzt"
  "02_sexual_dimorphism"
  "03_postmating"
  "04_maleness"
  "05_tissue_landscapes"
  "06_global_atlas"
  "07_tau_specificity"
  "08_upset_sex_bias"
  "10_crosstissue_correlation"
  "99_de_stats_master"
  "12_manuscript_plots"
)

failed=()
skipped=()
for stage in "${STAGES[@]}"; do
  script="03_analysis/scripts/${stage}.R"
  log="${LOG_DIR}/${stage}.log"
  ok="${LOG_DIR}/${stage}.ok"

  if [ -z "$FORCE" ] && [ -f "$ok" ]; then
    printf '=== %-32s skip (%s exists; FORCE=1 to redo)\n' "$stage" "${ok#03_analysis/}"
    skipped+=("$stage")
    continue
  fi

  rm -f "$ok"
  printf '=== %-32s ' "$stage"
  start=$(date +%s)
  if Rscript "$script" > "$log" 2>&1; then
    touch "$ok"
    echo "ok   ($(($(date +%s) - start))s)  -> $log"
  else
    echo "FAIL ($(($(date +%s) - start))s)  -> $log"
    failed+=("$stage")
  fi
done

echo
[ ${#skipped[@]} -gt 0 ] && echo "Skipped (already .ok): ${skipped[*]}"
if [ ${#failed[@]} -eq 0 ]; then
  echo "All stages up to date. Results in 04_results/."
  exit 0
else
  echo "Stages that FAILED (see their logs):"
  printf '  - %s\n' "${failed[@]}"
  exit 1
fi
