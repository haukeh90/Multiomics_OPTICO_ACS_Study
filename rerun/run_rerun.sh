#!/bin/bash
# ==============================================================================
# OPTICO-ACS | rerun of the code that exports figure panels and tables
#
# Works on a separate folder: the inputs are copied from the original project
# folder (SRC) into RERUN, all scripts run with OPTICO_DIR=RERUN, and nothing is
# written to SRC. Integration, clustering, annotation and pySCENIC are
# not recomputed; the steps start from data/optico_final.rds and
# T_reclustering_2025/data/T_reclustering_2025_annotated.rds.
#
#   bash run_rerun.sh                      # all default steps
#   STEPS="fig3de fig3f" bash run_rerun.sh # selected steps only
#   RUN_NMF=1 bash run_rerun.sh            # additionally rerun GeneNMF (Fig. 3b, 3c; long, not seeded)
#   RERUN_DRY=1 bash run_rerun.sh          # only list the sections each step would run
#
# Results:  $RERUN/_rerun_logs/summary.tsv, <step>.log, <step>_sections.tsv
#           $RERUN/_panels_check/<Figure>/  rerun and original panel side by side
#           $RERUN/_panels_check/table_comparison.csv
# ==============================================================================
set -u

SRC="${SRC:-/path/to/OPTICO_project}"
RERUN="${RERUN:-/path/to/OPTICO_rerun_check}"
HERE="$(cd "$(dirname "$0")" && pwd)"
DEFAULT_STEPS="fig2e fig3de fig3f fig3g fig4fgh fig5_6_T fig5g fig5i fig6"
STEPS="${STEPS:-$DEFAULT_STEPS}"
if [ "${RUN_NMF:-0}" = "1" ]; then STEPS="$STEPS nmf"; fi

RSCRIPT="$(command -v Rscript || true)"
for c in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript /Library/Frameworks/R.framework/Resources/bin/Rscript; do
  if [ -z "$RSCRIPT" ] && [ -x "$c" ]; then RSCRIPT="$c"; fi
done
[ -n "$RSCRIPT" ] || { echo "Rscript not found. Add R to PATH or set it in this script."; exit 1; }
[ -d "$SRC" ]     || { echo "Source folder not found: $SRC"; exit 1; }
case "$RERUN/" in "$SRC"/*) echo "RERUN must not be the source folder or lie inside it."; exit 1;; esac

LOG="$RERUN/_rerun_logs"
mkdir -p "$LOG" || exit 1
echo "Source : $SRC"
echo "Rerun  : $RERUN"
echo "Rscript: $RSCRIPT"
echo "Steps  : $STEPS"

# ---- 1) folder skeleton (directories only) -----------------------------------
echo; echo "== creating folder skeleton"
( cd "$SRC" && find . -type d \
    -not -path './github_repo*' -not -path './GEO_Optico*' -not -path './raw_data*' \
    -not -path './MOFA_export*' -not -path './Dropbox_Fogure_Backup*' \
    -not -path './data/TCR_*' -not -path '*/.*' ) | while IFS= read -r d; do
  mkdir -p "$RERUN/$d"
done

# ---- 2) inputs (copied once; existing copies are kept) -----------------------
echo "== copying inputs"
: > "$LOG/inputs.txt"
copy_in() {  # relative file or folder
  rel="$1"
  if [ ! -e "$SRC/$rel" ]; then echo "   MISSING INPUT: $rel"; return; fi
  mkdir -p "$(dirname "$RERUN/$rel")"
  if [ -d "$SRC/$rel" ]; then
    [ -d "$RERUN/$rel" ] || cp -R "$SRC/$rel" "$RERUN/$rel"
    ( cd "$SRC" && find "$rel" -type f -not -name '._*' ) >> "$LOG/inputs.txt"
  else
    [ -f "$RERUN/$rel" ] || { echo "   $rel"; cp "$SRC/$rel" "$RERUN/$rel"; }
    echo "$rel" >> "$LOG/inputs.txt"
  fi
}
copy_in "data/optico_final.rds"
copy_in "data/TCR_raw"
copy_in "T_reclustering_2025/data/T_reclustering_2025_annotated.rds"
copy_in "T_reclustering_2025/data/SCENIC/SCENIC_inspection_tables.xlsx"
copy_in "Figures/Figure_1/tables/TH_subtypes_for_spider_CR.xlsx"
copy_in "Figures/Figure_1/tables/T_Activation_for_spider_CR.xlsx"
copy_in "Optico_revisions/tables/FACS_stimulation/FACS_tidy.xlsx"
copy_in "Optico_revisions/data/CVD_TCR_Database.xlsx"

# ---- 3) steps ----------------------------------------------------------------
printf 'step\tstatus\tminutes\tlog\n' > "$LOG/summary.tsv"
for s in $STEPS; do
  echo; echo "== step $s  ($(date '+%H:%M:%S'))"
  t0=$(date +%s)
  OPTICO_DIR="$RERUN" OPTICO_SRC="$SRC" "$RSCRIPT" "$HERE/run_step.R" "$s" > "$LOG/$s.log" 2>&1
  rc=$?
  t1=$(date +%s)
  if [ $rc -eq 0 ]; then st="ok"; else st="FAILED"; fi
  mins=$(( (t1 - t0 + 30) / 60 ))
  printf '%s\t%s\t%s\t%s\n' "$s" "$st" "$mins" "$LOG/$s.log" >> "$LOG/summary.tsv"
  echo "   $st after ${mins} min"
  if [ -f "$LOG/${s}_sections.tsv" ]; then
    awk -F'\t' 'NR>1 && $3!="skipped" {printf "     %-10s %s%s\n", $3, $2, ($5!="" ? "  | " $5 : "")}' "$LOG/${s}_sections.tsv" | cut -c1-200
  else
    tail -n 5 "$LOG/$s.log" | sed 's/^/     /'
  fi
done

# ---- 4) collect panels, compare tables ---------------------------------------
if [ -z "${RERUN_DRY:-}" ]; then
  echo; echo "== collecting panels and comparing tables"
  OPTICO_DIR="$RERUN" OPTICO_SRC="$SRC" "$RSCRIPT" "$HERE/collect_and_compare.R" > "$LOG/collect_and_compare.log" 2>&1 \
    && sed -n '/===== table comparison/,$p' "$LOG/collect_and_compare.log" \
    || { echo "   collect_and_compare.R failed, see $LOG/collect_and_compare.log"; tail -n 5 "$LOG/collect_and_compare.log"; }
fi

echo; echo "== summary"
cat "$LOG/summary.tsv"
echo; echo "Logs:   $LOG"
echo "Panels: $RERUN/_panels_check"
