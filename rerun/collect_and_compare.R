#!/usr/bin/env Rscript
# ==============================================================================
# OPTICO-ACS | rerun check
# 1) copies the regenerated panel files and source tables of Fig. 2-6 into
#    <rerun folder>/_panels_check/<Figure>/, each next to the original version
# 2) compares every regenerated table (csv, xlsx) with the original table at the
#    same relative path and writes _panels_check/table_comparison.csv
#
# usage: OPTICO_DIR=<rerun folder> OPTICO_SRC=<original project folder> Rscript collect_and_compare.R
# ==============================================================================

NEW <- normalizePath(Sys.getenv("OPTICO_DIR"))
OLD <- normalizePath(Sys.getenv("OPTICO_SRC"))
OUT <- file.path(NEW, "_panels_check")
dir.create(OUT, showWarnings = FALSE)

panels <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, quote = "", text = "
figure   | panel | file
Figure_2 | 2e    | Figures/Figure_1/plots/radarchart_no_labels.pdf
Figure_3 | 3d    | Figures/Figure_2/plots/R2_umap.pdf
Figure_3 | 3e    | Figures/Figure_2/plots/Erosion_umap.pdf
Figure_3 | 3e    | Figures/Figure_2/plots/Rupture_umap.pdf
Figure_3 | 3f    | Optico_revisions/tables/cluster_abundance_Optico.xlsx
Figure_3 | 3g    | DE_Genes/Res1/CD4/CD4_markers.csv
Figure_3 | 3g    | DE_Genes/Res1/CD8/CD8_markers.csv
Figure_3 | 3b    | plots/jaccard_heatmap.pdf
Figure_3 | 3c    | tables/md_subset_R1.xlsx
Figure_4 | 4fgh  | Optico_revisions/tables/FACS_stimulation/Stats_summary.csv
Figure_4 | 4fgh  | Optico_revisions/tables/FACS_stimulation/FigA_*_Paired_Stim.pdf
Figure_4 | 4fgh  | Optico_revisions/tables/FACS_stimulation/FigB_FC_Pct_Ruptur_Erosion.pdf
Figure_5 | 5a    | T_reclustering_2025/plots/umap_annotated.pdf
Figure_5 | 5b    | T_reclustering_2025/plots/VlnPlot_*.pdf
Figure_5 | 5c    | T_reclustering_2025/tables/cluster_abundance_Optico_T.xlsx
Figure_5 | 5d    | T_reclustering_2025/tables/or_results.xlsx
Figure_5 | 5e    | T_reclustering_2025/plots/cloneSize_umap_magma_muted_black.pdf
Figure_5 | 5f    | T_reclustering_2025/tables/erosion_table_clonal_occupy_annotation.xlsx
Figure_5 | 5f    | T_reclustering_2025/tables/rupture_table_clonal_oocupy_annotation.xlsx
Figure_5 | 5g    | T_reclustering_2025/tables/expanded_DEG/DEG_expanded_vs_not.xlsx
Figure_5 | 5i    | T_reclustering_2025/plots/expanded_DEG/SCENIC_lollipop_ErosionVsRupture_CD*_labeled.pdf
Figure_6 | 6b    | Optico_revisions/plots/CVD_TCR/plots/circos_shared_TR*_Erosion_Rupture.pdf
Figure_6 | 6c    | T_reclustering_2025/plots/clonalSizeDistributuion.pdf
Figure_6 | 6d    | T_reclustering_2025/tables/clonalDiversity*.xlsx
Figure_6 | 6e    | Optico_revisions/plots/CVD_TCR/plots/OPTICO_CVD_shared_CDR3b_heatmap_landscape.pdf
Figure_6 | 6f    | Optico_revisions/plots/CVD_TCR/tables/jaccard_matrix_full.csv
Figure_6 | 6f    | Optico_revisions/plots/CVD_TCR/plots/jaccard_barplot.pdf
Figure_6 | 6g    | Optico_revisions/plots/CVD_TCR/tables/per_patient_freq_weighted_overlap.csv
Figure_6 | 6g    | Optico_revisions/plots/CVD_TCR/plots/freq_weighted_overlap_boxplot.pdf
")

# ---- 1) collect panels -------------------------------------------------------
rel <- function(p, root) substring(p, nchar(root) + 2)
rows <- list()
for (i in seq_len(nrow(panels))) {
  hits <- Sys.glob(file.path(NEW, panels$file[i]))
  if (!length(hits)) {
    rows[[length(rows) + 1]] <- data.frame(figure = panels$figure[i], panel = panels$panel[i], file = panels$file[i],
                                           rerun = "MISSING", original = NA, stringsAsFactors = FALSE)
    next
  }
  for (h in hits) {
    r <- rel(h, NEW)
    d <- file.path(OUT, panels$figure[i]); dir.create(d, showWarnings = FALSE)
    stem <- tools::file_path_sans_ext(basename(h)); ext <- tools::file_ext(h)
    file.copy(h, file.path(d, paste0("Fig", panels$panel[i], "__", stem, "__RERUN.", ext)), overwrite = TRUE)
    o <- file.path(OLD, r)
    has_o <- file.exists(o)
    if (has_o) file.copy(o, file.path(d, paste0("Fig", panels$panel[i], "__", stem, "__ORIGINAL.", ext)), overwrite = TRUE)
    rows[[length(rows) + 1]] <- data.frame(figure = panels$figure[i], panel = panels$panel[i], file = r,
                                           rerun = "ok", original = if (has_o) "found" else "not found", stringsAsFactors = FALSE)
  }
}
panel_report <- do.call(rbind, rows)
write.csv(panel_report, file.path(OUT, "panel_files.csv"), row.names = FALSE)

# ---- 2) compare tables -------------------------------------------------------
inputs_file <- file.path(NEW, "_rerun_logs", "inputs.txt")
inputs <- if (file.exists(inputs_file)) readLines(inputs_file) else character()

read_tabs <- function(p) {
  ext <- tolower(tools::file_ext(p))
  if (ext == "csv") {
    list(csv = utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE))
  } else {
    sh <- readxl::excel_sheets(p)
    stats::setNames(lapply(sh, function(s) as.data.frame(suppressMessages(readxl::read_excel(p, sheet = s)), check.names = FALSE)), sh)
  }
}
sort_df <- function(d) if (nrow(d) && ncol(d)) d[do.call(order, unname(lapply(d, function(x) as.character(x)))), , drop = FALSE] else d
cmp_df <- function(a, b) {
  if (!identical(dim(a), dim(b))) return(list(verdict = "dimensions differ", detail = sprintf("%dx%d vs %dx%d", nrow(a), ncol(a), nrow(b), ncol(b))))
  if (!identical(names(a), names(b))) return(list(verdict = "column names differ", detail = ""))
  one <- function(a, b) {
    maxd <- 0; ndiff <- 0
    for (j in seq_along(a)) {
      x <- a[[j]]; y <- b[[j]]
      if (is.numeric(x) && is.numeric(y)) {
        nas <- is.na(x) != is.na(y)
        d <- abs(x - y); d[is.na(d)] <- 0
        tol <- 1e-6 * pmax(1, abs(y), na.rm = TRUE); tol[is.na(tol)] <- 1e-6
        ndiff <- ndiff + sum(nas | d > tol); maxd <- max(maxd, d)
      } else {
        x <- as.character(x); y <- as.character(y)
        ndiff <- ndiff + sum(is.na(x) != is.na(y) | (!is.na(x) & !is.na(y) & x != y))
      }
    }
    c(ndiff = ndiff, maxd = maxd)
  }
  r <- one(a, b)
  if (r["ndiff"] == 0) return(list(verdict = "identical", detail = sprintf("max abs numeric diff %.3g", r["maxd"])))
  r2 <- one(sort_df(a), sort_df(b))
  if (r2["ndiff"] == 0) return(list(verdict = "identical after sorting rows", detail = ""))
  list(verdict = "values differ", detail = sprintf("%d of %d cells differ (as is), %d after sorting; max abs numeric diff %.3g",
                                                   r["ndiff"], nrow(a) * ncol(a), r2["ndiff"], r["maxd"]))
}

new_files <- list.files(NEW, pattern = "\\.(csv|xlsx)$", recursive = TRUE, ignore.case = TRUE)
new_files <- new_files[!grepl("^_panels_check/|^_rerun_logs/|(^|/)~\\$|(^|/)\\._", new_files)]
new_files <- setdiff(new_files, inputs)
res <- list()
for (f in new_files) {
  o <- file.path(OLD, f); n <- file.path(NEW, f)
  if (!file.exists(o)) { res[[length(res) + 1]] <- data.frame(file = f, sheet = "", verdict = "no original at this path", detail = ""); next }
  if (max(file.size(o), file.size(n)) > 60e6) { res[[length(res) + 1]] <- data.frame(file = f, sheet = "", verdict = "not compared (file > 60 MB)", detail = ""); next }
  out <- tryCatch({
    A <- read_tabs(n); B <- read_tabs(o)
    do.call(rbind, lapply(union(names(A), names(B)), function(s) {
      if (is.null(A[[s]]) || is.null(B[[s]])) return(data.frame(file = f, sheet = s, verdict = "sheet only in one file", detail = ""))
      r <- cmp_df(A[[s]], B[[s]])
      data.frame(file = f, sheet = s, verdict = r$verdict, detail = r$detail)
    }))
  }, error = function(e) data.frame(file = f, sheet = "", verdict = "could not be read", detail = conditionMessage(e)))
  res[[length(res) + 1]] <- out
}
cmp <- if (length(res)) do.call(rbind, res) else data.frame(file = character(), sheet = character(), verdict = character(), detail = character())
cmp$is_panel_table <- cmp$file %in% panel_report$file
cmp <- cmp[order(!cmp$is_panel_table, cmp$verdict, cmp$file), ]
write.csv(cmp, file.path(OUT, "table_comparison.csv"), row.names = FALSE)

cat("\n===== panel files =====\n")
print(panel_report, row.names = FALSE, right = FALSE)
cat("\n===== table comparison (rerun vs original) =====\n")
print(as.data.frame(table(verdict = cmp$verdict)), row.names = FALSE)
cat("\nTables behind main-figure panels:\n")
print(cmp[cmp$is_panel_table, c("file", "sheet", "verdict", "detail")], row.names = FALSE, right = FALSE)
cat("\nWritten: ", file.path(OUT, "panel_files.csv"), "\n         ", file.path(OUT, "table_comparison.csv"), "\n", sep = "")
