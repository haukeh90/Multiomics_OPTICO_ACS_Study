#!/usr/bin/env Rscript
# ==============================================================================
# OPTICO-ACS | rerun of the panel- and table-generating code
#
# Runs one step of the rerun (see run_rerun.sh). Each step executes the scripts
# in code/ section by section, starting from the stored objects
# (data/optico_final.rds, T_reclustering_2025_annotated.rds). Integration,
# clustering, NMF and SCENIC are not recomputed unless requested.
#
# usage: OPTICO_DIR=/path/to/rerun_folder Rscript run_step.R <step>
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
step <- if (length(args)) args[1] else stop("usage: run_step.R <step>")

cmd  <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
REPO <- normalizePath(file.path(dirname(self), ".."))
OPTICO_DIR <- Sys.getenv("OPTICO_DIR")
if (!nzchar(OPTICO_DIR) || !dir.exists(OPTICO_DIR)) stop("OPTICO_DIR is not set or does not exist")
OPTICO_DIR <- normalizePath(OPTICO_DIR)
Sys.setenv(OPTICO_DIR = OPTICO_DIR)
LOG_DIR <- file.path(OPTICO_DIR, "_rerun_logs")
dir.create(LOG_DIR, showWarnings = FALSE)

DRY <- nzchar(Sys.getenv("RERUN_DRY"))   # RERUN_DRY=1: only list the sections that would run
options(warn = 1, width = 160)
# plots that are only printed (not saved) go to a null device instead of Rplots.pdf
options(device = function(...) grDevices::pdf(file = NULL))

DASH_HEADER <- "^\\s*#.*-{4,}\\s*([0-9.]+)?\\s*$"   # "#Section title -------"
status <- data.frame(script = character(), section = character(), status = character(),
                     seconds = numeric(), message = character(), stringsAsFactors = FALSE)

# Evaluate a file expression by expression; errors are logged and evaluation continues.
eval_each <- function(tmp) {
  exprs <- tryCatch(parse(tmp, keep.source = FALSE, encoding = "UTF-8"), error = function(e) e)
  if (inherits(exprs, "error")) return(conditionMessage(exprs))
  fails <- character()
  for (e in exprs) {
    d <- deparse(e, width.cutoff = 120L)
    cat("\n> ", d[1], if (length(d) > 1) " ...", "\n", sep = "")
    r <- tryCatch({
      v <- withVisible(eval(e, envir = globalenv()))
      if (v$visible) print(v$value)
      NULL
    }, error = function(err) conditionMessage(err))
    if (!is.null(r)) {
      cat("[rerun] ERROR: ", r, "\n", sep = "")
      fails <- c(fails, paste0(substr(d[1], 1, 60), ": ", r))
    }
  }
  if (!length(fails)) "ok" else paste0(length(fails), " of ", length(exprs), " expressions failed; first: ", fails[1])
}

# Run a script section by section. A section starts at a header line and ends
# before the next one; lines before the first header are the preamble.
#   header   regex that identifies header lines
#   run      regex vector; only sections whose header matches one of them are run
#            (NULL = all). The preamble always runs unless preamble = FALSE.
#   skip     regex vector; sections whose header matches are not run
#   drop     regex; single lines matching it are removed from the sections that run
#   continue TRUE: evaluate a section expression by expression and carry on after
#            an error (for plotting code where the expressions are independent)
run_script <- function(file, header = DASH_HEADER, run = NULL, skip = NULL,
                       drop = NULL, preamble = TRUE, continue = FALSE, log_skips = TRUE,
                       uncomment = NULL) {
  path  <- file.path(REPO, file)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  # uncomment: regex; matching commented lines are activated for this run only
  if (!is.null(uncomment)) { u <- grepl(uncomment, lines); lines[u] <- sub("^\\s*#\\s?", "", lines[u]) }
  # the repo scripts set OPTICO_DIR themselves; keep the value of this run
  lines[grepl("^OPTICO_DIR <- Sys.getenv", lines)] <- ""
  h <- grep(header, lines, perl = TRUE)
  starts <- c(1, h); ends <- c(h - 1, length(lines))
  titles <- c("(preamble)", trimws(gsub("[^[:alnum:] _().,:/&+-]|-{3,}.*$", "", iconv(lines[h], "UTF-8", "ASCII", sub = ""))))
  keep <- ends >= starts
  starts <- starts[keep]; ends <- ends[keep]; titles <- titles[keep]
  sec <- lapply(seq_along(starts), function(i) lines[starts[i]:ends[i]])
  hdr <- vapply(seq_along(starts), function(i) if (titles[i] == "(preamble)") "" else lines[starts[i]], "")
  # if a header sits inside an expression the split is not valid: run as one block
  ok <- vapply(sec, function(s) !inherits(try(parse(text = s, keep.source = FALSE), silent = TRUE), "try-error"), TRUE)
  if (!all(ok)) {
    if (!is.null(run) || length(skip)) stop("sections of ", file, " cannot be split cleanly; nothing was run")
    message("[rerun] ", file, ": sections cannot be split cleanly, running the file as one block")
    sec <- list(lines); titles <- "(whole file)"; hdr <- ""; starts <- 1; ends <- length(lines)
  }
  for (i in seq_along(sec)) {
    is_pre <- titles[i] %in% c("(preamble)", "(whole file)")
    sel <- if (is_pre) preamble || titles[i] == "(whole file)" else
      (is.null(run) || any(vapply(run, grepl, TRUE, x = hdr[i], perl = TRUE))) &&
      !(length(skip) && any(vapply(skip, grepl, TRUE, x = hdr[i], perl = TRUE)))
    if (!sel) {
      if (log_skips) status[nrow(status) + 1, ] <<- list(file, titles[i], "skipped", 0, "")
      next
    }
    if (DRY) {
      status[nrow(status) + 1, ] <<- list(file, titles[i], "would run", 0, paste(starts[i], ends[i], sep = "-"))
      next
    }
    cat("\n\n######## [rerun] ", file, " :: ", titles[i], " ########\n", sep = "")
    txt <- sec[[i]]
    if (!is.null(drop)) txt[grepl(drop, txt)] <- ""
    tmp <- tempfile(fileext = ".R"); writeLines(txt, tmp, useBytes = TRUE)
    t0 <- Sys.time()
    res <- if (continue) eval_each(tmp) else tryCatch({
      source(tmp, echo = TRUE, print.eval = TRUE, max.deparse.length = 400,
             local = globalenv(), keep.source = FALSE, encoding = "UTF-8", spaced = FALSE)
      "ok"
    }, error = function(e) conditionMessage(e))
    while (!is.null(grDevices::dev.list())) grDevices::dev.off()
    secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    if (identical(res, "ok")) {
      status[nrow(status) + 1, ] <<- list(file, titles[i], "ok", secs, "")
    } else {
      cat("\n[rerun] ERROR in section '", titles[i], "': ", res, "\n", sep = "")
      status[nrow(status) + 1, ] <<- list(file, titles[i], "FAILED", secs, gsub("[\r\n\t]+", " ", res))
    }
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
steps <- list(

  # Fig. 2e: radar plot of culprit ratios
  fig2e = function() run_script("code/Figure_2/Fig2e_radar_culprit_ratios.R", continue = TRUE),

  # Fig. 3d, 3e: WNN UMAP, split by condition
  fig3de = function() run_script("code/Figure_3/Fig3de_WNN_UMAP.R", continue = TRUE),

  # Fig. 3f: cluster proportions per patient
  fig3f = function() run_script("code/Figure_3/Fig3f_cluster_abundance.R"),

  # Fig. 3g (+ gene lists for 3h): pseudobulk edgeR-LRT via Libra
  # The stored object has a Seurat v3 RNA assay; the commented conversion line in the script
  # is activated here so that the step runs with a current Libra installation.
  fig3g = function() run_script("code/Figure_3/Fig3g_pseudobulk_DE.R",
                                uncomment = "^# optico_final\\[\\[\"RNA\"\\]\\] <- as\\("),

  # Fig. 4f-h: stimulation assay plots and statistics
  fig4fgh = function() run_script("code/Figure_4/Fig4fgh_stimulation_assay.R", continue = TRUE),

  # Fig. 5a-f, 6c, 6d: plots and tables from the annotated T cell object.
  # Subclustering, annotation and the combineExpression step are not rerun.
  fig5_6_T = function() {
    f <- "code/Figure_5/Fig5a-f_Fig6cd_T_subclustering_TCR.R"
    run_script(f, run = character(0), log_skips = FALSE)   # libraries only
    if (!DRY) {
      suppressPackageStartupMessages({ library(tibble); library(tidyr); library(scRepertoire) })
      cat("[rerun] loading T_reclustering_2025_annotated.rds\n")
      assign("T_reclustering_2025",
             readRDS(file.path(OPTICO_DIR, "T_reclustering_2025/data/T_reclustering_2025_annotated.rds")),
             envir = globalenv())
    }
    run_script(f, preamble = FALSE, continue = TRUE, drop = "^T_reclustering_2025 <- readRDS\\(", run = c(
      "scTCR-seq: load contigs", "scTCR-seq: add condition labels",
      "Fig\\. 6c", "Fig\\. 6d", "Fig\\. 5a, 5e", "Fig\\. 5d", "Fig\\. 5b", "Fig\\. 5f", "Fig\\. 5c"
    ))
  },

  # Fig. 5g (+ gene lists for 5h): DE of expanded vs non-expanded T cells.
  # Section 6 (matrix export for pySCENIC) is not rerun.
  fig5g = function() run_script("code/Figure_5/Fig5g_oligoclonal_DEG.R",
                                header = "^# \\S+ [0-9]+b?\\) ", skip = "^# \\S+ 6\\) "),

  # Fig. 5i: lollipop plot from the stored pySCENIC result tables
  fig5i = function() run_script("code/Figure_5/Fig5i_2_regulon_lollipop.R"),

  # Fig. 6b, 6e, 6f, 6g: CVD-TCR mapping and panels (same session)
  fig6 = function() {
    hd <- "^# +[0-9]+[A-Z]?\\. .*-{4,}\\s*$"
    run_script("code/Figure_6/Fig6_1_CVD_TCR_mapping.R", header = hd)
    run_script("code/Figure_6/Fig6_2_CVD_TCR_panels.R", header = hd)
  },

  # optional, long and not seeded: Fig. 3b, 3c (GeneNMF with TCR genes excluded)
  nmf = function() run_script("code/Figure_3/Fig3bc_NMF_metaprograms.R")
)

if (!step %in% names(steps)) stop("unknown step '", step, "'. Available: ", paste(names(steps), collapse = ", "))
cat("[rerun] step ", step, " | OPTICO_DIR = ", OPTICO_DIR, " | ", format(Sys.time()), "\n", sep = "")
err <- tryCatch({ steps[[step]](); NULL }, error = function(e) conditionMessage(e))
if (!is.null(err)) {
  cat("\n[rerun] step aborted: ", err, "\n", sep = "")
  status[nrow(status) + 1, ] <- list("(step)", step, "FAILED", 0, gsub("[\r\n\t]+", " ", err))
}
write.table(status, file.path(LOG_DIR, paste0(step, "_sections.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n\n[rerun] section summary for step ", step, "\n", sep = "")
print(status[, c("section", "status", "seconds", "message")], right = FALSE, row.names = FALSE)
cat("\n[rerun] sessionInfo\n"); print(sessionInfo())
quit(status = if (any(status$status == "FAILED")) 1 else 0)
