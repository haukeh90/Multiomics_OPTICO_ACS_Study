# ==============================================================================
# OPTICO-ACS | Differential expression of oligoclonal (expanded) vs non-expanded T cells
#
# Panels : Fig. 5g (ranked DEGs), gene lists for Fig. 5h (Metascape), input matrices for Fig. 5i
# Input  : T_reclustering_2025/data/T_reclustering_2025_annotated.rds
# Output : T_reclustering_2025/tables/expanded_DEG/DEG_expanded_vs_not.xlsx (one sheet per cluster and condition)
#          T_reclustering_2025/data/SCENIC/*_RNA_*counts.csv, *_metadata.csv (input for pySCENIC)
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

# ============================================================
# DEG: Hyperexpanded vs. Not-expanded T cells
# per cluster × disease (Erosion / Rupture)
# → enrichGO (BP) on upregulated genes
# → Upset plot: shared DEGs between Erosion and Rupture
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)   # for UpSet — cleaner than UpSetR for publication
library(openxlsx)

set.seed(123)

# ── Paths ─────────────────────────────────────────────────────
base_dir  <- file.path(OPTICO_DIR, "T_reclustering_2025")
plot_dir  <- file.path(base_dir, "plots/expanded_DEG")
table_dir <- file.path(base_dir, "tables/expanded_DEG")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load object ───────────────────────────────────────────────
T_obj <- readRDS(file.path(base_dir, "data/T_reclustering_2025_annotated.rds"))
DefaultAssay(T_obj) <- "RNA"


# ── 1) DEFINE EXPANSION GROUPS ────────────────────────────────
# Hyperexpanded = "Large" (0.01 < X <= 0.1) — top tier in this dataset
# To also include "Medium", change the first case_when condition:
#    cloneSize %in% c("Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)")
#
# Not expanded  = "Rare" (< 0.01%) or NA (no TCR detected)
T_obj$expansion_group <- dplyr::case_when(
  T_obj$cloneSize %in% c(
    "Large (0.01 < X <= 0.1)",
    "Medium (0.001 < X <= 0.01)",
    "Small (1e-04 < X <= 0.001)"
  )                                ~ "Expanded",
  T_obj$cloneSize == "Rare (0 < X <= 1e-04)" ~ "Not_expanded",
  is.na(T_obj$cloneSize)           ~ "Not_expanded",
  TRUE                             ~ NA_character_
)

cat("Expansion group × cluster × disease:\n")
print(table(T_obj$expansion_group, T_obj$annotation, T_obj$disease, useNA = "ifany"))


# ── 1b) LINEAGE GROUPING for meta-cluster analysis ────────────
# Pools all CD8 or CD4 clusters to gain power where per-cluster n is low
cd8_clusters <- c("CD8_TEMRA", "CD69_CD8_TEMRA", "CD8_TCM", "CD8_TEM", "Naive_CD8")
cd4_clusters <- c("Naive_CD4", "CD4_TEM", "CD4_TCM", "Treg", "DPT")

T_obj$lineage <- dplyr::case_when(
  T_obj$annotation %in% cd8_clusters ~ "CD8_all",
  T_obj$annotation %in% cd4_clusters ~ "CD4_all",
  TRUE                                ~ NA_character_
)

cat("\nLineage × expansion × disease:\n")
print(table(T_obj$lineage, T_obj$expansion_group, T_obj$disease, useNA = "ifany"))


# ── 2) GENE FILTER: remove HLA and TCR genes ─────────────────
# Applied before DEG AND before enrichGO universe
all_genes   <- rownames(T_obj)
hla_genes   <- grep("^HLA-",       all_genes, value = TRUE)
tcr_genes   <- grep("^TR[ABGD][VDJC]", all_genes, value = TRUE)
genes_excl  <- unique(c(hla_genes, tcr_genes))
genes_keep  <- setdiff(all_genes, genes_excl)

cat("\nGenes excluded (HLA + TCR):", length(genes_excl),
    "\nGenes retained for analysis:", length(genes_keep), "\n")


# ── 3) DEG FUNCTION ───────────────────────────────────────────
# Within each cluster × disease: Hyperexpanded vs. Not_expanded
# Downsamples to equal n before FindMarkers

run_expanded_DEG <- function(obj, cluster, disease, min_n = 10,
                             group_col = "annotation") {

  # Subset to cluster × disease × defined expansion group only
  cells <- colnames(obj)[
    obj[[group_col]] == cluster &
    obj$disease      == disease &
    !is.na(obj$expansion_group)
  ]
  if (length(cells) == 0) return(NULL)
  sub <- obj[, cells]

  n_hyp <- sum(sub$expansion_group == "Expanded")
  n_not <- sum(sub$expansion_group == "Not_expanded")

  if (n_hyp < min_n | n_not < min_n) {
    message("Skipping ", cluster, " / ", disease,
            " — Expanded n=", n_hyp, ", Not_expanded n=", n_not)
    return(NULL)
  }

  # Equal downsampling to min(n_hyp, n_not)
  n_ds       <- min(n_hyp, n_not)
  cells_hyp  <- sample(colnames(sub)[sub$expansion_group == "Expanded"], n_ds)
  cells_not  <- sample(colnames(sub)[sub$expansion_group == "Not_expanded"],  n_ds)
  sub_ds     <- sub[, c(cells_hyp, cells_not)]

  cat("DEG:", cluster, "/", disease,
      "| n per group =", n_ds, "\n")

  Idents(sub_ds) <- "expansion_group"

  deg <- FindMarkers(
    sub_ds,
    ident.1         = "Expanded",
    ident.2         = "Not_expanded",
    test.use        = "MAST",
    features        = genes_keep,        # HLA/TCR already removed
    min.pct         = 0.10,
    logfc.threshold = 0.25,
    verbose         = FALSE
  )

  deg$gene    <- rownames(deg)
  deg$cluster <- cluster
  deg$disease <- disease
  deg$n_per_group <- n_ds
  deg
}


# ── 4) RUN ACROSS ALL CLUSTER × DISEASE ──────────────────────
clusters  <- levels(T_obj$annotation)
lineages  <- c("CD8_all", "CD4_all")
diseases  <- c("Erosion", "Rupture")

all_DEG <- list()

# Per-cluster analysis
for (dis in diseases) {
  for (cl in clusters) {
    key <- paste(cl, dis, sep = "__")
    res <- run_expanded_DEG(T_obj, cl, dis, group_col = "annotation")
    if (!is.null(res)) all_DEG[[key]] <- res
  }
}

# Meta-cluster analysis (CD8_all / CD4_all)
for (dis in diseases) {
  for (lin in lineages) {
    key <- paste(lin, dis, sep = "__")
    res <- run_expanded_DEG(T_obj, lin, dis, group_col = "lineage")
    if (!is.null(res)) all_DEG[[key]] <- res
  }
}

# Combine, apply global BH correction, filter significant upregulated DEGs
all_DEG_df       <- dplyr::bind_rows(all_DEG)
all_DEG_df$p_adj <- p.adjust(all_DEG_df$p_val, method = "BH")

sig_up <- all_DEG_df %>%
  dplyr::filter(avg_log2FC > 0, p_adj < 0.05)

cat("\nSignificant upregulated DEGs — per cluster:\n")
print(table(sig_up$cluster[!sig_up$cluster %in% lineages],
            sig_up$disease[!sig_up$cluster %in% lineages]))

cat("\nSignificant upregulated DEGs — meta-cluster (lineage):\n")
print(table(sig_up$cluster[sig_up$cluster %in% lineages],
            sig_up$disease[sig_up$cluster %in% lineages]))


# ── 5) EXPORT DEG TABLES ──────────────────────────────────────
wb_deg <- createWorkbook()
for (key in names(all_DEG)) {
  sn <- substr(gsub("__", "_", key), 1, 31)
  addWorksheet(wb_deg, sn)
  writeData(wb_deg, sn, all_DEG[[key]])
}
saveWorkbook(wb_deg,
             file.path(table_dir, "DEG_expanded_vs_not.xlsx"),
             overwrite = TRUE)


# ── 6) EXPORT FOR pySCENIC: CD4 and CD8 expanded vs not-expanded ──────────────
# Exports cells × genes CSV matrices (raw counts + normalised) per lineage,
# matching the established SCENIC workflow used previously (ApoB project).
#
# Gene filtering steps (applied in order):
#   1. Remove HLA / TCR genes (already defined as genes_keep above)
#   2. Remove genes expressed in < 1% of cells (>99% zeros)
#   3. Remove genes with very low total counts (colSums threshold)
#   4. Intersect with cisTarget database genes (feather CSV)
#
# Both Expanded and Not_expanded cells are exported TOGETHER per lineage so
# that pySCENIC infers one shared regulon set. The expansion_group label in
# the metadata CSV is used downstream in Python for AUCell group comparisons.
#
# Purpose: TFs regulating oligoclonal T cells

scenic_dir <- file.path(base_dir, "data/SCENIC")
dir.create(scenic_dir, recursive = TRUE, showWarnings = FALSE)

# ── Path to cisTarget feather database (converted to CSV with pandas) ─────────
# Update this path to the hg38 500bp or 10kb feather CSV for your OPTICO run
feather_db_path <- file.path(OPTICO_DIR, "external_data/SCENIC/feather_500_hg.csv")

# Load database gene list (column names only — avoids loading the full matrix)
cat("Loading cisTarget database gene list from:\n ", feather_db_path, "\n")
db           <- read.csv(feather_db_path, nrows = 1, check.names = FALSE)
db_genes     <- colnames(db)
rm(db); gc()
cat("  Genes in cisTarget DB:", length(db_genes), "\n")

# Metadata columns to embed in output CSV
meta_cols <- c("expansion_group", "disease", "annotation", "lineage")
patient_col <- intersect(
  c("patient", "sample", "orig.ident", "donor"),
  colnames(T_obj@meta.data)
)[1]
if (!is.na(patient_col)) {
  meta_cols <- c(meta_cols, patient_col)
  cat("Patient/sample column detected:", patient_col, "\n")
}

export_scenic_csv <- function(lineage_label) {

  cat(sprintf("\n── SCENIC export: %s ──────────────────────────────────\n",
              lineage_label))

  # ── Subset: defined lineage + defined expansion group only ───────────────────
  cells_use <- colnames(T_obj)[
    !is.na(T_obj$lineage) &
    T_obj$lineage == lineage_label &
    !is.na(T_obj$expansion_group)
  ]

  if (length(cells_use) == 0) {
    message("No cells found for lineage: ", lineage_label)
    return(invisible(NULL))
  }

  sub <- T_obj[, cells_use]
  cat(sprintf("Cells: %d total | Expanded: %d | Not_expanded: %d\n",
              ncol(sub),
              sum(sub$expansion_group == "Expanded"),
              sum(sub$expansion_group == "Not_expanded")))
  print(table(sub$expansion_group, sub$disease))

  # ── Seurat v5: join layers if split by sample ─────────────────────────────────
  if (inherits(sub[["RNA"]], "Assay5")) {
    sub <- JoinLayers(sub, assay = "RNA")
  }

  # ── Extract raw counts and normalised data (genes × cells) ───────────────────
  raw  <- GetAssayData(sub, assay = "RNA", layer = "counts")
  norm <- GetAssayData(sub, assay = "RNA", layer = "data")

  # ── Gene filter 1: Remove HLA / TCR genes (genes_keep defined in Section 2) ──
  raw  <- raw[intersect(genes_keep, rownames(raw)), ]
  norm <- norm[rownames(raw), ]
  cat(sprintf("After HLA/TCR removal:      %d genes\n", nrow(raw)))

  # ── Transpose to cells × genes (pySCENIC convention) ─────────────────────────
  counts      <- t(as.matrix(raw))
  norm_counts <- t(as.matrix(norm))

  # ── Gene filter 2: expressed in ≥ 1% of cells (< 99% zeros) ─────────────────
  keep_expr <- apply(counts, MARGIN = 2,
                     function(x) (length(which(x == 0)) / length(x)) < 0.99)
  counts      <- counts[, keep_expr, drop = FALSE]
  norm_counts <- norm_counts[, keep_expr, drop = FALSE]
  cat(sprintf("After sparsity filter:      %d genes\n", ncol(counts)))

  # ── Gene filter 3: minimum total count threshold ──────────────────────────────
  keep_sum <- colSums(counts) > 3 * 0.01 * nrow(counts)
  counts      <- counts[, keep_sum, drop = FALSE]
  norm_counts <- norm_counts[, keep_sum, drop = FALSE]
  cat(sprintf("After total-count filter:   %d genes\n", ncol(counts)))

  # ── Gene filter 4: intersect with cisTarget database ─────────────────────────
  shared_genes <- intersect(colnames(counts), db_genes)
  counts      <- counts[, shared_genes, drop = FALSE]
  norm_counts <- norm_counts[, shared_genes, drop = FALSE]
  cat(sprintf("After DB intersection:      %d genes\n", ncol(counts)))
  cat(sprintf("Final matrix: %d cells × %d genes\n", nrow(counts), ncol(counts)))

  # ── Write CSVs ────────────────────────────────────────────────────────────────
  counts_path <- file.path(scenic_dir,
                           paste0(lineage_label, "_RNA_counts.csv"))
  norm_path   <- file.path(scenic_dir,
                           paste0(lineage_label, "_RNA_norm_counts.csv"))

  write.csv(counts,      file = counts_path, row.names = TRUE, quote = FALSE)
  write.csv(norm_counts, file = norm_path,   row.names = TRUE, quote = FALSE)
  cat(sprintf("Raw counts written:         %s\n", counts_path))
  cat(sprintf("Normalised counts written:  %s\n", norm_path))

  # ── Write metadata CSV ────────────────────────────────────────────────────────
  meta_avail <- intersect(meta_cols, colnames(sub@meta.data))
  cell_meta  <- sub@meta.data[rownames(counts), meta_avail, drop = FALSE]

  meta_path  <- file.path(scenic_dir,
                          paste0(lineage_label, "_metadata.csv"))
  write.csv(cell_meta, file = meta_path, row.names = TRUE, quote = FALSE)
  cat(sprintf("Metadata written:           %s\n", meta_path))

  invisible(list(counts = counts_path, norm = norm_path, meta = meta_path))
}

# Export both lineages
export_scenic_csv("CD8_all")
export_scenic_csv("CD4_all")

cat("\n── pySCENIC output directory:", scenic_dir, "\n")
cat("Next steps in Python (pySCENIC):\n")
cat("  1. grn()    — infer co-expression modules (arboreto + GRNBoost2)\n")
cat("  2. ctx()    — prune to regulons using cisTarget motif databases\n")
cat("  3. aucell() — score each cell for regulon activity (AUC)\n")
cat("  4. Load metadata CSV; compare AUC: Expanded vs Not_expanded\n")
cat("             → Wilcoxon test per regulon, BH correction\n")
cat("  5. Stratify by disease (Erosion / Rupture) for differential regulon analysis\n")
