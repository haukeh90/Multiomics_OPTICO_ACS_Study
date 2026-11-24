# ==============================================================================
# OPTICO-ACS | Transcriptional similarity to plaque and blood reference cohorts
#
# Panels : Fig. 3a
# Input  : data/optico_final.rds
#          reference scRNA-seq cohorts (coronary and carotid plaque, healthy and STEMI PBMC; see README)
# Output : Optico_revisions/plots/similarity_analyses/tables/transcriptional_similarity_*.xlsx
#          (sheet Plaque_ref_similarity: plaque similarity index per donor, plotted in Prism)
# Notes  : Steps: (1) load and harmonise the reference cohorts, (2) merge with OPTICO,
#          (3) per-sample QC, (4) pseudobulk per donor and Spearman correlation,
#          (5) extension 1: Erosion and Rupture as separate groups,
#          (6) extension 2: combined plaque reference (carotid + coronary) = Fig. 3a.
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

##OPTICO REVISIOSN
####HAuke Hortsmann 


################################################################################
################################################################################

library(Seurat)
library(Matrix)
library(dplyr)
options(future.globals.maxSize = 20 * 1024^3)  # 20 GB
#load and preprocess Qian data ------------------------------------------------

# -----------------------------
# 1. Set parent directory
# -----------------------------
base_dir <- file.path(OPTICO_DIR, "Optico_revisions/data/GSE269269_RAW")
base_dir <- path.expand(base_dir)

sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

sample_dirs
length(sample_dirs)

# -----------------------------
# 2. Helper: check whether a folder contains 10x files
# -----------------------------
is_10x_dir <- function(x) {
  has_matrix <- file.exists(file.path(x, "matrix.mtx")) ||
    file.exists(file.path(x, "matrix.mtx.gz"))
  
  has_features <- file.exists(file.path(x, "features.tsv")) ||
    file.exists(file.path(x, "features.tsv.gz")) ||
    file.exists(file.path(x, "genes.tsv")) ||
    file.exists(file.path(x, "genes.tsv.gz"))
  
  has_barcodes <- file.exists(file.path(x, "barcodes.tsv")) ||
    file.exists(file.path(x, "barcodes.tsv.gz"))
  
  has_matrix && has_features && has_barcodes
}

# -----------------------------
# 3. More flexible finder
# -----------------------------
find_matrix_dir <- function(sample_path) {
  candidates <- c(
    sample_path,
    file.path(sample_path, "filtered_feature_bc_matrix"),
    file.path(sample_path, "outs", "filtered_feature_bc_matrix"),
    file.path(sample_path, "raw_feature_bc_matrix"),
    file.path(sample_path, "outs", "raw_feature_bc_matrix")
  )
  
  existing <- candidates[sapply(candidates, is_10x_dir)]
  
  if (length(existing) == 0) {
    cat("\nCould not find 10x files in:", sample_path, "\n")
    cat("Contents are:\n")
    print(list.files(sample_path, recursive = TRUE))
    stop(paste("No Cell Ranger / 10x matrix directory found for sample:", sample_path))
  }
  
  return(existing[1])
}

# -----------------------------
# 4. Read and preprocess one sample
# -----------------------------
read_and_filter_sample <- function(sample_path) {
  sample_name <- basename(sample_path)
  matrix_dir  <- find_matrix_dir(sample_path)
  
  cat("Reading:", sample_name, "\n")
  cat("From:", matrix_dir, "\n")
  
  counts <- Read10X(data.dir = matrix_dir)
  
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 0
  )
  
  obj$sample <- sample_name
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 4000 &
      nCount_RNA < 17500 &
      percent.mt < 8
  )
  
  return(obj)
}

# -----------------------------
# 5. Load all samples
# -----------------------------
obj_list <- lapply(sample_dirs, read_and_filter_sample)
names(obj_list) <- basename(sample_dirs)

# -----------------------------
# 6. Merge all samples
# -----------------------------
if (length(obj_list) == 1) {
  quin_data <- obj_list[[1]]
} else {
  quin_data <- merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = names(obj_list),
    project = "quin_data"
  )
}

quin_data
table(quin_data$sample)

# -----------------------------
# 7.Prepocess and save
# -----------------------------
quin_data <- NormalizeData(quin_data)
quin_data <- FindVariableFeatures(quin_data, selection.method = "vst", nfeatures = 2000)
quin_data <- ScaleData(quin_data)
quin_data <- JoinLayers(quin_data)
quin_data[["RNA3"]] <- as(object = quin_data[["RNA"]], Class = "Assay")

# -----------------------------
# 8.Clean metadata and assay
# -----------------------------

#remove V5 assay
DefaultAssay(quin_data) <- "RNA3"
quin_data[['RNA']] <- NULL

#harmonize medata 
quin_data$Pat.ID <- quin_data$orig.ident
Idents(quin_data) <- "Pat.ID"
levels(quin_data)
quin_data$Group <- "Periphery"

#save 
saveRDS(quin_data, 
        file.path(OPTICO_DIR, "Optico_revisions/data/quin_data.rds"))

#load Groningen data  ----------------------------------------------------------

# -----------------------------
# 1.Load
# -----------------------------
Groningen_only <- readRDS(file.path(OPTICO_DIR, "external_data/Groningen_only.rds"))            # healthy and STEMI PBMC

# -----------------------------
# 2. Subset on samples of interest
# -----------------------------
Idents(Groningen_only) <- "timepoint.final"
levels(Groningen_only)
Groningen_acute <- subset(Groningen_only, idents = c("Baseline", "UT"))

# -----------------------------
# 3. Subset as above 
# -----------------------------
table(Groningen_acute$Pat.ID)
Groningen_acute <- subset(
  Groningen_acute,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 4000 &
    nCount_RNA < 17500 &
    percent.mt < 8
)
table(Groningen_acute$Pat.ID)

# -----------------------------
# 4. Preprocess as above 
# -----------------------------
Groningen_acute <- NormalizeData(Groningen_acute)
Groningen_acute <- FindVariableFeatures(Groningen_acute, selection.method = "vst", nfeatures = 2000)
Groningen_acute <- ScaleData(Groningen_acute)
Groningen_acute <- JoinLayers(Groningen_acute)
Groningen_acute[["RNA3"]] <- as(object = Groningen_acute[["RNA"]], Class = "Assay")

# -----------------------------
# 5. Clean metadata and assay
# -----------------------------

#inspect 
Groningen_acute@meta.data

# Columns you want to keep
cols_to_keep <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "orig.ident",
  "Pat.ID",
  "timepoint.final"
)

# Keep only those columns that actually exist
cols_to_keep <- cols_to_keep[cols_to_keep %in% colnames(Groningen_acute@meta.data)]

# Replace metadata with cleaned version
Groningen_acute@meta.data <- Groningen_acute@meta.data[, cols_to_keep, drop = FALSE]

#remove V5 assay
DefaultAssay(Groningen_acute) <- "RNA3"
Groningen_acute[['RNA']] <- NULL

#harmonize metadata 
Groningen_acute$Group <- Groningen_acute$timepoint.final

# -----------------------------
# 6. Save
# -----------------------------
saveRDS(Groningen_acute, 
        file.path(OPTICO_DIR, "Optico_revisions/data/Groningen_acute.rds"))

#load Coronary data  ----------------------------------------------------------

# -----------------------------
# 1.Load
# -----------------------------
Coronary_data <- readRDS(file.path(OPTICO_DIR, "Optico_revisions/data/Coronary_data.rds"))
DefaultAssay(Coronary_data) <- "RNA"
Coronary_data <- DietSeurat(Coronary_data, assays = "RNA")

# -----------------------------
# 2. Subset on immune cells 
# -----------------------------
Idents(Coronary_data) <- "cell.type"
levels(Coronary_data)
Coronary_data <- subset(Coronary_data, idents = c("SMCPericyte", "ModSMC", "Fibroblast1", "Fibroblast2", "Endothelium", "Glia"), 
                                                  invert = T)
levels(Coronary_data)

# -----------------------------
# 3. Subset as above 
# -----------------------------
Coronary_data$Pat.ID <- Coronary_data$sampleID
table(Coronary_data$Pat.ID)
Coronary_data <- subset(
  Coronary_data,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 4000 &
    nCount_RNA < 17500 &
    percent.mt < 8
)
table(Coronary_data$Pat.ID)

# -----------------------------
# 4. Preprocess as above 
# -----------------------------
Coronary_data <- NormalizeData(Coronary_data)
Coronary_data <- FindVariableFeatures(Coronary_data, selection.method = "vst", nfeatures = 2000)
Coronary_data <- ScaleData(Coronary_data)

# -----------------------------
# 5. Clean metadata 
# -----------------------------

#inspect 
Coronary_data@meta.data
Coronary_data$Group <- "Coronary"

# Columns you want to keep
cols_to_keep <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "orig.ident",
  "Pat.ID",
  "Group"
)

# Keep only those columns that actually exist
cols_to_keep <- cols_to_keep[cols_to_keep %in% colnames(Coronary_data@meta.data)]

# Replace metadata with cleaned version
Coronary_data@meta.data <- Coronary_data@meta.data[, cols_to_keep, drop = FALSE]

# -----------------------------
# 6. Save
# -----------------------------
saveRDS(Coronary_data, 
        file.path(OPTICO_DIR, "Optico_revisions/data/Coronary_data.rds"))


#load Giannarelli Carotid data  -------------------------------------------------

# -----------------------------
# 1.Load
# -----------------------------
Giannareli_carotid <- readRDS(file.path(OPTICO_DIR, "external_data/Giannareli_carotid.rds"))        # carotid plaque
DimPlot(Giannareli_carotid)
DefaultAssay(Giannareli_carotid) <- "RNA"
Giannareli_carotid <- DietSeurat(Giannareli_carotid, assays = "RNA")

# -----------------------------
# 3. Subset as above 
# -----------------------------
Giannareli_carotid$Pat.ID <- Giannareli_carotid$new.ident
table(Giannareli_carotid$Pat.ID)
Giannareli_carotid <- subset(
  Giannareli_carotid,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 4000 &
    nCount_RNA < 17500 &
    percent.mt < 8
)
table(Giannareli_carotid$Pat.ID)

# -----------------------------
# 4. Preprocess as above 
# -----------------------------
Giannareli_carotid <- NormalizeData(Giannareli_carotid)
Giannareli_carotid <- FindVariableFeatures(Giannareli_carotid, selection.method = "vst", nfeatures = 2000)
Giannareli_carotid <- ScaleData(Giannareli_carotid)

# -----------------------------
# 5. Clean metadata 
# -----------------------------

#inspect 
Giannareli_carotid@meta.data
Giannareli_carotid$Group <- "Carotid"

# Columns you want to keep
cols_to_keep <- c(
  "nFeature_RNA",
  "nCount_RNA",
  "percent.mt",
  "orig.ident",
  "Pat.ID",
  "Group"
)

# Keep only those columns that actually exist
cols_to_keep <- cols_to_keep[cols_to_keep %in% colnames(Giannareli_carotid@meta.data)]

# Replace metadata with cleaned version
Giannareli_carotid@meta.data <- Giannareli_carotid@meta.data[, cols_to_keep, drop = FALSE]

# -----------------------------
# 6. Save
# -----------------------------
saveRDS(Giannareli_carotid, 
        file.path(OPTICO_DIR, "Optico_revisions/data/Giannareli_carotid.rds"))

# Align object assays and metadata ---------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)

#load and prepare plaque object
CVR_plaques_mapping <- readRDS(file.path(OPTICO_DIR, "external_data/CVR_plaques_mapping.rds"))       # carotid plaque reference
Idents(CVR_plaques_mapping) <- "species"
levels(CVR_plaques_mapping)
CVR_plaques_mapping <- subset(CVR_plaques_mapping, idents = "human")
CVR_plaques_mapping$Pat.ID <- CVR_plaques_mapping$orig.ident
CVR_plaques_mapping$Group <- "Carotid" 

#load other object
Giannareli_carotid <-  readRDS(file.path(OPTICO_DIR, "Optico_revisions/data/Giannareli_carotid.rds"))
Coronary_data <- readRDS(file.path(OPTICO_DIR, "Optico_revisions/data/Coronary_data.rds"))
Groningen_acute <- readRDS(file.path(OPTICO_DIR, "Optico_revisions/data/Groningen_acute.rds"))
quin_data <- readRDS(file.path(OPTICO_DIR, "Optico_revisions/data/quin_data.rds"))
optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))
optico_final$Pat.ID <- optico_final$patients
optico_final$Group <- "Culprit"

#Set default assay to RNA
DefaultAssay(CVR_plaques_mapping) <- "RNA"
DefaultAssay(Giannareli_carotid) <- "RNA"
DefaultAssay(Coronary_data) <- "RNA"
DefaultAssay(optico_final) <- "RNA"

#Fix RNA3 issue 
quin_data[['RNA']] = quin_data[['RNA3']]
Groningen_acute[['RNA']] = Groningen_acute[['RNA3']]
DefaultAssay(quin_data) <- "RNA"
DefaultAssay(Groningen_acute) <- "RNA"

#DietSeurat 
CVR_plaques_mapping <- DietSeurat(CVR_plaques_mapping, assay = "RNA", layer = c("counts", "data"))
Coronary_data <- DietSeurat(Coronary_data, assay = "RNA", layer = c("counts", "data"))
Groningen_acute <- DietSeurat(Groningen_acute, assay = "RNA", layer = c("counts", "data"))
quin_data <- DietSeurat(quin_data, assay = "RNA", layer = c("counts", "data"))
optico_final <- DietSeurat(optico_final, assay = "RNA", layer = c("counts", "data"))

#Add relevent metadata columns 
unique(CVR_plaques_mapping$Pat.ID)
unique(CVR_plaques_mapping$Group)
unique(Coronary_data$Pat.ID)
unique(Coronary_data$Group)
unique(Giannareli_carotid$Pat.ID)
unique(Giannareli_carotid$Group)
unique(Groningen_acute$Pat.ID)
unique(Groningen_acute$Group)
unique(quin_data$Pat.ID)
unique(quin_data$Group)
unique(optico_final$Pat.ID)
unique(optico_final$Group)

# Merge all objects
merged_1 <- merge(
  x = CVR_plaques_mapping,
  y = list(Coronary_data, Giannareli_carotid),
  project = "MergedPlaques"
)
rm(CVR_plaques_mapping, Coronary_data, Giannareli_carotid)
gc()

merged_2 <- merge(
  x = Groningen_acute,
  y = list(quin_data, optico_final),
  project = "MergedPlaques"
)
rm(Groningen_acute, quin_data, optico_final)
gc()

merged_obj <- merge(
  x = merged_1,
  y = merged_2,
  project = "MergedPlaques"
)

saveRDS(merged_obj, file.path(OPTICO_DIR, "Optico_revisions/data/merged_all_plaques.rds"))

rm(merged_1, merged_2)

#################################QC#############################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)

# ============================================================
# PATHS
# ============================================================
base_dir  <- file.path(OPTICO_DIR, "Optico_revisions/plots/similarity_analyses")
plot_dir  <- file.path(base_dir, "plots")
table_dir <- file.path(base_dir, "tables")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

#extract numbers for schematic
seurat_obj <- readRDS(file.path(OPTICO_DIR, "Optico_revisions/data/merged_all_plaques.rds"))

table(seurat_obj$Pat.ID, seurat_obj$Group)

# ============================================================
# 1) COMPUTE PER-CELL QC METRICS (if not already present)
# ============================================================
# Mitochondrial %
if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}
# Ribosomal % (large + small subunit)
if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
}
# Haemoglobin % (flags red blood cell contamination)
if (!"percent.hb" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^BP]")
}
# Transcriptional complexity: log10(genes) / log10(UMIs)
# Values well below ~0.8 indicate low-complexity / dying cells
seurat_obj[["log10GenesPerUMI"]] <- log10(seurat_obj$nFeature_RNA) /
  log10(seurat_obj$nCount_RNA)

# ============================================================
# 2) AGGREGATE TO PER-SAMPLE SUMMARY
# ============================================================
qc_per_sample <- seurat_obj@meta.data %>%
  group_by(Group, Pat.ID) %>%
  summarise(
    n_cells            = n(),
    mean_nCount        = mean(nCount_RNA),
    median_nCount      = median(nCount_RNA),
    mean_nFeature      = mean(nFeature_RNA),
    median_nFeature    = median(nFeature_RNA),
    mean_pct_mt        = mean(percent.mt),
    mean_pct_rb        = mean(percent.rb),
    mean_pct_hb        = mean(percent.hb),
    mean_complexity    = mean(log10GenesPerUMI),
    # Fraction of cells that look like low-quality / dying
    # (high mt is the classic marker; >25 % is a common threshold)
    frac_high_mt       = mean(percent.mt > 25),
    # Fraction of cells with very low counts (potential empty droplets)
    frac_low_counts    = mean(nCount_RNA < 500),
    .groups = "drop"
  ) %>%
  arrange(Group, Pat.ID)

cat("\nPer-sample QC summary:\n")
print(qc_per_sample, n = Inf)

# ============================================================
# 3) FLAG OUTLIERS (per group, using median ± 2 × MAD)
# ============================================================
flag_outliers <- function(x, k = 2) {
  m   <- median(x, na.rm = TRUE)
  mad <- mad(x, na.rm = TRUE)
  abs(x - m) > k * mad
}

qc_flagged <- qc_per_sample %>%
  group_by(Group) %>%
  mutate(
    flag_nCount     = flag_outliers(mean_nCount),
    flag_nFeature   = flag_outliers(mean_nFeature),
    flag_mt         = flag_outliers(mean_pct_mt),
    flag_rb         = flag_outliers(mean_pct_rb),
    flag_hb         = flag_outliers(mean_pct_hb),
    flag_complexity = flag_outliers(mean_complexity),
    flag_n_cells    = flag_outliers(n_cells),
    n_flags         = flag_nCount + flag_nFeature + flag_mt +
      flag_rb + flag_hb + flag_complexity + flag_n_cells
  ) %>%
  ungroup() %>%
  arrange(desc(n_flags), Group, Pat.ID)

cat("\nSamples with ≥1 outlier flag:\n")
print(filter(qc_flagged, n_flags >= 1), n = Inf)

# ============================================================
# 4) VISUALISE — dotplot of key metrics, coloured by group
# ============================================================
metrics_long <- qc_per_sample %>%
  pivot_longer(
    cols = c(mean_nCount, mean_nFeature, mean_pct_mt,
             mean_pct_rb, mean_pct_hb, mean_complexity, n_cells),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = recode(metric,
                         mean_nCount    = "Mean nCount_RNA",
                         mean_nFeature  = "Mean nFeature_RNA",
                         mean_pct_mt    = "Mean % MT",
                         mean_pct_rb    = "Mean % Ribosomal",
                         mean_pct_hb    = "Mean % Haemoglobin",
                         mean_complexity = "Mean Complexity\n(log10Gene/log10UMI)",
                         n_cells        = "N cells"
  ))

# Highlight the four suspected outlier Carotid donors
suspects <- c("CV7209", "CV7215", "CV7223", "CV7231A")

p_qc <- ggplot(metrics_long,
               aes(x = Group, y = value, colour = Group)) +
  geom_jitter(aes(shape = Pat.ID %in% suspects),
              width = 0.2, size = 2, alpha = 0.8) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8),
                     labels = c("Other", "Suspected outlier"),
                     name   = NULL) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 35, hjust = 1),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  ) +
  labs(x = NULL, y = NULL,
       title = "Per-sample QC metrics",
       colour = "Group")

ggsave(file.path(plot_dir, "per_sample_qc_dotplot.pdf"),
       plot = p_qc, width = 14, height = 18, units = "cm")

# ============================================================
# 5) EXPORT
# ============================================================
wb <- createWorkbook()
addWorksheet(wb, "QC_per_sample")
writeDataTable(wb, "QC_per_sample", qc_per_sample)
addWorksheet(wb, "QC_outlier_flags")
writeDataTable(wb, "QC_outlier_flags", qc_flagged)
saveWorkbook(wb, file.path(table_dir, "per_sample_qc.xlsx"), overwrite = TRUE)

cat("\nDone. QC table -> ", file.path(table_dir, "per_sample_qc.xlsx"),
    "\n      QC plot  -> ", file.path(plot_dir,  "per_sample_qc_dotplot.pdf"), "\n")


# ============================================================
# 6) Identify and remove low QC samples
# ============================================================

# --- Directional, within-group thresholds -------------------
# Only the low-quality tail: too few cells, too few genes, blood contamination
qc_per_sample <- qc_per_sample %>%
  group_by(Group) %>%
  mutate(
    # nFeature: flag if below within-group median - 2*MAD
    nFeature_threshold = median(mean_nFeature) - 2 * mad(mean_nFeature),
    flag_low_nFeature = mean_nFeature < nFeature_threshold & mean_nFeature < 1500,
    # Haemoglobin: absolute threshold (>0.15% is elevated regardless of group)
    flag_high_hb       = mean_pct_hb > 0.3,
    # Too few cells for a stable pseudobulk
    flag_low_cells     = n_cells < 200
  ) %>%
  ungroup() %>%
  mutate(
    exclude = flag_low_nFeature | flag_high_hb | flag_low_cells
  )

# --- Report what gets removed --------------------------------
excluded <- filter(qc_per_sample, exclude)
cat("\n=== Samples excluded by QC ===\n")
print(
  excluded %>%
    select(Group, Pat.ID, n_cells, mean_nFeature,
           mean_pct_hb, flag_low_cells, flag_low_nFeature, flag_high_hb),
  n = Inf
)
cat("\nTotal excluded:", nrow(excluded), "of", nrow(qc_per_sample), "donors\n")

# --- Remove their cells from the Seurat object --------------
Idents(seurat_obj) <- "Pat.ID"
levels(seurat_obj)
seurat_obj <- subset(seurat_obj, idents = c("S22", "S6", "CV7209", "CV7221", "CV7231A", "NPR2"), invert = T)

#######################SIMILARITY analysis#####################################
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(openxlsx)
library(pheatmap)

# ============================================================
# 0) PATHS  — only edit this block if you move the project
# ============================================================
base_dir  <- file.path(OPTICO_DIR, "Optico_revisions/plots/similarity_analyses")
plot_dir  <- file.path(base_dir, "plots")
table_dir <- file.path(base_dir, "tables")

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) LOAD & PREP
# ============================================================
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- FindVariableFeatures(seurat_obj)

# Create a unique donor ID combining Group + Pat.ID
# This avoids collisions if the same Pat.ID string appears in multiple datasets
seurat_obj$unique_donor_id <- paste(seurat_obj$Group, seurat_obj$Pat.ID, sep = "_")

# ============================================================
# 2) PSEUDOBULK (sum of normalised counts per donor)
# ============================================================
pb <- AggregateExpression(
  seurat_obj,
  group.by = "unique_donor_id",   # one pseudobulk profile per donor
  assays   = "RNA",
  slot     = "data"               # normalised (log1p) counts
)$RNA

# Restrict to variable genes present in the pseudobulk matrix
var_genes <- intersect(VariableFeatures(seurat_obj), rownames(pb))
cat("Variable genes used for similarity:", length(var_genes), "\n")

pb <- as.matrix(pb[var_genes, ])

# ============================================================
# 3) SPEARMAN CORRELATION MATRIX
# ============================================================
cor_mat <- cor(pb, method = "spearman")

# Map each column back to its Group.
# AggregateExpression can mangle the separator in unique_donor_id
# (e.g. "Baseline_S1" may become "Baseline-S1"), so separator-based
# matching is fragile.  Instead we use prefix matching: because
# unique_donor_id = paste(Group, Pat.ID, sep = <anything>), the column
# name will always START with the Group string regardless of separator.
cn                <- colnames(cor_mat)
all_group_values  <- unique(seurat_obj$Group)

group_lab <- sapply(cn, function(x) {
  hit <- all_group_values[startsWith(x, all_group_values)]
  if (length(hit) == 1L) hit else NA_character_
})

# Diagnostic: warn if any donors still unresolved
unmatched <- cn[is.na(group_lab)]
if (length(unmatched) > 0) {
  warning(length(unmatched), " pb column(s) could not be matched to a Group ",
          "— check that Group values do not share a common prefix:\n",
          paste(head(unmatched, 20), collapse = ", "))
}

cat("\nDonors per group:\n")
print(table(group_lab))

# ============================================================
# 4) OPTIONAL GLOBAL HEATMAP — all donors, all groups
# ============================================================
library(colorspace)

# Desaturated group annotation colours
all_groups   <- unique(group_lab)
base_cols    <- c("#E64B35", "#F39B7F", "#00A087", "#4DBBD5",
                  "#3C5488", "#a65628", "#f781bf", "#999999")
desat_cols   <- colorspace::desaturate(base_cols[seq_along(all_groups)], amount = 0.45)
group_colors <- setNames(desat_cols, all_groups)

ann_col_full <- data.frame(Group = group_lab, row.names = cn)

pheatmap::pheatmap(
  cor_mat,
  # reversed inferno: high correlation = dark (black/purple), low = cream
  color                    = rev(viridisLite::magma(100)),
  border_color             = "white",
  # draws black lines around the k main clusters — adjust k to taste
  cutree_rows              = 5,
  cutree_cols              = 5,
  annotation_col           = ann_col_full,
  annotation_colors        = list(Group = group_colors),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames            = TRUE,
  show_colnames            = TRUE,
  fontsize                 = 6,
  # tighten annotation bar and remove the legend gap
  annotation_legend        = TRUE,
  annotation_names_col     = FALSE,
  filename                 = file.path(plot_dir, "global_similarity_heatmap.pdf"),
  width = 16, height = 14
)

# ============================================================
# 5) CORE FUNCTION: run similarity analysis for one reference
# ============================================================
run_similarity <- function(cor_mat, group_lab, reference, output_prefix) {
  
  groups            <- unique(group_lab)
  comparison_groups <- setdiff(groups, reference)
  
  ref_cols  <- which(group_lab == reference)
  stopifnot(length(ref_cols) > 0)
  
  # Extract similarity submatrix for each comparison group vs. reference
  sim_list <- lapply(comparison_groups, function(g) {
    g_cols <- which(group_lab == g)
    cor_mat[g_cols, ref_cols, drop = FALSE]
  })
  names(sim_list) <- comparison_groups
  
  # ----------------------------------------------------------
  # 5a) Focused heatmap (reference + all comparison groups)
  # ----------------------------------------------------------
  order_idx <- c(ref_cols,
                 unlist(lapply(comparison_groups,
                               function(g) which(group_lab == g))))
  cor_sub       <- cor_mat[order_idx, order_idx]
  group_ordered <- group_lab[order_idx]
  
  ann_col_sub <- data.frame(Group = group_ordered,
                            row.names = colnames(cor_sub))
  
  pheatmap::pheatmap(
    cor_sub,
    annotation_col           = ann_col_sub,
    annotation_colors        = list(Group = group_colors),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames            = TRUE,
    show_colnames            = TRUE,
    fontsize                 = 7,
    filename                 = file.path(plot_dir, paste0(output_prefix, "_heatmap.pdf")),
    width = 14, height = 12
  )
  
  # ----------------------------------------------------------
  # 5b) Tidy data frame of similarity values
  # ----------------------------------------------------------
  df_sim <- bind_rows(lapply(comparison_groups, function(g) {
    data.frame(
      similarity = as.vector(sim_list[[g]]),
      group      = g,
      reference  = reference,
      contrast   = paste0(g, " vs ", reference)
    )
  }))
  
  # Order contrasts by median similarity (descending)
  group_medians <- df_sim %>%
    group_by(contrast) %>%
    summarise(med = median(similarity), .groups = "drop") %>%
    arrange(desc(med))
  df_sim$contrast <- factor(df_sim$contrast, levels = group_medians$contrast)
  
  # ----------------------------------------------------------
  # 5c) Pairwise Wilcoxon tests (two-sided, BH correction)
  # ----------------------------------------------------------
  pairs <- combn(comparison_groups, 2, simplify = FALSE)
  
  stats_list <- lapply(pairs, function(pair) {
    g1 <- pair[1]; g2 <- pair[2]
    wt <- wilcox.test(
      as.vector(sim_list[[g1]]),
      as.vector(sim_list[[g2]]),
      alternative = "two.sided"
    )
    data.frame(
      group_1    = g1,
      group_2    = g2,
      reference  = reference,
      wilcoxon_W = as.numeric(wt$statistic),
      p_value    = wt$p.value
    )
  })
  stats_tbl      <- bind_rows(stats_list)
  stats_tbl$p_adj <- p.adjust(stats_tbl$p_value, method = "BH")
  
  sig_label <- function(p) {
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  }
  stats_tbl$significance <- sig_label(stats_tbl$p_adj)
  
  cat("\n=== Reference:", reference,
      "— pairwise comparisons of similarity scores ===\n")
  print(stats_tbl)
  
  # ----------------------------------------------------------
  # 5d) Boxplot
  # ----------------------------------------------------------
  n_contrasts  <- length(comparison_groups)
  plot_palette <- setNames(
    base_cols[seq_along(comparison_groups)],
    levels(df_sim$contrast)
  )
  
  p_box <- ggplot(df_sim,
                  aes(x = contrast, y = similarity, fill = contrast)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha     = 0.75,
      width     = 0.55,
      color     = "black",
      linewidth = 0.4
    ) +
    geom_jitter(
      aes(color = contrast),
      width = 0.15, alpha = 0.55, size = 1.8
    ) +
    scale_fill_manual(values  = plot_palette) +
    scale_color_manual(values = plot_palette) +
    coord_flip() +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text       = element_text(size = 11),
      axis.title      = element_text(size = 12)
    ) +
    labs(
      x     = NULL,
      y     = paste0("Spearman r (similarity to ", reference, ")"),
      title = paste("Transcriptional similarity — reference:", reference)
    )
  
  plot_height <- max(8, 4 + n_contrasts * 2)
  ggsave(
    file.path(plot_dir, paste0(output_prefix, "_similarity_boxplot.pdf")),
    plot   = p_box,
    units  = "cm",
    height = plot_height,
    width  = 16
  )
  
  list(df_sim = df_sim, stats = stats_tbl, plot = p_box)
}

# ============================================================
# EXTENSION 1: Similarity analysis with Culprit split into
#            Erosion and Rupture subgroups
#
# Assumptions:
#   - seurat_obj already QC-filtered (as above)
#   - seurat_obj$Group == "Culprit" for all culprit samples
#   - seurat_obj$Condition (or equivalent) holds "Erosion" / "Rupture"
#   - run_similarity(), plot_dir, table_dir, base_cols, var_genes
#     are all already defined from the main script above
# ============================================================

condition_col <- "disease"   # metadata column holding Erosion / Rupture

# ============================================================
# A1) BUILD SPLIT GROUPING VARIABLE
#     Keep all non-Culprit groups unchanged;
#     replace "Culprit" with "Erosion" or "Rupture"
# ============================================================
seurat_obj$Group_split <- ifelse(
  seurat_obj$Group == "Culprit",
  seurat_obj@meta.data[[condition_col]],   # "Erosion" or "Rupture"
  seurat_obj$Group
)

cat("\nGroup_split distribution:\n")
print(table(seurat_obj$Group_split))

# ============================================================
# A2) NEW UNIQUE DONOR ID using Group_split
# ============================================================
seurat_obj$unique_donor_id_split <- paste(
  seurat_obj$Group_split,
  seurat_obj$Pat.ID,
  sep = "_"
)

# ============================================================
# A3) PSEUDOBULK with split grouping
# ============================================================
pb_split <- AggregateExpression(
  seurat_obj,
  group.by = "unique_donor_id_split",
  assays   = "RNA",
  slot     = "data"
)$RNA

# Restrict to the same variable genes as the main analysis
pb_split <- as.matrix(pb_split[intersect(var_genes, rownames(pb_split)), ])
cat("Donors in split pseudobulk:", ncol(pb_split), "\n")

# ============================================================
# A4) SPEARMAN CORRELATION MATRIX
# ============================================================
cor_mat_split <- cor(pb_split, method = "spearman")

# Map column names back to Group_split labels
cn_split          <- colnames(cor_mat_split)
all_split_groups  <- unique(seurat_obj$Group_split)

group_lab_split <- sapply(cn_split, function(x) {
  hit <- all_split_groups[startsWith(x, all_split_groups)]
  if (length(hit) == 1L) hit else NA_character_
})

unmatched_split <- cn_split[is.na(group_lab_split)]
if (length(unmatched_split) > 0) {
  warning(length(unmatched_split),
          " column(s) could not be matched to Group_split:\n",
          paste(head(unmatched_split, 20), collapse = ", "))
}

cat("\nDonors per split group:\n")
print(table(group_lab_split))

# ============================================================
# A5) COLOUR PALETTE — add Erosion and Rupture entries
# ============================================================
split_group_colors <- c(
  "Erosion"   = "#E64B35",   # red-orange  — erosion
  "Rupture"   = "#3C5488",   # navy        — rupture
  "Carotid"   = "#A08060",   # warm brown
  "Coronary"  = "#BDBDBD",   # neutral grey
  "Baseline"  = "#D4B483",
  "Periphery" = "#6495C8",
  "UT"        = "#444444"
)
split_group_colors <- split_group_colors[
  names(split_group_colors) %in% unique(group_lab_split)
]

# For run_similarity() the palette it uses internally is base_cols;
# override group_colors in its environment so the boxplots pick up
# the right colours for Erosion/Rupture
group_colors <- split_group_colors   # run_similarity() reads this from .GlobalEnv

# ============================================================
# A7) SIMILARITY vs CORONARY — Erosion and Rupture separately
# ============================================================
cat("\n\n########################################")
cat("\n# SPLIT ANALYSIS 1: Coronary as reference")
cat("\n########################################\n")

results_coronary_split <- run_similarity(
  cor_mat       = cor_mat_split,
  group_lab     = group_lab_split,
  reference     = "Coronary",
  output_prefix = "coronary_ref_erosion_rupture"
)

# ============================================================
# A8) SIMILARITY vs CAROTID — Erosion and Rupture separately
# ============================================================
cat("\n\n########################################")
cat("\n# SPLIT ANALYSIS 2: Carotid as reference")
cat("\n########################################\n")

results_carotid_split <- run_similarity(
  cor_mat       = cor_mat_split,
  group_lab     = group_lab_split,
  reference     = "Carotid",
  output_prefix = "carotid_ref_erosion_rupture"
)

# ============================================================
# A9) DIRECT EROSION vs RUPTURE SIMILARITY SCORES
#     (pairwise r values between the two subgroups)
# ============================================================
ero_idx <- which(group_lab_split == "Erosion")
rup_idx <- which(group_lab_split == "Rupture")

er_vs_rup <- cor_mat_split[ero_idx, rup_idx, drop = FALSE]
cat("\nMean Spearman r — Erosion vs Rupture:",
    round(mean(er_vs_rup), 4), "\n")

# ============================================================
# A10) CORRPLOT — group-level mean similarity (split version)
# ============================================================
groups_split_ordered <- c("Erosion", "Rupture",
                          setdiff(unique(group_lab_split),
                                  c("Erosion", "Rupture")))
n_grp_split  <- length(groups_split_ordered)

grp_sim_mat_split <- matrix(
  NA_real_,
  nrow     = n_grp_split, ncol = n_grp_split,
  dimnames = list(groups_split_ordered, groups_split_ordered)
)
for (g1 in groups_split_ordered) {
  for (g2 in groups_split_ordered) {
    i1 <- which(group_lab_split == g1)
    i2 <- which(group_lab_split == g2)
    grp_sim_mat_split[g1, g2] <- mean(cor_mat_split[i1, i2, drop = FALSE])
  }
}

col_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)

pdf(file.path(plot_dir, "erosion_rupture_split_corrplot.pdf"),
    width = 7, height = 6)
corrplot::corrplot(
  grp_sim_mat_split,
  method        = "color",
  type          = "upper",
  col           = col_pal,
  addCoef.col   = "black",
  number.cex    = 0.85,
  number.digits = 3,
  tl.col        = "black",
  tl.cex        = 0.95,
  tl.srt        = 45,
  cl.lim        = c(min(grp_sim_mat_split), 1),
  is.corr       = FALSE,
  mar           = c(0, 0, 2, 0),
  title         = "Group-level similarity — Erosion / Rupture split (mean Spearman r)"
)
dev.off()
cat("\nSplit corrplot saved ->",
    file.path(plot_dir, "erosion_rupture_split_corrplot.pdf"), "\n")

# ============================================================
# A11) EXPORT
# ============================================================
wb_split <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb_split, "Coronary_Erosion_Rupture_sim")
openxlsx::writeDataTable(wb_split, "Coronary_Erosion_Rupture_sim",
                         results_coronary_split$df_sim)

openxlsx::addWorksheet(wb_split, "Coronary_Erosion_Rupture_stats")
openxlsx::writeDataTable(wb_split, "Coronary_Erosion_Rupture_stats",
                         results_coronary_split$stats)

openxlsx::addWorksheet(wb_split, "Carotid_Erosion_Rupture_sim")
openxlsx::writeDataTable(wb_split, "Carotid_Erosion_Rupture_sim",
                         results_carotid_split$df_sim)

openxlsx::addWorksheet(wb_split, "Carotid_Erosion_Rupture_stats")
openxlsx::writeDataTable(wb_split, "Carotid_Erosion_Rupture_stats",
                         results_carotid_split$stats)

grp_sim_split_df <- as.data.frame(grp_sim_mat_split)
grp_sim_split_df <- cbind(Group = rownames(grp_sim_split_df), grp_sim_split_df)
openxlsx::addWorksheet(wb_split, "Group_mean_sim_split")
openxlsx::writeDataTable(wb_split, "Group_mean_sim_split", grp_sim_split_df)

openxlsx::saveWorkbook(
  wb_split,
  file.path(table_dir, "transcriptional_similarity_erosion_rupture_split.xlsx"),
  overwrite = TRUE
)

cat("\nDone — split analysis complete.\n",
    "  Plots  ->", plot_dir,  "\n",
    "  Tables ->", table_dir, "\n")

# ============================================================
# EXTENSION 2: Similarity analysis — Erosion & Rupture vs a
#              combined "Plaque" reference (Carotid + Coronary)
#
# Requires: Extension 1 already run, so the following exist:
#   cor_mat_split   — Spearman r matrix (donors × donors)
#   group_lab_split — named vector mapping donor IDs to groups
#   plot_dir, table_dir, run_similarity()
#
# What this does:
#   Relabels Carotid and Coronary as "Plaque" in the group
#   label vector, then runs run_similarity() with "Plaque" as
#   the reference. No pseudobulking needed — same cor matrix.
# ============================================================

# ============================================================
# B1) REMAP GROUP LABELS — merge Carotid + Coronary → "Plaque"
# ============================================================
group_lab_plaque <- group_lab_split
group_lab_plaque[group_lab_plaque %in% c("Carotid", "Coronary")] <- "Plaque"

cat("\nGroup_lab_plaque distribution:\n")
print(table(group_lab_plaque))

# ============================================================
# B2) COLOUR PALETTE for combined analysis
# ============================================================
plaque_group_colors <- c(
  "Erosion"   = "#E64B35",   # red-orange
  "Rupture"   = "#3C5488",   # navy
  "Plaque"    = "#7B6E8B",   # muted purple  — combined plaque reference
  "Baseline"  = "#D4B483",
  "Periphery" = "#6495C8",
  "UT"        = "#444444"
)
plaque_group_colors <- plaque_group_colors[
  names(plaque_group_colors) %in% unique(group_lab_plaque)
]

# Override group_colors so run_similarity() picks up new palette
group_colors <- plaque_group_colors

# ============================================================
# B4) SIMILARITY vs COMBINED PLAQUE REFERENCE
# ============================================================
cat("\n\n########################################")
cat("\n# COMBINED PLAQUE ANALYSIS: Plaque as reference")
cat("\n########################################\n")

results_plaque <- run_similarity(
  cor_mat       = cor_mat_split,
  group_lab     = group_lab_plaque,
  reference     = "Plaque",
  output_prefix = "combined_plaque_ref_erosion_rupture"
)

# ============================================================
# B5) CORRPLOT — group-level mean similarity (combined plaque)
# ============================================================
groups_plaque_ordered <- c("Erosion", "Rupture", "Plaque",
                           setdiff(unique(group_lab_plaque),
                                   c("Erosion", "Rupture", "Plaque")))
n_grp_plaque <- length(groups_plaque_ordered)

grp_sim_mat_plaque <- matrix(
  NA_real_,
  nrow     = n_grp_plaque, ncol = n_grp_plaque,
  dimnames = list(groups_plaque_ordered, groups_plaque_ordered)
)
for (g1 in groups_plaque_ordered) {
  for (g2 in groups_plaque_ordered) {
    i1 <- which(group_lab_plaque == g1)
    i2 <- which(group_lab_plaque == g2)
    grp_sim_mat_plaque[g1, g2] <- mean(cor_mat_split[i1, i2, drop = FALSE])
  }
}

cat("\nGroup-level mean similarity matrix (combined plaque):\n")
print(round(grp_sim_mat_plaque, 4))

col_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)

pdf(file.path(plot_dir, "combined_plaque_corrplot.pdf"),
    width = 7, height = 6)
corrplot::corrplot(
  grp_sim_mat_plaque,
  method        = "color",
  type          = "upper",
  col           = col_pal,
  addCoef.col   = "black",
  number.cex    = 0.85,
  number.digits = 3,
  tl.col        = "black",
  tl.cex        = 0.95,
  tl.srt        = 45,
  cl.lim        = c(min(grp_sim_mat_plaque), 1),
  is.corr       = FALSE,
  mar           = c(0, 0, 2, 0),
  title         = "Group similarity — Erosion / Rupture vs combined Plaque (mean Spearman r)"
)
dev.off()
cat("Corrplot saved ->",
    file.path(plot_dir, "combined_plaque_corrplot.pdf"), "\n")

# ============================================================
# B6) PRINT KEY COMPARISONS FOR REPORTING
# ============================================================
cat("\n=== Key similarity scores (mean Spearman r) ===\n")
cat("Erosion  vs Plaque :", round(grp_sim_mat_plaque["Erosion",  "Plaque"], 4), "\n")
cat("Rupture  vs Plaque :", round(grp_sim_mat_plaque["Rupture",  "Plaque"], 4), "\n")
cat("Erosion  vs Rupture:", round(grp_sim_mat_plaque["Erosion",  "Rupture"], 4), "\n")

# Print Wilcoxon result for Erosion vs Rupture similarity-to-plaque
stats <- results_plaque$stats
er_row <- stats[
  (stats$group_1 == "Erosion" & stats$group_2 == "Rupture") |
    (stats$group_1 == "Rupture" & stats$group_2 == "Erosion"), ]
if (nrow(er_row) > 0) {
  cat("\nWilcoxon test — Erosion vs Rupture similarity to Plaque:\n")
  print(er_row)
}

# ============================================================
# B7) EXPORT
# ============================================================
wb_plaque <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb_plaque, "Plaque_ref_similarity")
openxlsx::writeDataTable(wb_plaque, "Plaque_ref_similarity",
                         results_plaque$df_sim)

openxlsx::addWorksheet(wb_plaque, "Plaque_ref_stats")
openxlsx::writeDataTable(wb_plaque, "Plaque_ref_stats",
                         results_plaque$stats)

grp_sim_plaque_df <- as.data.frame(grp_sim_mat_plaque)
grp_sim_plaque_df <- cbind(Group = rownames(grp_sim_plaque_df), grp_sim_plaque_df)
openxlsx::addWorksheet(wb_plaque, "Group_mean_sim_plaque")
openxlsx::writeDataTable(wb_plaque, "Group_mean_sim_plaque", grp_sim_plaque_df)

openxlsx::saveWorkbook(
  wb_plaque,
  file.path(table_dir, "transcriptional_similarity_combined_plaque.xlsx"),
  overwrite = TRUE
)

cat("\nDone — combined plaque analysis complete.\n",
    "  Plots  ->", plot_dir,  "\n",
    "  Tables ->", table_dir, "\n")
