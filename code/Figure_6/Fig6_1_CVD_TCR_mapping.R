# ==============================================================================
# OPTICO-ACS | Mapping of OPTICO CDR3b sequences to the CVD-TCR database
#
# Panels : upstream of Fig. 6b, 6e, 6f, 6g (run first, same R session as Fig6_2)
# Input  : data/TCR_raw/Pool*/P.*_contigs.csv
#          Optico_revisions/data/CVD_TCR_Database.xlsx (curated CVD-TCR database)
# Output : Optico_revisions/plots/CVD_TCR/tables/OPTICO_TCR_formatted.xlsx, CVD_TCR_plus_OPTICO.xlsx
#          objects used by Fig6_2: combined_TCR_Optico, tcr_raw, df2, outdir, erosion_label, rupture_label
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readxl)
  library(writexl)
  library(tidyr)
  library(purrr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(scales)
  library(immGLIPH)
  library(scRepertoire)
  # scRepertoire must already be loaded (combined_TCR_Optico is in the workspace)
})

# ── 0. Settings ───────────────────────────────────────────────────────────────
outdir      <- file.path(OPTICO_DIR, "Optico_revisions/plots/CVD_TCR")
cvd_db_path <- file.path(OPTICO_DIR, "Optico_revisions/data/CVD_TCR_Database.xlsx")

dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"),  recursive = TRUE, showWarnings = FALSE)

# ── Disease mapping: one row per patient in combined_TCR_Optico ───────────────
# Sample names match the 'samples' argument used in your combineTCR() call.
disease_map <- data.frame(
  sample  = c("P1","P2","P3","P4",
              "P5","P6","P7","P8",
              "P9","P10","P11","P12",
              "P13","P14","P15","P16"),
  disease = c("Rupture","Rupture","Rupture","Rupture",
              "Rupture","Rupture","Rupture","Rupture",
              "Erosion","Erosion","Erosion","Erosion",
              "Erosion","Erosion","Erosion","Erosion"),
  stringsAsFactors = FALSE
)

# Filtration: minimum clone size (match CVD database: expanded clones > 1)
min_clone_size <- 2

# OPTICO category labels (will appear as columns in the heatmaps)
erosion_label  <- "ACS - Erosion"
rupture_label  <- "ACS - Rupture"

set.seed(42)

################################################################################


# 1. Extract and format OPTICO TCR data --------------------------------------
#
# combined_TCR_Optico is the list produced by scRepertoire::combineTCR().
# Each element is one patient's data frame. Expected columns (scRepertoire v1):
#   CTaa      — e.g. "CAVXXXXF_CASSXXXXF"  (TRA_TRB, underscore separator)
#   CTgene    — e.g. "TRAV12-1.TRAJ6_TRBV20-1.TRBJ2-7" (TRA_TRB, underscore)
#   Frequency — integer clone size
# The sample name (patient ID) lives in the list element name.
################################################################################
cat("\n── Step 1: Extracting OPTICO TCR data from combined_TCR_Optico ──\n")

#load data 
P1_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool1/P1_contigs.csv"))
P9_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool1/P9_contigs.csv"))
P2_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool2/P2_contigs.csv"))
P11_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool4/P11_contigs.csv"))
P10_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool4/P10_contigs.csv"))
P3_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool4/P3_contigs.csv"))
P4_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool5/P4_contigs.csv"))
P5_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool5/P5_contigs.csv"))
P12_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool5/P12_contigs.csv"))
P13_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool6/P13_contigs.csv"))
P7_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool6/P7_contigs.csv"))
P6_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool6/P6_contigs.csv"))
P15_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool7/P15_contigs.csv"))
P8_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool7/P8_contigs.csv"))
P16_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool7/P16_contigs.csv"))

P14_contigs <- read.csv(file.path(OPTICO_DIR, "data/TCR_raw/Pool8/P14_contigs.csv"))

contig_list <- list(P1_contigs, P9_contigs, 
                    P2_contigs,
                    P11_contigs, P10_contigs, P3_contigs,
                    P4_contigs, P5_contigs, P12_contigs,
                    P13_contigs, P7_contigs, P6_contigs,
                    P15_contigs, P8_contigs, P16_contigs,
                    P14_contigs)

#combine clones into contigs 
combined_TCR_Optico <- combineTCR(contig_list, 
                                  samples = c("P1", "P9", 
                                              "P2",
                                              "P11", "P10", "P3",
                                              "P4", "P5", "P12",
                                              "P13", "P7", "P6",
                                              "P15", "P8", "P16",
                                              "P14"),
                                  removeNA = TRUE, 
                                  removeMulti = FALSE, 
                                  filterMulti = FALSE)

head(combined_TCR_Optico[[1]])


# ── 1a. Bind list into a single data frame, add sample column ─────────────────
tcr_raw <- bind_rows(
  lapply(names(combined_TCR_Optico), function(s) {
    df <- combined_TCR_Optico[[s]]
    # scRepertoire v1 sometimes already embeds the sample name as a prefix in
    # the barcode, but we add it explicitly to be safe.
    df$sample <- s
    df
  })
)

cat("Total rows in combined_TCR_Optico:", nrow(tcr_raw), "\n")
cat("Samples found:", paste(sort(unique(tcr_raw$sample)), collapse = ", "), "\n")

# ── 1b. Attach disease label via lookup table ─────────────────────────────────
tcr_raw <- tcr_raw %>%
  left_join(disease_map, by = "sample")

missing_map <- tcr_raw$sample[is.na(tcr_raw$disease)]
if (length(missing_map) > 0)
  warning("These samples are not in disease_map and will be dropped: ",
          paste(unique(missing_map), collapse = ", "))

tcr_raw <- tcr_raw %>% filter(!is.na(disease))

# ── 1c. Compute clone frequency (scRepertoire v2 does not pre-compute it) ─────
# Frequency = number of cells sharing the same CTaa within a sample
tcr_raw <- tcr_raw %>%
  group_by(sample, CTaa) %>%
  mutate(Frequency = n()) %>%
  ungroup()

# Filter to expanded clones
tcr_raw <- tcr_raw %>%
  filter(!is.na(CTaa), Frequency >= min_clone_size)

cat("Rows after expansion filter (Frequency >=", min_clone_size, "):", nrow(tcr_raw), "\n")
cat("Disease breakdown:\n")
print(table(tcr_raw$disease))

# ── 1d. CDR3β is already in its own column ────────────────────────────────────
# scRepertoire v2 stores:
#   cdr3_aa2  — CDR3 beta amino acid sequence  (use directly)
#   cdr3_aa1  — CDR3 alpha amino acid sequence
#   TCR2      — beta gene string, e.g. "TRBV12-3..TRBJ2-2.TRBC2"
#   TCR1      — alpha gene string, e.g. "TRAV24.TRAJ53.TRAC"
tcr_raw$CDR3.beta.aa <- tcr_raw$cdr3_aa2
tcr_raw$CDR3.alpha.aa <- tcr_raw$cdr3_aa1

cat("CDR3β non-NA:", sum(!is.na(tcr_raw$CDR3.beta.aa)), "\n")

# ── 1e. Extract TRBV / TRAV from TCR gene columns ────────────────────────────
# TCR2 format: "TRBV12-3..TRBJ2-2.TRBC2"  → TRBV = token before first "."
# TCR1 format: "TRAV24.TRAJ53.TRAC"        → TRAV = token before first "."
tcr_raw$TRBV <- str_extract(tcr_raw$TCR2, "^TRBV[^.]+")
tcr_raw$TRAV <- str_extract(tcr_raw$TCR1, "^TRAV[^.]+")

# ── 1f. Build one row per unique CDR3β × condition ───────────────────────────
optico_tcr <- tcr_raw %>%
  filter(!is.na(CDR3.beta.aa), CDR3.beta.aa != "") %>%
  mutate(
    Category    = ifelse(disease == "Erosion", erosion_label, rupture_label),
    Species     = "Human",
    `PubMed ID` = "OPTICO-ACS",
    Patient     = sample
  ) %>%
  select(CDR3.beta.aa, TRBV, TRAV, Category, Species, `PubMed ID`,
         Patient, Frequency) %>%
  distinct(CDR3.beta.aa, Category, .keep_all = TRUE)   # one row per CDR3β × condition

cat("\nOPTICO CDR3β sequences per condition:\n")
print(table(optico_tcr$Category))

# Save OPTICO TCR table
write_xlsx(optico_tcr, file.path(outdir, "tables", "OPTICO_TCR_formatted.xlsx"))

################################################################################


# 2. Load and harmonise the CVD-TCR database ---------------------------------
################################################################################
cat("\n── Step 2: Loading CVD-TCR database ──\n")

df_cvd <- read_excel(cvd_db_path)

# Filter to human sequences only
df_cvd <- df_cvd[df_cvd$Species == "Human", ]
cat("CVD database (human): ", nrow(df_cvd), "rows\n")

# TRBV harmonization (from Richter et al. code)
harmonize_vgene <- function(v, prefix) {
  v %>%
    str_trim() %>%
    na_if("NA") %>%
    str_replace_all(paste0("^", prefix, "0+(?=\\d)"), prefix) %>%
    str_replace_all("-(0+)(?=\\d)", "-") %>%
    str_replace_all("/0+(?=\\d)", "/")
}

df_cvd <- df_cvd %>%
  mutate(
    TRBV_harmonized = harmonize_vgene(TRBV, "TRBV"),
    TRAV_harmonized = harmonize_vgene(TRAV, "TRAV") %>%
      str_replace_all("/DV0+(?=\\d)", "/DV"),

    TRBV_IMGT_nomenclature_compliant = case_when(
      is.na(TRBV_harmonized) ~ FALSE,
      str_detect(TRBV_harmonized, "/") ~ FALSE,
      str_detect(TRBV_harmonized, "^TRBV\\d+$") ~ FALSE,
      str_detect(TRBV_harmonized, "^TRBV\\d+-\\d+(\\*\\d+|\\*)?$") ~ TRUE,
      TRUE ~ FALSE
    ),
    TRAV_IMGT_nomenclature_compliant = case_when(
      is.na(TRAV_harmonized) ~ FALSE,
      str_detect(TRAV_harmonized, "/(?!DV)") ~ FALSE,
      str_detect(TRAV_harmonized, "^TRAV\\d+$") ~ FALSE,
      str_detect(TRAV_harmonized, "^TRAV\\d+(-\\d+)?(/DV\\d+)?$") ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  distinct()

################################################################################


# 3. Harmonise V genes and merge OPTICO with the database --------------------
################################################################################
cat("\n── Step 3: Merging OPTICO with CVD database ──\n")

optico_tcr <- optico_tcr %>%
  mutate(
    TRBV_harmonized = harmonize_vgene(TRBV, "TRBV"),
    TRAV_harmonized = harmonize_vgene(TRAV, "TRAV") %>%
      str_replace_all("/DV0+(?=\\d)", "/DV"),
    TRBV_IMGT_nomenclature_compliant = case_when(
      is.na(TRBV_harmonized) ~ FALSE,
      str_detect(TRBV_harmonized, "/") ~ FALSE,
      str_detect(TRBV_harmonized, "^TRBV\\d+$") ~ FALSE,
      str_detect(TRBV_harmonized, "^TRBV\\d+-\\d+(\\*\\d+|\\*)?$") ~ TRUE,
      TRUE ~ FALSE
    ),
    TRAV_IMGT_nomenclature_compliant = case_when(
      is.na(TRAV_harmonized) ~ FALSE,
      str_detect(TRAV_harmonized, "/(?!DV)") ~ FALSE,
      str_detect(TRAV_harmonized, "^TRAV\\d+$") ~ FALSE,
      str_detect(TRAV_harmonized, "^TRAV\\d+(-\\d+)?(/DV\\d+)?$") ~ TRUE,
      TRUE ~ FALSE
    )
  )

# Align columns for binding
shared_cols <- intersect(colnames(df_cvd), colnames(optico_tcr))
df2 <- bind_rows(
  df_cvd  %>% select(all_of(shared_cols)),
  optico_tcr %>% select(all_of(shared_cols))
)

cat("Combined dataset:", nrow(df2), "rows\n")
cat("Categories:\n")
print(table(df2$Category))

write_xlsx(df2, file.path(outdir, "tables", "CVD_TCR_plus_OPTICO.xlsx"))

################################################################################
