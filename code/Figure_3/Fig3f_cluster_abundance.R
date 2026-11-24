# ==============================================================================
# OPTICO-ACS | Cell-type proportions per patient
#
# Panels : Fig. 3f (group means shown as donut plots, drawn in Prism)
# Input  : data/optico_final.rds
# Output : Optico_revisions/tables/cluster_abundance_Optico.xlsx (% of cells per patient and Res2 cluster)
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(Seurat)
library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)


#load object 
optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))
table(optico_final$ID_disease, optico_final$disease)

#extract abudnance 
md <- optico_final@meta.data %>% as.data.table
cluster_counts_patient <- md[, .N, by = .(ID_disease, disease, Res2)]

cluster_counts_patient_wide <- dcast(
  cluster_counts_patient,
  ID_disease + disease ~ Res2,
  value.var = "N",
  fill = 0
)

cluster_counts_patient_wide

cluster_props_patient <- md[, .N, by = .(ID_disease, disease, Res2)] %>%
  group_by(ID_disease, disease) %>%
  mutate(
    total_cells = sum(N),
    proportion = N / total_cells,
    percent = proportion * 100
  ) %>%
  ungroup()

cluster_props_patient

cluster_percent_patient_wide <- cluster_props_patient %>%
  select(ID_disease, disease, Res2, percent) %>%
  pivot_wider(
    names_from = Res2,
    values_from = percent,
    values_fill = 0
  )

cluster_percent_patient_wide

cluster_abundance_long <- md[, .N, by = .(ID_disease, disease, Res2)] %>%
  group_by(ID_disease, disease) %>%
  mutate(
    total_cells = sum(N),
    proportion = N / total_cells,
    percent = proportion * 100
  ) %>%
  ungroup()

head(cluster_abundance_long)

#export abundance 
openxlsx::write.xlsx(cluster_percent_patient_wide,
                     file.path(OPTICO_DIR, "Optico_revisions/tables/cluster_abundance_Optico.xlsx"),
                     asTable =T, rowNames =T)
