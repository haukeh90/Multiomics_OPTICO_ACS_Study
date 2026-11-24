# ==============================================================================
# OPTICO-ACS | GeneNMF metaprograms (TCR genes excluded) and per-cell module scores
#
# Panels : Fig. 3b (Jaccard heatmap), Fig. 3c (rank plots of MP3 and MP6)
# Input  : data/optico_final.rds
# Output : plots/jaccard_heatmap.pdf
#          tables/metaprograms_genes_df_R1.xlsx (genes per metaprogram), tables/MP*_df_*_R1.xlsx (fgsea)
#          tables/md_subset_R1.xlsx (UCell score per cell and condition, ranked and plotted in Prism)
# Notes  : multiNMF is not seeded; metaprograms can differ slightly between runs.
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)

optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))
DefaultAssay(optico_final) <- "RNA"

# =============================================================================
# TCR gene exclusion
# TCR V/J/C/D genes should be excluded from DEG
# and module analyses to avoid confounding clonal-expansion signals with
# functional/phenotypic interpretation.
# Pattern covers all four chains (α, β, γ, δ): TRAV, TRAC, TRAJ, TRBV, TRBC,
# TRBJ, TRBD, TRGV, TRGC, TRGJ, TRDV, TRDC, TRDJ, TRDD.
# =============================================================================
TCR_GENE_PATTERN <- "^TR[ABGD][VCJD]"

all_genes     <- rownames(optico_final)
non_tcr_genes <- all_genes[!grepl(TCR_GENE_PATTERN, all_genes)]
n_tcr_removed <- length(all_genes) - length(non_tcr_genes)
message(sprintf("TCR genes removed before NMF: %d (%d remaining)",
                n_tcr_removed, length(non_tcr_genes)))

# TCR exclusion strategy — two options depending on available RAM:
#
# OPTION A (recommended, ~1x peak memory): subset the full object once BEFORE
# splitting. This creates only one intermediate copy rather than N_patients
# simultaneous copies that the lapply approach would require.
#
#optico_no_tcr      <- optico_final[non_tcr_genes, ]  # single subset in-place
#optico_final.list  <- SplitObject(optico_no_tcr, split.by = "patients")
#rm(optico_no_tcr); gc()   # release the intermediate object immediately

#geneNMF.programs <- multiNMF(optico_final.list, assay="RNA", slot="data", k=4:9, L1=c(0,0),
 #                            do_centering=TRUE, nfeatures = 2000)

# OPTION B (zero extra memory, use if Option A still crashes):
# Pass the non-TCR feature whitelist directly to multiNMF via features=.
# multiNMF will then select its HVGs only from within non_tcr_genes.
# Comment out Option A above and uncomment the three lines below:

# The blocklist = TCR genes to EXCLUDE (not the non-TCR genes)
tcr_genes_blocklist <- rownames(optico_final)[grepl(TCR_GENE_PATTERN, rownames(optico_final))]

optico_final.list <- SplitObject(optico_final, split.by = "patients")

rm(optico_final)
geneNMF.programs <- multiNMF(optico_final.list, assay="RNA", slot="data", k=4:9, L1=c(0,0),
                             center = TRUE, nfeatures = 2000,
                             hvg.blocklist = tcr_genes_blocklist)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=6,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.25)
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.5))
class(ph)
ggsave(plot = ph,
       file.path(OPTICO_DIR, "plots/jaccard_heatmap_R1.jpg"),
       unit = "cm", height = 25, width = 30)
ggsave(plot = ph,
       file.path(OPTICO_DIR, "plots/jaccard_heatmap_wide_R1.jpg"),
       unit = "cm", height = 25, width = 50)
ggsave(plot = ph,
       file.path(OPTICO_DIR, "plots/jaccard_heatmap_R1.pdf"),
       unit = "cm", height = 25, width = 30)
geneNMF.metaprograms$metaprograms.metrics

t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))

# Step 1: Find the maximum length
max_length <- max(sapply(geneNMF.metaprograms$metaprograms.genes, length))

# Step 2: Pad each list with NA values
padded_genes <- lapply(geneNMF.metaprograms$metaprograms.genes, function(x) {
  c(x, rep(NA, max_length - length(x)))
})

# Step 3: Convert to a data frame
metaprograms_genes_df <- as.data.frame(padded_genes)

openxlsx::write.xlsx(metaprograms_genes_df,
                     file.path(OPTICO_DIR, "tables/metaprograms_genes_df_R1.xlsx"))


#Interpretation of the results
library(msigdbr)
library(fgsea)
msigdbr::msigdbr_species()

# Use non_tcr_genes as the GSEA universe so that it reflects the actual gene
# space tested (i.e., consistent with what was input to NMF).
# Using the full rownames(optico_final) as universe would be anti-conservative
# for TCR-linked GO terms because those genes were never eligible for selection.
top_p_BP <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe = non_tcr_genes,   # TCR exclusion: universe restricted
          category = c("C5"),
          subcategory = c("BP"),
          species = "Homo sapiens")
})

top_p_KEGG <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe = non_tcr_genes,
          category    = c("C2"),
          subcategory = c("CP:KEGG_LEGACY"),   # updated from "KEGG" for msigdbr v10+
          species     = "Homo sapiens")
})


MP1_df_BP <- as.data.frame(top_p_BP$MP1)
openxlsx::write.xlsx(MP1_df_BP,
                     file.path(OPTICO_DIR, "tables/MP1_df_BP_R1.xlsx"))
MP2_df_BP <- as.data.frame(top_p_BP$MP2)
openxlsx::write.xlsx(MP2_df_BP,
                     file.path(OPTICO_DIR, "tables/MP2_df_BP_R1.xlsx"))
MP3_df_BP <- as.data.frame(top_p_BP$MP3)
openxlsx::write.xlsx(MP3_df_BP,
                     file.path(OPTICO_DIR, "tables/MP3_df_BP_R1.xlsx"))
MP4_df_BP <- as.data.frame(top_p_BP$MP4)
openxlsx::write.xlsx(MP4_df_BP,
                     file.path(OPTICO_DIR, "tables/MP4_df_BP_R1.xlsx"))
MP5_df_BP <- as.data.frame(top_p_BP$MP5)
openxlsx::write.xlsx(MP5_df_BP,
                     file.path(OPTICO_DIR, "tables/MP5_df_BP_R1.xlsx"))
MP6_df_BP <- as.data.frame(top_p_BP$MP6)
openxlsx::write.xlsx(MP6_df_BP,
                     file.path(OPTICO_DIR, "tables/MP6_df_BP_R1.xlsx"))

MP1_df_KEGG <- as.data.frame(top_p_KEGG$MP1)
openxlsx::write.xlsx(MP1_df_KEGG,
                     file.path(OPTICO_DIR, "tables/MP1_df_KEGG_R1.xlsx"))
MP2_df_KEGG <- as.data.frame(top_p_KEGG$MP2)
openxlsx::write.xlsx(MP2_df_KEGG,
                     file.path(OPTICO_DIR, "tables/MP2_df_KEGG_R1.xlsx"))
MP3_df_KEGG <- as.data.frame(top_p_KEGG$MP3)
openxlsx::write.xlsx(MP3_df_KEGG,
                     file.path(OPTICO_DIR, "tables/MP3_df_KEGG_R1.xlsx"))
MP4_df_KEGG <- as.data.frame(top_p_KEGG$MP4)
openxlsx::write.xlsx(MP4_df_KEGG,
                     file.path(OPTICO_DIR, "tables/MP4_df_KEGG_R1.xlsx"))
MP5_df_KEGG <- as.data.frame(top_p_KEGG$MP5)
openxlsx::write.xlsx(MP5_df_KEGG,
                     file.path(OPTICO_DIR, "tables/MP5_df_KEGG_R1.xlsx"))
MP6_df_KEGG <- as.data.frame(top_p_KEGG$MP6)
openxlsx::write.xlsx(MP6_df_KEGG,
                     file.path(OPTICO_DIR, "tables/MP6_df_KEGG_R1.xlsx"))

class(geneNMF.metaprograms)

#Add back to Seurat object
# Note: optico_final still contains all genes (including TCR) — we only excluded
# TCR genes from the NMF input; the UCell scoring step uses the MP gene lists
# derived from that TCR-free NMF and is unaffected.
optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))
mp.genes <- geneNMF.metaprograms$metaprograms.genes
rm(optico_final.list)
optico_final <- AddModuleScore_UCell(optico_final, features = mp.genes, assay="RNA",
                                         ncores=4, name = "")


#Export the mean per sample
library(data.table)
library(dplyr)
optico_final$ID_disease <- paste(optico_final$patients, optico_final$disease, sep = "_")
Idents(optico_final) <- "ID_disease"
md <- optico_final@meta.data %>% as.data.table
md <- md %>%
  select(ID_disease, MP1:MP8)
head(md)
module_means <- md %>%
  group_by(ID_disease) %>%
  summarize(across(starts_with("MP"), ~ mean(.x, na.rm = TRUE)))
openxlsx::write.xlsx(module_means,
                     file.path(OPTICO_DIR, "tables/module_means_R1.xlsx"))


Idents(optico_final) <- "Res2"
sampled_colors <- c("#BC3C2999", "#0072B599", "#FFDC9199", "#20854E99", "#6F99AD99", "#EE4C9799", "#E1872799", "#7876B199")
num_colors <- 20
color_palette <- colorRampPalette(sampled_colors)
created_colors <- color_palette(num_colors)
print(created_colors)
cell_populations <- levels(optico_final)
color_mapping <- as.list(setNames(created_colors, cell_populations))

scCustomize::VlnPlot_scCustom(optico_final, features = "MP1", group.by = "Res2", plot_boxplot = T, colors_use = color_mapping) &
  NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/MP1_vln_R1.pdf"), limitsize = F,
       units = "cm", height = 8, width = 45)
scCustomize::VlnPlot_scCustom(optico_final, features = "MP2", group.by = "Res2", plot_boxplot = T, colors_use = color_mapping) &
  NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/MP2_vln_R1.pdf"), limitsize = F,
       units = "cm", height = 8, width = 45)
scCustomize::VlnPlot_scCustom(optico_final, features = "MP3", group.by = "Res2", plot_boxplot = T, colors_use = color_mapping) &
  NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/MP3_vln_R1.pdf"), limitsize = F,
       units = "cm", height = 8, width = 45)
scCustomize::VlnPlot_scCustom(optico_final, features = "MP4", group.by = "Res2", plot_boxplot = T, colors_use = color_mapping) &
  NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/MP4_vln_R1.pdf"), limitsize = F,
       units = "cm", height = 8, width = 45)
scCustomize::VlnPlot_scCustom(optico_final, features = "MP5", group.by = "Res2", plot_boxplot = T, colors_use = color_mapping) &
  NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/MP5_vln_R1.pdf"), limitsize = F,
       units = "cm", height = 8, width = 45)
scCustomize::VlnPlot_scCustom(optico_final, features = "MP6", group.by = "Res2", plot_boxplot = T, colors_use = color_mapping) &
  NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/MP6_vln_R1.pdf"), limitsize = F,
       units = "cm", height = 8, width = 45)

#save optico object with added Metascores
#saveRDS(optico_final,
 #       file.path(OPTICO_DIR, "data/optico_final.rds"))


#Fig. 3c: export UCell scores per cell and condition -------------------------
library(data.table)
library(dplyr)
library(tidyr)

md <- optico_final@meta.data %>% as.data.table
head(md)
md_subset <- md[, .(disease, MP1, MP2, MP3, MP4, MP5, MP6)]
openxlsx::write.xlsx(md_subset,
                     file.path(OPTICO_DIR, "tables/md_subset_R1.xlsx"))
