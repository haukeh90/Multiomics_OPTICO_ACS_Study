# ==============================================================================
# OPTICO-ACS | T cell subclustering, annotation and scTCR-seq repertoire analysis
#
# Panels : Fig. 5a-f, Fig. 6c, 6d
# Input  : data/optico_final.rds
#          data/TCR_raw/Pool*/P.*_contigs.csv (filtered contig annotations per patient)
# Output : T_reclustering_2025/data/T_reclustering_2025_annotated.rds
#          T_reclustering_2025/plots/*.pdf (5a, 5b, 5e, 6c)
#          T_reclustering_2025/tables/*.xlsx (5c, 5d, 5f, 6d; plotted in Prism)
# Notes  : Run with scRepertoire 2.0.7; later versions return a different clonalDiversity() table.
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(Seurat)
library(scCustomize)
library(SeuratData)
library(cowplot)
library(dplyr)
library(ggplot2)
library(reticulate)
library(Signac)
library(UCell)

set.seed(123)


#Fig. 5a: ADT-based subclustering of CD4+ and CD8+ T cells (rLSI integration, Leiden) ----

#load object 
optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))
Idents(optico_final) <- "Res1"
levels(optico_final)
T_cells_optico <- subset(optico_final, idents = c("CD4", "CD8"))

#set default assay 
DefaultAssay(T_cells_optico) <- ("ADT")

# Extract current assay
old_assay <- GetAssay(T_cells_optico, assay = DefaultAssay(T_cells_optico))
data_mat <- GetAssayData(old_assay, slot = "data")
new_rownames <- gsub("[.]", "-", rownames(data_mat))
rownames(data_mat) <- new_rownames
new_assay <- CreateAssayObject(data = data_mat)
T_cells_optico[["ADT_integration"]] <- new_assay

#prepare object for integration 
DefaultAssay(T_cells_optico) <- ("ADT_integration")
features <- c(
  "Hu-CD3-UCHT1", "Hu-CD4-RPA-T4", "Hu-CD8", "Hu-CD5", "Hu-CD7", "Hu-CD27",
  "Hu-CD28", "Hu-CD45RA", "Hu-CD45RO", "Hu-CD127", "Hu-CD25", "Hu-CD69",
  "Hu-CD95", "Hu-CD279", "Hu-CD152", "Hu-CD134", "Hu-CD137", "Hu-CD223",
  "Hu-CD226-11A8", "Hu-CD244", "Hu-CD58", "Hu-TCR-AB", "Hu-TCR-Va7-2",
  "Hu-TCR-Vd2", "Hu-CD62L", "Hu-KLRG1", "Hu-CD57", "Hu-TIGIT", "Hu-CD39",
  "Hu-CD103", "Hu-CD49a", "Hu-CD49d", "Hu-CD2", "Hu-CD45-HI30"
)
features
VariableFeatures(T_cells_optico[["ADT_integration"]]) <- features
T_cells_optico <- ScaleData(T_cells_optico, features = features)

# split the dataset into a list of two seurat objects (stim and CTRL)
ADT.list <- SplitObject(T_cells_optico, split.by = "pools")
ADT.list

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
Features.adt <- SelectIntegrationFeatures(object.list = ADT.list)
ADT.list <- lapply(ADT.list, function(x) RunSVD(x, n = 15))

#integrate 
immune.anchors <- FindIntegrationAnchors(object.list = ADT.list, normalization.method = "LogNormalize",
                                         anchor.features = Features.adt, reduction = "rlsi", dims = 1:15)
immune.combined.adt <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", dims = 1:15)

#Downstream
# original unmodified data still resides in the 'RNA' assay
immune.combined.adt <- RenameAssays(immune.combined.adt, assay.name = "integrated", new.assay.name = "inetgrated.adt")
DefaultAssay(immune.combined.adt) <- "inetgrated.adt"

# Run the standard workflow for visualization and clustering
immune.combined.adt <- ScaleData(immune.combined.adt, verbose = T)
#VariableFeatures(T_cells_optico[["inetgrated.adt"]]) <- features
immune.combined.adt <- RunPCA(immune.combined.adt, verbose = T, reduction.name = "integrated.adt.svd", npcs = 30)
immune.combined.adt <- FindNeighbors(immune.combined.adt, reduction = "integrated.adt.svd")
immune.combined.adt <- FindClusters(immune.combined.adt, resolution = 0.6, cluster.name = "leiden_clusters", algorithm = 4)
immune.combined.adt <- RunUMAP(immune.combined.adt, reduction = 'integrated.adt.svd', dims = 1:30, assay = 'inetgrated.adt', 
                               reduction.name = 'integrated.adt.umap', reduction.key = 'adtUMAP_', a=1, b=1)

DimPlot(immune.combined.adt, label = T)

#siwtch axes
umap_coords <- Embeddings(immune.combined.adt, "integrated.adt.umap")
head(umap_coords)
umap_coords[, "adtUMAP_2"] <- -umap_coords[, "adtUMAP_2"]
head(umap_coords)
immune.combined.adt[["integrated.adt.umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "adtUMAP_", assay = DefaultAssay(immune.combined.adt))
DimPlot(immune.combined.adt, group.by = "seurat_clusters", label = T, reduction = "integrated.adt.umap")

ggsave("plots/adt_umap_all_clusters_T_specific.jpg", units = "cm", height = 14, width = 18)
DimPlot(immune.combined.adt, label = T, group.by = "Res1")
ggsave("plots/adt_umap_T_specific.jpg", units = "cm", height = 14, width = 18)

saveRDS(immune.combined.adt,
        file.path(OPTICO_DIR, "T_reclustering_2025/data/T_reclustering_2025.rds"))


#Cluster markers (ADT) used for annotation -----------------------------------

#load object 
T_reclustering_2025 <- readRDS(file.path(OPTICO_DIR, "T_reclustering_2025/data/T_reclustering_2025.rds"))
DimPlot(T_reclustering_2025, label = T)
rownames(T_reclustering_2025)

#FindAllMarkers for protein 
DefaultAssay(T_reclustering_2025) <- "ADT"
T_reclustering_2025 <- ScaleData(T_reclustering_2025, features = rownames(T_reclustering_2025))
leiden_markers <- FindAllMarkers(T_reclustering_2025, only.pos = T, min.pct = 0.2, logfc.threshold = 0.25)
significant_leiden_markers <- subset(leiden_markers, p_val_adj < 0.05)
significant_leiden_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
duplicated_genes_list <- top10$gene[duplicated(top10$gene)]
print(duplicated_genes_list)
unique_genes_mask <- !duplicated(top10$gene)
unique_top10 <- top10[unique_genes_mask, ]
unique_top10 %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3
DotPlot(T_reclustering_2025, 
        features = top3$gene, 
        group.by = "seurat_clusters", 
        cols = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

openxlsx::write.xlsx(leiden_markers, asTable = T, 
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/leiden_markers.xlsx"))


DotPlot(T_reclustering_2025, features = c("Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD45RA", "Hu.CD45RO",
                                          "Hu.CD62L", "Hu.CD27", "Hu.CD28", "Hu.CD127", "Hu.CD69", "Hu.KLRG1", "Hu.CD103", "Hu.CD95"), cols = "RdYlBu")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plots/adt_markers_for_annotation.pdf", units = "cm", height = 20, width = 32)


#scTCR-seq: load contigs and combine clonotypes (scRepertoire) ---------------

#install
#BiocManager::install("scRepertoire")

library(scRepertoire)
library(dplyr)

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


#scTCR-seq: add condition labels and order samples ---------------------------
#Add Erosion and Rupture Variable 
combined_TCR_Optico <- addVariable(combined_TCR_Optico, 
                                   variable.name = "disease", 
                                   variables = c("Rupture", "Erosion", 
                                                 "Rupture",
                                                 "Erosion", "Erosion", "Rupture",
                                                 "Rupture", "Rupture", "Erosion",
                                                 "Erosion", "Rupture", "Rupture",
                                                 "Erosion", "Rupture", "Erosion",
                                                 "Erosion"))

combined_TCR_Optico <- addVariable(combined_TCR_Optico, 
                                   variable.name = "disease_ID", 
                                   variables = c("P1_Rupture", "P9_Erosion", 
                                                 "P2_Rupture",
                                                 "P11_Erosion", "P10_Erosion", "P3_Rupture",
                                                 "P4_Rupture", "P5_Rupture", "P12_Erosion",
                                                 "P13_Erosion", "P7_Rupture", "P6_Rupture",
                                                 "P15_Erosion", "P8_Rupture", "P16_Erosion",
                                                 "P14_Erosion"))

head(combined_TCR_Optico[[1]])

#reorder so that Rupture and Erosion are seperated 
class(combined_TCR_Optico)
custom_order <- c("P1",
                  "P2",
                  "P3",
                  "P4",
                  "P5",
                  "P6",
                  "P7",
                  "P8",
                  "P9",
                  "P10",
                  "P11",
                  "P12",
                  "P13",
                  "P14",
                  "P15",
                  "P16")

# Create a factor with custom order
print(combined_TCR_Optico)
combined_TCR_Optico<- combined_TCR_Optico[c("P1",
                                            "P2",
                                            "P3",
                                            "P4",
                                            "P5",
                                            "P6",
                                            "P7",
                                            "P8",
                                            "P9",
                                            "P10",
                                            "P11",
                                            "P12",
                                            "P13",
                                            "P14",
                                            "P15",
                                            "P16")]


#Fig. 6c: clustering of patients by clonal-size distribution -----------------

clonalSizeDistribution(combined_TCR_Optico, 
                       cloneCall = "aa", 
                       method= "ward.D2", group.by = "disease_ID", threshold = 1) & NoAxes()

ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/clonalSizeDistributuion.pdf"),
       units = "cm", height = 12, width = 10)


#Fig. 6d: clonal diversity (Chao1; both chains, TRA, TRB) --------------------
# NOTE: run with scRepertoire 2.0.7. In that version clonalDiversity(exportTable = TRUE) returns
# one row per sample with all indices (shannon, inv.simpson, norm.entropy, gini.simpson, chao1,
# ACE); the chao1 column is plotted in Fig. 6d. Later versions (tested: 2.8.0) return a single
# metric in long format, so the tables written below are not reproduced with them.
set.seed(123)

clonalDiversity(combined_TCR_Optico, 
                chain = "both",
                cloneCall = "aa", 
                group.by = "disease_ID",
                x.axis = "disease", 
                n.boots = 100)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/clonal_diversity.jpg"),
       units = "cm", height = 16, width = 28)

clonalDiversity <- clonalDiversity(combined_TCR_Optico, 
                                   chain = "both",
                                   cloneCall = "aa", 
                                   group.by = "disease_ID",
                                   x.axis = "disease", 
                                   n.boots = 100,
                                   exportTable = T)

openxlsx::write.xlsx(clonalDiversity,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/clonalDiversity.xlsx"))

clonalDiversity(combined_TCR_Optico, 
                chain = "TRA",
                cloneCall = "aa", 
                group.by = "disease_ID",
                x.axis = "disease", 
                n.boots = 100)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/clonal_diversity_TRA.jpg"),
       units = "cm", height = 16, width = 28)

clonalDiversity_TRA <- clonalDiversity(combined_TCR_Optico, 
                                       chain = "TRA",
                                       cloneCall = "aa", 
                                       group.by = "disease_ID",
                                       x.axis = "disease", 
                                       n.boots = 100,
                                       exportTable = T)

openxlsx::write.xlsx(clonalDiversity_TRA,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/clonalDiversity_TRA.xlsx"))

clonalDiversity(combined_TCR_Optico, 
                chain = "TRB",
                cloneCall = "aa", 
                group.by = "disease_ID",
                x.axis = "disease", 
                n.boots = 100)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/clonal_diversity_TRB.jpg"),
       units = "cm", height = 16, width = 28)

clonalDiversity_TRB <- clonalDiversity(combined_TCR_Optico, 
                                       chain = "TRB",
                                       cloneCall = "aa", 
                                       group.by = "disease_ID",
                                       x.axis = "disease", 
                                       n.boots = 100,
                                       exportTable = T)

openxlsx::write.xlsx(clonalDiversity_TRB,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/clonalDiversity_TRB.xlsx"))


#Add clonotype information to the T cell object (clone size per cell) --------

colnames(T_reclustering_2025)

T_reclustering_2025 <- combineExpression(combined_TCR_Optico, 
                                            T_reclustering_2025, 
                                            cloneCall="aa", 
                                            group.by = "sample", 
                                            proportion = T)

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

unique(T_reclustering_2025$cloneSize)
Seurat::DimPlot(T_reclustering_2025, group.by = "cloneSize", reduction = "integrated.adt.umap"
                , raster = F, pt.size = 2) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))

Seurat::FeaturePlot(T_reclustering_2025, features = "clonalFrequency", order = T, reduction = "integrated.adt.umap")

clonalOccupy(T_reclustering_2025, 
             x.axis = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/clonalOccupy_all.jpg"),
       units = "cm", height = 16, width = 28)

clonalOccupy(T_reclustering_2025, 
             x.axis = "seurat_clusters",
             facet.by = "disease") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/clonalOccupy_split.jpg"),
       units = "cm", height = 16, width = 48)

clonal_occupy_combine<- clonalOccupy(T_reclustering_2025, 
             x.axis = "seurat_clusters", 
             proportion = F, 
             exportTable = T) 
openxlsx::write.xlsx(clonal_occupy_combine,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/clonal_occupy_combine.xlsx"))
                                    

clonal_occupy_split <- clonalOccupy(T_reclustering_2025, 
             x.axis = "seurat_clusters",
             facet.by = "disease", 
             proportion = F,
             exportTable = T) 

head(clonal_occupy_split)
openxlsx::write.xlsx(clonal_occupy_split,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/clonal_occupy_split.xlsx"))

library(tibble)
library(tidyr)
# Filter for Erosion
erosion_table <- clonal_occupy_split %>%
  filter(disease == "Erosion") %>%
  select(-disease) %>%
  pivot_wider(names_from = cloneSize, values_from = n, values_fill = 0) %>%
  arrange(seurat_clusters) %>%
  column_to_rownames("seurat_clusters")
erosion_table$seurat_clusters <- rownames(erosion_table)

# Filter for Rupture
rupture_table <- clonal_occupy_split %>%
  filter(disease == "Rupture") %>%
  select(-disease) %>%
  pivot_wider(names_from = cloneSize, values_from = n, values_fill = 0) %>%
  arrange(seurat_clusters) %>%
  column_to_rownames("seurat_clusters")
rupture_table$seurat_clusters <- rownames(rupture_table)

head(erosion_table)
openxlsx::write.xlsx(erosion_table, asTable = T,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/erosion_table_clonal_occup.xlsx"))
openxlsx::write.xlsx(rupture_table, asTable = T,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/rupture_table_clonal_oocupy.xlsx"))

#create nice umap
library(Seurat)
library(ggplot2)

# Set identity and define clone sizes
Idents(T_reclustering_2025) <- "cloneSize"
clone_levels <- c(
  "Large (0.01 < X <= 0.1)",
  "Medium (0.001 < X <= 0.01)",
  "Small (1e-04 < X <= 0.001)",
  "Rare (0 < X <= 1e-04)"
)

# Create Rest object by excluding defined clone sizes
Rest <- subset(T_reclustering_2025, idents = clone_levels, invert = TRUE)

# Combine all subsets into one Seurat object
T_reclustering_2025$cloneSize_group <- ifelse(
  Idents(T_reclustering_2025) %in% clone_levels,
  as.character(Idents(T_reclustering_2025)),
  "Rest"
)

# Get UMAP coordinates and metadata
umap_df <- as.data.frame(Embeddings(T_reclustering_2025, "integrated.adt.umap"))
umap_df$cloneSize <- T_reclustering_2025$cloneSize_group

# Desired plotting order (bottom to top)
plot_order <- c("Rest", "Rare (0 < X <= 1e-04)", "Small (1e-04 < X <= 0.001)", 
                "Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)")

# Set factor levels accordingly
umap_df$cloneSize <- factor(umap_df$cloneSize, levels = plot_order)

# Reorder the dataframe so higher-priority clone sizes appear later (on top)
umap_df <- umap_df[order(match(umap_df$cloneSize, plot_order)), ]

# Define color scheme
clone_colors <- c(
  "Rest" = "lightgrey",
  "Rare (0 < X <= 1e-04)" = "#fca5a5",
  "Small (1e-04 < X <= 0.001)" = "#ef4444",
  "Medium (0.001 < X <= 0.01)" = "#b91c1c",
  "Large (0.01 < X <= 0.1)" = "#7f1d1d"
)

# Plot
ggplot(umap_df, aes(x = adtUMAP_1, y = adtUMAP_2)) +
  geom_point(aes(fill = cloneSize), color = "white", shape = 21, size = 3, stroke = 0.3) +
  scale_fill_manual(values = clone_colors) +
  labs(title = "Clone Size Distribution", x = "UMAP 1", y = "UMAP 2", fill = "Clone Size") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) & NoLegend() & NoAxes() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/cloneSize_umap.pdf"),
       units = "cm", height = 12, width = 16)


#Fig. 5a: cluster annotation, removal of excluded clusters, final UMAP -------

#Rename
Idents(T_reclustering_2025) <- "seurat_clusters"
levels(T_reclustering_2025)
T_reclustering_2025 <- RenameIdents(T_reclustering_2025, 
                             `1` = "Naive_CD4",
                             `2` = "CD4_TEM",
                             `3` = "CD4_TCM",
                             `4` = "DPT", 
                             `5` = "CD8_TEMRA",
                             `6` = "exclude",
                             `7` = "CD69_CD8_TEMRA",
                             `8` = "CD8_TCM",
                             `9` = "CD8_TEM",
                             `10` = "Treg",
                             `11` = "exclude",
                             `12` = "Naive_CD8",
                             `13` = "exclude"
)
T_reclustering_2025[["annotation"]] <- Idents(object = T_reclustering_2025)  
DimPlot(T_reclustering_2025, group.by = "annotation", label = T, reduction = "integrated.adt.umap")

T_reclustering_2025 <- subset(T_reclustering_2025, idents = "exclude", invert = T)
DimPlot(T_reclustering_2025, group.by = "annotation", label = T, reduction = "integrated.adt.umap")

T_reclustering_2025 <- RunUMAP(T_reclustering_2025, reduction = 'integrated.adt.svd', dims = 1:30, assay = 'inetgrated.adt', 
                               reduction.name = 'integrated.adt.umap', reduction.key = 'adtUMAP_', a=1, b=1)
DimPlot(T_reclustering_2025, group.by = "annotation", label = T, reduction = "integrated.adt.umap")

# Access the UMAP coordinates
umap_coords <- Embeddings(T_reclustering_2025, reduction = "integrated.adt.umap")
umap_coords[,2] <- -1 * umap_coords[,2]
T_reclustering_2025[["integrated.adt.umap"]]@cell.embeddings <- umap_coords
DimPlot(T_reclustering_2025, group.by = "annotation", label = T, reduction = "integrated.adt.umap")
DimPlot(T_reclustering_2025, group.by = "annotation", label = T, reduction = "integrated.adt.umap", 
        cols = c("#a6cee3",
                 "#1f78b4",
                 "#b2df8a",
                 "#33a02c",
                 "#fb9a99",
                 "#e31a1c",
                 "#fdbf6f",
                 "#ff7f00",
                 "#cab2d6",
                 "#6a3d9a"))


#Fig. 5a, 5e: UMAP plots (annotation, clone size) ----------------------------

#Fig1a
DimPlot(T_reclustering_2025, group.by = "annotation", label = T, reduction = "integrated.adt.umap", 
        cols = c("#a6cee3",
                 "#1f78b4",
                 "#b2df8a",
                 "#33a02c",
                 "#fb9a99",
                 "#e31a1c",
                 "#fdbf6f",
                 "#ff7f00",
                 "#cab2d6",
                 "#6a3d9a"), alpha = 0.7, pt.size = 2) & NoLegend() & NoAxes()
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/umap_annotated_w_legend.jpg"),
       units = "cm", height = 12, width = 16)

# Extract UMAP coordinates and metadata
T_reclustering_2025@reductions$integrated.adt.umap
umap_df <- FetchData(T_reclustering_2025, vars = c("adtUMAP_1", "adtUMAP_2", "annotation"))
ggplot(umap_df, aes(x = adtUMAP_1, y = adtUMAP_2, fill = annotation)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.1, alpha = 0.7) +
  scale_fill_manual(values = c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
    "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"
  )) +
  theme_void() +
  theme(legend.position = "none")
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/umap_annotated.jpg"),
       units = "cm", height = 12, width = 16)

library(colorspace)

# Original palette
pal <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"
)

# Desaturate and slightly lighten
pal_nm <- lighten(desaturate(pal, amount = 0.35), amount = 0.10)

umap_df <- FetchData(
  T_reclustering_2025,
  vars = c("adtUMAP_1", "adtUMAP_2", "annotation")
)

ggplot(umap_df, aes(x = adtUMAP_1, y = adtUMAP_2, fill = annotation)) +
  geom_point(
    shape = 21,
    color = "black",
    size = 3,
    stroke = 0.1,
    alpha = 0.75
  ) +
  scale_fill_manual(values = pal_nm) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(OPTICO_DIR, "T_reclustering_2025/plots/umap_annotated.jpg"),
  units = "cm",
  height = 12,
  width = 16,
  dpi = 600
)
ggsave(
  file.path(OPTICO_DIR, "T_reclustering_2025/plots/umap_annotated.pdf"),
  units = "cm",
  height = 12,
  width = 16)

#Fig1d
library(Seurat)
library(ggplot2)
library(viridisLite)
library(colorspace)

# Extract UMAP coordinates
umap_df <- as.data.frame(Embeddings(T_reclustering_2025, "integrated.adt.umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

# Add clone size group
umap_df$cloneSize <- T_reclustering_2025$cloneSize_group

# Define plotting order
plot_order <- c(
  "Rest",
  "Rare (0 < X <= 1e-04)",
  "Small (1e-04 < X <= 0.001)",
  "Medium (0.001 < X <= 0.01)",
  "Large (0.01 < X <= 0.1)"
)

umap_df$cloneSize <- factor(umap_df$cloneSize, levels = plot_order)

# Plot background first, enriched clones last
umap_df <- umap_df[order(umap_df$cloneSize), ]

# Build a muted magma palette
magma_base <- viridisLite::magma(3, begin = 0.35, end = 0.82)
magma_muted <- colorspace::lighten(
  colorspace::desaturate(magma_base, amount = 0.25),
  amount = 0.4
)

# Final palette: black/charcoal + muted magma
clone_colors_magma_muted <- c(
  "Rest" = "#111111",
  "Rare (0 < X <= 1e-04)" = "#5a5a5a",
  "Small (1e-04 < X <= 0.001)" = magma_muted[1],
  "Medium (0.001 < X <= 0.01)" = magma_muted[2],
  "Large (0.01 < X <= 0.1)" = magma_muted[3]
)

# Optional: inspect colors
clone_colors_magma_muted
# colorspace::show_col(clone_colors_magma_muted)

p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(fill = cloneSize),
    shape = 21,
    color = "black",
    size = 3,
    stroke = 0.08,
    alpha = 0.80
  ) +
  scale_fill_manual(
    values = clone_colors_magma_muted,
    drop = FALSE
  ) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(fill = "Clone size")

p & NoLegend()

ggsave(
  file.path(OPTICO_DIR, "T_reclustering_2025/plots/cloneSize_umap_magma_muted_black.pdf"),
  units = "cm",
  height = 12,
  width = 16
)


#Fig1f
Idents(T_reclustering_2025) <- "disease"
levels(T_reclustering_2025)
Erosion <- subset(T_reclustering_2025, idents = "Erosion")
Rupture <- subset(T_reclustering_2025, idents = "Rupture")
Nebulosa::plot_density(Erosion,
                       features = "clonalFrequency_for_nebulosa",
                       method = "ks",
                       pal = "plasma",
                       size = 5)
Nebulosa::plot_density(Rupture,
                       features = "clonalFrequency_for_nebulosa",
                       method = "ks",
                       pal = "plasma",
                       size = 5)
plot_erosion <- Nebulosa::plot_density(Erosion,
                                       features = "clonalFrequency_for_nebulosa",
                                       method = "ks",
                                       pal = "plasma",
                                       size = 5)

plot_rupture <- Nebulosa::plot_density(Rupture,
                                       features = "clonalFrequency_for_nebulosa",
                                       method = "ks",
                                       pal = "plasma",
                                       size = 5)
plot_erosion_scaled <- plot_erosion + scale_fill_viridis_c(limits = c(0, 0.15), option = "plasma") #geht jetzs beides von 0-0.20
plot_rupture_scaled <- plot_rupture + scale_fill_viridis_c(limits = c(0, 0.15), option = "plasma")
plot_erosion_scaled & NoLegend() & NoAxes() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/nebulosa_erosion.jpg"),
       units = "cm", height = 12, width = 16)
plot_rupture_scaled & NoLegend() & NoAxes() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/nebulosa_rupture.jpg"),
       units = "cm", height = 12, width = 16)
plot_erosion_scaled & NoAxes() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/nebulosa_erosion_w_legend.jpg"),
       units = "cm", height = 12, width = 16)
plot_rupture_scaled & NoAxes() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "T_reclustering_2025/plots/nebulosa_rupture_w_legend.jpg"),
       units = "cm", height = 12, width = 16)


#Fig. 5d: odds ratios per cluster, Rupture vs Erosion ------------------------
library(epitools)

#preview 
table(T_reclustering_2025$annotation, T_reclustering_2025$disease)

# Extract table
cell_table <- table(T_reclustering_2025$annotation, T_reclustering_2025$disease)

# Create an empty dataframe to store results
or_results <- data.frame(
  CellType = rownames(cell_table),
  OR = NA,
  Lower_CI = NA,
  Upper_CI = NA,
  p_value = NA
)

# Loop through each cell type
for (i in 1:nrow(cell_table)) {
  # Get counts
  a <- cell_table[i, "Rupture"]    # Cell type in Rupture
  b <- cell_table[i, "Erosion"]    # Cell type in Erosion
  c <- sum(cell_table[-i, "Rupture"])  # All other cell types in Rupture
  d <- sum(cell_table[-i, "Erosion"])  # All other cell types in Erosion
  
  # 2x2 table
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  # Compute OR and CI
  or <- oddsratio(mat, method = "wald")
  
  # Save to results
  or_results$OR[i] <- or$measure[2, 1]
  or_results$Lower_CI[i] <- or$measure[2, 2]
  or_results$Upper_CI[i] <- or$measure[2, 3]
  or_results$p_value[i] <- or$p.value[2, 3]
}

#preview OR 
or_results <- or_results[order(or_results$OR, decreasing = TRUE), ]
print(or_results)

ggplot(or_results, aes(x = reorder(CellType, OR), y = OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  ylab("Odds Ratio (Rupture vs. Erosion)") +
  xlab("Cell Type") +
  theme_minimal()

head(or_results)
or_results$log2OR <- log2(or_results$OR)
or_results$log2_Lower_CI <- log2(or_results$Lower_CI)
or_results$log2_Upper_CI <- log2(or_results$Upper_CI)
ggplot(or_results, aes(x = reorder(CellType, log2OR), y = log2OR)) +
  geom_point() +
  geom_errorbar(aes(ymin = log2_Lower_CI, ymax = log2_Upper_CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  ylab("log2 Odds Ratio (Rupture vs. Erosion)") +
  xlab("Cell Type") +
  theme_minimal()
openxlsx::write.xlsx(or_results,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/or_results.xlsx"), 
                     asTable = T)


#Fig. 5b: ADT expression per cluster (violin plots) --------------------------
library(viridis)
DefaultAssay(T_reclustering_2025) <- "ADT"
rownames(T_reclustering_2025)
features <- c(
  "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD5", "Hu.CD7", "Hu.CD27",
  "Hu.CD28", "Hu.CD45RA", "Hu.CD45RO", "Hu.CD127", "Hu.CD25", "Hu.CD69",
  "Hu.CD95", "Hu.CD279", "Hu.CD152", "Hu.CD134", "Hu.CD137", "Hu.CD223",
  "Hu.CD226-11A8", "Hu.CD244", "Hu.CD58", "Hu.TCR.AB", "Hu.TCR.Va7.2",
  "Hu.TCR.Vd2", "Hu.CD62L", "Hu.KLRG1", "Hu.CD57", "Hu.TIGIT", "Hu.CD39",
  "Hu.CD103", "Hu.CD49a", "Hu.CD49d", "Hu.CD2", "Hu.CD45-HI30"
)

DotPlot(T_reclustering_2025, features = features, cols = "RdYlBu", group.by = "annotation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  

p <- DotPlot(T_reclustering_2025, 
             features = features, 
             group.by = "annotation") +
  scale_color_distiller(palette = "RdYlBu", direction = -1) +  # Red = high, Blue = low
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank()   # Remove y-axis title
  )
p

rownames(T_reclustering_2025)
features <- c("Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD45RA", "Hu.CD45RO",
              "Hu.CD62L", "Hu.CD28", "Hu.CD69", "Hu.CD279")


# Loop through features
for (feature in features) {
  p <- scCustomize::VlnPlot_scCustom(T_reclustering_2025, features = feature, pt.size = 0, group.by = "annotation", colors_use = pal_nm,
                                     plot_boxplot = T) +
    theme_void() +                        # Remove axes and grid
    theme(legend.position = "none",      # Remove legend
          plot.title = element_blank(),  # Remove title
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # Save as PDF
  ggsave(filename = paste0(file.path(OPTICO_DIR, "T_reclustering_2025/plots/VlnPlot_"), gsub("[:/ ]", "_", feature), ".pdf"),
         plot = p, width = 16, height = 2)
}


#Fig. 5f: clonal occupancy per cluster and condition -------------------------

clonal_occupy_split <- clonalOccupy(T_reclustering_2025, 
                                   x.axis = "annotation",
                                   facet.by = "disease", 
                                   proportion = F,
                                   exportTable = T) 

head(clonal_occupy_split)
openxlsx::write.xlsx(clonal_occupy_split,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/clonal_occupy_annotated_split.xlsx"))

# Filter for Erosion
erosion_table <- clonal_occupy_split %>%
  filter(disease == "Erosion") %>%
  select(-disease) %>%
  pivot_wider(names_from = cloneSize, values_from = n, values_fill = 0) %>%
  arrange(annotation) %>%
  column_to_rownames("annotation")
erosion_table$annotation <- rownames(erosion_table)

# Filter for Rupture
rupture_table <- clonal_occupy_split %>%
  filter(disease == "Rupture") %>%
  select(-disease) %>%
  pivot_wider(names_from = cloneSize, values_from = n, values_fill = 0) %>%
  arrange(annotation) %>%
  column_to_rownames("annotation")
rupture_table$annotation <- rownames(rupture_table)

head(erosion_table)
openxlsx::write.xlsx(erosion_table, asTable = T,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/erosion_table_clonal_occupy_annotation.xlsx"))
openxlsx::write.xlsx(rupture_table, asTable = T,
                     file.path(OPTICO_DIR, "T_reclustering_2025/tables/rupture_table_clonal_oocupy_annotation.xlsx"))


#Save T cell object with clonotype information -------------------------------
saveRDS(T_reclustering_2025,
        file.path(OPTICO_DIR, "T_reclustering_2025/data/T_reclustering_2025.rds"))


#Fig. 5c: cluster abundance per condition ------------------------------------
library(data.table)

#load 
T_reclustering_2025 <- readRDS(file.path(OPTICO_DIR, "T_reclustering_2025/data/T_reclustering_2025.rds"))

#export abundance per Group 
md <- T_reclustering_2025@meta.data %>% as.data.table
cluster_abundance <- md[, .N, by = c("disease", "annotation")] %>% 
  dcast(., disease ~ annotation, value.var = "N")
openxlsx::write.xlsx(cluster_abundance, file.path(OPTICO_DIR, "T_reclustering_2025/tables/cluster_abundance_Optico_T.xlsx"), asTable =F, rowNames =T)


#Save annotated object (input for Fig. 5g-i) ---------------------------------

saveRDS(T_reclustering_2025,
        file.path(OPTICO_DIR, "T_reclustering_2025/data/T_reclustering_2025_annotated.rds"))
