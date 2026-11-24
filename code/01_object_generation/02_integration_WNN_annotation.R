# ==============================================================================
# OPTICO-ACS | RPCA integration of RNA and ADT, WNN, clustering and annotation
#
# Panels : upstream of Fig. 3d-h, 4b-e, 5, 6 (builds the main object)
# Input  : data/optico_RNA.rds, data/optico_ADT.rds
# Output : data/optico_final.rds (annotation levels Res1 and Res2)
#          tables/optico.markers_*.csv (cluster markers, RNA and ADT)
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(Azimuth)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v3")

#Integrate RNA assay for variable test------------------------------------------ 
optico_RNA <- readRDS(file.path(OPTICO_DIR, "data/optico_RNA.rds"))
DefaultAssay(optico_RNA) <- "RNA"

#exclude patients P17 amd P18
Idents(optico_RNA) <- "test"
levels(optico_RNA)
optico_RNA <- subset(optico_RNA, idents = c("P17", "P18"), invert = T)

#rename for pools
Idents(optico_RNA) <- "test"
levels(optico_RNA)
optico_RNA <- RenameIdents(optico_RNA, 
                           `P1` = "Pool1", 
                           `P9` = "Pool1",
                           `P2` = "Pool2",
                           `P3.2` = "Pool4",
                           `P10.2` = "Pool4",
                           `P11` = "Pool4",
                           `P12` = "Pool5",
                           `P4` = "Pool5",
                           `P5` = "Pool5",
                           `P6` = "Pool6",
                           `P13` = "Pool6",
                           `P7` = "Pool6",
                           `P14` = "Pool7",
                           `P8` = "Pool7",
                           `P15` = "Pool7",
                           `P16` = "Pool8",
                           `P3` = "PoolP"
)

optico_RNA[["pools"]] <- Idents(object = optico_RNA)
levels(optico_RNA)

# split the dataset into a list of two seurat objects (stim and CTRL)
RNA.list <- SplitObject(optico_RNA, split.by = "pools")

# normalize and identify variable features for each dataset independently
RNA.list <- lapply(X = RNA.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = RNA.list)
RNA.list <- lapply(X = RNA.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})
rm(optico_RNA)

#integrate 
immune.anchors <- FindIntegrationAnchors(object.list = RNA.list, anchor.features = features, reduction = "rpca",
                                         k.anchor = 20, dims = 1:30, verbose = T)
rm(RNA.list)
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", dims = 1:30, verbose = T)

#Downstream
# original unmodified data still resides in the 'RNA' assay
immune.combined <- RenameAssays(immune.combined, assay.name = "integrated", new.assay.name = "inetgrated.rna")
DefaultAssay(immune.combined) <- "inetgrated.rna"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = T)
immune.combined <- RunPCA(immune.combined, verbose = T, reduction.name = "integrated.rna.pca")
immune.combined <- RunUMAP(immune.combined, reduction = 'integrated.rna.pca', dims = 1:30, assay = 'inetgrated.rna', 
                           reduction.name = 'integrated.rna.umap', reduction.key = 'rnaUMAP_', a=2, b=0.8)

#Export Plot for Batch effects 
DimPlot(immune.combined, reduction = 'integrated.rna.umap', group.by = 'pools', label = F, 
        repel = TRUE, label.size = 2.5) 
ggsave(file.path(OPTICO_DIR, "plots/integrated.rna.umap.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)

#save integrated RNA object 
saveRDS(immune.combined, file.path(OPTICO_DIR, "data/immune.combined_RNA.rds"))
rm(immune.anchors, immune.combined, features)

#Integrate ADT assay for variable test----------------------------------------
optico_ADT <- readRDS(file.path(OPTICO_DIR, "data/optico_ADT.rds"))
DefaultAssay(optico_ADT) <- ("ADT")

#exclude patients P17 amd P18
Idents(optico_ADT) <- "test"
levels(optico_ADT)
optico_ADT <- subset(optico_ADT, idents = c("P17", "P18"), invert = T)

Idents(optico_ADT) <- "test"
levels(optico_ADT)
optico_ADT <- RenameIdents(optico_ADT, 
                           `P1` = "Pool1", 
                           `P9` = "Pool1",
                           `P2` = "Pool2",
                           `P3.2` = "Pool4",
                           `P10.2` = "Pool4",
                           `P11` = "Pool4",
                           `P12` = "Pool5",
                           `P4` = "Pool5",
                           `P5` = "Pool5",
                           `P6` = "Pool6",
                           `P13` = "Pool6",
                           `P7` = "Pool6",
                           `P14` = "Pool7",
                           `P8` = "Pool7",
                           `P15` = "Pool7",
                           `P16` = "Pool8",
                           `P3` = "PoolP"
)

optico_ADT[["pools"]] <- Idents(object = optico_ADT)
levels(optico_ADT)

ADT.list <- SplitObject(optico_ADT, split.by = "pools")
ADT.list
Features.adt <- SelectIntegrationFeatures(object.list = ADT.list)
ADT.list <- lapply(X = ADT.list, FUN = RunPCA, features = Features.adt)
immune.anchors <- FindIntegrationAnchors(object.list = ADT.list, normalization.method = "LogNormalize",
                                         anchor.features = Features.adt, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined.adt <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", dims = 1:30)

#Downstream
# original unmodified data still resides in the 'RNA' assay
immune.combined.adt <- RenameAssays(immune.combined.adt, assay.name = "integrated", new.assay.name = "inetgrated.adt")
DefaultAssay(immune.combined.adt) <- "inetgrated.adt"

# Run the standard workflow for visualization and clustering
immune.combined.adt <- ScaleData(immune.combined.adt, verbose = T)
immune.combined.adt <- RunPCA(immune.combined.adt, verbose = T, reduction.name = "integrated.adt.pca")
immune.combined.adt <- RunUMAP(immune.combined.adt, reduction = 'integrated.adt.pca', dims = 1:30, assay = 'inetgrated.adt', 
                               reduction.name = 'integrated.adt.umap', reduction.key = 'adtUMAP_', a=2, b=0.8)

#Export Plot for Batch effects 
DimPlot(immune.combined.adt, reduction = 'integrated.adt.umap', group.by = 'pools', label = F, 
        repel = TRUE, label.size = 2.5) 
ggsave(file.path(OPTICO_DIR, "plots/integrated.adt.umap.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)

#save integrated ADT object 
saveRDS(immune.combined.adt, file.path(OPTICO_DIR, "data/immune.combined.adt.rds"))
rm(immune.anchors, immune.combined, features)


#Merge ADT Assay into RNA object and perform WNN 
immune.combined.adt <- readRDS(file.path(OPTICO_DIR, "data/immune.combined.adt.rds"))
immune.combined_RNA <- readRDS(file.path(OPTICO_DIR, "data/immune.combined_RNA.rds"))

immune.combined_RNA[["inetgrated.adt"]] <- immune.combined.adt[["inetgrated.adt"]]
immune.combined_RNA[["ADT"]] <- immune.combined.adt[["ADT"]]
immune.combined_RNA[["integrated.adt.pca"]] <- immune.combined.adt[["integrated.adt.pca"]]
immune.combined_RNA[["integrated.adt.umap"]] <- immune.combined.adt[["integrated.adt.umap"]]


optico_final <- immune.combined_RNA
rm(immune.combined.adt, immune.combined_RNA)
# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
optico_final <- FindMultiModalNeighbors(
  optico_final, reduction.list = list("integrated.rna.pca", "integrated.adt.pca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

#downstream 
optico_final <- RunUMAP(optico_final, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", a=0.8, b=0.8)
DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'pools', label = F, 
        repel = TRUE, label.size = 2.5)
optico_final <- FindClusters(optico_final, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = T)
DimPlot(optico_final, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
VlnPlot(optico_final, features = "RNA.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 0.1) +
  NoLegend()

#export wnn plots 
DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'seurat_clusters', label = T, 
        repel = TRUE, label.size = 2.5) & NoLegend()
ggsave(file.path(OPTICO_DIR, "plots/wnn_seurat_clusters.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)
DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'pools', label = F, 
        repel = TRUE, label.size = 2.5) 
ggsave(file.path(OPTICO_DIR, "plots/wnn_seurat_pools.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)
DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'seurat_clusters', label = F, 
        repel = TRUE, label.size = 2.5, split.by = "pools", ncol = 4) 
ggsave(file.path(OPTICO_DIR, "plots/wnn_seurat_pools_split.jpg"), limitsize = F, 
       units = "cm", height = 16, width = 28)

#save integrated ADT object 
saveRDS(optico_final, file.path(OPTICO_DIR, "data/optico_final.rds"))


#Calculate DE Genes in RNA and ADT slot for seurat clusters--------------------------------------
library(patchwork)
library(dplyr)
DefaultAssay(optico_final) <- "RNA" 
Idents(optico_final) <- "seurat_clusters"
optico_final <- NormalizeData(optico_final, verbose = T)
all.genes <- rownames(optico_final)
optico_final <- ScaleData(optico_final, features = all.genes, verbose = T)
optico.markers_RNA <- FindAllMarkers(optico_final,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25,
                                     densify = T)
write.csv(optico.markers_RNA, file.path(OPTICO_DIR, "tables/optico.markers_RNA_seurat_clusters.csv"))

optico.markers_RNA %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(optico_final, downsample = 100),
          features = top10$gene) & scale_fill_viridis_c(option = "magma", na.value = "white") & NoLegend()
ggsave(file.path(OPTICO_DIR, "plots/RNA_Heatmap.jpg"), limitsize = F, 
       units = "cm", height = 64, width = 48)

#Calculate DE Genes in RNA and ADT slot for seurat clusters
#check on slots first 
adt_counts <- as.data.frame(optico_final@assays$ADT@counts)
adt_data <- as.data.frame(optico_final@assays$ADT@data)

#now calculate upregulated markers withou threshold
DefaultAssay(optico_final) <- "ADT" 
Idents(optico_final) <- "seurat_clusters"
all.adts <- rownames(optico_final)
optico_final <- ScaleData(optico_final, features = all.adts)

optico.markers_ADT <- FindAllMarkers(optico_final,
                                     slot = "data",
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25,
                                     densify = T)
write.csv(optico.markers_ADT, file.path(OPTICO_DIR, "tables/optico.markers_ADT_seurat_clusters.csv"))

optico.markers_ADT %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(optico_final, downsample = 100),
          features = top10$gene)& scale_fill_viridis_c(option = "magma", na.value = "white") & NoLegend()
ggsave(file.path(OPTICO_DIR, "plots/ADT_Heatmap.jpg"), limitsize = F, 
       units = "cm", height = 64, width = 48)

#Rename Seurat clusters according to Azimuth and ADT Mapping------------------------
Idents(optico_final) <- "seurat_clusters"
levels(optico_final)

#Res1
optico_final <- RenameIdents(optico_final, 
                             `0` = "NK", 
                             `1` = "Mono",
                             `2` = "CD4",
                             `3` = "CD8",
                             `4` = "CD4",
                             `5` = "CD8",
                             `6` = "Mono",
                             `7` = "Mono",
                             `8` = "CD8",
                             `9` = "B",
                             `10` = "CD4",
                             `11` = "CD4",
                             `12` = "B",
                             `13` = "other_T",
                             `14` = "NK",
                             `15` = "other_T",
                             `16` = "B",
                             `17` = "B",
                             `18` = "CD8",
                             `19` = "Mono",
                             `20` = "CD8", 
                             `21` = "DC",
                             `22` = "Mono",
                             `23` = "Mono",
                             `24` = "NK",
                             `25` = "B",
                             `26` = "DC",
                             `27` = "CD8",
                             `28` = "other_T",
                             `29` = "CD4",
                             `30` = "other",
                             `31` = "B",
                             `32` = "NK",
                             `33` = "Mono",
                             `34` = "other",
                             `35` = "CD4",
                             `36` = "B",
                             `37` = "NK"
)

optico_final[["Res1"]] <- Idents(object = optico_final)
levels(optico_final)

sampled_colors <- c("#6F99AD99", "#0072B599", "#FFDC9199", "#20854E99", "#BC3C2999", "#EE4C9799", "#E1872799", "#7876B199")
levels(optico_final)
cols_named_list <- list(
  Mono = "#6F99AD99",
  NK = "#0072B599",
  CD4 = "#FFDC9199",
  CD8 = "#20854E99",
  B = "#BC3C2999",
  other_T = "#EE4C9799",
  DC ="#E1872799",
  other = "#7876B199")


DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'Res1', label = F, pt.size = 0.7,
        repel = TRUE, label.size = 2.5, cols = sampled_colors, order = F) 
ggsave(file.path(OPTICO_DIR, "plots/R1_umap.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)
DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'Res1', label = F, pt.size = 0.7,
        repel = TRUE, label.size = 2.5, cols = sampled_colors, order = F) & NoAxes() & NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/R1_umap_unlabeled.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)

#Res2
Idents(optico_final) <- "seurat_clusters"
levels(optico_final)
optico_final <- RenameIdents(optico_final, 
                             `0` = "NK", 
                             `1` = "CD14_Mono",
                             `2` = "CD4_TCM_TEM",
                             `3` = "CD8_TCM_TEM",
                             `4` = "CD4_Naive",
                             `5` = "CD8_TCM_TEM",
                             `6` = "CD16_Mono",
                             `7` = "CD14_Mono",
                             `8` = "CD8_TCM_TEM",
                             `9` = "B_Naive",
                             `10` = "Treg",
                             `11` = "CD4_TCM_TEM",
                             `12` = "B_Memory",
                             `13` = "MAIT",
                             `14` = "NK",
                             `15` = "gdT",
                             `16` = "B_Memory",
                             `17` = "B_intermediate",
                             `18` = "CD8_Naive",
                             `19` = "CD14_Mono",
                             `20` = "CD8_TCM_TEM", 
                             `21` = "cDC",
                             `22` = "CD14_Mono",
                             `23` = "CD14_Mono",
                             `24` = "NK",
                             `25` = "B_class_switch",
                             `26` = "ASDC",
                             `27` = "CD8_TCM_TEM",
                             `28` = "Dying_T",
                             `29` = "CD4_TCM_TEM",
                             `30` = "Proliferating",
                             `31` = "B_Naive",
                             `32` = "NK",
                             `33` = "CD16_Mono",
                             `34` = "HSPC",
                             `35` = "CD4_Naive",
                             `36` = "Plasmablast",
                             `37` = "NK"
)

optico_final[["Res2"]] <- Idents(object = optico_final)
levels(optico_final)

sampled_colors <- c("#BC3C2999", "#0072B599", "#FFDC9199", "#20854E99", "#6F99AD99", "#EE4C9799", "#E1872799", "#7876B199")

# Anzahl der gewünschten Farben (hier 21 Farben)
num_colors <- 21

# Farbpalette erstellen
color_palette <- colorRampPalette(sampled_colors)

# Farbabstufungen erstellen
created_colors <- color_palette(num_colors)

# Die erstellten Farben anzeigen
print(created_colors)

DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'Res2', label = T, pt.size = 0.7,
        repel = TRUE, label.size = 3, cols = created_colors, order = F) & NoLegend()
ggsave(file.path(OPTICO_DIR, "plots/R2_umap.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 14)
DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'Res2', label = F, pt.size = 0.7,
        repel = TRUE, label.size = 2.5, cols = created_colors, order = F) & NoAxes() & NoLegend() & labs(title = NULL)
ggsave(file.path(OPTICO_DIR, "plots/R2_umap_unlabeled.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 16)

#Calculate DE Genes for R2---------------------------------------------------------
library(patchwork)
library(dplyr)
DefaultAssay(optico_final) <- "RNA" 
Idents(optico_final) <- "Res2"
optico_final <- NormalizeData(optico_final, verbose = T)
all.genes <- rownames(optico_final)
optico_final <- ScaleData(optico_final, features = all.genes)
optico.markers_RNA_R2 <- FindAllMarkers(optico_final,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25,
                                        densify = T)
write.csv(optico.markers_RNA_R2, file.path(OPTICO_DIR, "tables/optico.markers_RNA_R2.csv"))

optico.markers_RNA_R2 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(optico_final, downsample = 50),
          features = top10$gene) & scale_fill_viridis_c(option = "magma", na.value = "white") & NoLegend()
ggsave(file.path(OPTICO_DIR, "plots/RNA_Heatmap_R2.jpg"), limitsize = F, 
       units = "cm", height = 64, width = 48)

#now calculate for ADT with threshold 
DefaultAssay(optico_final) <- "ADT" 
Idents(optico_final) <- "Res2"
all.adts <- rownames(optico_final)
optico_final <- ScaleData(optico_final, features = all.adts)

optico.markers_ADT_R2 <- FindAllMarkers(optico_final,
                                        slot = "data",
                                        only.pos = TRUE,
                                        densify = T, 
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)
write.csv(optico.markers_ADT_R2, file.path(OPTICO_DIR, "tables/optico.markers_ADT_R2.csv"))

optico.markers_ADT_R2 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(optico_final, downsample = 50), group.colors = created_colors,
          features = top10$gene)& scale_fill_viridis_c(option = "magma", na.value = "white") & NoLegend()
ggsave(file.path(OPTICO_DIR, "plots/ADT_Heatmap_R2.jpg"), limitsize = F, 
       units = "cm", height = 64, width = 48)

optico_final <- DietSeurat(optico_final, assays = c("RNA", "ADT"),
                           layers = c("counts", "data"), dimreducs = "wnn.umap")
saveRDS(optico_final, file.path(OPTICO_DIR, "data/optico_final.rds"))

#Rename Idents for Patients ---------------------------------------------------
Idents(optico_final) <- "test"
levels(optico_final)
optico_final <- RenameIdents(optico_final, 
                             `P1` = "P1", 
                             `P9` = "P9",
                             `P2` = "P2",
                             `P3.2` = "P3",
                             `P10.2` = "P10",
                             `P11` = "P11",
                             `P12` = "P12",
                             `P4` = "P4",
                             `P5` = "P5",
                             `P6` = "P6",
                             `P13` = "P13",
                             `P7` = "P7",
                             `P14` = "P14",
                             `P8` = "P8",
                             `P15` = "P15",
                             `P16` = "P16",
                             `P3` = "P3"
)
optico_final[["patients"]] <- Idents(object = optico_final)
levels(optico_final)

#Rename Idents for Disease 
Idents(optico_final) <- "patients"
levels(optico_final)
optico_final <- RenameIdents(optico_final, 
                             `P1` = "Rupture", 
                             `P9` = "Erosion",
                             `P2` = "Rupture",
                             `P3` = "Rupture",
                             `P10` = "Erosion",
                             `P11` = "Erosion",
                             `P12` = "Erosion",
                             `P4` = "Rupture",
                             `P5` = "Rupture",
                             `P6` = "Rupture",
                             `P13` = "Erosion",
                             `P7` = "Rupture",
                             `P14` = "Erosion",
                             `P8` = "Rupture",
                             `P15` = "Erosion",
                             `P16` = "Erosion"
)
optico_final[["disease"]] <- Idents(object = optico_final)
levels(optico_final)

saveRDS(optico_final, file.path(OPTICO_DIR, "data/optico_final.rds"))

#Export QC metrics--------------------------------------------------------------

##export cluster abundance per patient-----------------------------------------20.01.2023
library(data.table)
md <- optico_final@meta.data %>% as.data.table
seurat_abundance <- md[, .N, by = c("patients", "seurat_clusters")] %>% 
  dcast(., patients ~ seurat_clusters, value.var = "N")
write.csv(seurat_abundance,
          file.path(OPTICO_DIR, "tables/seurat_abundance.csv"),
          row.names=T)

R1_abundance <- md[, .N, by = c("patients", "Res1")] %>% 
  dcast(., patients ~ Res1, value.var = "N")
write.csv(R1_abundance,
          file.path(OPTICO_DIR, "tables/R1_abundance.csv"),
          row.names=T)

R2_abundance <- md[, .N, by = c("patients", "Res2")] %>% 
  dcast(., patients ~ Res2, value.var = "N")
write.csv(R2_abundance,
          file.path(OPTICO_DIR, "tables/R2_abundance.csv"),
          row.names=T)
