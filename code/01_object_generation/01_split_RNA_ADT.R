# ==============================================================================
# OPTICO-ACS | Split the processed CITE-seq object into RNA and ADT objects
#
# Panels : upstream of all CITE-seq panels
# Input  : data/ObjFull.rds (demultiplexed, QC-filtered CITE-seq object; RNA and DSB-normalised ADT assay)
# Output : data/optico_RNA.rds, data/optico_ADT.rds
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(Azimuth)
library(SeuratData)
library(patchwork)
library(SCpubr)
library(ggplot2)
library(harmony)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v3")
options(Seurat.object.assay.version = "v3")

#load object 
optico <- readRDS(file.path(OPTICO_DIR, "data/ObjFull.rds"))
optico@version
optico <- UpdateSeuratObject(optico)
optico@version

#first glance
DimPlot(optico, reduction = "wnn.umap", group.by = "test")

#Normalize RNA and CITE Assay and compute PCA for WNN 
DefaultAssay(optico) <- 'RNA'
optico <- NormalizeData(optico) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(optico) <- 'CITE'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(optico) <- rownames(optico[["CITE"]])
var(optico@assays$CITE@data %>% as.matrix()%>% t())
optico <-ScaleData(optico) %>% RunPCA()

#umaps withour adjustment 
optico <- RunUMAP(optico, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                  reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
optico <- RunUMAP(optico, reduction = 'apca', dims = 1:18, assay = 'ADT', 
                  reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
p1 <- DimPlot(optico, reduction = 'rna.umap', group.by = 'test', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(optico, reduction = 'adt.umap', group.by = 'test', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
ggsave(file.path(OPTICO_DIR, "plots/Umaps_wo_adjustment.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 28)

#umaps with adjustemnt 
optico <- RunUMAP(optico, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                  reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_', a=2, b=0.8)
optico <- RunUMAP(optico, reduction = 'apca', dims = 1:18, assay = 'ADT', 
                  reduction.name = 'adt.umap', reduction.key = 'adtUMAP_',  a=1, b=0.6)
p3 <- DimPlot(optico, reduction = 'rna.umap', group.by = 'test', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(optico, reduction = 'adt.umap', group.by = 'test', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

p1 + p3
ggsave(file.path(OPTICO_DIR, "plots/RNA_adjustment.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 28)

p2 + p4
ggsave(file.path(OPTICO_DIR, "plots/ADT_adjustment.jpg"), limitsize = F, 
       units = "cm", height = 12, width = 28)

#export metadaa from seurat object to make sure all objects are included 
Idents(optico) <- "test"
samples <- levels(optico)
write.csv2(samples, 
           file.path(OPTICO_DIR, "tables/samples.csv"))

#Adjust wnn 
#determine dimensionality for CITE 
ElbowPlot(optico, ndims = 30, reduction = "apca")
ElbowPlot(optico, ndims = 30, reduction = "pca")

# Cell-specific modality weights can be accessed at bm$RNA.weight
optico <- FindMultiModalNeighbors(
  optico, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

#umap without adjustement - 1:30, 1:18 
optico <- RunUMAP(optico, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
p5 <- DimPlot(optico, reduction = 'wnn.umap', group.by = 'test', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p5

#umap with adjustment 
optico <- RunUMAP(optico, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", a=25, b=0.75)
p6 <- DimPlot(optico, reduction = 'rna.umap', group.by = 'test', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()

#Delete unwanted assays
DefaultAssay(optico) <- "RNA"
optico[['SCT']] <- NULL
optico[['predicted_ADT']] <- NULL
optico[['prediction.score.celltype.l1']] <- NULL
optico[['prediction.score.celltype.l2']] <- NULL
optico[['prediction.score.celltype.l3']] <- NULL

optico <- RenameAssays(optico, assay.name = "CITE", new.assay.name = "ADT")
optico_both_assays <- DietSeurat(optico,assays = c("RNA", "ADT"))
optico_RNA <- optico_both_assays
optico_RNA[['ADT']] <- NULL
optico_ADT <- optico_both_assays
DefaultAssay(optico_ADT) <- "ADT"
optico_ADT[['RNA']] <- NULL
saveRDS(optico_RNA, file.path(OPTICO_DIR, "data/optico_RNA.rds"))
saveRDS(optico_ADT, file.path(OPTICO_DIR, "data/optico_ADT.rds"))
