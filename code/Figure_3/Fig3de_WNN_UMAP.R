# ==============================================================================
# OPTICO-ACS | WNN UMAP of all cells, coloured by cell type and split by condition
#
# Panels : Fig. 3d, 3e
# Input  : data/optico_final.rds
# Output : Figures/Figure_2/plots/R2_umap.pdf, Erosion_umap.pdf, Rupture_umap.pdf
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(Seurat)
library(SeuratObject)
library(dplyr)
library(data.table)
library(ggplot2)

#Fig. 3d: WNN UMAP coloured by annotation (Res2) ------------------------------
optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))

setwd(file.path(OPTICO_DIR, "Figures/Figure_2"))
getwd()

Idents(optico_final) <- "Res2"
sampled_colors <- c("#BC3C2999", "#0072B599", "#FFDC9199", "#20854E99", "#6F99AD99", "#EE4C9799", "#E1872799", "#7876B199")
num_colors <- 20
color_palette <- colorRampPalette(sampled_colors)
created_colors <- color_palette(num_colors)
print(created_colors)
cell_populations <- levels(optico_final)
color_mapping <- as.list(setNames(created_colors, cell_populations))

DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'Res2', label = F, pt.size = 0.7,
        repel = TRUE, label.size = 3, cols = color_mapping, order = F) & NoLegend() & NoAxes() & labs(title = NULL)
ggsave("plots/R2_umap.pdf", limitsize = F, 
       units = "cm", height = 14, width = 14)

DimPlot(optico_final, reduction = 'wnn.umap', group.by = 'Res2', label = F, pt.size = 2,
        repel = TRUE, label.size = 3, cols = color_mapping, order = F) & NoLegend() & NoAxes() & labs(title = NULL)
ggsave("plots/R2_umap_big.pdf", limitsize = F, 
       units = "cm", height = 14, width = 14)

#Fig. 3e: WNN UMAP split by condition ----------------------------------------
Idents(optico_final) <- "disease"
Erosion <- subset(optico_final, idents = "Erosion")
DimPlot(Erosion, reduction = 'wnn.umap', group.by = 'disease', label = F, pt.size = 0.7,
        repel = TRUE, label.size = 3, cols = "#904011", order = F) & NoLegend() & NoAxes() & labs(title = NULL)
ggsave("plots/Erosion_umap.pdf", limitsize = F, 
       units = "cm", height = 14, width = 14)

Rupture <- subset(optico_final, idents = "Rupture")
DimPlot(Rupture, reduction = 'wnn.umap', group.by = 'disease', label = F, pt.size = 0.7,
        repel = TRUE, label.size = 3, cols = "#0F1E5E", order = F) & NoLegend() & NoAxes() & labs(title = NULL)
ggsave("plots/Rupture_umap.pdf", limitsize = F, 
       units = "cm", height = 14, width = 14)
