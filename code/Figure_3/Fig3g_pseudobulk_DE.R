# ==============================================================================
# OPTICO-ACS | Pseudobulk differential expression, Erosion vs Rupture (Libra, edgeR-LRT)
#
# Panels : Fig. 3g (volcano plots CD4+ and CD8+ T cells), gene lists for Fig. 3h (Metascape)
# Input  : data/optico_final.rds
# Output : DE_Genes/Res1/<cluster>/<cluster>_markers.csv (plotted as volcano in Prism)
#          DE_Genes/Res*/<cluster>/*_Up_in_Rupture.csv, *_Up_in_Erosion.csv, metascape input (xlsx)
# Notes  : Requires neurorestore/Libra (GitHub), not the CRAN package of the same name.
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

#Pseudobulk DE per cluster (edgeR-LRT, patients as replicates) ----------------
library(Libra)
library(Seurat)
library(dplyr)

#load object 
optico_final <- readRDS(file.path(OPTICO_DIR, "data/optico_final.rds"))
DefaultAssay(optico_final) <- "RNA"
optico_final <- NormalizeData(optico_final, verbose = T)

#for seurat clusters 
Idents(optico_final) <- "seurat_clusters"

#check order of disease 
Idents(optico_final) <- "disease"
levels(optico_final)

# NOTE (added for the repository, not part of the original analysis):
# The object stores the RNA assay in the Seurat v3 "Assay" format and the analysis was run
# with the Libra version available at the time (neurorestore/Libra, 2023/24). Current Libra
# releases expect Seurat v5 assays and stop with
#   'no slot of name "layers" for this object of class "Assay"'.
# To rerun with a current Libra installation, convert the assay first by uncommenting the
# next line. This conversion was not used for the published results and was not tested;
# values may differ slightly from the published tables.
# optico_final[["RNA"]] <- as(optico_final[["RNA"]], "Assay5")
#
# Direction of avg_logFC: in the published tables positive values mean higher expression in
# Rupture (this is what "Up_in_Rupture" / "Up_in_Erosion" below rely on). A rerun with Libra
# 1.0.0 (November 2026) returned the same statistics with the opposite sign (Pearson r of
# avg_logFC = -0.9999 for CD4 and CD8). When using a current Libra version, check the sign
# against the Rupture.exp / Erosion.exp columns and swap the two gene lists if necessary.

DE_EdgR_seurat = run_de(optico_final, replicate_col = "patients", cell_type_col = "seurat_clusters", label_col = "disease",
                        de_family = 'pseudobulk', de_method = 'edgeR', de_type = 'LRT')
head(DE_EdgR_seurat)
write.csv(DE_EdgR_seurat, file = file.path(OPTICO_DIR, "DE_Genes/seurat_clusters/DE_EdgR_seurat.csv"), row.names = TRUE)

DE_EdgR_Res1 = run_de(optico_final, replicate_col = "patients", cell_type_col = "Res1", label_col = "disease",
                      de_family = 'pseudobulk', de_method = 'edgeR', de_type = 'LRT')
head(DE_EdgR_Res1)
write.csv(DE_EdgR_Res1, file = file.path(OPTICO_DIR, "DE_Genes/Res1/DE_EdgR_Res1.csv"), row.names = TRUE)

DE_EdgR_Res2 = run_de(optico_final, replicate_col = "patients", cell_type_col = "Res2", label_col = "disease",
                      de_family = 'pseudobulk', de_method = 'edgeR', de_type = 'LRT')
head(DE_EdgR_Res2)
write.csv(DE_EdgR_Res2, file = file.path(OPTICO_DIR, "DE_Genes/Res2/DE_EdgR_Res2.csv"), row.names = TRUE)

#Res2 (fine annotation): tables per cluster, gene lists for Metascape ---------

#For Loop over folders in Res2
Idents(optico_final) <- "Res2"
levels(optico_final)
folder_names <- c("NK", "CD14_Mono", "CD4_TCM_TEM", "CD8_TCM_TEM", "CD4_Naive", "CD16_Mono",
                  "B_Naive", "Treg", "B_Memory", "MAIT", "gdT", "B_intermediate", "CD8_Naive",
                  "cDC", "B_class_switch", "ASDC", "Dying_T", "Proliferating", "HSPC", "Plasmablast")

folder_names

# Specify the parent directory
parent_directory <- file.path(OPTICO_DIR, "DE_Genes/Res2")

# Create subfolders based on the vector
for (folder_name in folder_names) {
  subfolder_path <- file.path(parent_directory, folder_name)
  dir.create(subfolder_path, recursive = TRUE, showWarnings = T)
  if (file.exists(subfolder_path)) {
    cat("Subfolder", folder_name, "created successfully.\n")
  } else {
    warning("Failed to create subfolder", folder_name)
  }
}

#subset df for each folder 
for (cluster_name in folder_names) {
  cluster_markers <- DE_EdgR_Res2[DE_EdgR_Res2$cell_type == cluster_name, ]
  file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res2/"), cluster_name, "/", cluster_name, "_markers.csv")
  write.csv(cluster_markers, file_path, row.names = TRUE)
}

#load run_de output 
library(readr)
library(openxlsx)
DE_EdgR_Res2 <- read_csv(file.path(OPTICO_DIR, "DE_Genes/Res2/DE_EdgR_Res2.csv"))
DE_EdgR_Res2 <- select(DE_EdgR_Res2, -1)

# Loop through each folder
for (folder_name in folder_names) {
  # Construct the path to the current folder
  current_folder <- file.path(parent_directory, folder_name)
  
  # Create a subset of the data frame based on the current folder
  subset_dataframe <- DE_EdgR_Res2[DE_EdgR_Res2$cell_type == folder_name, ]
  
  # Perform operations on folder_name_markers
  folder_name_markers <- as_tibble(subset_dataframe)
  
  # Check column names and subset
  colnames(folder_name_markers)
  folder_name_markers <- subset(folder_name_markers, !grepl("^MT", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^RPS", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^RPL", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^HLA-", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^H2", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^ATP", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRA", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRB", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRC", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRD", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRG", gene))
  
  # Export Up_in_Rupture
  folder_name_Up_in_Rupture <- subset(folder_name_markers, p_val <= 0.05 & avg_logFC > 0)
  names(folder_name_Up_in_Rupture)[names(folder_name_Up_in_Rupture) == "name"] <- "Up_in_Rupture"
  write.csv(folder_name_Up_in_Rupture, row.names = TRUE, col.names = TRUE,
            file.path(current_folder, paste0(folder_name, "_Up_in_Rupture.csv")))
  
  # Export Up_in_Erosion
  folder_name_Up_in_Erosion <- subset(folder_name_markers, p_val <= 0.05 & avg_logFC < 0)
  names(folder_name_Up_in_Erosion)[names(folder_name_Up_in_Erosion) == "name"] <- "Up_in_Erosion"
  write.csv(folder_name_Up_in_Erosion, row.names = TRUE, col.names = TRUE,
            file.path(current_folder, paste0(folder_name, "_Up_in_Erosion.csv")))
}

#2nd loop for merged excel as metascpae input 
for (folder_name in folder_names) {
  # Construct the path to the current folder
  current_folder <- file.path(parent_directory, folder_name)
  
  # File paths
  Up_in_Rupture_path <- file.path(current_folder, paste0(folder_name, "_Up_in_Rupture.csv"))
  Up_in_Erosion_path <- file.path(current_folder, paste0(folder_name, "_Up_in_Erosion.csv"))
  
  # Read CSV files
  Up_in_Rupture <- read_csv(Up_in_Rupture_path)
  Up_in_Erosion <- read_csv(Up_in_Erosion_path)
  
  # Rename columns
  names(Up_in_Rupture)[names(Up_in_Rupture) == "gene"] <- "Up_in_Rupture"
  names(Up_in_Erosion)[names(Up_in_Erosion) == "gene"] <- "Up_in_Erosion"
  
  # Merge data frames
  mscape <- merge(Up_in_Rupture, Up_in_Erosion, all = TRUE)
  
  # Select columns
  mscape <- select(mscape, Up_in_Rupture, Up_in_Erosion)
  
  # Create the directory if it does not exist
  if (!file.exists(current_folder)) {
    dir.create(current_folder, recursive = TRUE)
  }
  
  # Write to Excel
  excel_path <- file.path(current_folder, paste0(folder_name, "_mscape.xlsx"))
  write.xlsx(mscape, excel_path)
}

#set up function to calculate genes 
calculate_gene_counts <- function(df) {
  # Exclude rows with NAs
  df <- df[complete.cases(df), ]
  
  # Filter genes with p_val < 0.05 and positive avg_log2FC (upregulated)
  upregulated <- df[df$p_val < 0.05 & df$avg_logFC > 0, ]
  
  # Filter genes with p_val < 0.05 and negative avg_log2FC (downregulated)
  downregulated <- df[df$p_val < 0.05 & df$avg_logFC < 0, ]
  
  # Filter genes with p_val > 0.05 (not regulated)
  not_regulated <- df[df$p_val > 0.05, ]
  
  # Calculate the number of 12_months genes in the dataframe
  total_genes <- nrow(df)
  
  # Calculate the number of significantly upregulated genes
  num_upregulated <- nrow(upregulated)
  
  # Calculate the number of significantly downregulated genes
  num_downregulated <- nrow(downregulated)
  
  # Calculate the number of genes that are not regulated
  num_not_regulated <- nrow(not_regulated)
  
  # Return a list with the results
  return(list(
    total_genes = total_genes,
    upregulated = num_upregulated,
    downregulated = num_downregulated,
    not_regulated = num_not_regulated
  ))
}


#check on whole df
head(DE_EdgR_Res2)
count_of_p_val_less_than_0.05 <- sum(DE_EdgR_Res2$p_val < 0.05, na.rm = TRUE)
cat("Count of p_val < 0.05:", count_of_p_val_less_than_0.05, "\n")

#for waterfall_months plot 
for (cluster_name in folder_names) {
  # Read the data for the current cluster
  file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res2/"), cluster_name, "/", cluster_name, "_markers.csv")
  cluster_data <- read_csv(file_path)
  
  # Calculate gene counts for the cluster_data
  cluster_data_counts <- as.data.frame(calculate_gene_counts(cluster_data))
  
  # Set row names for the resulting dataframe
  rownames(cluster_data_counts) <- cluster_name
  
  # Save the dataframe to a CSV file
  output_file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res2/"), cluster_name, "/", cluster_name, "_gene_counts.csv")
  write.csv(cluster_data_counts, output_file_path, row.names = TRUE)
}

#read gene counts
gene_counts_list <- list()

# Loop through each cluster to read the gene count CSV file
for (cluster_name in folder_names) {
  # Construct the file path
  file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res2/"), cluster_name, "/", cluster_name, "_gene_counts.csv")
  
  # Read the gene count CSV file and store it in the list
  gene_counts <- read.csv(file_path, row.names = 1)  # Assuming the row names are stored in the first column
  gene_counts_list[[cluster_name]] <- gene_counts
}

regulated_genes_R2 <- do.call(rbind, gene_counts_list)

openxlsx::write.xlsx(regulated_genes_R2, asTable =T, rowNames =T,
                     file.path(OPTICO_DIR, "DE_Genes/Res2/regulated_genes_R2.xlsx"))


#Fig. 3g: Res1 (CD4, CD8, ...): tables per cluster, gene lists for Metascape --

#For Loop over folders in Res1
Idents(optico_final) <- "Res1"
levels(optico_final)
folder_names <- c("NK", "Mono", "CD4", "CD8", "B", "other_T", "DC", "other")

folder_names

# Specify the parent directory
parent_directory <- file.path(OPTICO_DIR, "DE_Genes/Res1")

# Create subfolders based on the vector
for (folder_name in folder_names) {
  subfolder_path <- file.path(parent_directory, folder_name)
  dir.create(subfolder_path, recursive = TRUE, showWarnings = T)
  if (file.exists(subfolder_path)) {
    cat("Subfolder", folder_name, "created successfully.\n")
  } else {
    warning("Failed to create subfolder", folder_name)
  }
}

#subset df for each folder 
for (cluster_name in folder_names) {
  cluster_markers <- DE_EdgR_Res1[DE_EdgR_Res1$cell_type == cluster_name, ]
  file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res1/"), cluster_name, "/", cluster_name, "_markers.csv")
  write.csv(cluster_markers, file_path, row.names = TRUE)
}

#load run_de output 
library(readr)
library(openxlsx)
DE_EdgR_Res1 <- read_csv(file.path(OPTICO_DIR, "DE_Genes/Res1/DE_EdgR_Res1.csv"))
DE_EdgR_Res1 <- select(DE_EdgR_Res1, -1)

# Loop through each folder
for (folder_name in folder_names) {
  # Construct the path to the current folder
  current_folder <- file.path(parent_directory, folder_name)
  
  # Create a subset of the data frame based on the current folder
  subset_dataframe <- DE_EdgR_Res1[DE_EdgR_Res1$cell_type == folder_name, ]
  
  # Perform operations on folder_name_markers
  folder_name_markers <- as_tibble(subset_dataframe)
  
  # Check column names and subset
  colnames(folder_name_markers)
  folder_name_markers <- subset(folder_name_markers, !grepl("^MT", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^RPS", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^RPL", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^HLA-", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^H2", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^ATP", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRA", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRB", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRC", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRD", gene))
  folder_name_markers <- subset(folder_name_markers, !grepl("^TCRG", gene))
  
  # Export Up_in_Rupture
  folder_name_Up_in_Rupture <- subset(folder_name_markers, p_val <= 0.05 & avg_logFC > 0)
  names(folder_name_Up_in_Rupture)[names(folder_name_Up_in_Rupture) == "name"] <- "Up_in_Rupture"
  write.csv(folder_name_Up_in_Rupture, row.names = TRUE, col.names = TRUE,
            file.path(current_folder, paste0(folder_name, "_Up_in_Rupture.csv")))
  
  # Export Up_in_Erosion
  folder_name_Up_in_Erosion <- subset(folder_name_markers, p_val <= 0.05 & avg_logFC < 0)
  names(folder_name_Up_in_Erosion)[names(folder_name_Up_in_Erosion) == "name"] <- "Up_in_Erosion"
  write.csv(folder_name_Up_in_Erosion, row.names = TRUE, col.names = TRUE,
            file.path(current_folder, paste0(folder_name, "_Up_in_Erosion.csv")))
}

#2nd loop for merged excel as metascpae input 
for (folder_name in folder_names) {
  # Construct the path to the current folder
  current_folder <- file.path(parent_directory, folder_name)
  
  # File paths
  Up_in_Rupture_path <- file.path(current_folder, paste0(folder_name, "_Up_in_Rupture.csv"))
  Up_in_Erosion_path <- file.path(current_folder, paste0(folder_name, "_Up_in_Erosion.csv"))
  
  # Read CSV files
  Up_in_Rupture <- read_csv(Up_in_Rupture_path)
  Up_in_Erosion <- read_csv(Up_in_Erosion_path)
  
  # Rename columns
  names(Up_in_Rupture)[names(Up_in_Rupture) == "gene"] <- "Up_in_Rupture"
  names(Up_in_Erosion)[names(Up_in_Erosion) == "gene"] <- "Up_in_Erosion"
  
  # Merge data frames
  mscape <- merge(Up_in_Rupture, Up_in_Erosion, all = TRUE)
  
  # Select columns
  mscape <- select(mscape, Up_in_Rupture, Up_in_Erosion)
  
  # Create the directory if it does not exist
  if (!file.exists(current_folder)) {
    dir.create(current_folder, recursive = TRUE)
  }
  
  # Write to Excel
  excel_path <- file.path(current_folder, paste0(folder_name, "_mscape.xlsx"))
  write.xlsx(mscape, excel_path)
}

#set up function to calculate genes 
calculate_gene_counts <- function(df) {
  # Exclude rows with NAs
  df <- df[complete.cases(df), ]
  
  # Filter genes with p_val < 0.05 and positive avg_log2FC (upregulated)
  upregulated <- df[df$p_val < 0.05 & df$avg_logFC > 0, ]
  
  # Filter genes with p_val < 0.05 and negative avg_log2FC (downregulated)
  downregulated <- df[df$p_val < 0.05 & df$avg_logFC < 0, ]
  
  # Filter genes with p_val > 0.05 (not regulated)
  not_regulated <- df[df$p_val > 0.05, ]
  
  # Calculate the number of 12_months genes in the dataframe
  total_genes <- nrow(df)
  
  # Calculate the number of significantly upregulated genes
  num_upregulated <- nrow(upregulated)
  
  # Calculate the number of significantly downregulated genes
  num_downregulated <- nrow(downregulated)
  
  # Calculate the number of genes that are not regulated
  num_not_regulated <- nrow(not_regulated)
  
  # Return a list with the results
  return(list(
    total_genes = total_genes,
    upregulated = num_upregulated,
    downregulated = num_downregulated,
    not_regulated = num_not_regulated
  ))
}


#check on whole df
head(DE_EdgR_Res1)
count_of_p_val_less_than_0.05 <- sum(DE_EdgR_Res1$p_val < 0.05, na.rm = TRUE)
cat("Count of p_val < 0.05:", count_of_p_val_less_than_0.05, "\n")

#for waterfall_months plot 
for (cluster_name in folder_names) {
  # Read the data for the current cluster
  file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res1/"), cluster_name, "/", cluster_name, "_markers.csv")
  cluster_data <- read_csv(file_path)
  
  # Calculate gene counts for the cluster_data
  cluster_data_counts <- as.data.frame(calculate_gene_counts(cluster_data))
  
  # Set row names for the resulting dataframe
  rownames(cluster_data_counts) <- cluster_name
  
  # Save the dataframe to a CSV file
  output_file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res1/"), cluster_name, "/", cluster_name, "_gene_counts.csv")
  write.csv(cluster_data_counts, output_file_path, row.names = TRUE)
}

#read gene counts
gene_counts_list <- list()

# Loop through each cluster to read the gene count CSV file
for (cluster_name in folder_names) {
  # Construct the file path
  file_path <- paste0(file.path(OPTICO_DIR, "DE_Genes/Res1/"), cluster_name, "/", cluster_name, "_gene_counts.csv")
  
  # Read the gene count CSV file and store it in the list
  gene_counts <- read.csv(file_path, row.names = 1)  # Assuming the row names are stored in the first column
  gene_counts_list[[cluster_name]] <- gene_counts
}

regulated_genes_R1 <- do.call(rbind, gene_counts_list)

openxlsx::write.xlsx(regulated_genes_R1, asTable =T, rowNames =T,
                     file.path(OPTICO_DIR, "DE_Genes/Res1/regulated_genes_R1.xlsx"))
