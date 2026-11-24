# Multi-omics immunophenotyping of plaque erosion and rupture (OPTICO-ACS)

Analysis code for the main figures of Horstmann, Wolfram et al., *Cardiovascular Research* 2026
(flow cytometry, Olink proteomics, CITE-seq and scTCR-seq of culprit-proximal coronary blood in plaque
erosion vs. rupture).

Sequencing data: GEO [GSE344379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344379).

## Layout

```
code/
  01_object_generation/   builds the main Seurat object (run once, in order)
  Figure_2/ ... Figure_6/ one script per panel or panel group
rerun/                    optional: re-executes the panel-generating code on a copy of the inputs
```

Every script starts with a header (panels, input, output) and reads and writes relative to one project
folder:

```r
Sys.setenv(OPTICO_DIR = "/path/to/project")   # contains data/, tables/, plots/, T_reclustering_2025/, Optico_revisions/
```

Panels were either written directly by the script (R / Python) or drawn in GraphPad Prism from a table the
script exports. Schematics (Fig. 1, 2a, 2b, 4a, 6a), flow cytometry gating, Olink processing and Metascape
enrichment (Fig. 3h, 5h) were done outside R and are not part of this repository.

## Run order

1. `01_object_generation/01_split_RNA_ADT.R`, then `02_integration_WNN_annotation.R` -> `data/optico_final.rds`
2. Figure 3 and 4 scripts (independent of each other, all start from `data/optico_final.rds`)
3. `Figure_5/Fig5a-f_Fig6cd_T_subclustering_TCR.R` -> `T_reclustering_2025_annotated.rds`
4. `Figure_5/Fig5g_oligoclonal_DEG.R` -> `Fig5i_1_pySCENIC.py` -> `Fig5i_2_regulon_lollipop.R`
5. `Figure_6/Fig6_1_CVD_TCR_mapping.R`, then `Fig6_2_CVD_TCR_panels.R` in the same R session

## Figure-to-code map

| Panel | Content | Script | Table or plot written | Final panel |
|---|---|---|---|---|
| 2e | Radar plot of culprit ratios | `Figure_2/Fig2e_radar_culprit_ratios.R` | `radarchart_no_labels.pdf` | R |
| 3a | Plaque similarity index | `Figure_3/Fig3a_plaque_similarity.R` | `transcriptional_similarity_*.xlsx` | Prism |
| 3b | NMF metaprograms, Jaccard heatmap | `Figure_3/Fig3bc_NMF_metaprograms.R` | `jaccard_heatmap.pdf` | R |
| 3c | Rank plots MP3, MP6 | `Figure_3/Fig3bc_NMF_metaprograms.R` | `md_subset_R1.xlsx` | Prism |
| 3d, 3e | WNN UMAP, split by condition | `Figure_3/Fig3de_WNN_UMAP.R` | `R2_umap.pdf`, `Erosion_umap.pdf`, `Rupture_umap.pdf` | R |
| 3f | Cell-type proportions | `Figure_3/Fig3f_cluster_abundance.R` | `cluster_abundance_Optico.xlsx` | Prism |
| 3g | Volcano CD4+ / CD8+ T cells | `Figure_3/Fig3g_pseudobulk_DE.R` | `DE_Genes/Res1/CD4/CD4_markers.csv`, `.../CD8/CD8_markers.csv` | Prism |
| 3h | Pathway enrichment | gene lists from `Fig3g_pseudobulk_DE.R` | Metascape (web) | Prism |
| 4b | huCIRA enrichment heatmap | `Figure_4/Fig4bcde_huCIRA.ipynb` | `Heatmap_robust_landscape.pdf` | Python |
| 4c-e | huCIRA NES per cytokine | `Figure_4/Fig4bcde_huCIRA.ipynb` | `HuCIRA_all_enrichment_results.xlsx` | Prism |
| 4f-h | Stimulation assay | `Figure_4/Fig4fgh_stimulation_assay.R` | `Stats_summary.csv` | Prism |
| 5a | T cell subclustering UMAP | `Figure_5/Fig5a-f_Fig6cd_T_subclustering_TCR.R` | `umap_annotated.pdf` | R |
| 5b | ADT violin plots | same | `VlnPlot_*.pdf` | R |
| 5c | Cluster abundance | same | `cluster_abundance_Optico_T.xlsx` | Prism |
| 5d | Odds ratios per cluster | same | `or_results.xlsx` | Prism |
| 5e | Clone size on UMAP | same | `cloneSize_umap_magma_muted_black.pdf` | R |
| 5f | Clonal occupancy | same | `*_table_clonal_*_annotation.xlsx` | Prism |
| 5g | DEGs in oligoclonal T cells | `Figure_5/Fig5g_oligoclonal_DEG.R` | `DEG_expanded_vs_not.xlsx` | Prism |
| 5h | Pathway enrichment | gene lists from `Fig5g_oligoclonal_DEG.R` | Metascape (web) | Prism |
| 5i | Regulon activity | `Figure_5/Fig5i_1_pySCENIC.py`, `Fig5i_2_regulon_lollipop.R` | `SCENIC_lollipop_*.pdf` | R |
| 6b | Shared clonotypes (circos) | `Figure_6/Fig6_1_...R`, `Fig6_2_CVD_TCR_panels.R` (section 6) | `circos_shared_TR*_Erosion_Rupture.pdf` | R |
| 6c | Clustering by clonal-size distribution | `Figure_5/Fig5a-f_Fig6cd_T_subclustering_TCR.R` | `clonalSizeDistributuion.pdf` | R |
| 6d | Chao1 diversity | same | `clonalDiversity*.xlsx` | Prism |
| 6e | Shared CDR3b heatmap | `Fig6_2_CVD_TCR_panels.R` (section 5B) | `OPTICO_CVD_shared_CDR3b_heatmap_landscape.pdf` | R |
| 6f | Jaccard index per CVD category | `Fig6_2_CVD_TCR_panels.R` (section 4) | `jaccard_matrix_full.csv` | Prism |
| 6g | Frequency-weighted CDR3b overlap | `Fig6_2_CVD_TCR_panels.R` (section 4D) | `per_patient_freq_weighted_overlap.csv` | Prism |

## Sample identifiers

Patients are numbered P1-P16 in this repository (P1-P8 rupture, P9-P16 erosion). Sample names in the
metadata of the objects and the names of the per-patient contig files (`P<n>_contigs.csv`) have to follow
this numbering.

## External inputs

- `data/ObjFull.rds`: demultiplexed, QC-filtered CITE-seq object (Cell Ranger, hashtag demultiplexing, QC and
  DSB normalisation as described in the Methods).
- `external_data/`: reference cohorts for Fig. 3a (`Groningen_only.rds`, `Giannareli_carotid.rds`,
  `CVR_plaques_mapping.rds`; coronary plaque data GSE269269) and `external_data/SCENIC/` (cisTarget hg38 v10
  rankings, motif annotation, TF list) for Fig. 5i.
- `Optico_revisions/data/CVD_TCR_Database.xlsx`: curated CVD-TCR database (Richter et al. 2026).
- `Optico_revisions/data/optico_final.h5ad`: h5ad export of `data/optico_final.rds` for huCIRA.

## Software

R 4.x with Seurat 5, scRepertoire 2.0.7, Libra (`neurorestore/Libra` from GitHub, not the CRAN package of
the same name), GeneNMF, UCell, MAST, edgeR; Python with pySCENIC and hucira.

Notes for rerunning with current package versions:
- `Fig3g_pseudobulk_DE.R`: current Libra versions expect Seurat v5 assays (see the commented conversion
  line in the script) and report `avg_logFC` with the opposite sign; statistics are otherwise the same.
- Fig. 6d: `clonalDiversity()` returns a different table in scRepertoire versions later than 2.0.7.
- GeneNMF (Fig. 3b-c) is not seeded.

## Rerun check

```bash
bash rerun/run_rerun.sh            # edit SRC and RERUN at the top of the script, or pass them as variables
RERUN_DRY=1 bash rerun/run_rerun.sh
```

Copies the inputs into a separate folder, runs the panel-generating sections of the scripts there
(starting from the stored objects) and compares regenerated tables and panels with the originals.
