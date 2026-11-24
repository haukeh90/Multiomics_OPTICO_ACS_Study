# ==============================================================================
# OPTICO-ACS | CVD-TCR database overlap and clonotype sharing
#
# Panels : Fig. 6b, 6e, 6f, 6g (run after Fig6_1 in the same R session)
# Input  : objects from Fig6_1_CVD_TCR_mapping.R
# Output : Optico_revisions/plots/CVD_TCR/plots/circos_shared_TR*_Erosion_Rupture.pdf (6b)
#          Optico_revisions/plots/CVD_TCR/plots/OPTICO_CVD_shared_CDR3b_heatmap_landscape.pdf (6e)
#          Optico_revisions/plots/CVD_TCR/tables/jaccard_matrix_full.csv (6f, plotted in Prism)
#          Optico_revisions/plots/CVD_TCR/tables/per_patient_freq_weighted_overlap.csv (6g, plotted in Prism)
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

suppressPackageStartupMessages({
  library(ggpubr)
  library(rstatix)
  library(patchwork)
})

# ── Desaturated OPTICO colour palette ─────────────────────────────────────────
dis_pal <- c("Rupture" = "#0F1E5E", "Erosion" = "#904011")


################################################################################


# 4. Fig. 6f: Jaccard index between OPTICO repertoires and CVD-TCR categories ----
#
# Instead of replicating the full 9×9 pairwise matrix from the original review,
# we focus on the novel comparison: OPTICO (Erosion / Rupture) vs each CVD
# category.  Horizontal grouped bars ordered by descending Erosion similarity.
# The Erosion–Rupture Jaccard is annotated separately as a caption.
################################################################################
cat("\n── Section 4: Jaccard bar chart ──\n")

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

cdr3b_sets <- df2 %>%
  filter(!is.na(CDR3.beta.aa)) %>%
  group_by(Category) %>%
  summarise(CDR3b_set = list(unique(CDR3.beta.aa)), .groups = "drop")

pathologies <- cdr3b_sets$Category
jaccard_mat <- matrix(
  nrow = length(pathologies), ncol = length(pathologies),
  dimnames = list(pathologies, pathologies)
)
for (i in seq_along(pathologies))
  for (j in seq_along(pathologies))
    jaccard_mat[i, j] <- jaccard(cdr3b_sets$CDR3b_set[[i]], cdr3b_sets$CDR3b_set[[j]])

# Save full matrix as supplementary table
jac_df_full <- as.data.frame(jaccard_mat)
jac_df_full$Category <- rownames(jac_df_full)
write.csv(jac_df_full,
          file.path(outdir, "tables", "jaccard_matrix_full.csv"),
          row.names = FALSE)

# ── Focused comparison: OPTICO vs CVD database ────────────────────────────────
cvd_pathologies <- setdiff(pathologies, c(erosion_label, rupture_label))

optico_jac <- bind_rows(
  tibble(disease = "Erosion",
         cvd_cat = cvd_pathologies,
         jaccard = jaccard_mat[erosion_label, cvd_pathologies]),
  tibble(disease = "Rupture",
         cvd_cat = cvd_pathologies,
         jaccard = jaccard_mat[rupture_label, cvd_pathologies])
)

jac_er <- jaccard_mat[erosion_label, rupture_label]  # internal Erosion–Rupture

# Order by descending Erosion Jaccard
cat_ord_j <- optico_jac %>%
  filter(disease == "Erosion") %>%
  arrange(desc(jaccard)) %>%
  pull(cvd_cat)

optico_jac <- optico_jac %>%
  mutate(
    cvd_cat = factor(cvd_cat, levels = rev(cat_ord_j)),
    disease = factor(disease, levels = c("Erosion", "Rupture"))
  )

p_jac <- ggplot(optico_jac,
                aes(x = jaccard, y = cvd_cat, fill = disease)) +
  geom_col(position = position_dodge(0.68), width = 0.58, color = NA) +
  geom_text(
    aes(label = sprintf("%.4f", jaccard), group = disease),
    position  = position_dodge(0.68),
    hjust     = -0.12,
    size      = 3.2,
    color     = "grey35"
  ) +
  scale_fill_manual(values = dis_pal, name = "OPTICO phenotype") +
  scale_x_continuous(
    name   = "Jaccard similarity index (CDR3β overlap)",
    expand = expansion(mult = c(0, 0.28)),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y    = element_blank(),
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey93", linewidth = 0.4)
  ) +
  labs(
    title   = "CDR3β repertoire overlap: OPTICO vs CVD TCR database",
    subtitle = "Jaccard index  |  Ordered by Erosion overlap",
    caption = sprintf(
      "Internal Erosion–Rupture similarity: Jaccard = %.4f", jac_er
    )
  )

ggsave(file.path(outdir, "plots", "jaccard_barplot.pdf"),
       plot = p_jac, width = 8, height = 5)
cat("Saved: jaccard_barplot.pdf\n")


################################################################################


# 4B. Per-patient CDR3b overlap, unweighted (defines objects used in 4D) -----
#
# Overlap coefficient = unique CDR3β of patient shared with CVD category /
# total unique CDR3β of patient × 100.
# Raw (uncorrected) Wilcoxon rank-sum p-values annotated in each facet:
#   significant (p < 0.05) → * / ** / ***
#   non-significant        → exact p-value (e.g. "p = 0.31")
################################################################################
cat("\n── Section 4B: Per-patient CDR3β overlap boxplots ──\n")

cvd_cats <- setdiff(unique(df2$Category), c(erosion_label, rupture_label))

cvd_cat_sets <- df2 %>%
  filter(!is.na(CDR3.beta.aa), Category %in% cvd_cats) %>%
  group_by(Category) %>%
  summarise(cdr3b_set = list(unique(CDR3.beta.aa)), .groups = "drop")

patient_sets <- tcr_raw %>%
  filter(!is.na(CDR3.beta.aa), CDR3.beta.aa != "") %>%
  group_by(sample, disease) %>%
  summarise(
    pat_set  = list(unique(CDR3.beta.aa)),
    n_unique = n_distinct(CDR3.beta.aa),
    .groups  = "drop"
  )

sim_long <- expand_grid(
  sample  = patient_sets$sample,
  cvd_cat = cvd_cat_sets$Category
) %>%
  left_join(
    patient_sets %>% select(sample, disease, pat_set, n_unique),
    by = "sample"
  ) %>%
  left_join(
    cvd_cat_sets %>% rename(cvd_cat = Category, cvd_set = cdr3b_set),
    by = "cvd_cat"
  ) %>%
  rowwise() %>%
  mutate(
    n_shared    = length(intersect(pat_set, cvd_set)),
    overlap_pct = 100 * n_shared / max(n_unique, 1L)
  ) %>%
  ungroup() %>%
  select(sample, disease, cvd_cat, n_shared, overlap_pct)

write.csv(sim_long,
          file.path(outdir, "tables", "per_patient_overlap_pct.csv"),
          row.names = FALSE)
cat("Saved: per_patient_overlap_pct.csv\n")

# ── Raw Wilcoxon p-values (no multiple-testing correction) ────────────────────
pval_raw <- sim_long %>%
  group_by(cvd_cat) %>%
  rstatix::wilcox_test(overlap_pct ~ disease, paired = FALSE) %>%
  left_join(
    sim_long %>%
      group_by(cvd_cat) %>%
      summarise(y_top = max(overlap_pct, na.rm = TRUE), .groups = "drop"),
    by = "cvd_cat"
  ) %>%
  mutate(
    x     = 1.5,
    y_ann = y_top * 1.2,
    # Stars if significant, exact value if not
    label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ sprintf("p = %.2f", p)
    )
  )

p_box <- ggplot(sim_long,
                aes(x = disease, y = overlap_pct, fill = disease)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.85, color = "grey25"
  ) +
  geom_jitter(
    aes(color = disease),
    width = 0.14, height = 0, size = 2.2, alpha = 0.65, show.legend = FALSE
  ) +
  geom_text(
    data        = pval_raw,
    aes(x = x, y = y_ann, label = label),
    inherit.aes = FALSE,
    size        = 4,
    color       = "grey25"
  ) +
  facet_wrap(
    ~ cvd_cat, scales = "free_y", ncol = 4,
    labeller = labeller(cvd_cat = label_wrap_gen(18))
  ) +
  scale_fill_manual(values  = dis_pal, name = "") +
  scale_color_manual(values = dis_pal, guide = "none") +
  scale_x_discrete(
    labels = c("Erosion" = "Erosion\n(n = 8)", "Rupture" = "Rupture\n(n = 8)")
  ) +
  scale_y_continuous(
    name   = "% patient CDR3β shared with CVD category",
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.25)),
    limits = c(0, NA)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x     = element_blank(),
    legend.position  = "top",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(size = 9.5, face = "bold")
  ) +
  labs(
    title    = "Per-patient CDR3β overlap with CVD TCR database categories",
    subtitle = "Wilcoxon rank-sum (uncorrected)  |  * p<0.05 / ** p<0.01 / *** p<0.001"
  )

ggsave(file.path(outdir, "plots", "per_patient_overlap_boxplot.pdf"),
       plot = p_box, width = 12, height = 8)
cat("Saved: per_patient_overlap_boxplot.pdf\n")


################################################################################


# 4D. Fig. 6g: per-patient frequency-weighted CDR3b overlap ------------------
#
# Instead of binary presence/absence, each shared CDR3β contributes its clone
# frequency (number of cells carrying it) to the overlap score:
#
#   weighted_overlap = Σ(cells with shared CDR3β) / total patient cells × 100
#
# This rewards clonal expansion of cross-reactive sequences, not just their
# presence — a clone expanded to 50 cells counts 50× more than a singleton.
#
# Requires: tcr_raw, cvd_cat_sets (both from Section 4B)
################################################################################
cat("\n── Section 4D: Frequency-weighted CDR3β overlap boxplots ──\n")

# Per patient: cell count per unique CDR3β (sum across all CTaa sharing that CDR3β)
patient_clones_fw <- tcr_raw %>%
  filter(!is.na(CDR3.beta.aa), CDR3.beta.aa != "") %>%
  group_by(sample, disease, CDR3.beta.aa) %>%
  summarise(clone_freq = n(), .groups = "drop")

patient_totals_fw <- patient_clones_fw %>%
  group_by(sample, disease) %>%
  summarise(total_cells = sum(clone_freq), .groups = "drop")

# Unnest CVD sets to a long CDR3β lookup
cvd_cat_lookup_fw <- cvd_cat_sets %>%
  rename(cvd_cat = Category) %>%
  unnest(cdr3b_set) %>%
  rename(CDR3.beta.aa = cdr3b_set)

# Sum cell frequencies for CDR3β shared with each CVD category
shared_freq_fw <- patient_clones_fw %>%
  inner_join(cvd_cat_lookup_fw, by = "CDR3.beta.aa") %>%
  group_by(sample, cvd_cat) %>%
  summarise(shared_cells = sum(clone_freq), .groups = "drop")

# Full grid — zero for unobserved combos
fw_sim_long <- expand_grid(
  sample  = unique(patient_clones_fw$sample),
  cvd_cat = unique(cvd_cat_lookup_fw$cvd_cat)
) %>%
  left_join(patient_totals_fw, by = "sample") %>%
  left_join(shared_freq_fw, by = c("sample", "cvd_cat")) %>%
  mutate(
    shared_cells         = replace_na(shared_cells, 0L),
    weighted_overlap_pct = 100 * shared_cells / pmax(total_cells, 1L)
  )

write.csv(fw_sim_long,
          file.path(outdir, "tables", "per_patient_freq_weighted_overlap.csv"),
          row.names = FALSE)
cat("Saved: per_patient_freq_weighted_overlap.csv\n")

pval_fw <- fw_sim_long %>%
  group_by(cvd_cat) %>%
  rstatix::wilcox_test(weighted_overlap_pct ~ disease, paired = FALSE) %>%
  left_join(
    fw_sim_long %>%
      group_by(cvd_cat) %>%
      summarise(y_top = max(weighted_overlap_pct, na.rm = TRUE), .groups = "drop"),
    by = "cvd_cat"
  ) %>%
  mutate(
    x     = 1.5,
    y_ann = y_top * 1.2,
    label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ sprintf("p = %.2f", p)
    )
  )

p_fw_box <- ggplot(fw_sim_long,
                   aes(x = disease, y = weighted_overlap_pct, fill = disease)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.85, color = "grey25"
  ) +
  geom_jitter(
    aes(color = disease),
    width = 0.14, height = 0, size = 2.2, alpha = 0.65, show.legend = FALSE
  ) +
  geom_text(
    data        = pval_fw,
    aes(x = x, y = y_ann, label = label),
    inherit.aes = FALSE,
    size        = 4,
    color       = "grey25"
  ) +
  facet_wrap(
    ~ cvd_cat, scales = "free_y", ncol = 4,
    labeller = labeller(cvd_cat = label_wrap_gen(18))
  ) +
  scale_fill_manual(values  = dis_pal, name = "") +
  scale_color_manual(values = dis_pal, guide = "none") +
  scale_x_discrete(
    labels = c("Erosion" = "Erosion\n(n = 8)", "Rupture" = "Rupture\n(n = 8)")
  ) +
  scale_y_continuous(
    name   = "% patient cells carrying shared CDR3β",
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.25)),
    limits = c(0, NA)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x     = element_blank(),
    legend.position  = "top",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(size = 9.5, face = "bold")
  ) +
  labs(
    title    = "Frequency-weighted CDR3β overlap with CVD TCR database",
    subtitle = paste0(
      "Wilcoxon rank-sum (uncorrected)  |  * p<0.05 / ** p<0.01 / *** p<0.001\n",
      "Weight = cells carrying shared CDR3β / total patient cells"
    )
  )

ggsave(file.path(outdir, "plots", "freq_weighted_overlap_boxplot.pdf"),
       plot = p_fw_box, width = 12, height = 8)
cat("Saved: freq_weighted_overlap_boxplot.pdf\n")


################################################################################


# 5B. Fig. 6e: binary heatmap of CDR3b sequences shared with the CVD-TCR database ----
#
# Same structure as Section 5 (columns = disease categories, rows = CDR3β)
# but restricted to sequences that are:
#   (a) present in OPTICO (Erosion and/or Rupture), AND
#   (b) shared with ≥1 CVD database category (n_pathologies > 1 across all cats)
#
# This removes the CVD-only sequences that dominate Section 5, keeping the
# focus on OPTICO clonotypes that also appear in the published CVD literature.
#
# Row order: Both → Erosion only → Rupture only, then descending n_pathologies
# Column order: CVD categories → ACS-Erosion → ACS-Rupture (same as Section 5)
################################################################################
cat("\n── Section 5B: OPTICO-focused binary heatmap ──\n")

# ── Identify OPTICO sequences that overlap with ≥1 CVD category ───────────────
optico_cvd_meta <- df2 %>%
  filter(!is.na(CDR3.beta.aa)) %>%
  group_by(CDR3.beta.aa) %>%
  summarise(
    in_erosion   = any(Category == erosion_label),
    in_rupture   = any(Category == rupture_label),
    in_cvd       = any(!Category %in% c(erosion_label, rupture_label)),
    n_total_cats = n_distinct(Category),
    n_cvd_cats   = n_distinct(Category[!Category %in% c(erosion_label, rupture_label)]),
    n_pub_ids    = n_distinct(`PubMed ID`[!Category %in% c(erosion_label, rupture_label)]),
    .groups      = "drop"
  ) %>%
  filter((in_erosion | in_rupture) & in_cvd & n_total_cats > 1) %>%
  mutate(
    optico_group = case_when(
      in_erosion & in_rupture ~ "Both",
      in_erosion              ~ "Erosion only",
      TRUE                    ~ "Rupture only"
    )
  ) %>%
  arrange(
    factor(optico_group, levels = c("Both", "Erosion only", "Rupture only")),
    desc(n_total_cats),
    desc(n_cvd_cats)
  )

cat("OPTICO sequences overlapping CVD database:", nrow(optico_cvd_meta), "\n")
if (nrow(optico_cvd_meta) > 0) print(table(optico_cvd_meta$optico_group))

if (nrow(optico_cvd_meta) == 0) {

  cat("No overlapping sequences found — skipping Section 5B.\n")

} else {

  optico_cvd_seqs <- optico_cvd_meta$CDR3.beta.aa

  # ── Binary presence matrix (same as Section 5 but rows filtered) ─────────
  cdr3_pres5b <- df2 %>%
    filter(!is.na(CDR3.beta.aa), CDR3.beta.aa %in% optico_cvd_seqs) %>%
    distinct(Category, CDR3.beta.aa) %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = Category, values_from = present, values_fill = 0L)

  mat5b <- as.matrix(cdr3_pres5b[, -1])
  rownames(mat5b) <- cdr3_pres5b$CDR3.beta.aa

  # Column order: CVD first, then Erosion, Rupture
  pathology_order5b <- c(
    setdiff(colnames(mat5b), c(erosion_label, rupture_label)),
    erosion_label, rupture_label
  )
  pathology_order5b <- pathology_order5b[pathology_order5b %in% colnames(mat5b)]

  # Reorder rows to match optico_cvd_meta order
  mat5b <- mat5b[match(optico_cvd_seqs, rownames(mat5b)), , drop = FALSE]

  # Pre-apply column order so display indices align in annotations
  mat5b_disp <- mat5b[, match(pathology_order5b, colnames(mat5b)), drop = FALSE]

  # Metadata aligned to row order
  meta_ord5b <- optico_cvd_meta[match(rownames(mat5b_disp), optico_cvd_meta$CDR3.beta.aa), ]

  # ── Colour scheme (reuse Section 5 colours where possible) ────────────────
  n_cvd5b      <- length(pathology_order5b) - 2L
  cvd_cols5b   <- colorRampPalette(brewer.pal(8, "Dark2"))(max(n_cvd5b, 1L))
  optico_col5b <- unname(dis_pal[c("Erosion", "Rupture")])
  path_colors5b <- setNames(c(cvd_cols5b, optico_col5b), pathology_order5b)

  optico_grp_col5b <- c(
    "Both"         = "#5A1F6E",
    "Erosion only" = unname(dis_pal["Erosion"]),
    "Rupture only" = unname(dis_pal["Rupture"])
  )

  # ── Top annotation: pathology colour bar ──────────────────────────────────
  ha_top5b <- HeatmapAnnotation(
    Category = factor(pathology_order5b, levels = pathology_order5b),
    col = list(Category = path_colors5b),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Category = list(title = "Pathology", at = pathology_order5b,
                      labels_gp = gpar(fontsize = 9))
    )
  )

  # ── Right annotation: OPTICO group + N CVD categories + N studies ─────────
  max_cats5b <- max(meta_ord5b$n_total_cats, na.rm = TRUE)
  max_pubs5b <- max(meta_ord5b$n_pub_ids,    na.rm = TRUE)
  col_pub5b  <- colorRamp2(c(2, max(max_cats5b, 3L)), c("#FFEDA0", "#F03B20"))
  col_pap5b  <- colorRamp2(c(0, max(max_pubs5b, 1L)), c("white", "steelblue"))

  ha_right5b <- rowAnnotation(
    "OPTICO group"  = meta_ord5b$optico_group,
    col             = list("OPTICO group" = optico_grp_col5b),
    "# Pathologies" = anno_simple(
      meta_ord5b$n_total_cats,
      col   = col_pub5b,
      pch   = as.character(meta_ord5b$n_total_cats),
      pt_gp = gpar(fontsize = 9)
    ),
    "# Studies (CVD)" = anno_simple(
      meta_ord5b$n_pub_ids,
      col   = col_pap5b,
      pch   = as.character(meta_ord5b$n_pub_ids),
      pt_gp = gpar(fontsize = 9)
    ),
    gap = unit(2, "mm")
  )

  # Row label colours: match OPTICO group
  row_lbl5b <- unname(optico_grp_col5b[meta_ord5b$optico_group])

  # ── Draw ──────────────────────────────────────────────────────────────────
  n_rows5b <- nrow(mat5b_disp)
  n_cols5b <- ncol(mat5b_disp)

  tiff(file.path(outdir, "plots", "OPTICO_CVD_shared_CDR3b_heatmap.tiff"),
       units  = "in",
       width  = max(8, n_cols5b * 0.5 + 5),
       height = max(4, n_rows5b * 0.45 + 3),
       res    = 200)

  draw(Heatmap(
    mat5b_disp,
    name              = "Present",
    col               = c("0" = "white", "1" = "grey30"),
    width             = unit(n_cols5b * 0.5, "cm"),
    height            = unit(n_rows5b * 0.5, "cm"),
    rect_gp           = gpar(col = "grey90", lwd = 1),
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_column_names = FALSE,
    top_annotation    = ha_top5b,
    show_row_names    = TRUE,
    row_names_gp      = gpar(fontsize = 8, col = row_lbl5b),
    right_annotation  = ha_right5b,
    heatmap_legend_param = list(
      title     = "Present",
      at        = c(0, 1),
      labels    = c("absent", "present"),
      labels_gp = gpar(fontsize = 9)
    ),
    column_title = paste0(
      "OPTICO CDR3β shared with CVD database (n = ", n_rows5b, ")\n",
      "Label colour: purple = Both | sienna = Erosion only | navy = Rupture only"
    )
  ))

  dev.off()
  cat("Saved: OPTICO_CVD_shared_CDR3b_heatmap.tiff\n")

  write.csv(
    meta_ord5b %>%
      select(CDR3.beta.aa, optico_group, n_total_cats, n_cvd_cats, n_pub_ids),
    file.path(outdir, "tables", "OPTICO_CVD_shared_CDR3b.csv"),
    row.names = FALSE
  )
  cat("Saved: OPTICO_CVD_shared_CDR3b.csv\n")

  # ── 5B — Landscape (transposed) version ──────────────────────────────────
  # Transpose: rows = disease categories, columns = CDR3β sequences
  # Annotations rotate accordingly:
  #   • left row annotation  = pathology colour bar
  #   • right row annotation = disease category names
  #   • bottom col annotation = OPTICO group + # Pathologies + # Studies
  # ─────────────────────────────────────────────────────────────────────────
  cat("Saving landscape version...\n")

  mat5b_land <- t(mat5b_disp)   # [categories × CDR3β]

  n_r_land <- nrow(mat5b_land)  # number of disease categories
  n_c_land <- ncol(mat5b_land)  # number of CDR3β (19)

  # Left annotation — pathology colour bar (rows = categories)
  ha_left5b_land <- rowAnnotation(
    Pathology = factor(pathology_order5b, levels = pathology_order5b),
    col       = list(Pathology = path_colors5b),
    show_annotation_name = FALSE,
    width = unit(4, "mm")
  )

  # Right annotation — category name labels
  ha_right5b_land <- rowAnnotation(
    show_annotation_name = FALSE,
    annotation_name_side = "top",
    "Disease category" = anno_text(
      pathology_order5b,
      gp = gpar(fontsize = 9),
      just = "left"
    )
  )

  # Bottom annotation — OPTICO group + per-CDR3β stats (columns = CDR3β)
  ha_bottom5b_land <- HeatmapAnnotation(
    "OPTICO group"    = meta_ord5b$optico_group,
    col               = list("OPTICO group" = optico_grp_col5b),
    "# Pathologies"   = anno_simple(
      meta_ord5b$n_total_cats,
      col   = col_pub5b,
      pch   = as.character(meta_ord5b$n_total_cats),
      pt_gp = gpar(fontsize = 8)
    ),
    "# Studies (CVD)" = anno_simple(
      meta_ord5b$n_pub_ids,
      col   = col_pap5b,
      pch   = as.character(meta_ord5b$n_pub_ids),
      pt_gp = gpar(fontsize = 8)
    ),
    gap = unit(2, "mm"),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 8)
  )

  # Column (CDR3β) label colours matching OPTICO group
  col_lbl5b_land <- unname(optico_grp_col5b[meta_ord5b$optico_group])

  pdf(
    file.path(outdir, "plots", "OPTICO_CVD_shared_CDR3b_heatmap_landscape.pdf"),
    width  = max(10, n_c_land * 0.55 + 5),
    height = max(5,  n_r_land * 0.6  + 4)
  )

  draw(
    Heatmap(
      mat5b_land,
      name              = "Present",
      col               = c("0" = "white", "1" = "grey30"),
      width             = unit(n_c_land * 0.55, "cm"),
      height            = unit(n_r_land * 0.6,  "cm"),
      rect_gp           = gpar(col = "grey90", lwd = 3),
      cluster_rows      = FALSE,
      cluster_columns   = FALSE,
      show_row_names    = FALSE,          # shown via right annotation
      show_column_names = TRUE,
      column_names_rot  = 90,
      column_names_gp   = gpar(fontsize = 8, col = col_lbl5b_land),
      column_names_side = "bottom",
      left_annotation   = ha_left5b_land,
      right_annotation  = ha_right5b_land,
      bottom_annotation = ha_bottom5b_land,
      heatmap_legend_param = list(
        title     = "Present",
        at        = c(0, 1),
        labels    = c("absent", "present"),
        labels_gp = gpar(fontsize = 9)
      ),
      column_title = paste0(
        "OPTICO CDR3β shared with CVD database (n = ", n_c_land, ")   ",
        "Label colour: purple = Both | sienna = Erosion only | navy = Rupture only"
      ),
      column_title_gp = gpar(fontsize = 10)
    ),
    padding = unit(c(6, 2, 2, 2), "mm")
  )

  dev.off()
  cat("Saved: OPTICO_CVD_shared_CDR3b_heatmap_landscape.pdf\n")
}

################################################################################


# 6. Fig. 6b: circos plots of clonotypes shared between patients -------------
#
# Adapted from the CAD/NoCAD circos function.
#
# Key changes vs. original:
#   • No Seurat object required — disease group comes from disease_map
#   • CAD/NoCAD → Erosion/Rupture, colours from dis_pal
#   • min_count_any = 2 matches min_clone_size = 2 in the main pipeline
#     (i.e. only expanded clones — consistent with tcr_raw filtering)
#   • combined_TCR_Optico must be in the environment (created in Section 1
#     of OPTICO_CVD_TCR_pipeline.R)
#
# Outputs:
#   circos_shared_TRA_Erosion_Rupture.pdf
#   circos_shared_TRB_Erosion_Rupture.pdf
#   circos_shared_BOTH_CTaa_Erosion_Rupture.pdf
################################################################################
cat("\n── Section 6: Circos plots ──\n")

make_circos_optico <- function(
    combined_TCR,           # named list from combineTCR() — combined_TCR_Optico
    disease_map,            # data.frame: columns 'sample' and 'disease'
    chain_mode    = c("TRA", "TRB", "BOTH"),
    min_count_any = 2,      # ≥2 cells in ≥1 patient → aligns with min_clone_size
    top_n         = 60,
    out_pdf
) {
  chain_mode <- match.arg(chain_mode)

  clonotype_col <- switch(
    chain_mode,
    "TRA"  = "cdr3_aa1",
    "TRB"  = "cdr3_aa2",
    "BOTH" = "CTaa"
  )

  patients <- names(combined_TCR)

  # OPTICO desaturated palette (mirrors dis_pal)
  erosion_col <- "#904011"
  rupture_col <- "#0F1E5E"
  mixed_col   <- "grey75"

  # Named vector: sample → disease ("Erosion" / "Rupture")
  pat_status <- disease_map %>%
    mutate(sample = as.character(sample), disease = as.character(disease)) %>%
    tibble::deframe()

  # Erosion patients first, then Rupture; within each group alphabetical
  patients_ord <- patients[order(ifelse(pat_status[patients] == "Erosion", 0, 1),
                                 patients)]
  sector_cols <- ifelse(pat_status[patients_ord] == "Erosion",
                        erosion_col, rupture_col)
  names(sector_cols) <- patients_ord

  # ── Build patient × clonotype count table ─────────────────────────────────
  create_patient_clone_df <- function(pat) {
    df <- combined_TCR[[pat]]
    if (!(clonotype_col %in% colnames(df)))
      stop(paste("Missing column", clonotype_col, "in patient", pat))
    df %>%
      transmute(cdr3_aa = as.character(.data[[clonotype_col]])) %>%
      filter(!is.na(cdr3_aa), cdr3_aa != "") %>%
      count(cdr3_aa, name = paste0("replicate_count_", pat))
  }

  merged_pat <- map(patients, create_patient_clone_df) %>%
    reduce(full_join, by = "cdr3_aa")
  merged_pat[is.na(merged_pat)] <- 0

  # Expansion filter — consistent with min_clone_size in main pipeline
  merged_pat <- merged_pat %>%
    filter(if_any(starts_with("replicate_count_"), ~ .x >= min_count_any))

  cat(chain_mode, "— clonotypes after expansion filter (>=", min_count_any, "):",
      nrow(merged_pat), "\n")

  arc_sizes <- colSums(
    merged_pat[, paste0("replicate_count_", patients), drop = FALSE] > 0
  )
  names(arc_sizes) <- patients

  # ── Draw ──────────────────────────────────────────────────────────────────
  pdf(out_pdf, width = 12, height = 12)
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree   = 2,
    cell.padding = c(0.001, 0, 0.001, 0),
    track.margin = c(0.01, 0.01)
  )

  circos.initialize(
    factors = patients_ord,
    xlim    = cbind(rep(0, length(patients_ord)), arc_sizes[patients_ord])
  )

  # Outer track — patient labels
  circos.trackPlotRegion(
    ylim         = c(0, 2),
    bg.border    = NA,
    track.height = 0.08,
    track.margin = c(0.08, 0.25),
    panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim   <- get.cell.meta.data("xlim")
      grp    <- pat_status[sector]
      n_cl   <- arc_sizes[sector]
      circos.text(
        x = mean(xlim), y = 2.5,
        labels     = paste0(sector, "\n(", grp, ")\nn=", n_cl),
        facing     = "outside",
        niceFacing = TRUE,
        cex        = 0.6
      )
    }
  )

  # Outer arcs — coloured by disease group
  for (s in patients_ord)
    circos.rect(0, 0, arc_sizes[s], 1,
                sector.index = s,
                col = sector_cols[s], border = NA)

  # Sector coordinate lookup
  sector_ranges <- get.all.sector.index() %>%
    setNames(., .) %>%
    lapply(function(s) get.cell.meta.data("xlim", sector.index = s)) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  colnames(sector_ranges) <- c("xmin", "xmax")
  sector_ranges$sample    <- rownames(sector_ranges)

  # Inner blocks + link coordinates
  link_coords <- data.frame()

  for (s in patients_ord) {
    clone_counts <- merged_pat %>%
      select(cdr3_aa, !!sym(paste0("replicate_count_", s))) %>%
      rename(count = 2) %>%
      filter(count > 0) %>%
      arrange(desc(count)) %>%
      slice_head(n = top_n)

    total <- sum(clone_counts$count)
    if (total == 0) next

    width_scale <- sector_ranges[s, "xmax"] - sector_ranges[s, "xmin"]
    clone_counts <- clone_counts %>%
      mutate(
        width  = count / total * width_scale,
        xstart = cumsum(lag(width, default = 0)),
        xend   = xstart + width,
        sample = s
      )

    for (i in seq_len(nrow(clone_counts)))
      circos.rect(
        xleft = clone_counts$xstart[i], xright = clone_counts$xend[i],
        ybottom = -2, ytop = -0.7,
        sector.index = s,
        col = "white", border = "grey", lwd = 0.2
      )

    link_coords <- bind_rows(link_coords, clone_counts)
  }

  # Clonotypes shared across ≥2 patients
  shared <- link_coords %>%
    select(cdr3_aa, sample, xstart, xend) %>%
    group_by(cdr3_aa) %>%
    filter(n() > 1) %>%
    ungroup()

  pair_cols <- c(
    "Erosion_Erosion" = erosion_col,
    "Rupture_Rupture" = rupture_col,
    "Mixed"           = mixed_col
  )

  all_links <- list()
  for (cl in unique(shared$cdr3_aa)) {
    entries <- shared %>% filter(cdr3_aa == cl)
    if (nrow(entries) < 2) next
    combs <- combn(seq_len(nrow(entries)), 2)
    for (k in seq_len(ncol(combs))) {
      i1 <- combs[1, k]; i2 <- combs[2, k]
      s1 <- entries$sample[i1]; s2 <- entries$sample[i2]
      g1 <- pat_status[s1];     g2 <- pat_status[s2]
      color_group <- if      (g1 == "Erosion" & g2 == "Erosion") "Erosion_Erosion"
                     else if (g1 == "Rupture" & g2 == "Rupture") "Rupture_Rupture"
                     else "Mixed"
      all_links[[length(all_links) + 1]] <- data.frame(
        s1 = s1, s2 = s2,
        xstart1 = entries$xstart[i1], xend1 = entries$xend[i1],
        xstart2 = entries$xstart[i2], xend2 = entries$xend[i2],
        color_group = color_group, stringsAsFactors = FALSE
      )
    }
  }

  all_links_df <- bind_rows(all_links)

  if (nrow(all_links_df) > 0) {
    for (i in seq_len(nrow(all_links_df))) {
      is_mixed     <- all_links_df$color_group[i] == "Mixed"
      col_val      <- pair_cols[all_links_df$color_group[i]]
      link_alpha   <- ifelse(is_mixed, 0.12, 0.45)
      border_alpha <- ifelse(is_mixed, 0.08, 0.45)
      link_lwd     <- ifelse(is_mixed, 0.05, 0.18)
      circos.link(
        sector.index1 = all_links_df$s1[i],
        point1 = c(all_links_df$xstart1[i], all_links_df$xend1[i]),
        sector.index2 = all_links_df$s2[i],
        point2 = c(all_links_df$xstart2[i], all_links_df$xend2[i]),
        col    = adjustcolor(col_val, alpha.f = link_alpha),
        border = adjustcolor(col_val, alpha.f = border_alpha),
        lwd    = link_lwd
      )
    }
  }

  legend("topleft",
         legend = c("Erosion–Erosion", "Rupture–Rupture", "Erosion–Rupture (mixed)"),
         fill   = c(erosion_col, rupture_col, mixed_col),
         border = NA, bty = "n", cex = 0.9)

  title(main = paste0(
    "Shared clonotypes by patient  (OPTICO  |  ",
    chain_mode, "  |  min_count_any = ", min_count_any, ")"
  ))

  dev.off()
  circos.clear()
  cat("Saved:", out_pdf, "\n")
  invisible(list(merged_pat = merged_pat, arc_sizes = arc_sizes))
}

# ── Run all three chain modes ─────────────────────────────────────────────────
# min_count_any = 2  matches min_clone_size = 2 in the main pipeline
# (only expanded clones contribute; consistent with tcr_raw filtering)

make_circos_optico(
  combined_TCR  = combined_TCR_Optico,
  disease_map   = disease_map,
  chain_mode    = "TRA",
  min_count_any = 2,
  top_n         = 60,
  out_pdf = file.path(outdir, "plots", "circos_shared_TRA_Erosion_Rupture.pdf")
)

make_circos_optico(
  combined_TCR  = combined_TCR_Optico,
  disease_map   = disease_map,
  chain_mode    = "TRB",
  min_count_any = 2,
  top_n         = 60,
  out_pdf = file.path(outdir, "plots", "circos_shared_TRB_Erosion_Rupture.pdf")
)

make_circos_optico(
  combined_TCR  = combined_TCR_Optico,
  disease_map   = disease_map,
  chain_mode    = "BOTH",
  min_count_any = 2,
  top_n         = 60,
  out_pdf = file.path(outdir, "plots", "circos_shared_BOTH_CTaa_Erosion_Rupture.pdf")
)


cat("\n── Done ──\n")
cat("  Plots  ->", file.path(outdir, "plots"),  "\n")
cat("  Tables ->", file.path(outdir, "tables"), "\n")
