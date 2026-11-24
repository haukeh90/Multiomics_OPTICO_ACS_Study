# ==============================================================================
# OPTICO-ACS | Cytokine production after stimulation (flow cytometry)
#
# Panels : Fig. 4f-h (final panels drawn in Prism from the same table)
# Input  : Optico_revisions/tables/FACS_stimulation/FACS_tidy.xlsx (tidy table of gated frequencies)
# Output : Stats_summary.csv, FigA_*_Paired_Stim.pdf, FigB_FC_Pct_Ruptur_Erosion.pdf (same folder)
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

# ============================================================
#  OPTICO FACS — Cytokine analysis: Ruptur vs Erosion
#  Layout matching HuCIRA figure style
#
#  Figures produced:
#   FigA  — Paired stim/unstim dots per condition (Erosion | Ruptur)
#   FigB  — FC (Stim/NonStim) comparison Ruptur vs Erosion
#   FigC  — L/P ratio (StimL / StimP) Ruptur vs Erosion
#   FigD  — Summary heatmap: median log2 FC per group
#   Stats — Stats_summary.csv with clean Wilcoxon p-values
#
#  GzmB and TNF: sections present but commented out.
#  To activate: add those markers to WSP files, re-run
#  reformat_facs.py, then uncomment the marked sections below.
# ============================================================

# install.packages(c("readxl","dplyr","tidyr","ggplot2","ggpubr","rstatix","patchwork"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ── paths ─────────────────────────────────────────────────────────────────────

DATA_FILE <- file.path(OPTICO_DIR, "Optico_revisions/tables/FACS_stimulation/FACS_tidy.xlsx")   # set working dir to FACS_stimulation folder first
OUT_DIR   <- file.path(OPTICO_DIR, "Optico_revisions/tables/FACS_stimulation")

# ── load ──────────────────────────────────────────────────────────────────────

fc    <- read_excel(DATA_FILE, sheet = "FC_Stim_NonStim")
ratio <- read_excel(DATA_FILE, sheet = "Ratio_StimL_StimP")
tidy  <- read_excel(DATA_FILE, sheet = "Tidy_data")

tissue_labels <- c(L = "Lesional", P = "Peripheral")

# ── factor levels ──────────────────────────────────────────────────────────────
# When GzmB / TNF WSP data is available, add them to cytokine_levels below
# and uncomment the relevant plot sections.

cytokine_levels <- c("IFN-y", "IL-2", "IL-10")
# cytokine_levels <- c("IFN-y", "IL-2", "IL-10", "GzmB", "TNF")   # ← uncomment when ready

cellpop_levels  <- c("CD4+ TEM_TCM", "CD8+ TEM_TCM", "Treg")

factorize <- function(df) {
  df %>% mutate(
    Diagnose        = factor(Diagnose,        levels = c("Erosion", "Ruptur")),
    Cytokine        = factor(Cytokine,        levels = cytokine_levels),
    Cell_population = factor(Cell_population, levels = cellpop_levels)
  )
}

fc    <- factorize(fc)    %>% mutate(Tissue_label = factor(recode(Tissue, !!!tissue_labels), levels = c("Lesional","Peripheral")))
ratio <- factorize(ratio)
tidy  <- factorize(tidy)  %>% mutate(Tissue_label = factor(recode(Tissue, !!!tissue_labels), levels = c("Lesional","Peripheral")))

# ── aesthetics ────────────────────────────────────────────────────────────────
# Matching HuCIRA heatmap: blue = Ruptur, orange = Erosion

COL_DX   <- c(Erosion = "#D4763B", Ruptur = "#1F618D")
COL_STIM <- c("unstim." = "grey60", "stim." = "grey10")

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size),
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      panel.border     = element_rect(fill = NA, color = "grey70", linewidth = 0.4),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(color = "grey40", size = base_size - 1),
      axis.text.x      = element_text(hjust = 0.5)
    )
}

COMP <- list(c("Erosion", "Ruptur"))

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE A — Paired stim / unstim dots, split by Diagnose
#  Matches the existing flow validation panel style in HuCIRA figure
#  Rows = Diagnose (Erosion | Ruptur), Cols = Cytokine × Tissue
#  Separate page per Cell_population for legibility
# ══════════════════════════════════════════════════════════════════════════════

paired <- tidy %>%
  filter(!is.na(Pct)) %>%
  mutate(Stimulation = factor(Stimulation, levels = c("NonStim","Stim"),
                              labels = c("unstim.","stim.")))

make_figA <- function(cellpop, metric = "Pct", ylab = "% cytokine⁺ cells") {
  d <- paired %>% filter(Cell_population == cellpop)
  ggplot(d, aes(x = Stimulation, y = .data[[metric]])) +
    geom_line(aes(group = SampleID), color = "grey70", linewidth = 0.35, alpha = 0.8) +
    geom_point(aes(color = Stimulation), size = 1.6, alpha = 0.9) +
    stat_compare_means(comparisons = list(c("unstim.","stim.")),
                       paired = TRUE, method = "wilcox.test",
                       label = "p.format", label.y.npc = 0.90,
                       size = 2.8, tip.length = 0.01, bracket.size = 0.4) +
    scale_color_manual(values = COL_STIM) +
    facet_grid(Diagnose ~ Cytokine + Tissue_label, scales = "free_y") +
    labs(title    = paste0(cellpop, " — Stim vs Unstim"),
         subtitle = paste0(ylab, "  |  Wilcoxon paired, raw p"),
         x = NULL, y = ylab, color = NULL) +
    theme_pub() +
    theme(legend.position = "none")
}

figA_cd4  <- make_figA("CD4+ TEM_TCM")
figA_cd8  <- make_figA("CD8+ TEM_TCM")
figA_treg <- make_figA("Treg")

# MFI version (secondary) — swap metric = "MFI"
# figA_cd4_mfi  <- make_figA("CD4+ TEM_TCM",  metric = "MFI", ylab = "Geometric Mean (MFI)")
# figA_cd8_mfi  <- make_figA("CD8+ TEM_TCM",  metric = "MFI", ylab = "Geometric Mean (MFI)")
# figA_treg_mfi <- make_figA("Treg",           metric = "MFI", ylab = "Geometric Mean (MFI)")

# ── GzmB / TNF paired plots (uncomment when WSP data available) ───────────────
# figA_cd8_gzmb <- make_figA("CD8+ TEM_TCM") +
#   filter(Cytokine == "GzmB")    # restrict once GzmB is in cytokine_levels
# figA_all_tnf  <- make_figA("CD4+ TEM_TCM") +
#   filter(Cytokine == "TNF")

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE B — FC (Stim/NonStim) comparison: Ruptur vs Erosion
#  Box + individual dots; separate panels per tissue
# ══════════════════════════════════════════════════════════════════════════════

fc_pct <- fc %>% filter(!is.na(log2FC_Pct))

figB <- ggplot(fc_pct, aes(x = Diagnose, y = log2FC_Pct, fill = Diagnose, color = Diagnose)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.75) +
  stat_compare_means(comparisons = COMP, method = "wilcox.test",
                     label = "p.format", label.y.npc = 0.92,
                     size = 2.8, tip.length = 0.01, bracket.size = 0.4) +
  scale_fill_manual(values = COL_DX) +
  scale_color_manual(values = COL_DX) +
  facet_grid(Cell_population ~ Cytokine + Tissue_label, scales = "free_y") +
  labs(title    = "Stimulation fold change — Ruptur vs Erosion",
       subtitle = "log₂ FC (Stim / NonStim)  |  % cytokine⁺ cells  |  Wilcoxon, raw p",
       x = NULL, y = "log₂ FC (Stim / NonStim)",
       fill = "Diagnosis", color = "Diagnosis") +
  theme_pub()

# ── GzmB / TNF FC comparison (uncomment when data available) ─────────────────
# fc_extra <- fc %>% filter(!is.na(log2FC_Pct), Cytokine %in% c("GzmB","TNF"))
# figB_extra <- figB %+% fc_extra   # reuses same ggplot spec, swaps data

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE C — L/P ratio (Stim-Lesional / Stim-Peripheral): Ruptur vs Erosion
# ══════════════════════════════════════════════════════════════════════════════

ratio_pct <- ratio %>% filter(!is.na(log2Ratio_Pct))

figC <- ggplot(ratio_pct, aes(x = Diagnose, y = log2Ratio_Pct, fill = Diagnose, color = Diagnose)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.5,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.75) +
  stat_compare_means(comparisons = COMP, method = "wilcox.test",
                     label = "p.format", label.y.npc = 0.92,
                     size = 2.8, tip.length = 0.01, bracket.size = 0.4) +
  scale_fill_manual(values = COL_DX) +
  scale_color_manual(values = COL_DX) +
  facet_grid(Cell_population ~ Cytokine, scales = "free_y") +
  labs(title    = "Lesional enrichment — Ruptur vs Erosion",
       subtitle = "log₂ (Stim-Lesional / Stim-Peripheral)  |  % cytokine⁺ cells  |  Wilcoxon, raw p",
       x = NULL, y = "log₂ (Stim-Lesional / Stim-Peripheral)",
       fill = "Diagnosis", color = "Diagnosis") +
  theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE D — Summary heatmap: median log2 FC per Diagnose × cell pop × cytokine
#  Mirrors style of HuCIRA heatmap (panel b), closing the loop visually
# ══════════════════════════════════════════════════════════════════════════════

heatmap_data <- fc_pct %>%
  group_by(Diagnose, Cell_population, Cytokine, Tissue_label) %>%
  summarise(median_log2FC = median(log2FC_Pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = round(median_log2FC, 2))

figD <- ggplot(heatmap_data,
               aes(x = Cytokine, y = Cell_population, fill = median_log2FC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3, color = "grey20") +
  scale_fill_gradient2(low  = "#D4763B", mid = "white", high = "#1F618D",
                       midpoint = 0, name = "Median\nlog₂ FC") +
  facet_grid(Tissue_label ~ Diagnose) +
  labs(title    = "Summary: median log₂ FC (Stim/NonStim)",
       subtitle = "% cytokine⁺ cells",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "right")

# ── GzmB / TNF heatmap rows (uncomment when data available) ──────────────────
# Just adding GzmB/TNF to cytokine_levels above is sufficient; heatmap auto-expands.

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE
# ══════════════════════════════════════════════════════════════════════════════

ggsave(file.path(OUT_DIR, "FigA_CD4_Paired_Stim.pdf"),    figA_cd4,  width = 12, height = 7)
ggsave(file.path(OUT_DIR, "FigA_CD8_Paired_Stim.pdf"),    figA_cd8,  width = 12, height = 7)
ggsave(file.path(OUT_DIR, "FigA_Treg_Paired_Stim.pdf"),   figA_treg, width =  8, height = 7)
ggsave(file.path(OUT_DIR, "FigB_FC_Pct_Ruptur_Erosion.pdf"), figB,   width = 14, height = 9)
ggsave(file.path(OUT_DIR, "FigC_Ratio_LP.pdf"),            figC,     width = 10, height = 9)
ggsave(file.path(OUT_DIR, "FigD_Heatmap_FC.pdf"),          figD,     width =  8, height = 6)

message("All figures saved.")

# ══════════════════════════════════════════════════════════════════════════════
#  STATS — clean export, separate from plot objects
# ══════════════════════════════════════════════════════════════════════════════

stats_paired <- tidy %>%
  filter(!is.na(Pct)) %>%
  mutate(Stimulation = factor(Stimulation, levels = c("NonStim","Stim"))) %>%
  group_by(Cell_population, Cytokine, Tissue_label, Diagnose) %>%
  wilcox_test(Pct ~ Stimulation, paired = TRUE) %>%
  mutate(Metric = "Pct_paired_StimVsNonStim", Figure = "FigA")

stats_fc_pct <- fc_pct %>%
  group_by(Cell_population, Cytokine, Tissue_label) %>%
  wilcox_test(log2FC_Pct ~ Diagnose) %>%
  mutate(Tissue_label = as.character(Tissue_label),
         Metric = "log2FC_Pct", Figure = "FigB")

stats_ratio <- ratio_pct %>%
  group_by(Cell_population, Cytokine) %>%
  wilcox_test(log2Ratio_Pct ~ Diagnose) %>%
  mutate(Tissue_label = "L/P ratio", Metric = "log2Ratio_Pct", Figure = "FigC")

stats_all <- bind_rows(stats_paired, stats_fc_pct, stats_ratio) %>%
  select(Figure, Metric, Cell_population, Cytokine, Tissue_label,
         group1, group2, n1, n2, statistic, p) %>%
  mutate(p_round = round(p, 4),
         sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                         p < 0.05  ~ "*",   TRUE      ~ "ns")) %>%
  arrange(Figure, Cell_population, Cytokine)

write.csv(stats_all, file.path(OUT_DIR, "Stats_summary.csv"), row.names = FALSE)
message("Stats saved.")
print(stats_all %>% select(Figure, Cell_population, Cytokine, Tissue_label,
                            Diagnose = group1, p_round, sig), n = Inf)
