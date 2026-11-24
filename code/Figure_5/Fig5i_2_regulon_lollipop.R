# ==============================================================================
# OPTICO-ACS | Lollipop plot of regulon activity (delta AUC, Erosion - Rupture) in oligoclonal T cells
#
# Panels : Fig. 5i
# Input  : T_reclustering_2025/data/SCENIC/SCENIC_inspection_tables.xlsx
# Output : T_reclustering_2025/plots/expanded_DEG/SCENIC_lollipop_ErosionVsRupture_*.pdf
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# ── Paths ─────────────────────────────────────────────────────────────────────
scenic_dir <- file.path(OPTICO_DIR, "T_reclustering_2025/data/SCENIC")
plot_dir   <- file.path(OPTICO_DIR, "T_reclustering_2025/plots/expanded_DEG")
xlsx_path  <- file.path(scenic_dir, "SCENIC_inspection_tables.xlsx")

# ── Colours with transparency ─────────────────────────────────────────────────
col_erosion <- alpha("#904011", 0.75)
col_rupture <- alpha("#0F1E5E", 0.75)

# ── Read Erosion vs Rupture sheets ────────────────────────────────────────────
sheets         <- excel_sheets(xlsx_path)
ero_rup_sheets <- sheets[grepl("Erosion_vs_Rupture", sheets)]

dat <- lapply(ero_rup_sheets, function(sh) {
  read_excel(xlsx_path, sheet = sh) |>
    mutate(
      sheet   = sh,
      lineage = case_when(
        grepl("CD8", sh) ~ "CD8",
        grepl("CD4", sh) ~ "CD4"
      )
    )
}) |> bind_rows()

# ── Filter: FDR < 0.05 and known T cell TF only ───────────────────────────────
dat_filt <- dat |>
  filter(FDR < 0.05, known_Tcell_TF == TRUE) |>
  mutate(
    direction     = factor(direction, levels = c("Erosion", "Rupture")),
    regulon_label = gsub("\\(\\+\\)", "", regulon)
  )

# ── Plot function ─────────────────────────────────────────────────────────────
plot_lollipop <- function(df, lineage_label, show_labels = FALSE) {
  
  df <- df |>
    filter(lineage == lineage_label) |>
    arrange(delta_AUC) |>
    mutate(regulon_label = factor(regulon_label, levels = regulon_label))
  
  p <- ggplot(df, aes(x = delta_AUC, y = regulon_label, colour = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_segment(aes(x = 0, xend = delta_AUC,
                     y = regulon_label, yend = regulon_label),
                 linewidth = 0.5) +
    geom_point(size = 5, shape = 16) +
    scale_colour_manual(
      values = c("Erosion" = col_erosion, "Rupture" = col_rupture),
      name   = NULL
    ) + scale_x_continuous(
      limits = c(-0.01, 0.01),
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0))   # no extra padding beyond the fixed limits
    )
  
  if (show_labels) {
    p <- p +
      labs(
        title = paste0("Regulons in oligoclonal ", lineage_label, "\u207a T cells"),
        x     = expression(Delta * "AUC (Erosion \u2212 Rupture)"),
        y     = NULL
      ) +
      geom_text(
        aes(label = sprintf("q=%.3f  \u0394AUC=%.4f", FDR, delta_AUC),
            x     = delta_AUC + sign(delta_AUC) * 0.0005),
        hjust    = ifelse(df$delta_AUC > 0, 0, 1),
        size     = 2.2, colour = "grey35"
      ) +
      theme_classic(base_size = 9, base_family = "Helvetica") +
      theme(
        plot.title          = element_text(face = "plain", size = 8.5,
                                           hjust = 0.5, margin = margin(b = 4)),
        plot.title.position = "plot",
        panel.border        = element_blank(),
        plot.background     = element_blank(),
        axis.text.y         = element_text(face = "italic", size = 7.5, colour = "black"),
        axis.text.x         = element_text(size = 7.5, colour = "black"),
        axis.title.x        = element_text(size = 8),
        axis.line           = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks          = element_line(linewidth = 0.3),
        legend.position     = "top",
        legend.text         = element_text(size = 8),
        legend.key.size     = unit(0.35, "cm"),
        panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_line(colour = "grey92", linewidth = 0.25),
        panel.grid.minor    = element_blank()
      )
  } else {
    p <- p +
      labs(title = NULL, x = NULL, y = NULL) +
      theme_void(base_family = "Helvetica") +
      theme(
        legend.position = "none",
        plot.margin     = margin(4, 4, 4, 4)
      )
  }
  
  p
}

# ── Export ────────────────────────────────────────────────────────────────────
row_height_cm <- 0.55

# Calculate once — use the larger plot to set the shared canvas height
n_cd8 <- sum(dat_filt$lineage == "CD8")
n_cd4 <- sum(dat_filt$lineage == "CD4")
shared_height <- max(n_cd8, n_cd4) * row_height_cm  # both PDFs get this height
shared_width  <- 5   # clean version
labeled_width <- 14  # labeled version

for (lin in c("CD8", "CD4")) {
  
  p_clean <- plot_lollipop(dat_filt, lin, show_labels = FALSE)
  ggsave(
    file.path(plot_dir, paste0("SCENIC_lollipop_ErosionVsRupture_", lin, "_clean.pdf")),
    p_clean,
    width  = shared_width,
    height = shared_height,   # ← same for both
    units  = "cm",
    device = cairo_pdf
  )
  
  p_label <- plot_lollipop(dat_filt, lin, show_labels = TRUE)
  ggsave(
    file.path(plot_dir, paste0("SCENIC_lollipop_ErosionVsRupture_", lin, "_labeled.pdf")),
    p_label,
    width  = labeled_width,
    height = shared_height,   # ← same for both
    units  = "cm",
    device = cairo_pdf
  )
  
  cat("Saved:", lin, "\n")
}
