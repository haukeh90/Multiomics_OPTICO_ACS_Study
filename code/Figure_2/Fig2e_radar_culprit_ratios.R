# ==============================================================================
# OPTICO-ACS | Radar plot of culprit ratios of CD4+ T helper subsets
#
# Panels : Fig. 2e
# Input  : Figures/Figure_1/tables/TH_subtypes_for_spider_CR.xlsx (culprit ratios from flow cytometry)
# Output : Figures/Figure_1/plots/radarchart_no_labels.pdf
# ==============================================================================

OPTICO_DIR <- Sys.getenv("OPTICO_DIR", unset = "~/OPTICO")   # project folder (data/, tables/, plots/ ...)

library(readxl)     
library(tidyverse)  
library(fmsb)       

df <- read_excel(file.path(OPTICO_DIR, "Figures/Figure_1/tables/TH_subtypes_for_spider_CR.xlsx"))
head(df)

# Step 2: Reshape and calculate means
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    cell_type = sub("_.*", "", group),
    condition = sub(".*_", "", group)
  )

# Step 3: Calculate means for each group
df_mean <- df_long %>%
  group_by(cell_type, condition) %>%
  summarise(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_value)

# Step 4: Reformat for spider plot
# Add rows for max and min values
plot_data <- df_mean %>%
  column_to_rownames("cell_type") %>%
  t() %>%
  as.data.frame()

plot_data <- rbind(
  rep(max(plot_data, na.rm = TRUE), ncol(plot_data)),
  rep(min(plot_data, na.rm = TRUE), ncol(plot_data)),
  plot_data
)

# Step 5: Plot the radar chart
colors_border <- c("#904011", "#0F1E5E")
labels <- colnames(plot_data)

colors_fill <- c(
  rgb(144/255, 64/255, 17/255, alpha = 0.3),  # "#904011"
  rgb(15/255, 30/255, 94/255, alpha = 0.3)    # "#0F1E5E"
)

# Start PDF device
pdf(file.path(OPTICO_DIR, "Figures/Figure_1/plots/radarchart_no_labels.pdf"), width = 6, height = 6)

# Create radar chart without axis labels or variable labels
# Plot without axis and labels
radarchart(plot_data,
           axistype = 0,
           pcol = colors_border, 
           pfcol = colors_fill,
           plwd = 3,
           plty = 1,
           cglcol = "darkgrey", 
           cglty = 1,
           cglwd = 2,
           axislabcol = NA,
           vlcex = 0,         # Safe to keep
           title = NULL)

# Close PDF device
dev.off()
