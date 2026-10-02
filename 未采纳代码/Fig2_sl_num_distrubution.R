library(data.table)
library(tidyverse)

# ====================================================
# 1. Read the data
# ====================================================
df <- fread("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Synergy_Stats/all_gene_sl_counts.csv")

# Basic QC checks
print(colnames(df))
summary(df$SL_Count)

# ====================================================
# 2. Plot the main histogram (clean SCI-journal style)
# ====================================================
final_plot <- ggplot(df, aes(x = SL_Count)) +
  geom_histogram(bins = 120, fill = "#A6CEE3", color = "black", linewidth = 0.2) +
  scale_x_sqrt(breaks = c(1, 10, 50, 100, 500, 1000, 3000, 7000)) +
  labs(
    x = "Number of SL interactions",
    y = "Frequency"
  ) +
  # Set the base font size uniformly to 12
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# ====================================================
# 3. Print and export at high resolution
# ====================================================
print(final_plot)

# With the pie chart gone, narrow the width a bit (e.g. width = 6 or 6.5) for a tighter look
ggsave(
  "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Synergy_Stats/Global_SL_degree_distribution_optimized.pdf",
  final_plot,
  width = 6,
  height = 5,
  dpi = 600
)
