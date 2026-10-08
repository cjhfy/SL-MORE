#!/usr/bin/env Rscript
# =============================================================================
# Figure 2A: Feature-Method Provenance Matrix
#
# Panels:
#   1. Central matrix : 30 multi-omics features x 13 SL inference methods
#   2. Right barplot  : number of methods supporting each feature
#   3. Top barplot    : number of feature dimensions covered by each method
#
# Input : method_feature.csv (binary feature x method matrix)
# Output: Feature_method.pdf / Feature_use_count.pdf / Method_use_count.pdf
# =============================================================================

library(ggplot2)
library(reshape2)
library(grid)
library(gtable)
library(dplyr)
library(cowplot)

# ---- Configuration ----------------------------------------------------------

input_csv <- "/data/home/chenjiahao/nuaa/synlethDB/SL-MORE/data/Fig2/method_feature.csv"
out_dir   <- "/data/home/chenjiahao/nuaa/synlethDB/Feature_Method"

# ---- 1. Data preparation ----------------------------------------------------

# Six functional dimensions (one label per feature, in CSV row order)
data_types <- c(
  rep("Genomic Selection Signatures", 7),
  rep("Survival Benefit", 3),
  rep("Gene Dependency", 8),
  rep("Gene Function Similarity", 7),
  rep("Network Properties", 3),
  rep("Evolutionary Conservation", 2)
)

# Soft palette for the six dimensions
type_colors <- c(
  "Genomic Selection Signatures" = "#789ECF",
  "Survival Benefit"             = "#AAD4B9",
  "Gene Dependency"              = "#F2A471",
  "Gene Function Similarity"     = "#6AC3D3",
  "Network Properties"           = "#EDB9D2",
  "Evolutionary Conservation"    = "#A465A4"
)

# Read the binary feature x method matrix
mat <- read.csv(input_csv, row.names = 1, check.names = FALSE)

features_in_csv <- rownames(mat)
methods_in_csv  <- colnames(mat)

# Feature metadata: map each feature to its functional dimension
feature_meta <- data.frame(
  Feature  = features_in_csv,
  DataType = data_types,
  stringsAsFactors = FALSE
)

# Convert to long format and attach metadata
df <- melt(as.matrix(mat))
colnames(df) <- c("Feature", "Method", "Used")
df <- merge(df, feature_meta, by = "Feature")

# Lock factor levels to preserve the original ordering
df$Method   <- factor(df$Method, levels = methods_in_csv)
df$Feature  <- factor(df$Feature, levels = rev(features_in_csv))
df$DataType <- factor(df$DataType, levels = unique(data_types))

# ---- 2. Build the three panels ----------------------------------------------

# Panel 1: central provenance matrix
p <- ggplot(df, aes(x = Method, y = Feature)) +
  geom_tile(aes(fill = factor(Used)), color = "white", linewidth = 0.5) +
  facet_grid(DataType ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(position = "right") +
  scale_fill_manual(
    values = c("0" = "#F3F3F3", "1" = "#4A5A70"),
    labels = c("0" = "Not Derived", "1" = "Source"),
    name   = "Provenance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1,
                                      size = 12, face = "bold", color = "black"),
    axis.text.y        = element_text(size = 12, face = "bold", hjust = 0,
                                      color = "black"),
    strip.text.y.left  = element_text(angle = 0, face = "bold", size = 14,
                                      color = "white"),
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "grey", color = NA),
    panel.spacing      = unit(0.2, "lines"),
    panel.border       = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 12),
    legend.title       = element_text(size = 12, face = "bold"),
    plot.title         = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(
    title = "Feature-Method Provenance Matrix",
    x     = "Computational Methods",
    y     = ""
  )

# Panel 2: right barplot - number of methods supporting each feature
# Uses the same faceting and spacing as the matrix so rows align exactly
feat_counts <- df %>%
  group_by(Feature, DataType) %>%
  summarise(MethodCount = sum(Used), .groups = "drop")

p_feat_count <- ggplot(feat_counts, aes(y = Feature, x = MethodCount, fill = DataType)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9, show.legend = FALSE) +
  facet_grid(DataType ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = type_colors) +
  scale_x_continuous(
    breaks = 0:length(methods_in_csv),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal() +
  theme(
    axis.text.y      = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x      = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x     = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y     = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    panel.grid       = element_blank(),
    panel.spacing    = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, size = 12, face = "bold",
                                    margin = margin(b = 15))
  ) +
  labs(
    title = "Prevalence of Features",
    x     = "Method Count"
  )

# Panel 3: top barplot - number of feature dimensions covered by each method
method_cat_presence <- df %>%
  filter(Used == 1) %>%
  group_by(Method, DataType) %>%
  summarise(CategoryPresence = 1, .groups = "drop")

method_cat_presence$Method   <- factor(method_cat_presence$Method, levels = methods_in_csv)
method_cat_presence$DataType <- factor(method_cat_presence$DataType, levels = unique(data_types))

p_method_count <- ggplot(method_cat_presence,
                         aes(x = Method, y = CategoryPresence, fill = DataType)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",
           linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = type_colors) +
  scale_y_continuous(
    breaks = 0:length(unique(data_types)),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y      = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(hjust = 0.5, size = 12, face = "bold",
                                     margin = margin(b = 15))
  ) +
  labs(
    title = "Number of Feature Dimensions Used per Method",
    y     = "Category Count"
  )

# ---- 3. Align panels and color the facet strips -----------------------------

# Vertical alignment: top barplot edges match the central matrix
aligned_v <- cowplot::align_plots(p_method_count, p, align = "v", axis = "lr")
g_top  <- aligned_v[[1]]
g_main <- aligned_v[[2]]

# Horizontal alignment: right barplot panel height matches the matrix exactly
aligned_h     <- cowplot::align_plots(g_main, p_feat_count, align = "h", axis = "tb")
g_main_final  <- aligned_h[[1]]
g_right       <- aligned_h[[2]]

# Color the left facet strips by dimension (on the final aligned gtable)
strip_indices <- which(grepl("strip-l", g_main_final$layout$name))
ordered_types <- levels(df$DataType)

for (i in seq_along(strip_indices)) {
  index      <- strip_indices[i]
  fill_color <- type_colors[ordered_types[i]]
  g_main_final$grobs[[index]]$grobs[[1]]$children[[1]]$gp$fill <- fill_color
  g_main_final$grobs[[index]]$grobs[[1]]$children[[1]]$gp$col  <- fill_color
}

# ---- 4. Export --------------------------------------------------------------

ggsave(file.path(out_dir, "Feature_method.pdf"),
       plot = g_main_final, width = 9, height = 8)

ggsave(file.path(out_dir, "Feature_use_count.pdf"),
       plot = g_right, width = 4.5, height = 8, dpi = 300)

ggsave(file.path(out_dir, "Method_use_count.pdf"),
       plot = g_top, width = 9, height = 2, dpi = 300)
