# 1. Load required packages
library(ggplot2)
library(reshape2)
library(grid)
library(gtable)
library(ggtext)

methods <- c(
  "DAISY", "MiSL", "SiLi", "DiscoverSL", "ISLE", "Wang (2013)", "Wang (2019)",
  "underMutExSL", "De Kegel (2021)", "MetaSL", "Pandey (2010)", "SL2MF"
)

features <- c(
  "DiffExp", "DiffSof", "ExcluISLEmRNA", "ExcluISLECNV", "ExcluCNVAmp", "ExcluCNVdel", "ExcluMut", "ExcluAlt",
  "SurvISLEmRNA", "SurvISLECNV", "SurvSILI",
  "mRNAEssCrispr", "mRNaEssRNAi", "scnaEssCrispr", "scnaEssRNAi", "scnaEssCrisprDaisy",
  "scnaEssRNAiDaisy", "ppiEssScore", "ComplexEssScore",
  "coExp", "Complex_10", "GOscore", "paralog",
  "jacard_colocalisation", "Jacard_Domain", "pathway_p",
  "PPI_share", "PPI_union", "proteinInt",
  "geneConservationScore", "Mean_age"
)

data_types <- c(
  rep("Genomic Survival of the Fittest", 8),
  rep("Prognostic Relevance", 3),
  rep("Gene Dependency", 8),
  rep("Gene Function Similarity", 7),
  rep("Network Topological Attributes", 3),
  rep("Gene Evolutionary", 2)
)

# 3. Key change: match colors to data types
# Define 6 distinct colors
type_colors <- c(
  "Genomic Survival of the Fittest" = "#E64B35",
  "Prognostic Relevance" = "#4DBBD5",
  "Gene Dependency" = "#00A087",
  "Gene Function Similarity" = "#3C5488",
  "Network Topological Attributes" = "#F39B7F",
  "Gene Evolutionary" = "#8491B4"
)

feature_meta <- data.frame(Feature = features, DataType = data_types, stringsAsFactors = FALSE)

mat <- read.csv("E:\\合成致死\\论文绘图代码\\方法_特征对应表.csv", row.names = 1)
rownames(mat) <- features
colnames(mat) <- methods
df <- melt(as.matrix(mat))
colnames(df) <- c("Feature", "Method", "Used")

# Create HTML-rendered Y-axis labels
feature_meta$StyledLabel <- paste0(
  "<span style='color:",
  type_colors[feature_meta$DataType],
  "'>", feature_meta$Feature, "</span>"
)

# Merge the data
df <- merge(df, feature_meta, by = "Feature")

# Fix the order: Y-axis follows the original features list (reversed so it reads top-down)
df$StyledLabel <- factor(df$StyledLabel, levels = rev(feature_meta$StyledLabel))
# Fix the facet order to the order in which data_types appear
df$DataType <- factor(df$DataType, levels = unique(data_types))

# 4. Plot
p <- ggplot(df, aes(x = Method, y = StyledLabel)) +
  geom_tile(aes(fill = factor(Used)), color = "white", linewidth = 0.5) +
  # space = "free_y" lets facets with different feature counts size themselves
  facet_grid(DataType ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c("0" = "#f0f0f0", "1" = "#34495e"),
    labels = c("Not Used", "Used"), name = "Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    # Use element_markdown to render HTML colors
    axis.text.y = element_markdown(size = 9, face = "bold"),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8, color = "white"),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "top",
    strip.background = element_rect(fill = "grey", color = NA), # Initial colors
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "Feature-Method Utilization Matrix",
    x = "Computational Methods",
    y = "Multi-omics Features"
  )

# 5. Customize the background color of the right-side facet strips
g <- ggplot_gtable(ggplot_build(p))
# Find the indices of all right-side strips
strip_indices <- which(grepl("strip-r", g$layout$name))
# Assign colors in facet-level order
ordered_types <- levels(df$DataType)

for (i in seq_along(strip_indices)) {
  index <- strip_indices[i]
  fill_color <- type_colors[ordered_types[i]]
  # Modify the fill of the strip rectangles
  g$grobs[[index]]$grobs[[1]]$children[[1]]$gp$fill <- fill_color
  g$grobs[[index]]$grobs[[1]]$children[[1]]$gp$col <- fill_color # Match the border color too
}

# 6. Plot output
grid.newpage()
grid.draw(g)

ggsave("E:\\合成致死\\论文绘图代码\\result.pdf", g, width = 10, height = 12, device = cairo_pdf)
