# =========================================================
# Cancer-specific SL landscape analysis
# High-quality SCI-style visualization (Optimized version)
# =========================================================

library(data.table)
library(arrow)
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(RColorBrewer)
library(corrplot)

# ─────────────────────────────────────────────
# 1. Path & Configuration
# ─────────────────────────────────────────────
base_out_dir <- "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/"
gt4_dir <- file.path(base_out_dir, "Cancer_HighConf_GT4_Features")
gene_stat_dir <- file.path(base_out_dir, "Synergy_Stats/Gene_SL_Stats")
dir.create(gene_stat_dir, recursive = TRUE, showWarnings = FALSE)

cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

# ─────────────────────────────────────────────
# 2. Load protein-coding genes
# ─────────────────────────────────────────────
gene_anno <- fread("/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE
)
gene_anno <- gene_anno[gene_anno$type_of_gene == "protein-coding", c(2, 3, 5)]
total_protein_coding <- nrow(gene_anno)
cat("Total protein-coding genes:", total_protein_coding, "\n")

# ─────────────────────────────────────────────
# 3. Read and process data (Optimized)
# ─────────────────────────────────────────────
cancer_genes <- list()
cancer_pairs <- list()
sl_data_cache <- list() # Cache the raw data to avoid re-reading
stats <- data.frame(Cancer = cancers, Gene_Count = 0, Pair_Count = 0)

for (cancer in cancers) {
  file_path <- file.path(gt4_dir, paste0(cancer, "_GT4_Features.parquet"))
  if (file.exists(file_path)) {
    sl_data <- read_parquet(file_path)

    # Standardize: create the Pair column and deduplicate
    sl_processed <- sl_data %>%
      mutate(Pair = paste(Symbol1, Symbol2, sep = "_")) %>%
      distinct(Symbol1, Symbol2, .keep_all = TRUE) %>%
      select(Symbol1, Symbol2, Pair)

    # Extract genes and interaction pairs
    cancer_genes[[cancer]] <- unique(c(sl_processed$Symbol1, sl_processed$Symbol2))
    cancer_pairs[[cancer]] <- sl_processed$Pair

    # Cache the standardized data (for the later Jaccard computation)
    sl_data_cache[[cancer]] <- sl_processed

    # Update the summary statistics
    stats[stats$Cancer == cancer, c("Gene_Count", "Pair_Count")] <-
      c(length(cancer_genes[[cancer]]), length(cancer_pairs[[cancer]]))

    cat(
      cancer, "- Genes:", length(cancer_genes[[cancer]]),
      "Pairs:", length(cancer_pairs[[cancer]]), "\n"
    )
  }
}

# Remove cancer types that failed to load
valid_cancers <- names(cancer_genes)
n <- length(valid_cancers)
cat("\nSuccessfully loaded", n, "cancer types\n")

# ─────────────────────────────────────────────
# 4. Order cancer types
# ─────────────────────────────────────────────
stats <- stats %>%
  filter(Cancer %in% valid_cancers) %>%
  arrange(desc(Gene_Count))
stats$Cancer <- factor(stats$Cancer, levels = stats$Cancer)

# =========================================================
# 5. Gene count barplot
# =========================================================
p1 <- ggplot(stats, aes(x = Cancer, y = Gene_Count, fill = Gene_Count)) +
  geom_bar(stat = "identity", width = 0.78, color = "black", linewidth = 0.28) +
  geom_hline(yintercept = total_protein_coding, linetype = "dashed", color = "#8B1E3F", linewidth = 0.9) +
  annotate("text",
    x = length(valid_cancers) - 1, y = total_protein_coding,
    label = paste0("Protein-coding genes: ", comma(total_protein_coding)),
    color = "#8B1E3F", hjust = 1, vjust = -0.6, size = 3.8, fontface = "italic"
  ) +
  # geom_text(aes(label = comma(Gene_Count)), vjust = -0.45, size = 3) +
  scale_fill_gradient(low = "#DCEAF7", high = "#0B3C5D") +
  labs(
    title = "Genes involved in SL interactions across cancer types",
    subtitle = paste0("Background protein-coding genes: ", comma(total_protein_coding)),
    x = "Cancer type", y = "Gene count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none"
  )

ggsave(file.path(gene_stat_dir, "Gene_Count_Barplot.pdf"), p1, width = 7, height = 6)

# =========================================================
# 6. SL pair count barplot
# =========================================================
p2 <- ggplot(stats, aes(x = Cancer, y = Pair_Count, fill = Pair_Count)) +
  geom_bar(stat = "identity", width = 0.78, color = "black", linewidth = 0.28) +
  geom_text(aes(label = comma(Pair_Count)), vjust = -0.45, size = 3) +
  scale_fill_gradient(low = "#D8F0F0", high = "#145C75") +
  labs(
    title = "SL interaction pairs across cancer types",
    x = "Cancer type", y = "SL pair count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none"
  )

ggsave(file.path(gene_stat_dir, "Pair_Count_Barplot.pdf"), p2, width = 14, height = 7)
ggsave(file.path(gene_stat_dir, "Pair_Count_Barplot.png"), p2, width = 14, height = 7, dpi = 600)

# =========================================================
# 7. Jaccard similarity — gene sets
# =========================================================
jaccard_gene <- matrix(0, n, n, dimnames = list(valid_cancers, valid_cancers))

for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      jaccard_gene[i, j] <- 1
    } else {
      inter <- length(intersect(cancer_genes[[i]], cancer_genes[[j]]))
      union_sz <- length(union(cancer_genes[[i]], cancer_genes[[j]]))
      jaccard_gene[i, j] <- inter / union_sz
    }
  }
}

pdf(file.path(gene_stat_dir, "CancerGeneIntersect.pdf"), width = 12, height = 10)
corrplot(jaccard_gene,
  type = "upper", order = "original",
  col = rev(brewer.pal(n = 10, name = "RdBu")),
  tl.col = "black", tl.cex = 0.8, tl.srt = 90,
  is.corr = FALSE, diag = TRUE, cl.cex = 0.8,
  col.lim = c(0, 1), font = 3,
  title = "Jaccard similarity of SL gene sets",
  mar = c(0, 0, 2, 0)
)
dev.off()

# =========================================================
# 8. Conditional Jaccard similarity — SL pairs
# (only consider shared genes)
# =========================================================
jaccard_pair <- matrix(0, n, n, dimnames = list(valid_cancers, valid_cancers))

for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      jaccard_pair[i, j] <- 1
    } else {
      cancer_i <- valid_cancers[i]
      cancer_j <- valid_cancers[j]

      # Get shared genes
      common_genes <- intersect(cancer_genes[[cancer_i]], cancer_genes[[cancer_j]])
      # Subset the cached data (note the columns are Symbol1 and Symbol2)
      sl_i <- sl_data_cache[[cancer_i]]
      sl_j <- sl_data_cache[[cancer_j]]

      # Keep only pairs whose genes are in the shared set (column names corrected)
      sl_i_filtered <- sl_i %>%
        filter(Symbol1 %in% common_genes & Symbol2 %in% common_genes)

      sl_j_filtered <- sl_j %>%
        filter(Symbol1 %in% common_genes & Symbol2 %in% common_genes)

      # Set Jaccard to 0 when no pairs remain after filtering
      if (nrow(sl_i_filtered) == 0 || nrow(sl_j_filtered) == 0) {
        jaccard_pair[i, j] <- 0
        next
      }

      # Use the Pair column directly
      pairs_i <- sl_i_filtered$Pair
      pairs_j <- sl_j_filtered$Pair

      # Compute conditional Jaccard
      inter <- length(intersect(pairs_i, pairs_j))
      union_sz <- length(union(pairs_i, pairs_j))
      jaccard_pair[i, j] <- inter / union_sz
    }
  }
}

# Save the conditional Jaccard matrix
write.csv(jaccard_pair, file.path(gene_stat_dir, "Conditional_Jaccard_SLpairs_Matrix.csv"))

# Plot the conditional Jaccard heatmap
pdf(file.path(gene_stat_dir, "CancerSlIntersect.pdf"), width = 12, height = 10)
corrplot(jaccard_pair,
  type = "upper", order = "original",
  col = rev(brewer.pal(n = 10, name = "RdBu")),
  tl.col = "black", tl.cex = 0.8, tl.srt = 90,
  is.corr = FALSE, diag = TRUE, cl.cex = 0.8,
  col.lim = c(0, 1), font = 3,
  title = "Conditional Jaccard similarity of SL pairs (shared genes only)",
  mar = c(0, 0, 2, 0)
)
dev.off()

# =========================================================
# 9. Combined half heatmap (Nature style)
# =========================================================
# Prepare the data
gene_df <- melt(jaccard_gene)
pair_df <- melt(jaccard_pair)
colnames(gene_df) <- c("Cancer1", "Cancer2", "Value")
colnames(pair_df) <- c("Cancer1", "Cancer2", "Value")

gene_df$i <- match(gene_df$Cancer1, valid_cancers)
gene_df$j <- match(gene_df$Cancer2, valid_cancers)
pair_df$i <- match(pair_df$Cancer1, valid_cancers)
pair_df$j <- match(pair_df$Cancer2, valid_cancers)

# Upper triangle: gene-set Jaccard; lower triangle: conditional SL-pair Jaccard
upper_df <- gene_df %>% filter(i <= j)
lower_df <- pair_df %>% filter(i > j)
plot_df <- bind_rows(upper_df, lower_df)

# Set factor levels
plot_df$Cancer1 <- factor(plot_df$Cancer1, levels = valid_cancers)
plot_df$Cancer2 <- factor(plot_df$Cancer2, levels = rev(valid_cancers))

# Color scheme
cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)

# Plot
p <- ggplot(plot_df, aes(x = Cancer1, y = Cancer2, fill = Value)) +
  geom_tile(color = alpha("white", 0.7), linewidth = 0.25) +
  scale_fill_gradientn(
    colours = cols,
    limits = c(0, 1),
    values = rescale(c(0, 0.05, 0.15, 0.35, 0.6, 1)),
    name = "Jaccard\nsimilarity"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      face = "italic", size = 10, color = "black"
    ),
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey30", linewidth = 0.6),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    legend.key.height = unit(2, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
  ) +
  labs(
    title = "Jaccard similarity across cancer types",
    subtitle = "Upper triangle: Gene set similarity | Lower triangle: Conditional SL pair similarity"
  )

ggsave(file.path(gene_stat_dir, "Combined_Jaccard_RdBu_Style.pdf"), p, width = 10, height = 9)
print(p)

# =========================================================
# 10. Additional analysis: Jaccard comparison plot
# =========================================================
# Run the correlation analysis
compare_df <- data.frame(
  Gene_Jaccard = jaccard_gene[upper.tri(jaccard_gene)],
  Pair_Jaccard = jaccard_pair[upper.tri(jaccard_pair)]
)

cor_value <- cor(compare_df$Gene_Jaccard, compare_df$Pair_Jaccard, use = "complete.obs")

p3 <- ggplot(compare_df, aes(x = Gene_Jaccard, y = Pair_Jaccard)) +
  geom_point(alpha = 0.5, size = 2, color = "#145C75") +
  geom_smooth(method = "lm", se = TRUE, color = "#8B1E3F", fill = "#8B1E3F", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text",
    x = 0.2, y = 0.9,
    label = paste("Pearson r =", round(cor_value, 3)),
    size = 4.5, fontface = "bold"
  ) +
  labs(
    title = "Correlation between gene set and SL pair Jaccard similarity",
    subtitle = "Each point represents a pair of cancer types",
    x = "Gene set Jaccard similarity",
    y = "Conditional SL pair Jaccard similarity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30")
  )

ggsave(file.path(gene_stat_dir, "Jaccard_Correlation.pdf"), p3, width = 7, height = 6)

# =========================================================
# 11. Save all results
# =========================================================
write.csv(stats, file.path(gene_stat_dir, "Cancer_Stats.csv"), row.names = FALSE)
write.csv(jaccard_gene, file.path(gene_stat_dir, "Jaccard_GeneSet_Matrix.csv"))
write.csv(jaccard_pair, file.path(gene_stat_dir, "Conditional_Jaccard_SLpairs_Matrix.csv"))

# Save session info
sink(file.path(gene_stat_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n========================================\n")
cat("All analyses completed successfully!\n")
cat("Results saved in:", gene_stat_dir, "\n")
cat("Valid cancer types:", n, "\n")
cat("========================================\n")
