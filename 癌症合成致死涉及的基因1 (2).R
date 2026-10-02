# =========================================================
# Cancer-specific SL landscape analysis
# High-quality SCI-style visualization
# =========================================================

library(data.table)
library(arrow)
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)

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
# 3. Read data
# ─────────────────────────────────────────────
cancer_genes <- list()
cancer_pairs <- list()
stats <- data.frame(Cancer = cancers, Gene_Count = 0, Pair_Count = 0)

for (cancer in cancers) {
  file_path <- file.path(gt4_dir, paste0(cancer, "_GT4_Features.parquet"))
  if (file.exists(file_path)) {
    sl_data <- read_parquet(file_path)
    genes <- unique(c(sl_data$Symbol1, sl_data$Symbol2))
    pairs <- apply(
      sl_data[, c("Symbol1", "Symbol2")], 1,
      function(x) paste(sort(x), collapse = "_")
    )
    cancer_genes[[cancer]] <- genes
    cancer_pairs[[cancer]] <- unique(pairs)
    stats[stats$Cancer == cancer, c("Gene_Count", "Pair_Count")] <- c(length(genes), length(unique(pairs)))
  }
}

# ─────────────────────────────────────────────
# 4. Order cancer types
# ─────────────────────────────────────────────
stats <- stats %>% arrange(desc(Gene_Count))
stats$Cancer <- factor(stats$Cancer, levels = stats$Cancer)

# =========================================================
# 5. Gene count barplot
# =========================================================
p1 <- ggplot(stats, aes(x = Cancer, y = Gene_Count, fill = Gene_Count)) +
  geom_bar(stat = "identity", width = 0.78, color = "black", linewidth = 0.28) +
  geom_hline(yintercept = total_protein_coding, linetype = "dashed", color = "#8B1E3F", linewidth = 0.9) +
  annotate("text",
    x = length(cancers) - 1, y = total_protein_coding,
    label = paste0("Protein-coding genes: ", comma(total_protein_coding)),
    color = "#8B1E3F", hjust = 1, vjust = -0.6, size = 3.8, fontface = "italic"
  ) +
  geom_text(aes(label = comma(Gene_Count)), vjust = -0.45, size = 3) +
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

ggsave(file.path(gene_stat_dir, "Gene_Count_Barplot.pdf"), p1, width = 14, height = 7)
ggsave(file.path(gene_stat_dir, "Gene_Count_Barplot.png"), p1, width = 14, height = 7, dpi = 600)

# =========================================================
# 6. SL pair count barplot
# =========================================================
p2 <- ggplot(stats, aes(x = Cancer, y = Pair_Count, fill = Pair_Count)) +
  geom_bar(stat = "identity", width = 0.78, color = "black", linewidth = 0.28) +
  geom_text(aes(label = comma(Pair_Count)), vjust = -0.45, size = 3) +
  scale_fill_gradient(low = "#D8F0F0", high = "#145C75") +
  labs(title = "SL interaction pairs across cancer types", x = "Cancer type", y = "SL pair count") +
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
valid_cancers <- names(cancer_genes)
n <- length(valid_cancers)

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

jaccard_gene_melt <- melt(jaccard_gene)
colnames(jaccard_gene_melt) <- c("Cancer1", "Cancer2", "Jaccard")

p3 <- ggplot(jaccard_gene_melt, aes(x = Cancer1, y = Cancer2, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)),
    size = 2.4,
    color = ifelse(jaccard_gene_melt$Jaccard > 0.5, "white", "black")
  ) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"), limits = c(0, 1)) +
  labs(title = "Jaccard similarity of SL gene sets", x = "", y = "") +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 8),
    panel.grid = element_blank(), legend.title = element_blank()
  )

ggsave(file.path(gene_stat_dir, "Jaccard_GeneSet_Heatmap.pdf"), p3, width = 14, height = 12)
ggsave(file.path(gene_stat_dir, "Jaccard_GeneSet_Heatmap.png"), p3, width = 14, height = 12, dpi = 600)

# =========================================================
# 8. Jaccard similarity — SL pairs
# =========================================================
jaccard_pair <- matrix(0, n, n, dimnames = list(valid_cancers, valid_cancers))
for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      jaccard_pair[i, j] <- 1
    } else {
      inter <- length(intersect(cancer_pairs[[i]], cancer_pairs[[j]]))
      union_sz <- length(union(cancer_pairs[[i]], cancer_pairs[[j]]))
      jaccard_pair[i, j] <- inter / union_sz
    }
  }
}

jaccard_pair_melt <- melt(jaccard_pair)
colnames(jaccard_pair_melt) <- c("Cancer1", "Cancer2", "Jaccard")

p4 <- ggplot(jaccard_pair_melt, aes(x = Cancer1, y = Cancer2, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Jaccard)),
    size = 2.4,
    color = ifelse(jaccard_pair_melt$Jaccard > 0.5, "white", "black")
  ) +
  scale_fill_gradientn(colours = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"), limits = c(0, 1)) +
  labs(title = "Jaccard similarity of SL interaction pairs", x = "", y = "") +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 8),
    panel.grid = element_blank(), legend.title = element_blank()
  )

ggsave(file.path(gene_stat_dir, "Jaccard_SLpairs_Heatmap.pdf"), p4, width = 14, height = 12)
ggsave(file.path(gene_stat_dir, "Jaccard_SLpairs_Heatmap.png"), p4, width = 14, height = 12, dpi = 600)

# =========================================================
# 9. Save results
# =========================================================
write.csv(stats, file.path(gene_stat_dir, "Cancer_Stats.csv"), row.names = FALSE)
write.csv(jaccard_gene, file.path(gene_stat_dir, "Jaccard_GeneSet_Matrix.csv"))
write.csv(jaccard_pair, file.path(gene_stat_dir, "Jaccard_SLpairs_Matrix.csv"))

cat("All analyses completed successfully.\nResults saved in:", gene_stat_dir, "\n")
