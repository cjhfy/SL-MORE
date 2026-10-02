############################################################
# Synthetic Lethality Feature Analysis (Big Data Edition)
# Core approach: disk-level filtering (Arrow) plus in-memory analysis (data.table)
############################################################

library(arrow)
library(dplyr)
library(data.table)
library(tidyverse)
library(pheatmap)
library(ggplot2)

### === Step 0: Load libraries ===
cancer <- "COAD"
load("/data/home/chenjiahao/nuaa/SL_norm/sl.golden.set.RData")
sl_ISLE <- gd$sr0
# sl_ISLE<-sl_ISLE[!duplicated(paste(sl_ISLE[,1],sl_ISLE[,2],sep = "\t")),]
type <- gd$dat
sl_ISLE <- sl_ISLE[which(type[, 2] == cancer), ]
gene_anno <- data.table::fread(file = "/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info", header = T, sep = "\t", stringsAsFactors = F, data.table = F)
gene_anno <- gene_anno[which(gene_anno$type_of_gene == "protein-coding"), ]
gene_anno <- gene_anno[, c(2, 3, 5)]
gene_anno[, 3] <- gsub("\\|", ";", gene_anno[, 3])
genelist <- unique(c(sl_ISLE[, 1], sl_ISLE[, 2]))
gene_anno1 <- gene_anno
gene_anno1 <- gene_anno1[gene_anno1[, 2] %in% genelist, ]
rownames(gene_anno1) <- gene_anno1[, 2]
####
index1 <- c()
gene_anno2 <- data.frame(GeneID = NULL, Symbol = NULL)
genelist1 <- genelist[!genelist %in% gene_anno1[, 2]]
for (i in 1:length(genelist1)) {
  index <- unique(c(
    grep(paste(";", genelist1[i], ";", sep = ""), gene_anno[, 3], ignore.case = T),
    grep(paste("^", genelist1[i], ";", sep = ""), gene_anno[, 3], ignore.case = T),
    grep(paste(";", genelist1[i], "$", sep = ""), gene_anno[, 3], ignore.case = T),
    grep(paste("^", genelist1[i], "$", sep = ""), gene_anno[, 3], ignore.case = T)
  ))
  if (length(index) == 1) {
    gene_anno2[i, 1] <- gene_anno[index, 1]
    gene_anno2[i, 2] <- genelist1[i]
  }
  # Fix: check length(index) != 1, not length(res) != 1
  if (length(index) != 1) {
    index1 <- c(index1, i)
  }
}
gene_anno2 <- na.omit(gene_anno2)
rownames(gene_anno2) <- gene_anno2[, 2]
colnames(gene_anno2) <- c("GeneID", "Symbol")
gene_anno1 <- rbind(gene_anno1[, 1:2], gene_anno2)
###
sl_ISLE <- sl_ISLE[intersect(which(sl_ISLE[, 1] %in% gene_anno1[, 2] == T), which(sl_ISLE[, 2] %in% gene_anno1[, 2] == T)), ]
sl_ISLE_entrz <- cbind(gene_anno1[sl_ISLE[, 1], 1], gene_anno1[sl_ISLE[, 2], 1])
inx <- which(sl_ISLE_entrz[, 1] > sl_ISLE_entrz[, 2])
sl_ISLE_entrz[inx, ] <- sl_ISLE_entrz[inx, c(2, 1)]
sl_ISLE_entrz1 <- unique(paste(sl_ISLE_entrz[, 1], sl_ISLE_entrz[, 2], sep = ";"))
sl_list <- sl_ISLE_entrz1
# ─────────────────────────────────────────────
# 1. Basic configuration
# ─────────────────────────────────────────────
cancer_type <- "COAD"
base_dir <- "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/"
pan_path <- file.path(base_dir, "PanCancer_Shared/PanCancer_Shared_Features.parquet")
spec_path <- file.path(base_dir, "Cancer_Specific", paste0(cancer_type, "_Specific_Features.parquet"))

# Assume sl_list already exists in the environment
# Example format: sl_list <- c("9700;80198", "7189;80198", ...)

# ─────────────────────────────────────────────
# 2. Disk-level streaming filter (avoids out-of-memory on the 200M-row merge)
# ─────────────────────────────────────────────
cat(">>> 正在执行磁盘级流式过滤...\n")

# A. Prepare the filter list
sl_pairs <- data.table(pair = sl_list) %>%
  separate(pair, into = c("gene1", "gene2"), sep = ";", convert = TRUE) %>%
  mutate(gene1 = as.integer(gene1), gene2 = as.integer(gene2))

# B. Create an Arrow scan object (no memory cost)
pan_ds <- open_dataset(pan_path)
spec_ds <- open_dataset(spec_path)

# C. Run the semi-join (keep only gene pairs present in sl_list)
# Arrow performs this during the scan, so memory overhead is minimal
matched_pan <- pan_ds %>%
  semi_join(sl_pairs, by = c("gene1", "gene2")) %>%
  collect() %>%
  as.data.table()

matched_spec <- spec_ds %>%
  semi_join(sl_pairs, by = c("gene1", "gene2")) %>%
  collect() %>%
  as.data.table()

cat(">>> 正在进行 Key-based 内存合并...\n")

# A. Coerce types to avoid int-vs-double match failures
matched_pan[, `:=`(gene1 = as.integer(gene1), gene2 = as.integer(gene2))]
matched_spec[, `:=`(gene1 = as.integer(gene1), gene2 = as.integer(gene2))]

# B. Set the key: this sorts and builds a binary index, fixing mismatched joins
setkey(matched_pan, gene1, gene2)
setkey(matched_spec, gene1, gene2)

# C. Run the inner join
# This joins on the key and keeps only the 3068 rows present on both sides
FeatureMat_raw <- merge(matched_pan, matched_spec, by = c("gene1", "gene2"), all = FALSE)

cat(sprintf(
  ">>> 合并完成。最终特征矩阵维度: %d 行 x %d 列\n",
  nrow(FeatureMat_raw), ncol(FeatureMat_raw)
))

cat(sprintf(">>> 过滤完成。原始行数: ~2亿 | 提取后行数: %d\n", nrow(FeatureMat_raw)))

# ─────────────────────────────────────────────
# 3. Clean and preprocess the feature matrix
# ─────────────────────────────────────────────

# Select feature columns, excluding IDs and the survival HR columns
all_cols <- names(FeatureMat_raw)
feature_cols <- setdiff(all_cols, c("gene1", "gene2", grep("_HR$", all_cols, value = TRUE)))

FeatureMat <- as.data.frame(FeatureMat_raw[, ..feature_cols])

# A. Handle P-value features (convert to -log10)
# Match the naming used in the streaming script
p_val_cols <- c(
  "Exp_Comp", "CNV_Comp", "Exp_Excl", "CNV_Excl_ISLE", "CNV_Excl_DSL",
  "Mut_Excl", "MutCNV_Excl", "Exp_SB_SiLi", "CNV_SB", "Exp_SB_ISLE", "CoExp"
)
p_val_cols <- intersect(p_val_cols, names(FeatureMat))

for (col in p_val_cols) {
  FeatureMat[[col]] <- -log10(as.numeric(FeatureMat[[col]]) + 1e-10)
}

# B. Impute missing values (median)
for (i in 1:ncol(FeatureMat)) {
  if (any(is.na(FeatureMat[, i]))) {
    FeatureMat[is.na(FeatureMat[, i]), i] <- median(FeatureMat[, i], na.rm = TRUE)
  }
}

# C. Normalize (and drop all-zero columns)
keep_cols <- colSums(abs(FeatureMat), na.rm = TRUE) > 0
FeatureMat_clean <- FeatureMat[, keep_cols]
FeatureMat_scaled <- as.data.frame(scale(FeatureMat_clean))

# ─────────────────────────────────────────────
# 4. Correlation analysis and visualization (heatmap)
# ─────────────────────────────────────────────
cor_mat <- cor(FeatureMat_scaled, method = "spearman")
plot_mat <- cor_mat
diag(plot_mat) <- NA # Hide the diagonal for better contrast

# Define the feature categories (mapping)
featureCategories <- list(
  Mutation_Exclusion = c("Exp_Excl", "CNV_Excl_ISLE", "CNV_Excl_DSL", "Mut_Excl", "MutCNV_Excl"),
  Dependency = c(
    "Exp_Inact_CR_Dep", "Exp_Inact_RNAi_Dep", "scnaEssCrispr",
    "scnaEssCrisprDaisy", "scnaEssRNAi", "scnaEssRNAiDaisy", "PPI_Agg_Dep", "PCom_Agg_Dep"
  ),
  Functional_Sim = c("CoExp", "GO_Sim", "Paral_Rel", "SubLoc", "Domain", "Path_CoM"),
  Interaction = c("Direct_Int", "PCom_CoM", "PPI_Ovlp", "PPI_Union"),
  Evolution = c("Conservation_score", "Gene_Age")
)

# Build the annotation data frame
feature_group <- setNames(rep("Other", ncol(cor_mat)), colnames(cor_mat))
for (cat_name in names(featureCategories)) {
  common_feats <- intersect(featureCategories[[cat_name]], names(feature_group))
  feature_group[common_feats] <- cat_name
}

annotation_df <- data.frame(Category = feature_group, row.names = colnames(cor_mat))

# Plot the heatmap
pheatmap(
  plot_mat,
  clustering_distance_rows = as.dist(1 - cor_mat),
  clustering_distance_cols = as.dist(1 - cor_mat),
  clustering_method = "ward.D2",
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  na_col = "white",
  border_color = NA,
  main = paste(cancer_type, "SL Feature Synergy Matrix")
)

# ─────────────────────────────────────────────
# 5. Dimensionality reduction of the feature space (PCA)
# ─────────────────────────────────────────────
# We care about correlations between features, so run PCA on the transposed matrix
pca_res <- prcomp(t(FeatureMat_scaled), scale. = TRUE)
pca_var <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)

pca_plot_df <- data.frame(
  Feature = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Category = feature_group[rownames(pca_res$x)]
)

ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Category)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Feature), vjust = -1, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "PCA of Synthetic Lethality Features",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  )
