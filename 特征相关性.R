############################################################
# Synthetic Lethality Feature Synergy Analysis Pipeline
# Author: cjh
############################################################

library(tidyverse)
library(data.table)
library(pheatmap)

############################################################
# Step 1 Load Feature Matrix
############################################################

load("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/COAD_SL_Matrix.RData")
FeatureMat <- SlFeatures


cat("Number of SL pairs:", nrow(FeatureMat), "\n")
cat("Number of features:", ncol(FeatureMat), "\n")


############################################################
# Step 2 Feature Normalization
############################################################

### p-value features → information content

p_cols <- c(
  "DiffExp", "DiffSof",
  "ExcluISLEmRNA", "ExcluISLECNV", "ExcluCNVdel", "ExcluMut", "ExcluAlt",
  "mRNAEssCrispr", "mRNaEssRNAi",
  "scnaEssCrispr", "scnaEssCrisprDaisy",
  "scnaEssRNAi", "scnaEssRNAiDaisy",
  "pathway_p", "PPI_share"
)

FeatureMat[, p_cols] <- -log10(FeatureMat[, p_cols] + 1e-10)


############################################################
# Survival features
############################################################

hr_cols <- grep("_HR$", colnames(FeatureMat), value = TRUE)

for (hr_col in hr_cols) {
  p_col <- sub("_HR$", "_p", hr_col)

  FeatureMat[, p_col] <-
    -log10(FeatureMat[, p_col] + 1e-10) *
      as.numeric(FeatureMat[, hr_col] < 1)
}

FeatureMat <- FeatureMat[, !colnames(FeatureMat) %in%
  grep("_HR$", colnames(FeatureMat), value = TRUE)]


############################################################
# Co-expression feature
############################################################

FeatureMat$coExp <- -log10(FeatureMat$coExp_p + 1e-10) *
  pmax(FeatureMat$coExp_r, 0)

FeatureMat <- FeatureMat %>% select(-coExp_p, -coExp_r)


############################################################
# Step 3 Score and Binary Features
############################################################

binary_cols <- c(
  "Complex_10",
  "paralog",
  "proteinInt"
)

score_cols <- c(
  "ComplexEssScore",
  "geneConservationScore",
  "GOscore",
  "jacard_colocalisation",
  "Jacard_Domain",
  "Mean_age",
  "PPI_union",
  "ppiEssScore"
)

FeatureMat[, score_cols] <- scale(FeatureMat[, score_cols])

FeatureMat[, binary_cols] <- lapply(FeatureMat[, binary_cols], as.numeric)


############################################################
# Step 4 Missing Value Imputation
############################################################

for (i in 1:ncol(FeatureMat)) {
  FeatureMat[is.na(FeatureMat[, i]), i] <-
    median(FeatureMat[, i], na.rm = TRUE)
}


############################################################
# Step 5 Correlation Matrix
############################################################

# Spearman is more robust
cor_mat <- cor(FeatureMat, method = "spearman")


############################################################
# Step 6 Feature Categories
############################################################

featureCategories <- list(
  GeneSoF = c(
    "DiffExp", "DiffSof", "ExcluISLEmRNA", "ExcluISLECNV", "ExcluCNVdel", "ExcluMut", "ExcluAlt"
  ),
  GeneSurv = c(
    "SurvISLEmRNA_p", "SurvISLECNV_p", "SurvSILI_p"
  ),
  GeneDep = c(
    "mRNAEssCrispr", "mRNaEssRNAi", "scnaEssCrispr", "scnaEssCrisprDaisy",
    "scnaEssRNAi", "scnaEssRNAiDaisy", "ppiEssScore", "ComplexEssScore"
  ),
  GeneFunSim = c(
    "coExp", "Complex_10", "GOscore", "jacard_colocalisation",
    "Jacard_Domain", "paralog", "pathway_p"
  ),
  GeneEvolution = c(
    "geneConservationScore", "Mean_age"
  ),
  GenePPI = c(
    "PPI_share", "PPI_union", "proteinInt"
  )
)

############################################################
# Feature → Category mapping
############################################################

feature_group <- c()

for (cat in names(featureCategories)) {
  features <- featureCategories[[cat]]

  feature_group[features] <- cat
}

# Keep the same feature order as the heatmap
feature_group <- feature_group[colnames(cor_mat)]

annotation_df <- data.frame(
  Category = feature_group
)

rownames(annotation_df) <- colnames(cor_mat)


############################################################
# Category colors
############################################################

ann_colors <- list(
  Category = c(
    GeneSoF = "#E41A1C",
    GeneSurv = "#377EB8",
    GeneDep = "#4DAF4A",
    GeneFunSim = "#984EA3",
    GeneEvolution = "#FF7F00",
    GenePPI = "#A65628"
  )
)


############################################################
# Step 7 Heatmap
############################################################

pheatmap(
  cor_mat,
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-0.4, 1, length = 101),
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  border_color = NA,
  main = "Correlation among Synthetic Lethality Features"
)


### Hide the diagonal
############################################################
# Step 5 Correlation Matrix
############################################################

# Spearman is more robust
cor_mat <- cor(FeatureMat, method = "spearman")

# Matrix used for clustering (complete)
cluster_mat <- cor_mat

# Matrix used for the heatmap (diagonal hidden)
plot_mat <- cor_mat
diag(plot_mat) <- NA


############################################################
# Step 6 Feature Categories
############################################################

featureCategories <- list(
  GeneSoF = c(
    "DiffExp", "DiffSof", "ExcluISLEmRNA", "ExcluISLECNV",
    "ExcluCNVdel", "ExcluMut", "ExcluAlt"
  ),
  GeneSurv = c(
    "SurvISLEmRNA_p", "SurvISLECNV_p", "SurvSILI_p"
  ),
  GeneDep = c(
    "mRNAEssCrispr", "mRNaEssRNAi", "scnaEssCrispr",
    "scnaEssCrisprDaisy", "scnaEssRNAi", "scnaEssRNAiDaisy",
    "ppiEssScore", "ComplexEssScore"
  ),
  GeneFunSim = c(
    "coExp", "Complex_10", "GOscore",
    "jacard_colocalisation", "Jacard_Domain",
    "paralog", "pathway_p"
  ),
  GeneEvolution = c(
    "geneConservationScore", "Mean_age"
  ),
  GenePPI = c(
    "PPI_share", "PPI_union", "proteinInt"
  )
)


############################################################
# Feature → Category mapping
############################################################

feature_group <- c()

for (cat in names(featureCategories)) {
  features <- featureCategories[[cat]]

  feature_group[features] <- cat
}

feature_group <- feature_group[colnames(cluster_mat)]

annotation_df <- data.frame(
  Category = feature_group
)

rownames(annotation_df) <- colnames(cluster_mat)


############################################################
# Category colors
############################################################

ann_colors <- list(
  Category = c(
    GeneSoF = "#E41A1C",
    GeneSurv = "#377EB8",
    GeneDep = "#4DAF4A",
    GeneFunSim = "#984EA3",
    GeneEvolution = "#FF7F00",
    GenePPI = "#A65628"
  )
)


############################################################
# Step 7 Heatmap
############################################################

pheatmap(
  plot_mat, # Display using the matrix with the diagonal hidden
  clustering_distance_rows = as.dist(1 - cluster_mat),
  clustering_distance_cols = as.dist(1 - cluster_mat),
  clustering_method = "ward.D2",
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  na_col = "white",
  main = "Feature Correlation"
)


############################################################
# Step 8 PCA of Features
############################################################

############################################################
# Step 8 PCA of Features
############################################################

library(ggplot2)

# PCA input: feature x SL pair
feature_mat <- t(FeatureMat)

# PCA
pca_res <- prcomp(feature_mat, scale. = TRUE)

# PCA variance
pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# PCA dataframe
pca_df <- data.frame(
  Feature = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Category = feature_group
)

############################################################
# PCA plot
############################################################

ggplot(pca_df, aes(PC1, PC2, color = Category)) +
  geom_point(size = 4) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PC1 (", round(pca_var[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2] * 100, 1), "%)"),
    title = "PCA of Synthetic Lethality Features"
  ) +
  scale_color_manual(values = c(
    GeneSoF = "#E41A1C",
    GeneSurv = "#377EB8",
    GeneDep = "#4DAF4A",
    GeneFunSim = "#984EA3",
    GeneEvolution = "#FF7F00",
    GenePPI = "#A65628"
  ))

# PCA
pca_res <- prcomp(feature_mat, scale. = TRUE)

# PCA variance
pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# PCA dataframe
pca_df <- data.frame(
  Feature = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Category = feature_group
)

############################################################
# PCA plot
############################################################

ggplot(pca_df, aes(PC1, PC2, color = Category)) +
  geom_point(size = 4) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PC1 (", round(pca_var[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2] * 100, 1), "%)"),
    title = "PCA of Synthetic Lethality Features"
  ) +
  scale_color_manual(values = c(
    GeneSoF = "#E41A1C",
    GeneSurv = "#377EB8",
    GeneDep = "#4DAF4A",
    GeneFunSim = "#984EA3",
    GeneEvolution = "#FF7F00",
    GenePPI = "#A65628"
  ))
