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
FeatureBin <- SlFeatures
colSums(is.na(FeatureBin))
FeatureBin[is.na(FeatureBin[, "jacard_colocalisation"]), "jacard_colocalisation"] <- 0
FeatureBin[is.na(FeatureBin[, "Jacard_Domain"]), "Jacard_Domain"] <- 0
############################################################
# Step 1 Binarize features
############################################################
### p-value evidence
p_cols <- c(
  "DiffExp", "DiffSof",
  "ExcluISLEmRNA", "ExcluISLECNV", "ExcluCNVdel", "ExcluMut", "ExcluAlt",
  "mRNAEssCrispr", "mRNaEssRNAi",
  "scnaEssCrispr", "scnaEssCrisprDaisy",
  "scnaEssRNAi", "scnaEssRNAiDaisy",
  "pathway_p", "PPI_share"
)
FeatureBin[, p_cols] <- FeatureBin[, p_cols] < 0.01
### co-expression
FeatureBin <- FeatureBin %>%
  mutate(coExp = coExp_p < 0.01 & coExp_r > 0.5)
# Step 2: drop the coExp_r column
FeatureBin <- FeatureBin %>%
  select(-c(coExp_r, coExp_p))
### sur
FeatureBin <- FeatureBin %>%
  mutate(SurvISLEmRNA = SurvISLEmRNA_p < 0.01 & SurvISLEmRNA_HR < 1) %>%
  select(-c(SurvISLEmRNA_HR, SurvISLEmRNA_p))
FeatureBin <- FeatureBin %>%
  mutate(SurvISLECNV = SurvISLECNV_p < 0.01 & SurvISLECNV_HR < 1) %>%
  select(-c(SurvISLECNV_HR, SurvISLECNV_p))
FeatureBin <- FeatureBin %>%
  mutate(SurvSILI = SurvSILI_p < 0.01 & SurvSILI_HR < 1) %>%
  select(-c(SurvSILI_HR, SurvSILI_p))
### score features
FeatureBin[, "ComplexEssScore"] <- FeatureBin[, "ComplexEssScore"] > 0.9983349
FeatureBin[, "ppiEssScore"] <- FeatureBin[, "ppiEssScore"] > 0.9990749
FeatureBin[, "GOscore"] <- FeatureBin[, "GOscore"] > 0.9165151
FeatureBin[, "jacard_colocalisation"] <- FeatureBin[, "jacard_colocalisation"] >= 1
FeatureBin[, "Jacard_Domain"] <- FeatureBin[, "Jacard_Domain"] > 0
FeatureBin[, "geneConservationScore"] <- FeatureBin[, "geneConservationScore"] >= 199
FeatureBin[, "Mean_age"] <- FeatureBin[, "Mean_age"] > 2555.0
FeatureBin[, "PPI_union"] <- FeatureBin[, "PPI_union"] > 812.0
# Convert all logical columns to 0/1
FeatureBin <- FeatureBin %>%
  mutate(across(where(is.logical), as.integer))
#### Integrate the six feature classes
featureCategories <- list(
  GeneSoF = c(
    "DiffExp", "DiffSof", "ExcluISLEmRNA", "ExcluISLECNV", "ExcluCNVdel", "ExcluMut", "ExcluAlt"
  ),
  GeneSurv = c(
    "SurvISLEmRNA", "SurvISLECNV", "SurvSILI_p"
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
# Step 1 Convert feature → category
############################################################

CategoryMat <- data.frame(row.names = rownames(FeatureBin))

for (cat in names(featureCategories)) {
  features <- featureCategories[[cat]]

  # Keep only features present in FeatureBin
  features <- intersect(features, colnames(FeatureBin))

  if (length(features) > 0) {
    CategoryMat[[cat]] <- as.integer(rowSums(FeatureBin[, features, drop = FALSE], na.rm = TRUE) > 0)
  }
}

############################################################
# Step 2 Count supporting categories
############################################################

category_count <- rowSums(CategoryMat)
library(dplyr)

df <- data.frame(
  CategoryCount = category_count
)

df2 <- df %>%
  group_by(CategoryCount) %>%
  summarise(Count = n())

mean_cat <- mean(category_count)

ggplot(df2, aes(CategoryCount, Count)) +
  geom_bar(
    stat = "identity",
    fill = "#377EB8",
    color = "black"
  ) +
  geom_vline(
    xintercept = mean_cat,
    linetype = "dashed",
    color = "red"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of Supporting Evidence Categories",
    y = "Number of SL Gene Pairs",
    title = "Evidence Category Support for Validated SL Interactions"
  )

#
############################################################
# Step 1 Feature → Category
############################################################
library(UpSetR)

CategoryMat <- data.frame(row.names = rownames(FeatureBin))
# Convert NA to 0
for (cat in names(featureCategories)) {
  features <- featureCategories[[cat]]

  features <- intersect(features, colnames(FeatureBin))

  if (length(features) > 0) {
    CategoryMat[[cat]] <- as.integer(rowSums(FeatureBin[, features]) > 0)
  }
}
CategoryMat[is.na(CategoryMat)] <- 0
# Adjust graphical parameters
par(mar = c(5, 5, 5, 5)) # Increase margins
upset(
  CategoryMat,
  order.by = "freq",
  nintersects = 25,
  point.size = 3.5,
  line.size = 1,
  text.scale = 1.5,
  set_size.show = TRUE # Show set sizes
)


jaccard <- function(x, y) {
  inter <- sum(x & y)
  uni <- sum(x | y)

  if (uni == 0) {
    return(0)
  } else {
    return(inter / uni)
  }
}

n <- ncol(FeatureBin)

jac_mat <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    jac_mat[i, j] <- jaccard(
      FeatureBin[, i],
      FeatureBin[, j]
    )
  }
}

colnames(jac_mat) <- colnames(FeatureBin)
rownames(jac_mat) <- colnames(FeatureBin)

jac_mat[is.na(jac_mat)] <- 0
diag(jac_mat) <- 1

# Step 3 Heatmap
pheatmap(
  jac_mat,
  clustering_method = "ward.D2",
  color = colorRampPalette(
    c("white", "#FDB863", "#E66101")
  )(100),
  main = "Feature Evidence Overlap (Jaccard)"
)
