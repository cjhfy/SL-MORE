# ─────────────────────────────────────────────
# 0. Load packages
# ─────────────────────────────────────────────
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table)
library(arrow)
library(dplyr)
library(ggplot2)

# ─────────────────────────────────────────────
# 1. Basic configuration
# ─────────────────────────────────────────────
cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

feature_groups <- list(
  GeneSoF = c("Exp_Comp", "CNV_Comp", "Exp_Excl", "CNV_Excl_ISLE", "CNV_Excl_DSL", "Mut_Excl", "MutCNV_Excl"),
  GeneSurv = c("Exp_SB_ISLE", "CNV_SB", "Exp_SB_SiLi"),
  GeneDep = c(
    "Exp_Low_CR_Dep", "Exp_Low_RNAi_Dep", "CNV_Loss_CR_Dep",
    "MutDelLow_CR_Dep", "CNV_Loss_RNAi_Dep", "MutDelLow_RNAi_Dep",
    "PPI_Agg_Dep", "PCom_Agg_Dep"
  ),
  GeneFunSim = c("CoExp", "PCom_CoM", "GO_Sim", "SubLoc", "Domain", "Paral_Rel", "Path_CoM"),
  GenePPI = c("PPI_Ovlp", "PPI_Union", "Direct_Int"),
  GeneEvolution = c("Conservation_score", "Gene_Age")
)


all_features <- unlist(feature_groups)

detail_out_dir <- "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Cancer_Merged_Features"

# ─────────────────────────────────────────────
# 2. Initialize the matrix
# ─────────────────────────────────────────────
mat_raw <- matrix(0, nrow = length(cancers), ncol = length(all_features))
rownames(mat_raw) <- cancers
colnames(mat_raw) <- all_features

high_conf_sl_counts <- numeric(length(cancers))
names(high_conf_sl_counts) <- cancers

cat("\n>>> Starting data aggregation...\n")
# ─────────────────────────────────────────────
# 3. Read the parquet file and compute statistics
# ─────────────────────────────────────────────
for (cancer in cancers) {
  file_path <- file.path(detail_out_dir, paste0(cancer, "_Full_Features.parquet"))

  if (file.exists(file_path)) {
    cols_to_read <- c("CategoryCount", all_features)

    dt <- read_parquet(
      file_path,
      col_select = all_of(cols_to_read),
      as_data_frame = FALSE
    ) %>% as.data.table()

    missing_cols <- setdiff(all_features, names(dt))
    if (length(missing_cols) > 0) {
      for (col in missing_cols) dt[, (col) := 0]
    }

    mat_raw[cancer, ] <- colSums(dt[, ..all_features], na.rm = TRUE)
    high_conf_sl_counts[cancer] <- sum(dt$CategoryCount > 3, na.rm = TRUE)

    rm(dt)
    gc()
  }
}

# ─────────────────────────────────────────────
# 4. Filter and sort
# ─────────────────────────────────────────────
order_idx <- order(high_conf_sl_counts, decreasing = TRUE)
mat_plot <- mat_raw[order_idx, , drop = FALSE]
high_conf_sl_counts <- high_conf_sl_counts[order_idx]

### ------------5. fig2A-Feature Heatmap--------------------###
group_cols <- c(
  "GeneSoF" = "#E41A1C",
  "GeneSurv" = "#377EB8",
  "GeneDep" = "#4DAF4A",
  "GeneFunSim" = "#984EA3",
  "GenePPI" = "#FFFF33",
  "GeneEvolution" = "#FF7F00"
)

column_split_factor <- factor(
  rep(names(feature_groups), lengths(feature_groups)),
  levels = names(feature_groups)
)

top_anno <- HeatmapAnnotation(
  Category = column_split_factor,
  col = list(Category = group_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Color mapping
mid_val <- quantile(mat_plot[mat_plot > 0], 0.5, na.rm = TRUE)
max_val <- quantile(mat_plot[mat_plot > 0], 0.95, na.rm = TRUE)

col_fun <- colorRamp2(
  c(0, mid_val, max_val),
  c("#f7fbff", "#6baed6", "#08306b")
)

ht <- Heatmap(
  mat_plot,
  name = "Raw Count",
  column_split = column_split_factor,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  top_annotation = top_anno,
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12, angle = 45, hjust = 1),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Significant\nCounts",
    direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  row_gap = unit(1, "mm"),
  column_gap = unit(2, "mm")
)

pdf("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Synergy_Stats/Cancer_Feature_Heatmap_ONLY.pdf",
  width = 12, height = 10
)

draw(
  ht,
  column_title = "30-Feature Ensemble Atlas (Significant Counts)",
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  padding = unit(c(2, 2, 2, 10), "mm")
)

dev.off()

cat("\n>>> Heatmap saved\n")
# ─────────────────────────────────────────────
# print sl number of per cancer
# ─────────────────────────────────────────────
sl_df <- data.frame(
  Cancer = names(high_conf_sl_counts),
  Count = high_conf_sl_counts
)
fwrite(sl_df, "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Synergy_Stats/cancer_Sl_number.csv")
