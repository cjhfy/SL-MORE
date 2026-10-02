library(data.table)
library(arrow)
library(dplyr)
library(tidyr)

# ─────────────────────────────────────────────
# 1. Basic configuration
# ─────────────────────────────────────────────
cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

cols_main_cat <- c("GeneSoF", "GeneSurv", "GeneDep", "GeneFunSim", "GeneEvolution", "GenePPI")

cols_sub_features <- c(
  "Exp_Comp", "CNV_Comp", "Exp_Excl", "CNV_Excl_ISLE", "CNV_Excl_DSL",
  "Mut_Excl", "MutCNV_Excl", "Exp_SB_ISLE", "Exp_SB_SiLi", "CNV_SB",
  "Exp_Low_CR_Dep", "Exp_Low_RNAi_Dep", "CNV_Loss_CR_Dep",
  "CNV_Loss_RNAi_Dep", "MutDelLow_CR_Dep", "MutDelLow_RNAi_Dep",
  "PPI_Agg_Dep", "PCom_Agg_Dep", "CoExp", "GO_Sim", "Paral_Rel",
  "SubLoc", "Domain", "Path_CoM", "PCom_CoM", "PPI_Ovlp", "PPI_Union",
  "Direct_Int", "Conservation_score", "Gene_Age"
)

base_out_dir <- "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/"
pan_out_dir <- file.path(base_out_dir, "PanCancer_Shared")
spec_out_dir <- file.path(base_out_dir, "Cancer_Specific")
detail_out_dir <- file.path(base_out_dir, "Cancer_Merged_Features")
stat_out_dir <- file.path(base_out_dir, "Synergy_Stats")

if (!dir.exists(detail_out_dir)) dir.create(detail_out_dir, recursive = TRUE)
if (!dir.exists(stat_out_dir)) dir.create(stat_out_dir, recursive = TRUE)

# ─────────────────────────────────────────────
# 2. Preprocess the pan-cancer features (expensive; run only once)
# ─────────────────────────────────────────────
cat("\n>>> Loading and preprocessing PanCancer shared features (200M rows)...\n")
pan_file <- file.path(pan_out_dir, "PanCancer_Shared_Features.parquet")
# Read only the necessary columns to save memory
pan_dt <- read_parquet(pan_file, as_data_frame = FALSE) %>% as.data.table()

# Pre-fill missing values with the updated column names
pan_dt[is.na(SubLoc), SubLoc := 0]
pan_dt[is.na(Domain), Domain := 0]

# Binarize in place using the new column names
pan_dt[, `:=`(
  Exp_Low_CR_Dep = as.integer(!is.na(Exp_Low_CR_Dep) & Exp_Low_CR_Dep < 0.01),
  Exp_Low_RNAi_Dep = as.integer(!is.na(Exp_Low_RNAi_Dep) & Exp_Low_RNAi_Dep < 0.01),
  CNV_Loss_CR_Dep = as.integer(!is.na(CNV_Loss_CR_Dep) & CNV_Loss_CR_Dep < 0.01),
  MutDelLow_CR_Dep = as.integer(!is.na(MutDelLow_CR_Dep) & MutDelLow_CR_Dep < 0.01),
  CNV_Loss_RNAi_Dep = as.integer(!is.na(CNV_Loss_RNAi_Dep) & CNV_Loss_RNAi_Dep < 0.01),
  MutDelLow_RNAi_Dep = as.integer(!is.na(MutDelLow_RNAi_Dep) & MutDelLow_RNAi_Dep < 0.01),
  Path_CoM = as.integer(!is.na(Path_CoM) & Path_CoM < 0.01),
  PPI_Ovlp = as.integer(!is.na(PPI_Ovlp) & PPI_Ovlp < 0.01),
  PCom_Agg_Dep = as.integer(!is.na(PCom_Agg_Dep) & PCom_Agg_Dep > 0.998),
  PPI_Agg_Dep = as.integer(!is.na(PPI_Agg_Dep) & PPI_Agg_Dep > 0.999),
  GO_Sim = as.integer(!is.na(GO_Sim) & GO_Sim > 0.917),
  SubLoc = as.integer(SubLoc >= 1.0),
  Domain = as.integer(Domain > 0.0),
  Conservation_score = as.integer(!is.na(Conservation_score) & Conservation_score >= 199),
  Gene_Age = as.integer(!is.na(Gene_Age) & Gene_Age > 2555.0),
  PPI_Union = as.integer(!is.na(PPI_Union) & PPI_Union > 842)
)]

setkey(pan_dt, gene1, gene2)

# ─────────────────────────────────────────────
# 3. Loop over cancer types (memory-optimized strategy)
# ─────────────────────────────────────────────
results_list <- list()

for (cancer in cancers) {
  spec_file <- file.path(spec_out_dir, paste0(cancer, "_Specific_Features.parquet"))
  if (!file.exists(spec_file)) next

  cat(sprintf("\n>>> Processing cancer type: %s\n", cancer))
  spec_dt <- read_parquet(spec_file, as_data_frame = FALSE) %>% as.data.table()

  # 3.1 Binarize the per-cancer features
  spec_dt[, `:=`(
    Exp_Comp      = as.integer(!is.na(Exp_Comp) & Exp_Comp < 0.01),
    CNV_Comp      = as.integer(!is.na(CNV_Comp) & CNV_Comp < 0.01),
    Exp_Excl      = as.integer(!is.na(Exp_Excl) & Exp_Excl < 0.01),
    CNV_Excl_ISLE = as.integer(!is.na(CNV_Excl_ISLE) & CNV_Excl_ISLE < 0.01),
    CNV_Excl_DSL  = as.integer(!is.na(CNV_Excl_DSL) & CNV_Excl_DSL < 0.01),
    Mut_Excl      = as.integer(!is.na(Mut_Excl) & Mut_Excl < 0.01),
    MutCNV_Excl   = as.integer(!is.na(MutCNV_Excl) & MutCNV_Excl < 0.01),
    CoExp         = as.integer(!is.na(CoExp) & CoExp < 0.01 & !is.na(CoExp_r) & CoExp_r > 0.5),
    Exp_SB_ISLE   = as.integer(!is.na(Exp_SB_ISLE) & Exp_SB_ISLE < 0.01 & !is.na(Exp_SB_ISLE_HR) & Exp_SB_ISLE_HR < 1),
    CNV_SB        = as.integer(!is.na(CNV_SB) & CNV_SB < 0.01 & !is.na(CNV_SB_HR) & CNV_SB_HR < 1),
    Exp_SB_SiLi   = as.integer(!is.na(Exp_SB_SiLi) & Exp_SB_SiLi < 0.01 & !is.na(Exp_SB_SiLi_HR) & Exp_SB_SiLi_HR < 1)
  )]

  spec_cols <- c(
    "gene1", "gene2", "Exp_Comp", "CNV_Comp", "Exp_Excl",
    "CNV_Excl_ISLE", "CNV_Excl_DSL", "Mut_Excl", "MutCNV_Excl", "CoExp",
    "Exp_SB_ISLE", "CNV_SB", "Exp_SB_SiLi"
  )
  spec_dt <- spec_dt[, ..spec_cols]
  setkey(spec_dt, gene1, gene2)

  # 3.2 Memory-friendly merge
  all_pairs <- unique(rbind(pan_dt[, .(gene1, gene2)], spec_dt[, .(gene1, gene2)]))
  setkey(all_pairs, gene1, gene2)

  dt_merged <- pan_dt[all_pairs, on = .(gene1, gene2)]

  dt_merged[spec_dt, `:=`(
    Exp_Comp      = i.Exp_Comp,
    CNV_Comp      = i.CNV_Comp,
    Exp_Excl      = i.Exp_Excl,
    CNV_Excl_ISLE = i.CNV_Excl_ISLE,
    CNV_Excl_DSL  = i.CNV_Excl_DSL,
    Mut_Excl      = i.Mut_Excl,
    MutCNV_Excl   = i.MutCNV_Excl,
    CoExp         = i.CoExp,
    Exp_SB_ISLE   = i.Exp_SB_ISLE,
    CNV_SB        = i.CNV_SB,
    Exp_SB_SiLi   = i.Exp_SB_SiLi
  ), on = .(gene1, gene2)]

  # 3.3 Fill NAs introduced by the merge
  target_cols <- setdiff(names(dt_merged), c("gene1", "gene2"))
  for (col in target_cols) {
    set(dt_merged, i = which(is.na(dt_merged[[col]])), j = col, value = 0)
  }

  # 3.4 Compute the category scores
  dt_merged[, `:=`(
    GeneSoF = as.integer((Exp_Comp + CNV_Comp + Exp_Excl + CNV_Excl_ISLE + CNV_Excl_DSL + Mut_Excl + MutCNV_Excl) > 0),
    GeneSurv = as.integer((Exp_SB_ISLE + CNV_SB + Exp_SB_SiLi) > 0),
    GeneDep = as.integer((Exp_Low_CR_Dep + Exp_Low_RNAi_Dep + CNV_Loss_CR_Dep +
      MutDelLow_CR_Dep + CNV_Loss_RNAi_Dep + MutDelLow_RNAi_Dep +
      PPI_Agg_Dep + PCom_Agg_Dep) > 0),
    GeneFunSim = as.integer((CoExp + GO_Sim + SubLoc + Domain + Paral_Rel + Path_CoM + PCom_CoM) > 0),
    GeneEvolution = as.integer((Conservation_score + Gene_Age) > 0),
    GenePPI = as.integer((PPI_Ovlp + PPI_Union + Direct_Int) > 0)
  )]

  dt_merged[, CategoryCount := GeneSoF + GeneSurv + GeneDep + GeneFunSim + GeneEvolution + GenePPI]
  final_col_order <- c("gene1", "gene2", "CategoryCount", cols_main_cat, cols_sub_features)
  setcolorder(dt_merged, final_col_order)

  # 3.5 Write to file
  write_parquet(dt_merged, file.path(detail_out_dir, paste0(cancer, "_Full_Features.parquet")))

  # 3.6 Summary statistics
  summary_stats <- dt_merged[, .(Count = .N), by = CategoryCount][order(CategoryCount)]
  summary_stats[, Cancer := cancer]
  results_list[[cancer]] <- summary_stats

  cat(sprintf("[%s] Done. Total rows: %d\n", cancer, nrow(dt_merged)))

  rm(spec_dt, dt_merged, all_pairs)
  gc()
}

# ─────────────────────────────────────────────
# 4. Save the summary statistics table
# ─────────────────────────────────────────────
final_summary <- rbindlist(results_list)
wide_summary <- dcast(final_summary, Cancer ~ CategoryCount, value.var = "Count", fill = 0)
fwrite(wide_summary, file.path(stat_out_dir, "All_Cancer_CategoryCount_Summary.csv"))
cat("\n>>> All processing finished!\n")

# ─────────────────────────────────────────────
# 5. Examples of follow-up analysis and filtering
# ─────────────────────────────────────────────
# Example 1: filter a specific cancer type (e.g. BRCA) for analysis
brca_file <- file.path(detail_out_dir, "BRCA_Full_Features.parquet")
brca_dt <- read_parquet(brca_file, as_data_frame = FALSE) %>% as.data.table()

# Filter pairs with CategoryCount > 3 that contain a given gene (e.g. gene ID 1956)
filtered_dt1 <- brca_dt[
  CategoryCount > 3 &
    (gene1 == 1956 | gene2 == 1956)
]

# Filter pairs with GeneSoF == 1
filtered_dt2 <- brca_dt[GeneSoF == 1]
