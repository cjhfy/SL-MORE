library(data.table)
library(qs)
library(arrow)
library(parallel) # Load the parallel computing library

# ─────────────────────────────────────────────
# 1. Basic configuration and paths
# ─────────────────────────────────────────────
cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

base_out_dir <- "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/"
pan_out_dir <- file.path(base_out_dir, "PanCancer_Shared")
spec_out_dir <- file.path(base_out_dir, "Cancer_Specific")
if (!dir.exists(pan_out_dir)) dir.create(pan_out_dir, recursive = TRUE)
if (!dir.exists(spec_out_dir)) dir.create(spec_out_dir, recursive = TRUE)

base_paths <- list(
  sof_dir = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/DiffSOF/",
  DiffExp_dir = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/DiffExp/",
  dir_ExcluISLE = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/mutexISLE/",
  dir_ExcluCNVdel = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/ExclusiveScnaDel/",
  dir_ExcluMut = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/ExclusiveMut/",
  dir_ExcluAlt = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/ExclusiveAlt/",
  dir_sili = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/surv_R/SiLi_mRNA/",
  dir_ISLESCNA = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/surv_R/ISLE_CNV/",
  dir_ISLEmRNA = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/surv_R/ISLE_mRNA/",
  dir_coExp = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/coExp/"
)

crispr_path_map <- list(
  Exp_Inact_CR_Dep = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/CRISPR/ISLE_mRNA1/Pancancer.RData",
  Exp_Inact_RNAi_Dep = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/RNAi/ISLE_mRNA/pancancer.RData",
  CNV_Loss_CR_Dep_ISLE = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/CRISPR/ISLE_SCNA1/Pancancer.RData",
  CNV_Loss_CR_Dep_DAISY = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/CRISPR/DAISY/pancancer.RData",
  CNV_Loss_RNAi_Dep_ISLE = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/RNAi/ISLE_SCNA/pancancer.RData",
  CNV_Loss_RNAi_Dep_DAISY = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/RNAi/DAISY/pancancer.RData"
)

# ─────────────────────────────────────────────
# 2. Memory-optimized helper functions
# ─────────────────────────────────────────────

# [Key optimization 1]: drop string concatenation entirely; generate integer gene1/gene2 directly
split_rownames_to_genes <- function(rn) {
  parts <- tstrsplit(rn, ";", fixed = TRUE)
  list(gene1 = as.integer(parts[[1]]), gene2 = as.integer(parts[[2]]))
}

replace_sentinel <- function(vec, sentinel = 1000) {
  vec[vec == sentinel] <- NA
  vec
}

# Optimized CRISPR parsing: downcast early to reduce memory use
parse_crispr <- function(obj, feat_name, drop_col4 = FALSE, is_pval = TRUE) {
  if (drop_col4 && ncol(obj) >= 4) obj <- obj[, -4]
  vals <- replace_sentinel(obj[, 3])

  # Trim precision up front
  if (is_pval) vals <- signif(vals, 3) else vals <- round(vals, 3)

  dt <- data.table(
    gene1 = as.integer(obj[, 1]),
    gene2 = as.integer(obj[, 2]),
    value = vals
  )
  setnames(dt, "value", feat_name)
  return(dt)
}

# Generic function wrapping loading and parsing
safe_load_and_parse <- function(path, var_name, parser) {
  if (!file.exists(path)) {
    return(NULL)
  }
  env <- new.env(parent = emptyenv())
  tryCatch(
    {
      load(path, envir = env)
      if (!exists(var_name, envir = env)) {
        return(NULL)
      }
      obj <- get(var_name, envir = env)
      dt <- parser(obj)
      rm(list = ls(envir = env), envir = env)
      return(dt)
    },
    error = function(e) {
      NULL
    }
  )
}

# Preprocessing: build the protein-coding gene whitelist
gene_anno <- fread("/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info")
protein_coding_ids <- sort(unique(gene_anno[
  `#tax_id` == 9606 & type_of_gene == "protein-coding",
  as.integer(GeneID)
]))
cat(sprintf("提取到人类蛋白质编码基因: %d 个\n", length(protein_coding_ids)))
rm(gene_anno)
gc()

# [Key optimization 2]: streaming merge to reduce peak memory
build_features_streaming <- function(feature_configs, valid_genes) {
  master_dt <- data.table(gene1 = integer(), gene2 = integer())

  for (i in seq_along(feature_configs)) {
    cfg <- feature_configs[[i]]

    # 1. Read each table and immediately downcast types and parse
    dt_current <- safe_load_and_parse(cfg$path, cfg$var, cfg$parser)
    if (is.null(dt_current) || nrow(dt_current) == 0) next

    # 2. Immediately filter by whitelist to discard unused rows
    dt_current <- dt_current[gene1 %in% valid_genes & gene2 %in% valid_genes]
    if (nrow(dt_current) == 0) {
      rm(dt_current)
      gc()
      next
    }

    # 3. Full outer join with the main table
    if (nrow(master_dt) == 0) {
      master_dt <- dt_current
    } else {
      master_dt <- merge(master_dt, dt_current, by = c("gene1", "gene2"), all = TRUE)
    }

    # 4. Free each table immediately and trigger GC to keep memory in check
    rm(dt_current)
    gc()
  }

  if (nrow(master_dt) == 0) {
    return(NULL)
  }

  # Zero-fill specific columns
  zero_fill_cols <- intersect(c("Direct_Int", "Paral_Rel", "PCom_CoM"), names(master_dt))
  for (col in zero_fill_cols) {
    master_dt[is.na(get(col)), (col) := 0L]
  }

  setkey(master_dt, gene1, gene2)
  return(master_dt)
}

optimized_write_parquet <- function(dt, file_path) {
  write_parquet(
    dt, file_path,
    compression = "zstd", compression_level = 8,
    chunk_size = 500000, write_statistics = TRUE, use_dictionary = TRUE
  )
}

# ─────────────────────────────────────────────
# 3. Process the pan-cancer features (streaming version)
# ─────────────────────────────────────────────
cat("\n>>> [步骤 1] 构建泛癌种共享特征库 (流式计算中...) <<<\n")

pan_configs <- list(
  list(path = crispr_path_map$Exp_Inact_CR_Dep, var = "mRNAEssCrispr", parser = function(obj) parse_crispr(obj, "Exp_Low_CR_Dep", is_pval = TRUE)),
  list(path = crispr_path_map$Exp_Inact_RNAi_Dep, var = "mRNaEssRNAi", parser = function(obj) parse_crispr(obj, "Exp_Low_RNAi_Dep", is_pval = TRUE)),
  list(path = crispr_path_map$CNV_Loss_CR_Dep_ISLE, var = "scnaEssCrispr", parser = function(obj) parse_crispr(obj, "CNV_Loss_CR_Dep", is_pval = TRUE)),
  list(path = crispr_path_map$CNV_Loss_CR_Dep_DAISY, var = "scnaEssCrisprDaisy", parser = function(obj) parse_crispr(obj, "MutDelLow_CR_Dep", drop_col4 = TRUE, is_pval = TRUE)),
  list(path = crispr_path_map$CNV_Loss_RNAi_Dep_ISLE, var = "scnaEssRNAi", parser = function(obj) parse_crispr(obj, "CNV_Loss_RNAi_Dep", is_pval = TRUE)),
  list(path = crispr_path_map$CNV_Loss_RNAi_Dep_DAISY, var = "scnaEssRNAiDaisy", parser = function(obj) parse_crispr(obj, "MutDelLow_RNAi_Dep", drop_col4 = TRUE, is_pval = TRUE)),

  # Downcast early: integers via as.integer, doubles via round/signif
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/Complex_10.RData", var = "Complex_10", parser = function(obj) {
    g <- split_rownames_to_genes(rownames(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, PCom_CoM = as.integer(obj[, 3]))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/ComplexEssScore.RData", var = "ComplexEssScore", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, PCom_Agg_Dep = round(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/geneConservationScore.RData", var = "geneConservationScore", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, Conservation_score = as.integer(obj))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/GO_score/GOscore.RData", var = "GOscore", parser = function(obj) {
    g <- split_rownames_to_genes(obj[, 1])
    data.table(gene1 = g$gene1, gene2 = g$gene2, GO_Sim = round(as.numeric(obj[, 4]), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/jacard_colocalisation.RData", var = "jacard_colocalisation", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, SubLoc = round(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/Jacard_Domain.RData", var = "Jacard_Domain", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, Domain = round(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/Mean_age.RData", var = "Mean_age", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, Gene_Age = round(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/sixFeatures0.01/GeneFunSim/paralog.RData", var = "paralog", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    unique(data.table(gene1 = g$gene1, gene2 = g$gene2, Paral_Rel = as.integer(obj)), by = c("gene1", "gene2"))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/pathway_p.RData", var = "pathway_p", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, Path_CoM = signif(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/PPI_share.RData", var = "PPI_share", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, PPI_Ovlp = signif(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/PPI_union.RData", var = "PPI_union", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, PPI_Union = as.integer(obj))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/ppiEssScore.RData", var = "ppiEssScore", parser = function(obj) {
    g <- split_rownames_to_genes(names(obj))
    data.table(gene1 = g$gene1, gene2 = g$gene2, PPI_Agg_Dep = round(as.numeric(obj), 3))
  }),
  list(path = "/data/home/chenjiahao/nuaa/synlethDB/project_source/PPI/proteinInt.RData", var = "proteinInt", parser = function(obj) {
    df <- obj[[2]]
    swap_idx <- df[, 1] > df[, 2]
    df[swap_idx, c(1, 2)] <- df[swap_idx, c(2, 1)]
    unique(data.table(gene1 = as.integer(df[, 1]), gene2 = as.integer(df[, 2]), Direct_Int = 1L), by = c("gene1", "gene2"))
  })
)

PanCancer_Shared_Matrix <- build_features_streaming(pan_configs, protein_coding_ids)
optimized_write_parquet(PanCancer_Shared_Matrix, file.path(pan_out_dir, "PanCancer_Shared_Features.parquet"))
cat(sprintf("泛癌种特征保存成功，共 %d 行。\n", nrow(PanCancer_Shared_Matrix)))
rm(PanCancer_Shared_Matrix, pan_configs)
gc()

# ─────────────────────────────────────────────────────────────
# 4. Loop over the 32 cancer types sequentially
# ─────────────────────────────────────────────────────────────
cat("\n>>> [步骤 2] 批量构建癌种特异性特征库 (顺序处理中...) <<<\n")

for (cancer in cancers) {
  cat(sprintf("\n[开始处理] %s ...\n", cancer))

  spec_configs <- list(
    list(path = paste0(base_paths$DiffExp_dir, "TCGA-", cancer, ".RData"), var = "diffExp_p", parser = function(obj) {
      g <- split_rownames_to_genes(names(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, Exp_Comp = signif(as.numeric(obj), 3))
    }),
    list(path = paste0(base_paths$dir_ExcluISLE, "TCGA-", cancer, ".RData"), var = "mutexISLE", parser = function(obj) {
      g <- split_rownames_to_genes(rownames(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, Exp_Excl = signif(obj[, 1], 3), CNV_Excl_ISLE = signif(obj[, 2], 3))
    }),
    list(path = paste0(base_paths$sof_dir, "TCGA-", cancer, ".RData"), var = "SCNV_diff_p", parser = function(obj) {
      g <- split_rownames_to_genes(names(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, CNV_Comp = signif(as.numeric(obj), 3))
    }),
    list(path = paste0(base_paths$dir_ExcluCNVdel, "TCGA-", cancer, ".RData"), var = "scnaDel_res", parser = function(obj) {
      g <- split_rownames_to_genes(names(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, CNV_Excl_DSL = signif(as.numeric(obj), 3))
    }),
    list(path = paste0(base_paths$dir_ExcluMut, "TCGA-", cancer, ".RData"), var = "mut_res", parser = function(obj) {
      g <- split_rownames_to_genes(names(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, Mut_Excl = signif(as.numeric(obj), 3))
    }),
    list(path = paste0(base_paths$dir_ExcluAlt, "TCGA-", cancer, ".RData"), var = "scnaAlt", parser = function(obj) {
      g <- split_rownames_to_genes(names(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, MutCNV_Excl = signif(as.numeric(obj), 3))
    }),
    list(path = paste0(base_paths$dir_sili, cancer, "/", cancer, "_final.RData"), var = "final_results", parser = function(obj) {
      g <- split_rownames_to_genes(rownames(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, Exp_SB_SiLi = signif(obj[, "P_value"], 3), Exp_SB_SiLi_HR = round(obj[, "HR"], 3))
    }),
    list(path = paste0(base_paths$dir_ISLESCNA, cancer, "/", cancer, "_final.RData"), var = "final_results", parser = function(obj) {
      g <- split_rownames_to_genes(rownames(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, CNV_SB = signif(obj[, "P_value"], 3), CNV_SB_HR = round(obj[, "HR"], 3))
    }),
    list(path = paste0(base_paths$dir_ISLEmRNA, cancer, "/", cancer, "_final.RData"), var = "final_results", parser = function(obj) {
      g <- split_rownames_to_genes(rownames(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, Exp_SB_ISLE = signif(obj[, "P_value"], 3), Exp_SB_ISLE_HR = round(obj[, "HR"], 3))
    }),
    list(path = paste0(base_paths$dir_coExp, "TCGA-", cancer, ".RData"), var = "co_EXP_p", parser = function(obj) {
      g <- split_rownames_to_genes(rownames(obj))
      data.table(gene1 = g$gene1, gene2 = g$gene2, CoExp = signif(obj[, 1], 3), CoExp_r = round(obj[, 2], 3))
    })
  )

  Cancer_Spec_Matrix <- build_features_streaming(spec_configs, protein_coding_ids)

  if (!is.null(Cancer_Spec_Matrix)) {
    optimized_write_parquet(
      Cancer_Spec_Matrix,
      file.path(spec_out_dir, paste0(cancer, "_Specific_Features.parquet"))
    )
    cat(sprintf("[完成] %s | 行数: %d\n", cancer, nrow(Cancer_Spec_Matrix)))
  } else {
    cat(sprintf("[警告] %s 未能生成特征矩阵\n", cancer))
  }

  rm(spec_configs, Cancer_Spec_Matrix)
  gc()
}

cat("\n所有特征顺序处理与保存完毕！\n")
