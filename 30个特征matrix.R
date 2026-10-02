library(tidyverse)
library(data.table)
library(infotheo)
library(pheatmap)
library(igraph)

### === Step 0: Load libraries ===
cancer <- "LAML"
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

# Define feature names
FeatureNames <- c(
  "DiffExp", "DiffSof", "ExcluISLEmRNA", "ExcluISLECNV", "ExcluCNVdel", "ExcluMut", "ExcluAlt",
  "SurvISLEmRNA_p", "SurvISLEmRNA_HR", "SurvISLECNV_p", "SurvISLECNV_HR", "SurvSILI_p", "SurvSILI_HR",
  "coExp_p", "coExp_r", "mRNAEssCrispr", "mRNaEssRNAi", "scnaEssCrispr", "scnaEssCrisprDaisy", "scnaEssRNAi",
  "scnaEssRNAiDaisy", "Complex_10", "ComplexEssScore", "geneConservationScore", "GOscore", "jacard_colocalisation",
  "Jacard_Domain", "Mean_age", "paralog", "pathway_p", "PPI_share", "PPI_union", "ppiEssScore", "proteinInt"
)

# Define base directory paths
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

# Initialize feature matrix
SlFeatures <- matrix(NA, nrow = length(sl_list), ncol = length(FeatureNames)) %>% as.data.frame()
rownames(SlFeatures) <- sl_list
colnames(SlFeatures) <- FeatureNames

# Load and process each feature type
# 1. DiffSof
cat("    -> 正在加载 DiffSof 特征数据...\n")
load(paste0(base_paths$sof_dir, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(names(SCNV_diff_p), sl_list)
SlFeatures[sl_intersect, "DiffSof"] <- SCNV_diff_p[sl_intersect]
cat(sprintf("       已匹配 %d 个基因对的 DiffSof 数据\n", length(sl_intersect)))
# Free memory
rm(SCNV_diff_p)
gc()

# 2. DiffExp
cat("    -> 正在加载 DiffExp 特征数据...\n")
load(paste0(base_paths$DiffExp_dir, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(names(diffExp_p), sl_list)
SlFeatures[sl_intersect, "DiffExp"] <- diffExp_p[sl_intersect]
cat(sprintf("       已匹配 %d 个基因对的 DiffExp 数据\n", length(sl_intersect)))
# Free memory
rm(diffExp_p)
gc()

# 3. ExcluISLE
cat("    -> 正在加载 ExcluISLE 特征数据...\n")
load(paste0(base_paths$dir_ExcluISLE, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(rownames(mutexISLE), sl_list)
SlFeatures[sl_intersect, c("ExcluISLEmRNA", "ExcluISLECNV")] <- mutexISLE[sl_intersect, ]
cat(sprintf("       已匹配 %d 个基因对的 ExcluISLE 数据\n", length(sl_intersect)))
# Free memory
rm(mutexISLE)
gc()

# 4. ExcluCNVdel
cat("    -> 正在加载 ExcluCNVdel 特征数据...\n")
load(paste0(base_paths$dir_ExcluCNVdel, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(names(scnaDel_res), sl_list)
SlFeatures[sl_intersect, "ExcluCNVdel"] <- scnaDel_res[sl_intersect]
cat(sprintf("       已匹配 %d 个基因对的 ExcluCNVdel 数据\n", length(sl_intersect)))
# Free memory
rm(scnaDel_res)
gc()

# 5. ExcluMut
cat("    -> 正在加载 ExcluMut 特征数据...\n")
load(paste0(base_paths$dir_ExcluMut, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(names(mut_res), sl_list)
SlFeatures[sl_intersect, "ExcluMut"] <- mut_res[sl_intersect]
cat(sprintf("       已匹配 %d 个基因对的 ExcluMut 数据\n", length(sl_intersect)))
# Free memory
rm(mut_res)
gc()

# 6. ExcluAlt
cat("    -> 正在加载 ExcluAlt 特征数据...\n")
load(paste0(base_paths$dir_ExcluAlt, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(names(scnaAlt), sl_list)
SlFeatures[sl_intersect, "ExcluAlt"] <- scnaAlt[sl_intersect]
cat(sprintf("       已匹配 %d 个基因对的 ExcluAlt 数据\n", length(sl_intersect)))
# Free memory
rm(scnaAlt)
gc()

# 7. Survival features (SiLi mRNA)
cat("    -> 正在加载 SILI 特征数据...\n")
load(paste0(base_paths$dir_sili, cancer, "/", cancer, "_final.RData"))
sl_intersect <- intersect(rownames(final_results), sl_list)
SlFeatures[sl_intersect, c("SurvSILI_p", "SurvSILI_HR")] <- final_results[sl_intersect, c("P_value", "HR")]
cat(sprintf("       已匹配 %d 个基因对的 SILI 数据\n", length(sl_intersect)))
# Free memory
rm(final_results)
gc()

# 8. Survival features (ISLE SCNA)
cat("    -> 正在加载 ISLESCNA 特征数据...\n")
load(paste0(base_paths$dir_ISLESCNA, cancer, "/", cancer, "_final.RData"))
sl_intersect <- intersect(rownames(final_results), sl_list)
SlFeatures[sl_intersect, c("SurvISLECNV_p", "SurvISLECNV_HR")] <- final_results[sl_intersect, c("P_value", "HR")]
cat(sprintf("       已匹配 %d 个基因对的 ISLESCNA 数据\n", length(sl_intersect)))
# Free memory
rm(final_results)
gc()

# 9. Survival features (ISLE mRNA)
cat("    -> 正在加载 ISLEmRNA 特征数据...\n")
load(paste0(base_paths$dir_ISLEmRNA, cancer, "/", cancer, "_final.RData"))
sl_intersect <- intersect(rownames(final_results), sl_list)
SlFeatures[sl_intersect, c("SurvISLEmRNA_p", "SurvISLEmRNA_HR")] <- final_results[sl_intersect, c("P_value", "HR")]
cat(sprintf("       已匹配 %d 个基因对的 ISLEmRNA 数据\n", length(sl_intersect)))
# Free memory
rm(final_results)
gc()

# 10. Co-expression
cat("    -> 正在加载 Co-expression 特征数据...\n")
load(paste0(base_paths$dir_coExp, "TCGA-", cancer, ".RData"))
sl_intersect <- intersect(rownames(co_EXP_p), sl_list)
SlFeatures[sl_intersect, c("coExp_p", "coExp_r")] <- co_EXP_p[sl_intersect, ]
cat(sprintf("       已匹配 %d 个基因对的 Co-expression 数据\n", length(sl_intersect)))
# Free memory
rm(co_EXP_p)
gc()

### Gene knockout features (also pan-cancer; no per-cancer edits needed)
cat("    -> 正在加载 CRISPR/RNAi 特征数据...\n")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/CRISPR/DAISY/pancancer.RData")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/CRISPR/ISLE_mRNA1/Pancancer.RData")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/CRISPR/ISLE_SCNA1/Pancancer.RData")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/RNAi/DAISY/pancancer.RData")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/RNAi/ISLE_mRNA/pancancer.RData")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/RNAi/ISLE_SCNA/pancancer.RData")

# 1. Remove the 4th column of scnaEssRNAiDaisy (already present)
scnaEssRNAiDaisy <- scnaEssRNAiDaisy[, -4]
mRNAEssCrispr$pval_1[which(mRNAEssCrispr$pval_1 == 1000)] <- NA
mRNaEssRNAi$pval_1[which(mRNaEssRNAi$pval_1 == 1000)] <- NA
scnaEssCrispr$pval_1[which(scnaEssCrispr$pval_1 == 1000)] <- NA
scnaEssCrisprDaisy$pval_1[which(scnaEssCrisprDaisy$pval_1 == 1000)] <- NA
scnaEssRNAi$pval_1[which(scnaEssRNAi$pval_1 == 1000)] <- NA
scnaEssRNAiDaisy$pval_1[which(scnaEssRNAiDaisy$pval_1 == 1000)] <- NA

# mRNAEssCrispr
length(which(mRNAEssCrispr[, 1] > mRNAEssCrispr[, 2]))
rownames(mRNAEssCrispr) <- paste(mRNAEssCrispr[, 1], mRNAEssCrispr[, 2], sep = ";")
sl_intersect <- intersect(rownames(mRNAEssCrispr), sl_list)
SlFeatures[sl_intersect, "mRNAEssCrispr"] <- mRNAEssCrispr[sl_intersect, 3]

# mRNaEssRNAi
length(which(mRNaEssRNAi[, 1] > mRNaEssRNAi[, 2]))
rownames(mRNaEssRNAi) <- paste(mRNaEssRNAi[, 1], mRNaEssRNAi[, 2], sep = ";")
sl_intersect <- intersect(rownames(mRNaEssRNAi), sl_list)
SlFeatures[sl_intersect, "mRNaEssRNAi"] <- mRNaEssRNAi[sl_intersect, 3]

# scnaEssCrispr
length(which(scnaEssCrispr[, 1] > scnaEssCrispr[, 2]))
rownames(scnaEssCrispr) <- paste(scnaEssCrispr[, 1], scnaEssCrispr[, 2], sep = ";")
sl_intersect <- intersect(rownames(scnaEssCrispr), sl_list)
SlFeatures[sl_intersect, "scnaEssCrispr"] <- scnaEssCrispr[sl_intersect, 3]

# scnaEssCrisprDaisy
length(which(scnaEssCrisprDaisy[, 1] > scnaEssCrisprDaisy[, 2]))
rownames(scnaEssCrisprDaisy) <- paste(scnaEssCrisprDaisy[, 1], scnaEssCrisprDaisy[, 2], sep = ";")
sl_intersect <- intersect(rownames(scnaEssCrisprDaisy), sl_list)
SlFeatures[sl_intersect, "scnaEssCrisprDaisy"] <- scnaEssCrisprDaisy[sl_intersect, 3]

# scnaEssRNAi
length(which(scnaEssRNAi[, 1] > scnaEssRNAi[, 2]))
rownames(scnaEssRNAi) <- paste(scnaEssRNAi[, 1], scnaEssRNAi[, 2], sep = ";")
sl_intersect <- intersect(rownames(scnaEssRNAi), sl_list)
SlFeatures[sl_intersect, "scnaEssRNAi"] <- scnaEssRNAi[sl_intersect, 3]

# scnaEssRNAiDaisy
length(which(scnaEssRNAiDaisy[, 1] > scnaEssRNAiDaisy[, 2]))
rownames(scnaEssRNAiDaisy) <- paste(scnaEssRNAiDaisy[, 1], scnaEssRNAiDaisy[, 2], sep = ";")
sl_intersect <- intersect(rownames(scnaEssRNAiDaisy), sl_list)
SlFeatures[sl_intersect, "scnaEssRNAiDaisy"] <- scnaEssRNAiDaisy[sl_intersect, 3]

### Pancancer_feature (also pan-cancer; no per-cancer edits needed)
# 1. Complex_10
cat("    -> 正在加载 Pancancer 特征数据...\n")
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/Complex_10.RData")
sl_intersect <- intersect(rownames(Complex_10), sl_list)
SlFeatures[sl_intersect, "Complex_10"] <- Complex_10[sl_intersect, 3]
inx <- setdiff(sl_list, rownames(Complex_10))
SlFeatures[inx, "Complex_10"] <- 0

# 2. ComplexEssScore
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/ComplexEssScore.RData")
sl_intersect <- intersect(names(ComplexEssScore), sl_list)
SlFeatures[sl_intersect, "ComplexEssScore"] <- ComplexEssScore[sl_intersect]
rm(ComplexEssScore)
gc()

# 3. geneConservationScore
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/geneConservationScore.RData")
sl_intersect <- intersect(names(geneConservationScore), sl_list)
SlFeatures[sl_intersect, "geneConservationScore"] <- geneConservationScore[sl_intersect]
rm(geneConservationScore)
gc()

# 4. GoScore
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/GO_score/GOscore.RData")
rownames(GOscore) <- GOscore[, 1]
sl_intersect <- intersect(rownames(GOscore), sl_list)
SlFeatures[sl_intersect, "GOscore"] <- GOscore[sl_intersect, 4]
rm(GOscore)
gc()

# 5. jacard_colocalisation
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/jacard_colocalisation.RData")
sl_intersect <- intersect(names(jacard_colocalisation), sl_list)
SlFeatures[sl_intersect, "jacard_colocalisation"] <- jacard_colocalisation[sl_intersect]
rm(jacard_colocalisation)
gc()

# 6. Jacard_Domain
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/Jacard_Domain.RData")
sl_intersect <- intersect(names(Jacard_Domain), sl_list)
SlFeatures[sl_intersect, "Jacard_Domain"] <- Jacard_Domain[sl_intersect]
rm(Jacard_Domain)
gc()

# 7. Mean_age
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/Mean_age.RData")
sl_intersect <- intersect(names(Mean_age), sl_list)
SlFeatures[sl_intersect, "Mean_age"] <- Mean_age[sl_intersect]
rm(Mean_age)
gc()

# 8. paralog
load(file = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/sixFeatures0.01/GeneFunSim/paralog.RData")
sl_intersect <- intersect(names(paralog), sl_list)
SlFeatures[sl_intersect, "paralog"] <- paralog[sl_intersect]
# Set gene pairs absent from paralog to 0
sl_not_in_paralog <- setdiff(sl_list, names(paralog))
SlFeatures[sl_not_in_paralog, "paralog"] <- 0
rm(paralog)
gc()

# 9. pathway_p
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/pathway_p.RData")
sl_intersect <- intersect(names(pathway_p), sl_list)
SlFeatures[sl_intersect, "pathway_p"] <- pathway_p[sl_intersect]
rm(pathway_p)
gc()

# 10. PPI_share
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/PPI_share.RData")
sl_intersect <- intersect(names(PPI_share), sl_list)
SlFeatures[sl_intersect, "PPI_share"] <- PPI_share[sl_intersect]
rm(PPI_share)
gc()

# 11. PPI_union
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/PPI_union.RData")
sl_intersect <- intersect(names(PPI_union), sl_list)
SlFeatures[sl_intersect, "PPI_union"] <- PPI_union[sl_intersect]
rm(PPI_union)
gc()

# 12. ppiEssScore
load("/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA01/ppiEssScore.RData")
sl_intersect <- intersect(names(ppiEssScore), sl_list)
SlFeatures[sl_intersect, "ppiEssScore"] <- ppiEssScore[sl_intersect]
rm(ppiEssScore)
gc()

# 13. proteinInt
load(file = "/data/home/chenjiahao/nuaa/synlethDB/project_source/PPI/proteinInt.RData")
proteinInt <- proteinInt[[2]]
# Keep gene pairs ordered consistently (smaller gene ID first)
inx <- which(proteinInt[, 1] > proteinInt[, 2])
proteinInt[inx, ] <- proteinInt[inx, c(2, 1)]
# Deduplicate, keeping unique gene pairs
proteinInt <- proteinInt[!duplicated(paste(proteinInt[, 1], proteinInt[, 2], sep = ";")), ]
slName <- paste(proteinInt[, 1], proteinInt[, 2], sep = ";")
proteinInt <- rep(1, times = length(slName))
names(proteinInt) <- slName
# Intersect with sl_list
sl_intersect <- intersect(names(proteinInt), sl_list)
SlFeatures[sl_intersect, "proteinInt"] <- proteinInt[sl_intersect]
# Set gene pairs absent from proteinInt to 0
sl_not_in_proteinInt <- setdiff(sl_list, names(proteinInt))
SlFeatures[sl_not_in_proteinInt, "proteinInt"] <- 0
rm(proteinInt)
gc()

cat("    -> 正在对特征数据进行数值格式化...\n")
SlFeatures <- round(SlFeatures, digits = 4)
save(SlFeatures, file = paste0("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/", cancer, "_SL_Matrix.RData"))
