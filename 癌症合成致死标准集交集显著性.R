# ====================================================
# Statistical overlap between predicted SLs and
# gold-standard SL datasets across cancers
# ====================================================

library(tidyverse)
library(data.table)
library(arrow)
library(stringr)

#----------------------------------------------------
# Cancer types
#----------------------------------------------------
cancers <- c("COAD", "BRCA", "OV", "CESC", "LUAD", "SKCM", "KIRC")

#----------------------------------------------------
# Load gold-standard SL data
#----------------------------------------------------
load("/data/home/chenjiahao/nuaa/SL_norm/sl.golden.set.RData")

sl_ISLE_all <- gd$sr0
type <- gd$dat

#----------------------------------------------------
# Load gene annotation
#----------------------------------------------------
gene_anno <- fread("/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE
)

# keep protein-coding genes
gene_anno <- gene_anno[gene_anno$type_of_gene == "protein-coding", c(2, 3, 5)]

# alias symbols
gene_anno[, 3] <- gsub("\\|", ";", gene_anno[, 3])

# ====================================================
# Function: Convert gene symbols -> Entrez pairs
# ====================================================
convert_to_entrez_pairs <- function(sl_df, gene_anno, cancer_type = NULL, type_info = NULL) {
  if (!is.null(cancer_type) && !is.null(type_info)) {
    sl_df <- sl_df[type_info[, 2] == cancer_type, ]
  }

  if (nrow(sl_df) == 0) {
    return(list(pairs = character(0), mapping_df = data.frame()))
  }

  genelist <- unique(c(sl_df[, 1], sl_df[, 2]))

  # Direct matching
  gene_anno1 <- gene_anno[gene_anno[, 2] %in% genelist, ]
  rownames(gene_anno1) <- gene_anno1[, 2]

  # Alias matching
  genelist1 <- genelist[!genelist %in% gene_anno1[, 2]]
  gene_anno2 <- data.frame(GeneID = NULL, Symbol = NULL) # Restore NULL initialization to avoid coercion to character

  if (length(genelist1) > 0) {
    for (i in seq_along(genelist1)) {
      index <- unique(c(
        grep(paste0(";", genelist1[i], ";"), gene_anno[, 3], ignore.case = TRUE),
        grep(paste0("^", genelist1[i], ";"), gene_anno[, 3], ignore.case = TRUE),
        grep(paste0(";", genelist1[i], "$"), gene_anno[, 3], ignore.case = TRUE),
        grep(paste0("^", genelist1[i], "$"), gene_anno[, 3], ignore.case = TRUE)
      ))

      if (length(index) == 1) {
        gene_anno2[i, 1] <- gene_anno[index, 1]
        gene_anno2[i, 2] <- genelist1[i]
      }
    }
    gene_anno2 <- na.omit(gene_anno2)

    if (nrow(gene_anno2) > 0) {
      rownames(gene_anno2) <- gene_anno2[, 2]
      colnames(gene_anno2) <- c("GeneID", "Symbol")
      gene_anno_final <- rbind(gene_anno1[, 1:2], gene_anno2)
    } else {
      gene_anno_final <- gene_anno1[, 1:2]
    }
  } else {
    gene_anno_final <- gene_anno1[, 1:2]
  }

  rownames(gene_anno_final) <- gene_anno_final[, 2]

  # Keep gene pairs that can be mapped
  sl_df <- sl_df[sl_df[, 1] %in% gene_anno_final[, 2] & sl_df[, 2] %in% gene_anno_final[, 2], ]

  if (nrow(sl_df) == 0) {
    return(list(pairs = character(0), mapping_df = data.frame()))
  }

  # Extract Entrez IDs from row names and coerce to numeric for exact comparison
  id1 <- as.numeric(gene_anno_final[sl_df[, 1], 1])
  id2 <- as.numeric(gene_anno_final[sl_df[, 2], 1])

  # Sort and join with the same pmin/pmax logic used for pred_pairs
  # Avoid the "9" > "10" bug caused by lexicographic sorting
  sl_pairs_sorted <- paste(pmin(id1, id2), pmax(id1, id2), sep = ";")
  sl_pairs <- unique(sl_pairs_sorted)

  # Build the mapping record matrix
  mapping_df <- data.frame(
    Cancer = ifelse(is.null(cancer_type), "All", cancer_type),
    Symbol_1 = sl_df[, 1],
    Symbol_2 = sl_df[, 2],
    GeneID_1 = id1,
    GeneID_2 = id2,
    Sorted_ID_Pair = sl_pairs_sorted,
    stringsAsFactors = FALSE
  )

  mapping_df <- unique(mapping_df)

  return(list(pairs = sl_pairs, mapping_df = mapping_df))
}

# ====================================================
# Main loop
# ====================================================
result_list <- list()
mapping_list <- list()

for (cancer in cancers) {
  cat("\n==============================\nProcessing:", cancer, "\n==============================\n")

  # Gold-standard SLs
  conversion_res <- convert_to_entrez_pairs(sl_ISLE_all, gene_anno, cancer_type = cancer, type_info = type)
  sl_gold_pairs <- conversion_res$pairs
  mapping_list[[cancer]] <- conversion_res$mapping_df

  cat("Mapped gold-standard pairs:", length(sl_gold_pairs), "\n")

  if (length(sl_gold_pairs) == 0) {
    cat("Warning: No gold-standard pairs mapped for", cancer, "\n")
    next
  }

  # Predicted SLs (HighConf_GT4)
  pred_file <- paste0("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Cancer_HighConf_GT4_Features/", cancer, "_GT4_Features.parquet")
  pred_df <- read_parquet(pred_file) %>% as.data.table()

  # Ensure comparison is numeric
  pred_id1 <- as.numeric(pred_df$gene1)
  pred_id2 <- as.numeric(pred_df$gene2)
  pred_pairs <- unique(paste(pmin(pred_id1, pred_id2), pmax(pred_id1, pred_id2), sep = ";"))

  cat("Predicted SL pairs:", length(pred_pairs), "\n")

  # Intersection
  intersection <- intersect(pred_pairs, sl_gold_pairs)
  cat("Overlap:", length(intersection), "\n")

  # Hypergeometric enrichment (using Full_Features for background)
  bg_file <- paste0("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Cancer_Merged_Features/", cancer, "_Full_Features.parquet")
  bg_df <- read_parquet(bg_file) %>% as.data.table()
  bg_size <- nrow(bg_df)

  q <- length(intersection)
  m <- length(sl_gold_pairs)
  n <- bg_size - m
  k <- length(pred_pairs)

  if (q > 0) {
    pvalue <- phyper(q = q - 1, m = m, n = n, k = k, lower.tail = FALSE)
  } else {
    pvalue <- 1
  }

  result_list[[cancer]] <- data.frame(
    Source = "ISLE",
    Cancer = cancer,
    Predicted = length(pred_pairs),
    Gold = length(sl_gold_pairs),
    Overlap = length(intersection),
    Pvalue = pvalue
  )
}

# ====================================================
# Save Gene Mapping Result
# ====================================================
mapping_final_df <- bind_rows(mapping_list)
write.csv(mapping_final_df, "SL_gold_standard_gene_conversion_mapping.csv", row.names = FALSE)
cat("\nSaved gene conversion mapping to 'SL_gold_standard_gene_conversion_mapping.csv'\n")

# ====================================================
# Final result table
# ====================================================
result_df <- bind_rows(result_list)

result_df <- result_df %>%
  mutate(
    logP = -log10(Pvalue),
    OverlapRatio = Overlap / Gold,
    Recall = Overlap / Gold,
    Precision = Overlap / Predicted,
    F1 = 2 * (Precision * Recall) / (Precision + Recall),
    Significant = ifelse(Pvalue < 0.05, "Yes", "No") # Update the 0.05 threshold for the radar plot as well
  )

result_df$logP[is.infinite(result_df$logP)] <- max(result_df$logP[is.finite(result_df$logP)])

print(result_df)
write.csv(result_df, "SL_gold_standard_overlap_statistics.csv", row.names = FALSE)

# ====================================================
# Publication-quality figure: enrichment significance of SL predictions
# Data: result_df (requires Cancer, logP, Overlap, and similar columns)
# ====================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)

# Assume result_df exists; uncomment the next line to read from CSV:
# result_df <- read.csv("SL_gold_standard_overlap_statistics.csv", stringsAsFactors = FALSE)

# Ensure the data types are correct
result_df <- result_df %>%
  mutate(
    Cancer = as.character(Cancer),
    logP = as.numeric(logP),
    Overlap = as.integer(Overlap),
    Pvalue = as.numeric(Pvalue)
  )

# Optional: sort by -log10(P-value) descending (most significant on the left)
cancer_order <- result_df %>%
  arrange(desc(logP)) %>%
  pull(Cancer)
result_df$Cancer <- factor(result_df$Cancer, levels = cancer_order)

# Significance threshold: -log10(0.05)
sig_threshold <- -log10(0.05)

# Load required packages
library(ggplot2)
library(ggsci) # Provides palettes commonly used in SCI journals
library(ggrepel) # If text labels are needed

# Plot
p <- ggplot(result_df, aes(x = Cancer, y = logP)) +
  geom_col(width = 0.7, fill = "#3C5488FF", color = "black", size = 0.3) +
  geom_hline(
    yintercept = sig_threshold, linetype = "dashed",
    color = "red", size = 0.8
  ) +
  labs(
    x = "Cancer type",
    y = expression(-log[10](italic(P - value)))
  ) +
  theme_classic(base_size = 12) + # Key change: classic theme with axis lines only
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

# Save as high-resolution figures (vector PDF plus PNG backup)
ggsave("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Synergy_Stats/intersect/SL_enrichment_significance.pdf", p, width = 4, height = 4, dpi = 600)

# Display the plot
print(p)
