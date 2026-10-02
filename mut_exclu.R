library(data.table)
library(ggplot2)
library(stringi)

### ----------------- 1. Environment setup and paths -----------------###
base_paths <- list(
  gene_info    = "/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info",
  mut_dir      = "/data/home/chenjiahao/nuaa/synlethDB/project_source/TCGA-mut/",
  dir_ExcluMut = "/data/home/chenjiahao/nuaa/synlethDB/project_source/features/TCGA/ExclusiveMut/" # Update to the new path you provided
)

# List of 32 cancer types
cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "READ", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO",
  "OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
  "UCEC", "UCS"
)

# Extract protein-coding genes (for the left-axis mutated-gene count)
gene_anno <- fread(file = base_paths$gene_info, header = TRUE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
gene_anno <- gene_anno[which(gene_anno$type_of_gene == "protein-coding"), ]
genesProtein <- as.character(gene_anno[, 2])

### ----------------- 2. Loop over precomputed features and summarize -----------------###
plot_data <- data.table(
  Cancer      = character(),
  Mut_Genes   = numeric(), # Left axis: total mutated genes
  Sig_Pairs   = numeric() # Right axis: significant mutually exclusive pairs
)

for (cancer in cancers) {
  cat(sprintf("[正在提取] %s ...\n", cancer))

  file_mut <- paste0(base_paths$mut_dir, "TCGA-", cancer, ".RData")
  file_ExcluMut <- paste0(base_paths$dir_ExcluMut, "TCGA-", cancer, ".RData")

  if (!file.exists(file_mut) || !file.exists(file_ExcluMut)) {
    cat(sprintf("  [跳过] 找不到配对的文件。\n"))
    next
  }

  # ---- Left axis: count mutated genes from the raw mutation data ----
  load(file_mut) # Loaded variable is COAD_mut
  mut_code <- colnames(COAD_mut)
  mut_code <- substr(mut_code, 1, 12)
  inx <- duplicated(mut_code)
  COAD_mut <- COAD_mut[, !inx, drop = FALSE]
  colnames(COAD_mut) <- substr(colnames(COAD_mut), 1, 12)
  COAD_mut <- COAD_mut[rownames(COAD_mut) %in% genesProtein, , drop = FALSE]
  COAD_mut <- COAD_mut[rowSums(COAD_mut > 0) > 0, , drop = FALSE]
  # Mutation frequency per gene across samples, averaged and converted to percent
  gene_frequencies <- rowSums(COAD_mut > 0) / ncol(COAD_mut)
  mean_mutation_frequency <- mean(gene_frequencies) * 100

  # ---- Right axis: count significant mutually exclusive gene pairs via the parser ----
  load(file_ExcluMut) # Loaded variable is mut_res

  # Count significant mutually exclusive pairs (threshold assumed at P < 0.05; adjust as needed)
  num_sig_pairs <- length(which(mut_res < 0.01))

  # Store in the plot data frame
  plot_data <- rbind(plot_data, data.table(
    Cancer    = cancer,
    Mean_Freq = mean_mutation_frequency,
    Sig_Pairs = num_sig_pairs
  ))
}

### ----------------- 3. Publication-quality dual-axis bar chart -----------------###
cat("\n[作图] 正在绘制双轴图...\n")

# Sort the x-axis by total mutated genes, descending
plot_data <- plot_data[order(-Mut_Genes)]
plot_data[, Cancer := factor(Cancer, levels = Cancer)]

# Compute the dual-axis scaling factor automatically
max_left <- max(plot_data$Mut_Genes, na.rm = TRUE)
max_right <- max(plot_data$Sig_Pairs, na.rm = TRUE)
if (max_right == 0) max_right <- 1
scale_factor <- max_left / max_right

# Start plotting
p <- ggplot(plot_data, aes(x = Cancer)) +
  # Bar 1: total genes (left axis, dark teal)
  geom_bar(aes(y = Mut_Genes, fill = "Mutated Genes"), stat = "identity", width = 0.35, position = position_nudge(x = -0.2)) +
  # Bar 2: significant pair count (right axis, light orange)
  geom_bar(aes(y = Sig_Pairs * scale_factor, fill = "Significant Pairs (P < 0.05)"), stat = "identity", width = 0.35, position = position_nudge(x = 0.2)) +

  # Dual Y-axis scaling
  scale_y_continuous(
    name = "Total Number of Mutated Genes",
    limits = c(0, max_left * 1.1),
    # Use a standard anonymous function instead of a formula
    sec_axis = sec_axis(trans = function(x) x / scale_factor, name = "Number of Significant Mut_Excl Pairs")
  ) +

  # Academic palette following journal conventions
  scale_fill_manual(values = c("Mutated Genes" = "#2c7fb8", "Significant Pairs (P < 0.05)" = "#fdae61")) +
  labs(x = "Cancer Types (TCGA)", fill = "Metrics") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.title.y = element_text(color = "#2c7fb8", face = "bold"),
    axis.text.y = element_text(color = "#2c7fb8"),
    axis.title.y.right = element_text(color = "#e6550d", face = "bold"),
    axis.text.y.right = element_text(color = "#e6550d"),
    panel.grid.major = element_line(color = "#f5f5f5"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

# Save high-resolution figures
fig_output_path <- paste0(base_paths$dir_ExcluMut, "Supplementary_Fig_S1.png")
ggsave(fig_output_path, plot = p, width = 11, height = 5.5, dpi = 300)

cat(sprintf("[完成] 图表已生成并保存至: %s\n", fig_output_path))
