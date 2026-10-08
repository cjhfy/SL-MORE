library(data.table)
library(arrow)
library(dplyr)
library(ggplot2)
library(ggridges)
library(viridis)
# ─────────────────────────────────────────────
# 1. Paths and basic configuration
# ─────────────────────────────────────────────
base_out_dir <- "/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/"
detail_out_dir <- file.path(base_out_dir, "Cancer_Merged_Features")
stat_out_dir <- file.path(base_out_dir, "Synergy_Stats")

gt4_out_dir <- file.path(base_out_dir, "Cancer_HighConf_GT4_Features")
if (!dir.exists(gt4_out_dir)) dir.create(gt4_out_dir, recursive = TRUE)

gene_stat_dir <- file.path(stat_out_dir, "Gene_SL_Stats")
if (!dir.exists(gene_stat_dir)) dir.create(gene_stat_dir, recursive = TRUE)

cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

# ─────────────────────────────────────────────
# 2. Load and filter the gene annotation data (matches the original logic exactly)
# ─────────────────────────────────────────────
cat("\n>>> Loading and filtering protein-coding gene annotations...\n")
gene_info_path <- "/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info"

gene_anno <- fread(gene_info_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Strict filter: keep protein-coding genes only
gene_anno <- gene_anno[type_of_gene == "protein-coding", .(GeneID, Symbol, Synonyms)]
# Handle separators in aliases
gene_anno[, Synonyms := gsub("\\|", ";", Synonyms)]

# Build the mapping dictionary (GeneID as character for easier matching)
id_to_symbol_map <- setNames(gene_anno$Symbol, as.character(gene_anno$GeneID))

cat(sprintf("[+] Filtering done: kept %d protein-coding genes for mapping.\n", nrow(gene_anno)))

# ─────────────────────────────────────────────
# 3. Loop over the cancer data
# ─────────────────────────────────────────────
counts_list <- list()
gene_sl_counts_list <- list()

for (cancer in cancers) {
  input_file <- file.path(detail_out_dir, paste0(cancer, "_Full_Features.parquet"))
  if (!file.exists(input_file)) next

  # Read and convert to data.table
  dt <- read_parquet(input_file, as_data_frame = FALSE) %>% as.data.table()

  # Keep CategoryCount > 3
  dt_filtered <- dt[CategoryCount > 3]

  if (nrow(dt_filtered) > 0) {
    # Map Symbol
    dt_filtered[, Symbol1 := id_to_symbol_map[as.character(gene1)]]
    dt_filtered[, Symbol2 := id_to_symbol_map[as.character(gene2)]]

    # Key step: keep only pairs where both genes are protein-coding and mapped
    dt_filtered <- dt_filtered[!is.na(Symbol1) & !is.na(Symbol2)]

    pair_count <- nrow(dt_filtered)
    if (pair_count == 0) {
      cat(sprintf("[!] %s has no protein-coding gene pairs after filtering; skipping\n", cancer))
      next
    }
    # Count per-gene SL frequency (based on filtered Symbols)
    gene_sl_counts <- rbind(
      dt_filtered[, .(Symbol = Symbol1)],
      dt_filtered[, .(Symbol = Symbol2)]
    )[, .(SL_Count = .N), by = Symbol][order(-SL_Count)]

    gene_sl_counts[, Cancer := cancer]
    gene_sl_counts_list[[cancer]] <- gene_sl_counts
    counts_list[[cancer]] <- data.table(Cancer = cancer, PairCount = pair_count)

    # Export the results
    write_parquet(dt_filtered, file.path(gt4_out_dir, paste0(cancer, "_GT4_Features.parquet")))
    fwrite(gene_sl_counts, file.path(gene_stat_dir, paste0(cancer, "_Gene_SL_Counts.csv")))

    cat(sprintf("[+] %s done | coding gene pairs: %d\n", cancer, pair_count))
  }
  rm(dt, dt_filtered)
  gc()
}
all_gene_sl_counts <- rbindlist(gene_sl_counts_list)
df_counts <- rbindlist(counts_list)
# Export the summary CSV
fwrite(df_counts, file.path(stat_out_dir, "Cancer_CategoryCount_GT4_Summary.csv"))
fwrite(all_gene_sl_counts, file.path(stat_out_dir, "all_gene_sl_counts.csv"))

### Heatmap of the Top 5 genes per cancer type
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)
all_gene_sl_counts <- fread(file.path(stat_out_dir, "all_gene_sl_counts.csv"))
# 1. Extract the Top 5 gene data
top5_summary <- all_gene_sl_counts %>%
  group_by(Cancer) %>%
  slice_max(order_by = SL_Count, n = 5, with_ties = TRUE) %>%
  ungroup()

# 2. Count and sort the number of cancers per gene (key step)
# Order by cancer count descending, then by Symbol alphabetically
gene_order_df <- top5_summary %>%
  group_by(Symbol) %>%
  summarise(n_cancer = n()) %>%
  arrange(desc(n_cancer), Symbol)

top5_union_genes <- gene_order_df$Symbol

# 3. Build the main matrix in the sorted order
plot_data <- all_gene_sl_counts[Symbol %in% top5_union_genes]
mat <- dcast(plot_data, Symbol ~ Cancer, value.var = "SL_Count", fill = 0)
matrix_mat <- as.matrix(mat[, -1])
rownames(matrix_mat) <- mat$Symbol

# Ensure the matrix row order matches the breadth order computed above
matrix_mat <- matrix_mat[top5_union_genes, ]

# 4. Build the marker matrix (TRUE/FALSE)
mark_mat <- matrix(FALSE, nrow = nrow(matrix_mat), ncol = ncol(matrix_mat))
rownames(mark_mat) <- rownames(matrix_mat)
colnames(mark_mat) <- colnames(matrix_mat)
for (i in 1:nrow(top5_summary)) {
  g <- top5_summary$Symbol[i]
  c <- top5_summary$Cancer[i]
  if (g %in% rownames(mark_mat)) mark_mat[g, c] <- TRUE
}

# 5. Colors and annotations
max_val <- max(matrix_mat)
col_fun <- colorRamp2(c(0, max_val * 0.2, max_val), c("#F7FBFF", "#41B6C4", "#08306B"))

# Right-side bar annotation (must also match the sorted order)
right_anno <- rowAnnotation(
  "Cancers Count" = anno_barplot(
    gene_order_df$n_cancer,
    baseline = 0,
    gp = gpar(fill = "#41B6C4", col = "white"),
    border = FALSE,
    width = unit(3, "cm"),
    axis_param = list(side = "bottom", at = c(0, 10, 20, 30))
  ),
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

# 6. Plot
pdf(file.path(stat_out_dir, "Top_Genes_Heatmap_Ranked_By_Count.pdf"), width = 13, height = 16)
h <- Heatmap(matrix_mat,
  name = "SL Count",
  col = col_fun,
  cluster_rows = FALSE, # Already sorted manually, so turn off clustering
  cluster_columns = FALSE,
  width = unit(13, "cm"),
  height = unit(22, "cm"),

  # Cell marker: solid white diamond
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (mark_mat[i, j]) {
      grid.points(x, y, pch = 18, size = unit(3, "mm"), gp = gpar(col = "white"))
    }
  },

  # Annotations and layout
  right_annotation = right_anno,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "plain"),
  column_names_gp = gpar(fontsize = 12, fontface = "plain"),
  column_title = "TOP 5 SL Hub Genes of per cancer",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
    title = "Interaction Count",
    direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    legend_width = unit(8, "cm")
  ),
  border = TRUE
)

draw(h, heatmap_legend_side = "bottom")
dev.off()

### Redraw the transposed version
# 1. Extract the Top 5 gene data
top5_summary <- all_gene_sl_counts %>%
  group_by(Cancer) %>%
  slice_max(order_by = SL_Count, n = 5, with_ties = TRUE) %>%
  ungroup()

# 2. Count and sort the number of cancers per gene
gene_order_df <- top5_summary %>%
  group_by(Symbol) %>%
  summarise(n_cancer = n()) %>%
  arrange(desc(n_cancer), Symbol)

top5_union_genes <- gene_order_df$Symbol

# 3. Build the main matrix in the sorted order
plot_data <- all_gene_sl_counts[Symbol %in% top5_union_genes]
mat <- dcast(plot_data, Symbol ~ Cancer, value.var = "SL_Count", fill = 0)
matrix_mat <- as.matrix(mat[, -1])
rownames(matrix_mat) <- mat$Symbol

# Ensure the matrix row order matches the computed breadth order
matrix_mat <- matrix_mat[top5_union_genes, ]

# [Key change 1]: transpose the matrix so rows are cancers and columns are genes
t_matrix_mat <- t(matrix_mat)

# 4. Build the marker matrix (TRUE/FALSE)
mark_mat <- matrix(FALSE, nrow = nrow(matrix_mat), ncol = ncol(matrix_mat))
rownames(mark_mat) <- rownames(matrix_mat)
colnames(mark_mat) <- colnames(matrix_mat)
for (i in 1:nrow(top5_summary)) {
  g <- top5_summary$Symbol[i]
  c <- top5_summary$Cancer[i]
  if (g %in% rownames(mark_mat)) mark_mat[g, c] <- TRUE
}

# [Key change 2]: transpose the marker matrix accordingly
t_mark_mat <- t(mark_mat)

# 5. Colors and annotations
max_val <- max(t_matrix_mat)
col_fun <- colorRamp2(c(0, max_val * 0.2, max_val), c("#F7FBFF", "#41B6C4", "#08306B"))

# [Key change 3]: move the right-side row annotation to a top column annotation (HeatmapAnnotation)
top_anno <- HeatmapAnnotation(
  "Cancers Count" = anno_barplot(
    gene_order_df$n_cancer, # The order here still matches the column (gene) order
    baseline = 0,
    gp = gpar(fill = "#41B6C4", col = "white"),
    border = FALSE,
    height = unit(3, "cm"), # Former width becomes height
    axis_param = list(side = "left", at = c(0, 10, 20, 30)) # Move the axis to the left
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

# 6. Plot
# [Key change 4]: swap PDF width and height for a landscape layout
pdf(file.path(stat_out_dir, "Top_Genes_Heatmap_Ranked_By_Count_Rotated.pdf"), width = 16, height = 13)

h <- Heatmap(t_matrix_mat,
  name = "SL Count",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  # Swap width and height of the heatmap body
  width = unit(22, "cm"),
  height = unit(13, "cm"),

  # Cell marker: solid white diamond
  cell_fun = function(j, i, x, y, width, height, fill) {
    # i indexes rows (Cancer) and j indexes columns (Gene); index the transposed mark_mat directly
    if (t_mark_mat[i, j]) {
      grid.points(x, y, pch = 18, size = unit(3, "mm"), gp = gpar(col = "white"))
    }
  },

  # Annotation and layout adjustments
  top_annotation = top_anno,

  # Rows are now cancers: left side, bold
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "plain"),

  # Columns are now genes: bottom side, rotated 45 degrees to avoid overlap
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 12, fontface = "plain"),
  column_names_rot = 45,

  # Put the main title on top
  column_title = "TOP 5 SL Hub Genes of per cancer",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
    title = "Interaction Count",
    direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    legend_width = unit(8, "cm")
  ),
  border = TRUE
)

draw(h, heatmap_legend_side = "bottom")
dev.off()
## Possibly a supplementary figure
### Top 10 genes per cancer type
# Assume all_gene_sl_counts already has Symbol, SL_Count, and Cancer
top10_genes <- all_gene_sl_counts %>%
  group_by(Cancer) %>%
  slice_max(order_by = SL_Count, n = 10, with_ties = FALSE) %>%
  ungroup()

# 2. Plot with facets
# reorder_within() from tidytext keeps genes ordered within each facet
# If tidytext is not installed: install.packages("tidytext")
library(tidytext)

p_top10 <- ggplot(top10_genes, aes(
  x = reorder_within(Symbol, SL_Count, Cancer),
  y = SL_Count, fill = Cancer
)) +
  geom_col(show.legend = FALSE, width = 0.8) +
  scale_x_reordered() + # Used with reorder_within to strip the suffix
  coord_flip() + # Horizontal layout so gene names are readable
  facet_wrap(~Cancer, scales = "free", ncol = 4) + # 4-column layout; adjust as needed
  theme_bw(base_size = 12) +
  labs(
    x = "Gene Symbol",
    y = "Number of SL Interactions",
    title = "Top 10 SL Hub Genes across Cancers"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"), # Facet strip boxes
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "plain"),
    axis.text.x = element_text(size = 12),
    panel.spacing = unit(0.5, "lines"), # Adjust facet spacing
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

# 3. Save (many cancer types, so the PDF needs enough height)
ggsave(file.path(stat_out_dir, "PanCancer_Top10_Genes.pdf"), p_top10, width = 12, height = 18)
