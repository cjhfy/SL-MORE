library(data.table)
library(arrow)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(viridis)
library(tidytext)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ReactomePA)
library(igraph)
library(ggraph)

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
# 2. Load and filter the gene annotation data
# ─────────────────────────────────────────────
cat("\n>>> Loading and filtering protein-coding gene annotations...\n")
gene_info_path <- "/data/home/chenjiahao/nuaa/synlethDB/gene_info/Homo_sapiens.gene_info"
gene_anno <- fread(gene_info_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_anno <- gene_anno[type_of_gene == "protein-coding", .(GeneID, Symbol, Synonyms)]
gene_anno[, Synonyms := gsub("\\|", ";", Synonyms)]
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

  dt <- read_parquet(input_file, as_data_frame = FALSE) %>% as.data.table()
  dt_filtered <- dt[CategoryCount > 3]

  if (nrow(dt_filtered) > 0) {
    dt_filtered[, Symbol1 := id_to_symbol_map[as.character(gene1)]]
    dt_filtered[, Symbol2 := id_to_symbol_map[as.character(gene2)]]
    dt_filtered <- dt_filtered[!is.na(Symbol1) & !is.na(Symbol2)]

    pair_count <- nrow(dt_filtered)
    if (pair_count == 0) next

    gene_sl_counts <- rbind(
      dt_filtered[, .(Symbol = Symbol1)],
      dt_filtered[, .(Symbol = Symbol2)]
    )[, .(SL_Count = .N), by = Symbol][order(-SL_Count)]

    gene_sl_counts[, Cancer := cancer]
    gene_sl_counts_list[[cancer]] <- gene_sl_counts
    counts_list[[cancer]] <- data.table(Cancer = cancer, PairCount = pair_count)
    write_parquet(dt_filtered, file.path(gt4_out_dir, paste0(cancer, "_GT4_Features.parquet")))
    fwrite(gene_sl_counts, file.path(gene_stat_dir, paste0(cancer, "_Gene_SL_Counts.csv")))
    cat(sprintf("[+] %s done | coding gene pairs: %d\n", cancer, pair_count))
  }
  rm(dt, dt_filtered)
  gc()
}

all_gene_sl_counts <- rbindlist(gene_sl_counts_list)
df_counts <- rbindlist(counts_list)
fwrite(df_counts, file.path(stat_out_dir, "Cancer_CategoryCount_GT4_Summary.csv"))
fwrite(all_gene_sl_counts, file.path(stat_out_dir, "all_gene_sl_counts.csv"))


all_gene_sl_counts <- fread("/data/home/chenjiahao/nuaa/synlethDB/FeatureMatrix/Synergy_Stats/all_gene_sl_counts.csv", data.table = F)

# ─────────────────────────────────────────────
# 4. Recurrent SL hub genes
# ─────────────────────────────────────────────
top100_hubs <- all_gene_sl_counts %>%
  group_by(Cancer) %>%
  slice_max(order_by = SL_Count, n = 100, with_ties = FALSE) %>%
  ungroup()

recurrent_hubs <- top100_hubs %>%
  group_by(Symbol) %>%
  summarise(CancerCount = n()) %>%
  filter(CancerCount >= 5) %>%
  arrange(desc(CancerCount))

hub_genes <- recurrent_hubs$Symbol

# ─────────────────────────────────────────────
# 5. Gene ID conversion
# ─────────────────────────────────────────────
gene_df <- bitr(
  hub_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
entrez_genes <- unique(gene_df$ENTREZID)

# ─────────────────────────────────────────────
# 6. Hallmark enrichment
# ─────────────────────────────────────────────
# Background set: all genes in the GT4 consensus SL pairs (potential hubs)
background_symbols <- unique(all_gene_sl_counts$Symbol)
background_df <- bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
background_entrez <- unique(background_df$ENTREZID)
length(background_entrez)
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_res <- enricher(
  gene = entrez_genes,
  TERM2GENE = msig_hallmark[, c("gs_name", "entrez_gene")],
  universe = background_entrez,
  pAdjustMethod = "BH"
)
hallmark_df <- as.data.frame(hallmark_res)
write.csv(hallmark_df, file.path(stat_out_dir, "Hallmark_Enrichment.csv"), row.names = FALSE)


library(ggplot2)
library(dplyr)
top_pathways <- hallmark_df %>%
  arrange(p.adjust)
heat_df <- top_pathways %>%
  dplyr::select(Description, GeneRatio, p.adjust)
# Convert to numeric
heat_df$GeneRatioNum <- sapply(heat_df$GeneRatio, function(x) eval(parse(text = x)))
# Tidy up pathway names
heat_df$Description <- gsub("HALLMARK_", "", heat_df$Description)
heat_df$Description <- gsub("_", " ", heat_df$Description)
heat_df$Description <- factor(heat_df$Description, levels = rev(heat_df$Description))

# Plot
# Force barplot font to exactly 12 pt
p <- ggplot(heat_df, aes(x = GeneRatioNum, y = Description, fill = -log10(p.adjust))) +
  geom_col(width = 0.6) +
  # Force label font size to exactly 12 pt
  geom_text(aes(label = round(GeneRatioNum, 2)), color = "white", size = 12 / 2.84, fontface = "bold", hjust = 1.1) +
  scale_fill_gradient(low = "#DCE6F1", high = "#08306B", name = expression(-log[10](FDR))) +
  labs(x = "Gene Ratio", y = NULL) + # Drop the title; it belongs in the figure legend
  theme_classic(base_size = 12) + # Global base font size 12
  theme(
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.title.x = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 0.5), # size replaced by linewidth (newer ggplot2 convention)
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank()
  )

ggsave(file.path(stat_out_dir, "Hallmark_Enrichment_Heatmap_SCIstyle.pdf"),
  p,
  width = 6, height = 4, device = "pdf"
)
# ─────────────────────────────────────────────
# 9. Pathway-Gene network (4 Hallmarks)
# ─────────────────────────────────────────────
top_pathway_names <- hallmark_df$Description
# Extract the pathway-gene mapping
net_df <- msig_hallmark %>%
  dplyr::filter(gs_name %in% top_pathway_names)
id_map <- bitr(
  unique(as.character(net_df$entrez_gene)),
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)
# Unify types
net_df$entrez_gene <- as.character(net_df$entrez_gene)
# merge
net_df <- dplyr::left_join(net_df, id_map, by = c("entrez_gene" = "ENTREZID"))
# Keep recurrent hubs
net_df <- net_df %>% dplyr::filter(SYMBOL %in% hub_genes)
# Tidy up pathway names
net_df$gs_name <- gsub("HALLMARK_", "", net_df$gs_name)
net_df$gs_name <- gsub("_", " ", net_df$gs_name)
# edges
edges <- net_df %>%
  dplyr::select(gs_name, SYMBOL) %>%
  unique()
colnames(edges) <- c("Pathway", "Gene")
# ─────────────────────────────────────────────
# 9. Pathway-Gene network (GB Style - Compact Version)
# ─────────────────────────────────────────────
# --- [1. Graph object and topological attributes] ---
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$degree <- degree(g)
V(g)$type <- ifelse(V(g)$name %in% edges$Gene, "Gene", "Pathway")

# --- [2. Shorten edges: rebuild a compact two-ring layout] ---
pathway_nodes <- which(V(g)$type == "Pathway")
gene_nodes <- which(V(g)$type == "Gene")
coords <- matrix(0, nrow = vcount(g), ncol = 2)

# Inner ring (pathways): unchanged or slightly adjusted
theta_p <- seq(0, 2 * pi, length.out = length(pathway_nodes) + 1)[1:length(pathway_nodes)]
coords[pathway_nodes, 1] <- cos(theta_p) * 1.0
coords[pathway_nodes, 2] <- sin(theta_p) * 1.0

# Outer ring (genes): shrink the radius from 3.5 to 2.2 to shorten the edges
theta_g <- seq(0, 2 * pi, length.out = length(gene_nodes) + 1)[1:length(gene_nodes)]
coords[gene_nodes, 1] <- cos(theta_g) * 2.2
coords[gene_nodes, 2] <- sin(theta_g) * 2.2

# --- [3. Fine-tune labels for the compact layout] ---
# With the tighter layout, reduce gene label size to 2.8 to avoid overlap
V(g)$label_size <- ifelse(V(g)$type == "Pathway", 4.8, 2.8)
V(g)$label_font <- ifelse(V(g)$type == "Pathway", "bold", "plain")
V(g)$label_color <- ifelse(V(g)$type == "Pathway", "#111111", "#333333")

sci_colors <- c("Pathway" = "#E64B35", "Gene" = "#4DBBD5")

# --- [4. High-resolution rendering] ---
pdf(
  file.path(stat_out_dir, "Pathway_Gene_Network_Compact.pdf"),
  width = 7.0, # Slightly smaller canvas for a tighter composition
  height = 7.0
)

# --- [3. Fine-tune labels for the compact layout] ---
# Pathway labels at 12 and gene labels at 10 to avoid overlap
V(g)$label_size <- ifelse(V(g)$type == "Pathway", 12 / 2.84, 10 / 2.84)

p_net <- ggraph(g, layout = coords) +
  geom_edge_arc(aes(alpha = after_stat(index)), colour = "grey85", strength = 0.05, linewidth = 0.45) +
  geom_node_point(aes(size = ifelse(type == "Pathway", degree * 1.2 + 8, 4.5)), color = "white", stroke = 0) +
  geom_node_point(aes(color = type, size = ifelse(type == "Pathway", degree * 1.2 + 6, 3.5)), alpha = 0.95) +
  geom_node_text(
    aes(label = name, size = label_size, fontface = label_font, color = label_color),
    repel = TRUE, max.overlaps = 100, box.padding = 0.1, point.padding = 0.15,
    segment.color = "grey75", segment.linewidth = 0.1
  ) +
  scale_color_manual(name = "Entities", values = sci_colors, labels = c("Hub Gene", "Hallmark Pathway")) +
  scale_size_identity() +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"), # Legend title 12
    legend.text = element_text(size = 12), # Legend text 12
    plot.margin = margin(20, 20, 20, 20)
  )
print(p_net)

dev.off()
