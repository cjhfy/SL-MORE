# 1. 加载必要的包
library(ggplot2)
library(reshape2)
library(grid)
library(gtable)
library(dplyr) 
library(cowplot) # 用于完美对齐绘图面板

# ==============================================================================
# --- 1. 数据准备 ---
# ==============================================================================

# 1.1 从提供的 Cell 文章插图中提取 6 种配色
data_types <- c(rep("Genomic Selection Signatures", 7), 
                rep("Survival Benefit", 3), 
                rep("Gene Dependency", 8), 
                rep("Gene Function Similarity", 7),
                rep("Network Properties", 3),
                rep("Evolutionary Conservation", 2))

# 采用高级柔和配色方案
type_colors <- c(
  "Genomic Selection Signatures" = "#789ECF", # 柔和蓝
  "Survival Benefit"             = "#AAD4B9", # 豆沙绿
  "Gene Dependency"              = "#F2A471", # 柔和橘
  "Gene Function Similarity"     = "#6AC3D3", # 湖水蓝
  "Network Properties"           = "#EDB9D2", # 藕粉色
  "Evolutionary Conservation"    = "#A465A4"  # 灰紫色
)

# 1.2 读取数据
mat <- read.csv("/data/home/chenjiahao/nuaa/synlethDB/SL-MORE/data/Fig2/method_feature.csv", row.names = 1, check.names = FALSE)

features_in_csv <- rownames(mat)
methods_in_csv  <- colnames(mat)

# 1.3 构建特征元数据
feature_meta <- data.frame(Feature = features_in_csv, DataType = data_types, stringsAsFactors = FALSE)

# 1.4 转换为长表并融合
df <- melt(as.matrix(mat))
colnames(df) <- c("Feature", "Method", "Used")
df <- merge(df, feature_meta, by = "Feature")

# 1.5 锁定因子级别，保证顺序不乱
df$Method <- factor(df$Method, levels = methods_in_csv) 
df$Feature <- factor(df$Feature, levels = rev(features_in_csv))
df$DataType <- factor(df$DataType, levels = unique(data_types))


# ==============================================================================
# --- 2. 构建基础图表对象 ---
# ==============================================================================

# 【图一】矩阵图 (中心热图)
p <- ggplot(df, aes(x = Method, y = Feature)) + 
  geom_tile(aes(fill = factor(Used)), color = "white", linewidth = 0.5) +
  facet_grid(DataType ~ ., scales = "free_y", space = "free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  scale_fill_manual(
      values = c("0" = "#F3F3F3", "1" = "#4A5A70"), 
      labels = c("0" = "Not Derived", "1" = "Source"), 
      name = "Provenance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", hjust = 0, color = "black"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 14, color = "white"),
    strip.placement = "outside",
    panel.spacing = unit(0.2, "lines"),
    panel.border = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey", color = NA),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  labs(
      title = "Feature-Method Provenance Matrix", 
      x = "Computational Methods", 
      y = ""
  )

# 【图二】右侧柱状图：统计每个特征被使用的方法数
# 【核心修复】：引入与主图完全一致的分面（facet_grid）并克隆间距，隐藏标签
feat_counts <- df %>%
  group_by(Feature, DataType) %>%
  summarise(MethodCount = sum(Used), .groups = "drop")

p_feat_count <- ggplot(feat_counts, aes(y = Feature, x = MethodCount, fill = DataType)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9, show.legend = FALSE) +
  facet_grid(DataType ~ ., scales = "free_y", space = "free_y") + # 同步点1：建立相同的分面块
  scale_fill_manual(values = type_colors) +
  scale_x_continuous(breaks = 0:length(methods_in_csv), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15)),
    strip.background = element_blank(), # 同步点2：隐藏右侧分面的灰色背景条
    strip.text = element_blank(),       # 同步点3：隐藏右侧分面的文字
    panel.spacing = unit(0.2, "lines")  # 同步点4：保持与主图绝对一致的间块间距
  ) +
  labs(
    title = "Prevalence of Features",
    x = "Method Count"
  )

# 【图三】上方堆叠柱状图：统计每个方法涵盖的特征类别数量
method_cat_presence <- df %>%
  filter(Used == 1) %>%
  group_by(Method, DataType) %>%
  summarise(CategoryPresence = 1, .groups = "drop")

method_cat_presence$Method <- factor(method_cat_presence$Method, levels = methods_in_csv)
method_cat_presence$DataType <- factor(method_cat_presence$DataType, levels = unique(data_types))

p_method_count <- ggplot(method_cat_presence, aes(x = Method, y = CategoryPresence, fill = DataType)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = type_colors) +
  scale_y_continuous(breaks = 0:length(unique(data_types)), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", margin = margin(b = 15))
  ) +
  labs(
    title = "Number of Feature dimensions Used per Method",
    y = "Category Count"
  )


# ==============================================================================
# --- 3. 核心修复：多维面板对齐与分面着色 ---
# ==============================================================================

# 3.1 垂直对齐：强制对齐 上方柱状图 和 中心主图 的左右边缘
aligned_v <- cowplot::align_plots(p_method_count, p, align = "v", axis = "lr")
g_top  <- aligned_v[[1]]
g_main <- aligned_v[[2]]

# 3.2 水平对齐：强制对齐 中心主图 和 右侧柱状图 的上下边缘
# 这一步会动态微调右图的 padding，使其核心绘图区（Panel）高度与主图像素级契合
aligned_h <- cowplot::align_plots(g_main, p_feat_count, align = "h", axis = "tb")
g_main_final <- aligned_h[[1]]
g_right      <- aligned_h[[2]]

# 3.3 动态着色侧边分面条带 (操作最终对齐后的 g_main_final)
strip_indices <- which(grepl("strip-l", g_main_final$layout$name))
ordered_types <- levels(df$DataType)

for (i in seq_along(strip_indices)) {
  index <- strip_indices[i]
  fill_color <- type_colors[ordered_types[i]]
  g_main_final$grobs[[index]]$grobs[[1]]$children[[1]]$gp$fill <- fill_color
  g_main_final$grobs[[index]]$grobs[[1]]$children[[1]]$gp$col  <- fill_color
}


# ==============================================================================
# --- 4. 最终图片导出 ---
# ==============================================================================

# 保存图一：完美对齐且修改好颜色的中心主热图 (宽9，高8)
ggsave("/data/home/chenjiahao/nuaa/synlethDB/Feature_Method/Feature_method.pdf", 
       plot = g_main_final, width = 9, height = 8)

# 保存图二：【已修复】高度和内部特征行均完美对齐的右侧统计图 (宽4.5，高8)
ggsave("/data/home/chenjiahao/nuaa/synlethDB/Feature_Method/Feature_use_count.pdf", 
       plot = g_right, width = 4.5, height = 8, dpi = 300) 

# 保存图三：完美匹配主图宽度的上方统计图 (宽9，高3.5)
ggsave("/data/home/chenjiahao/nuaa/synlethDB/Feature_Method/Method_use_count.pdf", 
       plot = g_top, width = 9, height = 2, dpi = 300)