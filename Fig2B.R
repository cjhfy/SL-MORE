#!/usr/bin/env Rscript
# =============================================================================
# Figure 2B: Sankey Flow of Multi-Omics Features
#
# Three-level Sankey diagram linking each of the 30 curated features to its
# underlying data type and source database
# (Feature -> Data Type -> Database).
#
# Input : feature_data.csv (feature, data, database; ";" separated multiples)
# Output: Feature_Sankey.pdf
# =============================================================================

library(tidyverse)
library(ggalluvial)
library(RColorBrewer)

# ---- Configuration ----------------------------------------------------------

input_csv <- "/data/home/chenjiahao/nuaa/synlethDB/Feature_Method/feature_data.csv"
out_pdf   <- "/data/home/chenjiahao/nuaa/synlethDB/Feature_Method/Feature_Sankey.pdf"

# ---- 1. Load data -----------------------------------------------------------

raw_data <- read.csv(input_csv, sep = ",", header = FALSE)
colnames(raw_data) <- c("feature", "data", "database")

# ---- 2. Clean data (avoid Cartesian products on multi-value fields) ---------

df_clean <- raw_data %>%
  mutate(
    data_list = strsplit(data, ";"),
    db_list   = strsplit(database, ";")
  ) %>%
  pmap_dfr(function(feature, data, database, data_list, db_list, ...) {
    d   <- trimws(unlist(data_list))
    b   <- trimws(unlist(db_list))
    n_d <- length(d)
    n_b <- length(b)

    if (n_d == n_b) {
      # One-to-one mapping, e.g. "PPI;CRISPR | BioGRID;DepMap"
      tibble(feature = feature, data = d, database = b)
    } else if (n_b == 1) {
      # Single database shared by all data types
      tibble(feature = feature, data = d, database = b[1])
    } else if (n_d == 1) {
      # Single data type linked to multiple databases
      tibble(feature = feature, data = d[1], database = b)
    } else {
      # Fallback for complex cases (not expected in current data)
      expand.grid(data = d, database = b, stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        mutate(feature = feature)
    }
  })

# ---- 3. Define orderings ----------------------------------------------------

features_order <- c(
  "Exp_Comp", "CNV_Comp", "Exp_Excl", "CNV_Excl_ISLE", "CNV_Excl_DSL",
  "Mut_Excl", "MutCNV_Excl", "Exp_SB_ISLE", "Exp_SB_SiLi", "CNV_SB",
  "Exp_Low_CR_Dep", "Exp_Low_RNAi_Dep", "CNV_Loss_CR_Dep",
  "CNV_Loss_RNAi_Dep", "MutDelLow_CR_Dep", "MutDelLow_RNAi_Dep",
  "PPI_Agg_Dep", "PCom_Agg_Dep", "CoExp", "GO_Sim", "Paral_Rel",
  "SubLoc", "Domain", "Path_CoM", "PCom_CoM", "PPI_Ovlp", "PPI_Union",
  "Direct_Int", "conservation_score", "Gene_Age"
)

data_levels <- df_clean %>% count(data) %>% arrange(desc(n)) %>% pull(data)
db_levels   <- df_clean %>% count(database) %>% arrange(desc(n)) %>% pull(database)

# ---- 4. Summarize and normalize flow widths ---------------------------------
# Normalize counts within each feature so that all Feature-column strata have
# equal height and every feature's flows sum to total_pairs / n_features.

df_plot <- df_clean %>%
  group_by(feature, data, database) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    feature  = factor(feature, levels = features_order),
    data     = factor(data, levels = data_levels),
    database = factor(database, levels = db_levels)
  )

total_pairs <- sum(df_plot$count)
n_features  <- length(features_order)

df_plot <- df_plot %>%
  group_by(feature) %>%
  mutate(Freq = count / sum(count) * (total_pairs / n_features)) %>%
  ungroup()

# ---- 5. Convert to lodes format ---------------------------------------------

df_lodes <- to_lodes_form(
  df_plot,
  key   = "x",
  value = "stratum",
  id    = "alluvium",
  axes  = 1:3
)

df_lodes <- df_lodes %>%
  mutate(
    x = case_when(
      x == "feature"  ~ "Feature",
      x == "data"     ~ "Data Type",
      x == "database" ~ "Database",
      TRUE            ~ as.character(x)
    ),
    x = factor(x, levels = c("Feature", "Data Type", "Database"))
  )

# Color each flow by its originating feature
flow_labels <- df_lodes %>%
  filter(x == "Feature") %>%
  select(alluvium, flow_fill = stratum)

df_final <- df_lodes %>%
  left_join(flow_labels, by = "alluvium")

# ---- 6. Fixed-height strata for the Feature column --------------------------

feature_equal <- data.frame(
  stratum = factor(features_order, levels = features_order),
  x       = factor("Feature", levels = c("Feature", "Data Type", "Database")),
  Freq    = total_pairs / n_features
)

# ---- 7. Colors --------------------------------------------------------------

feature_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(features_order))
names(feature_colors) <- features_order

# ---- 8. Plot ----------------------------------------------------------------

p_refined <- ggplot() +
  geom_alluvium(
    data = df_final,
    aes(x = x, y = Freq, stratum = stratum, alluvium = alluvium, fill = flow_fill),
    width = 1 / 4, alpha = 0.4, knot.pos = 0.4
  ) +
  geom_stratum(
    data = feature_equal,
    aes(x = x, y = Freq, stratum = stratum, fill = stratum),
    width = 1 / 4, color = "white", linewidth = 0.3
  ) +
  geom_stratum(
    data = filter(df_final, x != "Feature"),
    aes(x = x, y = Freq, stratum = stratum),
    fill = "#F0F0F0", width = 1 / 4, color = "#999999", linewidth = 0.2
  ) +
  # Stratum labels at 12 pt (size = 12 / .pt converts pt to ggplot units)
  geom_text(
    data = feature_equal,
    aes(x = x, y = Freq, stratum = stratum, label = as.character(stratum)),
    stat = "stratum", size = 12 / .pt, fontface = "bold"
  ) +
  geom_text(
    data = filter(df_final, x != "Feature"),
    aes(x = x, y = Freq, stratum = stratum, label = as.character(stratum)),
    stat = "stratum", size = 12 / .pt, fontface = "bold"
  ) +
  scale_fill_manual(values = feature_colors) +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  labs(title = "Sankey Flow of Multi-Omics Features", y = "") +
  theme_minimal() +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(size = 12, face = "plain", color = "black"),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# ---- 9. Export --------------------------------------------------------------

ggsave(out_pdf, p_refined, width = 15, height = 10)
