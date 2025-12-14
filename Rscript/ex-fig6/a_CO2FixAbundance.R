rm(list=ls())

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# ========== 1. load data ==========
path.MAG.KEGG.CO2 <- './data/ex-fig6/a_kegg_level3_pathway_allMAG_CO2Fixation.xlsx'
MAG.KEGG.CO2 <- as.data.frame(read_excel(path.MAG.KEGG.CO2))

path.MAG.info <- './data/ex-fig6/a_MAG_Genus_Merged.xlsx'
MAG.info <- as.data.frame(read_excel(path.MAG.info))

# ========== 2. Calculation of CO2 fixation genes ==========
MAG_CO2_sum <- MAG.KEGG.CO2 %>%
  group_by(MAG, Level3) %>%
  summarise(GeneNum = sum(`Gene num`, na.rm=TRUE), .groups = "drop")

MAG_total <- MAG_CO2_sum %>%
  group_by(MAG) %>%
  summarise(TotalGeneNum = sum(GeneNum), .groups = "drop")

top20_MAG <- MAG_total %>%
  arrange(desc(TotalGeneNum)) %>%
  slice_head(n = 20) %>%
  pull(MAG)

MAG_CO2_top20 <- MAG_CO2_sum %>%
  filter(MAG %in% top20_MAG)

MAG_order <- MAG_total %>%
  filter(MAG %in% top20_MAG) %>%
  arrange(TotalGeneNum) %>%
  pull(MAG)

MAG_CO2_top20$MAG <- factor(MAG_CO2_top20$MAG, levels = MAG_order)

# ========== 3. Merge information ==========
MAG.info.top20 <- MAG.info %>%
  filter(MAG %in% top20_MAG)

MAG.info.top20$MAG <- factor(MAG.info.top20$MAG, levels = MAG_order)

MAG2Genus <- MAG.info.top20 %>%
  select(MAG, Genus) %>%
  distinct()
MAG2Genus <- MAG2Genus[match(MAG_order, MAG2Genus$MAG), ]  

if(any(duplicated(MAG2Genus$Genus))) {
  MAG2Genus$Genus_unique <- paste0(MAG2Genus$Genus, " (", MAG2Genus$MAG, ")")
} else {
  MAG2Genus$Genus_unique <- MAG2Genus$Genus
}

# ========== 4. plotdata preparation ==========
# 1. heatmap data
MAG.abund.top20 <- MAG.info.top20 %>%
  select(MAG, Genus, Initial, Dark, IR) %>%
  pivot_longer(cols = c("Initial", "Dark", "IR"), names_to = "Condition", values_to = "Abundance")

MAG.abund.top20$MAG <- factor(MAG.abund.top20$MAG, levels = MAG_order)

MAG.abund.top20 <- MAG.abund.top20 %>%
  left_join(MAG2Genus[,c("MAG","Genus_unique")], by = "MAG")

MAG.abund.top20$Genus_unique <- factor(MAG.abund.top20$Genus_unique, levels = MAG2Genus$Genus_unique)
MAG.abund.top20$Condition <- factor(MAG.abund.top20$Condition, levels = c("Initial", "Dark", "IR"))

# 2. barplot data
MAG_CO2_top20 <- MAG_CO2_top20 %>%
  left_join(MAG2Genus, by = "MAG")
MAG_CO2_top20$Genus_unique <- factor(MAG_CO2_top20$Genus_unique, levels = MAG2Genus$Genus_unique)

# ========== 4.5. barplot annotation ==========
MAG_total_top20 <- MAG_total %>%
  filter(MAG %in% top20_MAG) %>%
  left_join(MAG2Genus[,c("MAG","Genus_unique")], by = "MAG")

bar_label_df <- MAG_CO2_top20 %>%
  group_by(Genus_unique) %>%
  summarise(
    max_x = sum(GeneNum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(MAG_total_top20[,c("Genus_unique","TotalGeneNum")], by = "Genus_unique")

# ========== 5. Visualization theme ==========

cns_theme <- theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(size = 0.8, color = 'black'),
    axis.ticks = element_line(size = 0.8, color = 'black'),
    axis.text.x = element_text(color = 'black', size = 15, face = "plain"),
    axis.text.y = element_text(color = 'black', size = 12, face = "plain", family = "sans", margin = margin(r = 2), hjust = 1),
    axis.title.x = element_text(face = 'plain', size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(face = 'plain', size = 16, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5,  size = 18, margin = margin(b = 10), face = "plain"),
    legend.position = "right", # 先右侧，后面patchwork合并
    legend.title = element_text(size = 15, face = "plain"),
    legend.text = element_text(size = 13),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    strip.background = element_blank(),
    strip.text = element_text(face = 'plain', size = 15)
  )

cns_co2_palette <- c(
  "#FF9800", # 黄色
  "#673AB7", # 紫色
  "#1B2A41", # 深蓝
  "#D7263D", # 高对比红
  "#21A179", # 绿色
  "#F4A259", # 橙色
  "#3F88C5", # 蓝色
  "#F49D6E" # 浅橙
)
n_pathway <- length(unique(MAG_CO2_top20$Level3))
if (n_pathway > length(cns_co2_palette)) {
  cns_co2_palette <- colorRampPalette(cns_co2_palette)(n_pathway)
} else {
  cns_co2_palette <- cns_co2_palette[1:n_pathway]
}
names(cns_co2_palette) <- unique(MAG_CO2_top20$Level3)

abund_palette <- colorRampPalette(c("#F5F5F5", "#FFD491", "#FF4066"))(100) # "#4575b4",

library(patchwork)
library(ggnewscale)

MAG.abund.top20$Abundance_log <- log10(MAG.abund.top20$Abundance + 1)

p_heatmap <- ggplot(MAG.abund.top20, aes(x = Condition, y = Genus_unique, fill = Abundance_log)) +
  geom_tile(color = "grey80", width = 0.98, height = 0.95) +
  scale_fill_gradientn(
    colors = abund_palette,
    name = expression(log[10]*"(Abundance+1)")
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void(base_size = 16) +
  theme(
    axis.text.y = element_text(color = 'black', size = 12, face = "plain", family = "sans", hjust = 1), # 右对齐，略小
    axis.text.x = element_text(color = 'black', size = 12, face = "plain", angle = 45, hjust = 0.5),
    legend.position = "right", # 先右侧，后面patchwork合并
    legend.title = element_text(size = 15, face = "plain"),
    legend.text = element_text(size = 13),
    plot.margin = margin(15, 0, 15, 0) # 右侧margin为0，紧贴bar
  )


max_gene_sum <- max(bar_label_df$max_x, na.rm = TRUE)

x_axis_max <- ceiling(max_gene_sum * 1.12)
if (x_axis_max < 90) x_axis_max <- 90

p_bar <- ggplot(MAG_CO2_top20, aes(x = GeneNum, y = Genus_unique, fill = Level3)) +
  geom_bar(
    stat = "identity",
    width = 0.7,  
    color = "white",
    size = 0.5,
    alpha = 0.7
  ) +
  scale_fill_manual(values = cns_co2_palette) +
  labs(
    x = "Gene counts",
    y = NULL,
    fill = "CO2 fixation pathway",
    title = ""
  ) +
  cns_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(15, 15, 15, 0) 
  ) +

  scale_x_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, x_axis_max)) +
  
  geom_text(
    data = bar_label_df,
    aes(x = max_x + max(x_axis_max, na.rm = TRUE)*0.03, y = Genus_unique, label = TotalGeneNum),
    inherit.aes = FALSE,
    hjust = 0, 
    vjust = 0.5,
    size = 4.2,
    fontface = "plain",
    color = "black"
  )

# 6. merge
library(patchwork)

p_combined <- p_heatmap + p_bar + 
  plot_layout(widths = c(0.35, 0.65), guides = "collect") & 
  theme(legend.position = "right")

print(p_combined)
