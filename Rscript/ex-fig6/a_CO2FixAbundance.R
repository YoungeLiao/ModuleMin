rm(list=ls())
# 加载所需R包
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# ========== 1. 读取CO2固定通路基因注释数据 ==========
path.MAG.KEGG.CO2 <- './data/ex-fig6/a_kegg_level3_pathway_allMAG_CO2Fixation.xlsx'
MAG.KEGG.CO2 <- as.data.frame(read_excel(path.MAG.KEGG.CO2))

# 加载MAG.info数据（含分类、丰度等信息）
path.MAG.info <- './data/ex-fig6/a_MAG_Genus_Merged.xlsx'
MAG.info <- as.data.frame(read_excel(path.MAG.info))

# ========== 2. 统计每个MAG在不同CO2固定通路的基因数 ==========
MAG_CO2_sum <- MAG.KEGG.CO2 %>%
  group_by(MAG, Level3) %>%
  summarise(GeneNum = sum(`Gene num`, na.rm=TRUE), .groups = "drop")

# 计算每个MAG的CO2固定相关基因总数
MAG_total <- MAG_CO2_sum %>%
  group_by(MAG) %>%
  summarise(TotalGeneNum = sum(GeneNum), .groups = "drop")

# 选取CO2固定基因数最多的前20个MAG
top20_MAG <- MAG_total %>%
  arrange(desc(TotalGeneNum)) %>%
  slice_head(n = 20) %>%
  pull(MAG)

# 过滤出top20 MAG的数据
MAG_CO2_top20 <- MAG_CO2_sum %>%
  filter(MAG %in% top20_MAG)

# 为MAG按总基因数排序（用于y轴顺序美化）
MAG_order <- MAG_total %>%
  filter(MAG %in% top20_MAG) %>%
  arrange(TotalGeneNum) %>%
  pull(MAG)

MAG_CO2_top20$MAG <- factor(MAG_CO2_top20$MAG, levels = MAG_order)

# ========== 3. 合并MAG.info信息 ==========
# 只保留top20 MAG的MAG.info
MAG.info.top20 <- MAG.info %>%
  filter(MAG %in% top20_MAG)

# 保证MAG顺序与柱状图一致
MAG.info.top20$MAG <- factor(MAG.info.top20$MAG, levels = MAG_order)

# 生成MAG到Genus的映射（MAG_order顺序）
MAG2Genus <- MAG.info.top20 %>%
  select(MAG, Genus) %>%
  distinct()
MAG2Genus <- MAG2Genus[match(MAG_order, MAG2Genus$MAG), ]  # 保证顺序

# 检查Genus是否有重复，若有则为其添加MAG信息以保证唯一性
if(any(duplicated(MAG2Genus$Genus))) {
  MAG2Genus$Genus_unique <- paste0(MAG2Genus$Genus, " (", MAG2Genus$MAG, ")")
} else {
  MAG2Genus$Genus_unique <- MAG2Genus$Genus
}

# ========== 4. 数据准备：为绘图数据添加Genus列用于y轴 ==========
# 1. 热图数据
MAG.abund.top20 <- MAG.info.top20 %>%
  select(MAG, Genus, Initial, Dark, IR) %>%
  pivot_longer(cols = c("Initial", "Dark", "IR"), names_to = "Condition", values_to = "Abundance")

MAG.abund.top20$MAG <- factor(MAG.abund.top20$MAG, levels = MAG_order)

# 生成Genus_unique列
MAG.abund.top20 <- MAG.abund.top20 %>%
  left_join(MAG2Genus[,c("MAG","Genus_unique")], by = "MAG")

MAG.abund.top20$Genus_unique <- factor(MAG.abund.top20$Genus_unique, levels = MAG2Genus$Genus_unique)
MAG.abund.top20$Condition <- factor(MAG.abund.top20$Condition, levels = c("Initial", "Dark", "IR"))

# 2. 堆叠柱状图数据
MAG_CO2_top20 <- MAG_CO2_top20 %>%
  left_join(MAG2Genus, by = "MAG")
MAG_CO2_top20$Genus_unique <- factor(MAG_CO2_top20$Genus_unique, levels = MAG2Genus$Genus_unique)

# ========== 4.5. 为柱状图添加总CO2固定基因数信息 ==========
# 生成每个MAG的总CO2固定基因数和Genus_unique的映射
MAG_total_top20 <- MAG_total %>%
  filter(MAG %in% top20_MAG) %>%
  left_join(MAG2Genus[,c("MAG","Genus_unique")], by = "MAG")

# 计算每个Genus_unique的最大x值（用于标注位置）
bar_label_df <- MAG_CO2_top20 %>%
  group_by(Genus_unique) %>%
  summarise(
    max_x = sum(GeneNum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(MAG_total_top20[,c("Genus_unique","TotalGeneNum")], by = "Genus_unique")

# ========== 5. CNS正刊风格美化堆叠柱状图 ==========

# 1. 优化CNS正刊风格主题（去除加粗，纵坐标右对齐，字体略小，优化图例位置）
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

# 2. 自定义色盲友好高对比色盘（可根据通路数量调整）
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

# 3. 丰度热图色盘
abund_palette <- colorRampPalette(c("#F5F5F5", "#FFD491", "#FF4066"))(100) # "#4575b4",

# 4. 绘制丰度热图（左侧，丰度取对数，避免极大值掩盖其他MAG变化）
library(patchwork)
library(ggnewscale)

# 先对丰度取对数（加1防止log(0)）
MAG.abund.top20$Abundance_log <- log10(MAG.abund.top20$Abundance + 1)

# 用Genus_unique作为y轴
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

# 5. 堆叠柱状图（用Genus_unique作为y轴，宽度更窄），并标注总CO2固定基因数
# 计算横坐标最大值，留出足够空间用于标注
max_gene_sum <- max(bar_label_df$max_x, na.rm = TRUE)
# 让x轴最大值比最大堆叠值多10（或可根据实际情况调整）
x_axis_max <- ceiling(max_gene_sum * 1.12)
if (x_axis_max < 90) x_axis_max <- 90

p_bar <- ggplot(MAG_CO2_top20, aes(x = GeneNum, y = Genus_unique, fill = Level3)) +
  geom_bar(
    stat = "identity",
    width = 0.7,  # 更窄
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
    plot.margin = margin(15, 15, 15, 0) # 左侧margin为0，紧贴热图
  ) +
  # 设置x轴范围，留出空间给右侧标注
  scale_x_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, x_axis_max)) +
  # 添加总CO2固定基因数标注
  geom_text(
    data = bar_label_df,
    aes(x = max_x + max(x_axis_max, na.rm = TRUE)*0.03, y = Genus_unique, label = TotalGeneNum),
    inherit.aes = FALSE,
    hjust = 0, # 左对齐
    vjust = 0.5,
    size = 4.2,
    fontface = "plain",
    color = "black"
  )

# 6. 拼图（热图更宽，bar更窄，图例合并到一起）
library(patchwork)
# 使用guides = "collect"合并图例，图例放在右侧
p_combined <- p_heatmap + p_bar + 
  plot_layout(widths = c(0.35, 0.65), guides = "collect") & 
  theme(legend.position = "right")

print(p_combined)
