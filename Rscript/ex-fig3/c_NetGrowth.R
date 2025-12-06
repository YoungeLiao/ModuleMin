rm(list=ls())
# load library
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# ========== 1. Load basic data and metadata ==========

# 加载数据
path.MAG.info <- './data/ex-fig3/c_MAG_Genus_Merged.xlsx'
MAG.info <- as.data.frame(read_excel(path.MAG.info))

# 计算各样本（Dark1/2/3, IR1/2/3）相较于Initial的净增长（NetGrowth）
# 保留原始abundance列和netgrowth列
netgrowth_long <- MAG.info %>%
  select(Genus, MAG, Initial, Dark1, Dark2, Dark3, IR1, IR2, IR3) %>%
  pivot_longer(
    cols = c(Initial, Dark1, Dark2, Dark3, IR1, IR2, IR3),
    names_to = "Condition",
    values_to = "Abundance"
  ) %>%
  # 标记分组和重复
  mutate(
    Group = case_when(
      grepl("^Dark", Condition) ~ "Dark",
      grepl("^IR", Condition) ~ "IR",
      Condition == "Initial" ~ "Initial"
    ),
    Replicate = str_extract(Condition, "[123]")
  )

# 计算NetGrowth（每个MAG、Genus、样本的Abundance - Initial）
# 先获得Initial abundance
initial_df <- netgrowth_long %>%
  filter(Condition == "Initial") %>%
  select(Genus, MAG, Initial_Abundance = Abundance)

# 合并Initial abundance
netgrowth_long <- netgrowth_long %>%
  left_join(initial_df, by = c("Genus", "MAG")) %>%
  mutate(
    NetGrowth = Abundance - Initial_Abundance
  )

# 只保留非Initial的行用于后续统计
netgrowth_plotdf <- netgrowth_long %>%
  filter(Group %in% c("Dark", "IR"))

# 计算每个Genus-Group的净增长均值和标准差
summary_df <- netgrowth_plotdf %>%
  group_by(Genus, Group) %>%
  summarise(
    Mean_NetGrowth = mean(NetGrowth, na.rm = TRUE),
    SD_NetGrowth = sd(NetGrowth, na.rm = TRUE),
    n = sum(!is.na(NetGrowth)),
    SE_NetGrowth = SD_NetGrowth / sqrt(n)
  ) %>%
  ungroup()

# 仅可视化净增长丰度top20的genus（按净增长均值排序，负值表示减少）
genus_netgrowth_mean <- summary_df %>%
  group_by(Genus) %>%
  summarise(NetGrowthMean = mean(Mean_NetGrowth, na.rm = TRUE)) %>%
  arrange(desc(NetGrowthMean)) %>%
  slice_head(n = 25)

top20_genus <- genus_netgrowth_mean$Genus

summary_df <- summary_df %>%
  filter(Genus %in% top20_genus)

# 设定Genus顺序（按top20净增长均值排序，负值表示减少）
summary_df$Genus <- factor(summary_df$Genus, levels = top20_genus)

# 对NetGrowth取对数（加1防止负值和0，保留正负号信息）
# 这里采用对称log变换：sign(x) * log10(abs(x) + 1)
summary_df <- summary_df %>%
  mutate(
    Mean_NetGrowth_log = sign(Mean_NetGrowth) * log10(abs(Mean_NetGrowth) + 1),
    SE_NetGrowth_log = SE_NetGrowth / (abs(Mean_NetGrowth) + 1) / log(10)  # 近似变换后的误差
  )

# 计算Dark与IR之间的显著性差异（每个Genus做t检验）
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggsignif)

# 只保留top genus的原始数据
netgrowth_points <- netgrowth_plotdf %>%
  filter(Genus %in% levels(summary_df$Genus)) %>%
  mutate(
    NetGrowth_log = sign(NetGrowth) * log10(abs(NetGrowth) + 1)
  )

# 保证Genus因子顺序与summary_df一致
netgrowth_points$Genus <- factor(netgrowth_points$Genus, levels = levels(summary_df$Genus))

# 计算每个Genus在Dark和IR之间的t检验p值
stat_test <- netgrowth_points %>%
  filter(Group %in% c("Dark", "IR")) %>%
  group_by(Genus) %>%
  summarise(
    p.value = tryCatch(
      t.test(NetGrowth_log ~ Group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p.signif = case_when(
      is.na(p.value) ~ "",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 计算显著性标记的y轴位置（在两组bar的最大值上方加一定偏移）
# 先获得每个Genus的最大y值
y_max_df <- summary_df %>%
  group_by(Genus) %>%
  summarise(y_max = max(Mean_NetGrowth_log + SE_NetGrowth_log, na.rm = TRUE), .groups = "drop")

# 合并p值和y轴位置
stat_test <- stat_test %>%
  left_join(y_max_df, by = "Genus") %>%
  mutate(
    y.position = y_max + 0.18  # 可根据需要调整偏移
  )
# ===== CNS正刊风格美化版可视化 =====

# 1. CNS正刊风格主题（更精致，去除多余元素，增强对比，优化字体）
cns_theme <- theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(size = 0.8, color = 'black'),
    axis.ticks = element_line(size = 0.8, color = 'black'),
    axis.text.x = element_text(color = 'black', size = 15, angle = 45, hjust = 1, vjust = 1, face = "italic"),
    axis.text.y = element_text(color = 'black', size = 15),
    axis.title.x = element_text(face = 'bold', size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(face = 'bold', size = 18, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5,  size = 20, margin = margin(b = 10)),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold', size = 16)
  )

# 2. CNS色盲友好高对比分组配色
# 注意：这里顺序要和Group的实际因子顺序一致，且legend显示顺序也要一致
cns_colors <- c(
  "Dark" = "#1B2A41",   # 高对比深蓝
  "IR" = "#D7263D"      # 高对比红
)

# 3. 绘图
p <- ggplot(summary_df, aes(x = Genus, y = Mean_NetGrowth_log, group = Group)) +
  # 空心柱状图（CNS风格，线条加粗，去除填充）
  geom_bar(
    aes(color = Group),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.8,
    fill = NA,
    size = 1
  ) +
  # 误差线（加粗，黑色）
  geom_errorbar(
    aes(ymin = Mean_NetGrowth_log - SE_NetGrowth_log, ymax = Mean_NetGrowth_log + SE_NetGrowth_log, group = Group),
    position = position_dodge(width = 0.7),
    width = 0.3,
    size = 1,
    color = "black"
  ) +
  # 原始数据点（实心圆，分组色，轻微抖动，带黑色描边）
  geom_point(
    data = netgrowth_points,
    aes(x = Genus, y = NetGrowth_log, color = Group, fill = Group),
    position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.7),
    size = 2.6,
    shape = 21,
    stroke = 0.7,
    alpha = 0.6
  ) +
  # 显著性标记（加粗，居中，略上移）
  geom_text(
    data = stat_test,
    aes(x = Genus, y = y.position, label = p.signif),
    inherit.aes = FALSE,
    size = 6.5,
    vjust = 0,
    fontface = "bold",
    color = "#D7263D"
  ) +
  scale_color_manual(
    values = cns_colors,
    name = "Condition",
    breaks = c("Dark", "IR"),
    labels = c("Dark", "IR")
  ) +
  scale_fill_manual(
    values = cns_colors,
    name = "Condition",
    breaks = c("Dark", "IR"),
    labels = c("Dark", "IR")
  ) +
  labs(
    x = "",
    y = expression("Net Abundance Growth (log)"),
    title = NULL
  ) +
  cns_theme +
  coord_cartesian(clip = "off") +
  guides(
    color = guide_legend(
      override.aes = list(size = 4, shape = 21), # , fill = cns_colors[c("Dark", "IR")]
      order = 1
    ),
    fill = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.22)),
    breaks = pretty_breaks(n = 7)
  )

# 4. CNS风格导出（高分辨率，白色背景，边距充足）
ggsave("Fig2d_NetGrowth_Barplot.pdf", p, width = 9, height = 5, dpi = 400, units = "in")
# 保存可视化用的数据 summary_df 到指定目录
write.csv(summary_df, file = "results_data/Fig2/Fig2e_NetGrowth_summary_df.csv", row.names = FALSE)

print(p)
