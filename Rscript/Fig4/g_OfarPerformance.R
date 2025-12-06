rm(list=ls())
# source('./scr/functions.R')
# ===== default setting =====
library(readxl)
library(ggplot2)
library(ggpubr)

# ===== 1. load data =====
path <- './data/Fig4/g_Ofar.xlsx' # NH4CO3-MF.xlsx # Optogenetics.xlsx # './data/Hydraulic.xlsx' # './data/MultiSludge.xlsx' # './data/Optogenetics.xlsx' # './data/Sum-Batch.xlsx'
sheet.name <- 'SBR' # MF-1020 #'MD&C'  # 'MAB_MDC' # 'Bottle_5Microbes' # Bottle_MAB12-0505-0525 #'Bottle'# 'Bottle_MAB12' # 'MAB'
rawdata <- data.frame(read_excel(path, sheet = sheet.name))
head(rawdata)


library(scales)
library(tidyr)
library(ggplot2)

# 假设 rawdata 已经读取
# 选取需要的列
df <- rawdata[, c("Time", "NO3", "NO2", "NH4", "HCO3")]
# 转换为长格式
df_long <- pivot_longer(df, cols = -Time, names_to = "Indicator", values_to = "Concentration")

# 计算缩放因子
scale_factor <- max(df$HCO3, na.rm = TRUE) / max(df$NO3, df$NO2, df$NH4, na.rm = TRUE)

# 添加一个新列，缩放后的 HCO3
df_long$Concentration_plot <- ifelse(df_long$Indicator == "HCO3",
       df_long$Concentration / scale_factor,
       df_long$Concentration)

# 自定义颜色
custom_colors <- c(
  "NO3" = "#E1D7B7",
  "NO2" = "#7C93C3",
  "NH4" = "#1E2A5E",
  "HCO3" = "#FF204E"
)

# 绘图
ggplot() +
  # HCO3 用柱状图
  geom_col(
  data = subset(df_long, Indicator == "HCO3"),
  aes(x = Time, y = Concentration_plot, fill = Indicator),
  width = 0.9,
  alpha = 0.3
  ) +
  # 其他指标用点线图
  geom_point(
  data = subset(df_long, Indicator != "HCO3"),
  aes(x = Time, y = Concentration_plot, color = Indicator),
  size = 4, alpha = 0.7
  ) +
  geom_line(
  data = subset(df_long, Indicator != "HCO3"),
  aes(x = Time, y = Concentration_plot, color = Indicator),
  size = 0.7
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
  name = "N (mg L-1)",
  sec.axis = sec_axis(~ . * scale_factor, name = "C (mg L-1)")
  ) +
  labs(x = "Time (day)", color = "Indicator", fill = "Indicator") +
  # theme_classic() +
  # theme(
  # plot.margin = unit(c(0.5, 0.1, 0.1, 0.5), "cm"),
  # plot.background = element_rect(fill = 'transparent', color = NA),
  # panel.background = element_rect(fill = 'transparent'),
  # panel.grid = element_blank(),
  # axis.title = element_text(size = 16),
  # axis.text = element_text(size = 16),
  # # axis.text.x = element_text(angle = 30, hjust = 1),
  # legend.position = c(0.35, 0.8),
  # legend.title = element_blank(),
  # legend.text = element_text(size = 16),
  # legend.background = element_blank()
  # ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
    axis.line = element_line(size = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.background = element_rect(fill = 'transparent'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 调整边距，使y轴第一个图距离x轴更近
    plot.margin = margin(t = 0.2, r = 1, b = 0.2, l = 0.5, unit = "cm"),
    # 调整y轴位置
    axis.ticks.length = unit(0.1, "cm"),
    # # 优化图例布局
    # legend.box = "horizontal",
    # legend.spacing.x = unit(-0.1, "cm"),
    legend.background = element_blank(),
    # 调整legend位置
    legend.position = c(0.91, 0.75)
  ) 
