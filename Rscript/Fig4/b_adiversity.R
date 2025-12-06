rm(list=ls())
# source('./Rscript/functions.R')
# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

# load data
path <- './data/Fig4/b_adiversity.xlsx' # ex.fig2e.Shannon.xlsx
sheet.name <- 'sheet1' 
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

# ensure types
rawdata$Day <- as.numeric(rawdata$Day)
rawdata$Light <- as.factor(rawdata$Light)
rawdata$Group <- as.factor(rawdata$Group)

# summary by Light and Day (用于绘制每个 Light 条件下的趋势线)
summary_df <- rawdata %>%
  group_by(Light, Day) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),
    sd_Shannon = sd(Shannon, na.rm = TRUE),
    n = sum(!is.na(Shannon)),
    .groups = "drop"
  )

# publication-quality palette mapped to Group levels (auto-expand if needed)
# 修改调色板设置
pal <- setNames(c("#305291", "#CF545C"), c("Dark", "IR"))  # 灰黑色和浅红色

# groups <- levels(rawdata$Group)
# n_groups <- length(groups)
# if (n_groups > length(pal)) {
#   palette_full <- grDevices::colorRampPalette(pal)(n_groups)
# } else {
#   palette_full <- pal[seq_len(n_groups)]
# }
# pal <- setNames(palette_full, groups)

# 绘图：美化以符合高质量期刊（清晰的线条、统一刻度、均值与误差条、可读的文字）
ggplot(rawdata, aes(x = factor(Day), y = Shannon, fill = Light)) +
  # boxplots with black outlines
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black",
               width = 0.65,
               alpha = 0.9,
               linewidth = 0.5) +
  # points as filled circles with black border
  geom_jitter(aes(fill = Light),
              shape = 21,
              color = "black",
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.8,
              stroke = 0.25,
              alpha = 0.85,
              show.legend = TRUE) +
  # mean +/- sd error bars per Light x Day
  geom_errorbar(data = summary_df, inherit.aes = FALSE,
                aes(x = factor(Day), ymin = mean_Shannon - sd_Shannon, ymax = mean_Shannon + sd_Shannon),
                width = 0.15, color = "black", size = 0.5) +
  # mean trend line and points (aggregated within each Light; facets will subset automatically)
  geom_line(data = summary_df, inherit.aes = FALSE,
            aes(x = factor(Day), y = mean_Shannon, group = 1),
            color = "black", size = 0.5) +
  geom_point(data = summary_df, inherit.aes = FALSE,
             aes(x = factor(Day), y = mean_Shannon),
             color = "black", size = 2) +
  # facet by Light with consistent y-scale for easier comparison across panels
  facet_wrap(~ Light, scales = "fixed", ncol = 2) +
  scale_fill_manual(values = pal) +
  labs(x = "Time (day)", y = "Shannon index", fill = "") +
  theme_classic(base_size = 12) +
  theme(
    # text = element_text(family = "Arial"),  # 添加Arial字体设置
    # text = element_text(family = "sans"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 0.5),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.ticks = element_line(size = 0.5),
    strip.text = element_text(size = 14),
    panel.spacing = grid::unit(1, "lines"),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = grid::unit(0.8, "lines")
  )
