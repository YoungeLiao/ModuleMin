rm(list=ls())
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

# 加载数据
# path <- './data/Fig1/ex.fig1/ex-fig1b.CO2.MultiSludgeGroup.xlsx' 
path <- './data/ex-fig1/ac_CO2.MultiSludgeGroup.xlsx'
sheet.name <- 'MultiSludge' # MultiSludge # CN
rawdata <- data.frame(read_excel(path, sheet = sheet.name))
head(rawdata, n=10)

# 计算每个组在每个时间点的均值和标准差
summary_data <- rawdata %>%
  group_by(Group, Time) %>%
  summarise(
    Mean_HCO3 = mean(HCO3),
    SD_HCO3 = sd(HCO3),
    n = n(),
    SE_HCO3 = SD_HCO3 / sqrt(n)
  ) %>%
  ungroup()

# 设置CNS风格的图形参数
theme_set(theme_classic() +
            theme(
              # axis.line = element_line(size = 0.8),
              # axis.ticks = element_line(size = 0.8),
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 14),
              legend.position = "right",
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              plot.margin = unit(c(1, 1, 1, 1), "lines"),
              plot.background = element_rect(fill = 'transparent', color = NA),
              panel.background = element_rect(fill = 'transparent')
            ))

# 创建点线图
p <- ggplot() +
  # 添加原始数据点（使用较小的点）
  geom_point(data = rawdata, 
             aes(x = Time, y = HCO3, color = Group), 
             size = 2, alpha = 0.5) +
  # 添加均值点（较大的点）
  geom_point(data = summary_data, 
             aes(x = Time, y = Mean_HCO3, color = Group), 
             size = 4) +
  # 添加误差棒（可选）
  geom_errorbar(data = summary_data,
                aes(x = Time, y = Mean_HCO3, ymin = Mean_HCO3 - SE_HCO3, 
                    ymax = Mean_HCO3 + SE_HCO3, color = Group),
                width = 0.2, size = 0.8) +
  # 添加连接均值的线
  geom_line(data = summary_data, 
            aes(x = Time, y = Mean_HCO3, color = Group), 
            size = 1) +
  # # 添加置信区间阴影
  # geom_ribbon(data = summary_data,
  #             aes(x = Time, 
  #                 ymin = Mean_HCO3 - 1.96 * SE_HCO3, 
  #                 ymax = Mean_HCO3 + 1.96 * SE_HCO3,
  #                 fill = Group),
  #             alpha = 0.2) +
  # 设置颜色
  scale_color_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  # 添加坐标轴标签
  # labs(x = "Time (d)", y = expression(HCO[3]^-~""~"(mg/L)")) +
  labs(x = "Time (d)", y = expression(HCO[3]^-~""~"(mg·"~L^"-1"~")")) +
  # # 添加网格线
  # theme(panel.grid.major = element_line(color = "grey90", size = 0.3),
  #       panel.grid.minor = element_blank()) +
  ylim(0, max(rawdata$HCO3) * 1.1) +
  # 调整图例
  guides(color = guide_legend(override.aes = list(size = 3)))

# 显示图形
print(p)

# 如果需要保存图片
# ggsave("HCO3_Time_Series.png", p, width = 8, height = 6, dpi = 300)

# ===== for sheet 'CN' ======
# 创建点线图
p <- ggplot() +
  # 添加原始数据点（使用较小的点）
  geom_point(data = rawdata, 
             aes(x = Time, y = HCO3, color = Group), 
             size = 2, alpha = 0.5) +
  # 添加均值点（较大的点）
  geom_point(data = summary_data, 
             aes(x = Time, y = Mean_HCO3, color = Group), 
             size = 4) +
  # 添加误差棒（可选）
  geom_errorbar(data = summary_data,
                aes(x = Time, y = Mean_HCO3, ymin = Mean_HCO3 - SE_HCO3, 
                    ymax = Mean_HCO3 + SE_HCO3, color = Group),
                width = 0.2, size = 0.8) +
  # 添加置信区间阴影
  geom_ribbon(data = summary_data,
              aes(x = Time, 
                  ymin = Mean_HCO3 - 1.96 * SE_HCO3, 
                  ymax = Mean_HCO3 + 1.96 * SE_HCO3,
                  fill = Group),
              alpha = 0.2) +
  # 添加连接均值的线
  geom_line(data = summary_data, 
            aes(x = Time, y = Mean_HCO3, color = Group), 
            size = 1) +
  # 设置颜色
  scale_color_manual(values = c("#363062","#B31312", "#67729D", "#E7BCDE", "#BBE2EC",  "#F6EACB")) +
  scale_fill_manual(values = c("#363062","#B31312", "#67729D", "#E7BCDE", "#BBE2EC",  "#F6EACB")) +
  # 添加坐标轴标签
  labs(x = "Time (d)", y = expression(HCO[3]^-~""~"(mg·"~L^"-1"~")")) +
  # 设置Y轴范围
  ylim(0, max(rawdata$HCO3) * 1.1) +
  # 调整图例
  guides(color = guide_legend(override.aes = list(size = 3)))
p
