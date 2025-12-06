# 加载必要的包
library(readxl)
library(dplyr)
library(broom)
library(ggpmisc) # 用于添加拟合方程和R²

# 加载数据
path <- './data/ex.fig1/d_CO2Rate.LinearFitting.xlsx'
sheet.name <- 'Sheet1' # 'MultiSludge' #
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

# 筛选Time >=4的数据
filtered_data <- rawdata %>% 
  filter(Time >= 4)


# 按Group分组进行HCO3随Time的线性拟合，并获取详细统计信息
fit_results <- filtered_data %>%
  group_by(Group) %>%
  do({
    model <- lm(HCO3 ~ Time, data = .)
    # 获取系数和统计信息
    coef_data <- tidy(model)
    glance_data <- glance(model)
    # 整合结果
    data.frame(
      Group = first(.$Group),
      Intercept = coef_data$estimate[coef_data$term == "(Intercept)"],
      Slope = coef_data$estimate[coef_data$term == "Time"],
      Intercept_SE = coef_data$std.error[coef_data$term == "(Intercept)"],
      Slope_SE = coef_data$std.error[coef_data$term == "Time"],
      Intercept_p = coef_data$p.value[coef_data$term == "(Intercept)"],
      Slope_p = coef_data$p.value[coef_data$term == "Time"],
      R_squared = glance_data$r.squared,
      Adj_R_squared = glance_data$adj.r.squared,
      p_value = glance_data$p.value
    )
  })

# 打印拟合结果表格
print("线性拟合结果:")
print(fit_results)

# 定义自定义色盘
palete <- c("#363062","#B31312", "#67729D", "#E7BCDE", "#BBE2EC",  "#F6EACB")

# 设置CNS风格的图形参数
theme_set(theme_classic() +
            theme(
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              legend.position = "right",
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              plot.margin = unit(c(1, 1, 1, 1), "lines"),
              plot.background = element_rect(fill = 'transparent', color = NA),
              panel.background = element_rect(fill = 'transparent'),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text = element_text(size = 12, face = "bold")
            ))

# 创建可视化图表
plot <- ggplot(filtered_data, aes(x = Time, y = HCO3, color = Group, fill = Group)) +
  geom_point(size = 3, shape = 21, stroke = 1, alpha = 0.8) +  # 添加数据点，使用带边框的圆点
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, size = 1.2) +  # 添加拟合线和置信区间
  facet_wrap(~ Group, scales = "free") +  # 按组分面
  labs(
    title = "",
    x = "Time",
    y = "HCO3"
  ) +
  scale_color_manual(values = palete) +  # 应用自定义色盘
  scale_fill_manual(values = palete) +   # 应用自定义色盘
  # 添加拟合方程和R²
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    parse = TRUE,
    size = 4,
    label.x = 0.05,
    label.y = 0.95
  ) +
  guides(color = guide_legend(override.aes = list(fill = palete, alpha = 0.5)))  # 调整图例

# 显示图表
print(plot)

# 保存图表（如果需要）
# ggsave("HCO3_vs_Time_LinearFit.png", plot = plot, width = 10, height = 8, dpi = 300, bg = "transparent")
