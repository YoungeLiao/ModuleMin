rm(list=ls())

# ===== default setting =====
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales) # 用于格式化坐标轴
library(RColorBrewer) # 提供更专业的配色方案
library(showtext)
# font_add("Arial", "arial.ttf")  # Windows系统
font_add("Arial", "/System/Library/Fonts/Arial.ttf")  # Mac系统
showtext_auto()

# library(extrafont)
# font_import()
# ===== 1. load data =====
path.rawdata <- './data/rawdata.shared/asv_taxon.xlsx'
path.metadata <- './data/Fig3/e_metadata.xlsx'
output_dir <- "./output/Fig3/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 读取数据
rawdata <- data.frame(read_excel(path.rawdata))
metadata <- data.frame(read_excel(path.metadata))
# ===== 2. 数据处理 =====
# 设定目标属 - 可以根据需要修改
target_genus <- "g__Marimonas" 
# ROUND 1: # Aliidiomarina # Vicingus # Halomonas # Kangiella  # Aequorivita
# Enhanced: # Owenweeksia # Aliidiomarina # Vicingus # Idiomarina
# Inhibited: # Nitrosomonas # g__norank_o__Actinomarinales # g__norank_f__Paracoccaceae # g__Acidiferrimicrobium # g__Muriiphilus # g__unclassified_f__Flavobacteriaceae
# Neutral: # Candidatus_Scalindua # Marinobacter # unclassified_f__Paracoccaceae # Wenzhouxiangella # Marinicella # Nitrospira # Membranihabitans # unclassified_o__Aggregatilineales  # Aggregatilinea # Rhodanobacter # norank_f__A4b # IheB3-7


# 打印目标属的数据统计信息
cat("处理属:", target_genus, "\n")

# 筛选目标属的数据
rawdata.sub <- rawdata[rawdata$Genus == target_genus, ]

if(nrow(rawdata.sub) == 0) {
  stop("未找到指定的属: ", target_genus)
}

cat("找到", nrow(rawdata.sub), "条关于", target_genus, "的记录\n")

# 提取样本丰度数据（去除分类信息列）
taxonomy_cols <- c("ASV_ID", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
abundance_cols <- setdiff(colnames(rawdata.sub), taxonomy_cols)
abundance_data <- rawdata.sub[, abundance_cols]

# 按样本合并（同属下多个ASV求和）
abundance_by_sample <- colSums(abundance_data)

# 转为长表
df_long <- data.frame(Sample = names(abundance_by_sample), Abundance = as.numeric(abundance_by_sample))

# 合并分组和时间信息
df_long <- merge(df_long, metadata, by.x = "Sample", by.y = "Sample", all.x = TRUE)

# 转换Time为数值型
df_long$Time <- as.numeric(df_long$Time)

# 只保留有时间点和分组的样本
df_plot <- df_long[!is.na(df_long$Time) & !is.na(df_long$Group), ]

# 查看处理后的样本数据分布情况
sample_summary <- df_plot %>%
  group_by(Group, Time) %>%
  summarise(Sample_count = n())

print(sample_summary)

# 重命名处理组名为更科学的标签
df_plot$Group <- factor(df_plot$Group, 
                        levels = c("IR", "Dark"),
                        labels = c("Infrared", "Dark"))
df_plot <- df_plot[!is.na(df_plot$Group), ]

# ===== 3. 数据统计 =====
# 计算每组每个时间点的统计量
df_stat <- df_plot %>%
  group_by(Group, Time) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sd_abundance = sd(Abundance, na.rm = TRUE),
    se_abundance = sd_abundance / sqrt(n()), # 标准误
    n = n(),
    .groups = "drop"
  )
# ===== 4. 绘图 =====
# 设置CNS期刊风格的颜色方案
main_colors <- c("Infrared" = "#d55e00", "Dark" = "#0072b2")
point_shapes <- c("Infrared" = 16, "Dark" = 17)  # 实心圆和三角形
line_types <- c("Infrared" = "solid", "Dark" = "dashed")

# 创建包含原始数据点和阴影区域的时序图
p <- ggplot() +
  # 添加标准误差阴影区域
  geom_ribbon(data = df_stat, 
              aes(x = Time, 
                  ymin = mean_abundance - se_abundance, 
                  ymax = mean_abundance + se_abundance,
                  fill = Group),
              alpha = 0.2) +

  # 添加原始数据点 - 使用较小的点并增加透明度
  geom_point(data = df_plot, 
             aes(x = Time, y = Abundance, color = Group, shape = Group),
             size = 2, alpha = 0.4, position = position_dodge(width = 0.5)) +

  # 添加均值线
  geom_line(data = df_stat, 
            aes(x = Time, y = mean_abundance, color = Group, 
                group = Group, linetype = Group),
            size = 1.2) +

  # 添加均值点 - 使用较大的点
  geom_point(data = df_stat, 
             aes(x = Time, y = mean_abundance, color = Group, shape = Group, fill = Group),
             size = 3.5, stroke = 0.8) +

  # 设置颜色、形状和线型
  scale_color_manual(values = main_colors, name = '') +
  scale_shape_manual(values = point_shapes, name = '') +
  scale_linetype_manual(values = line_types, name = '') +
  scale_fill_manual(values = main_colors, name = '') +  # 新增：为阴影区域设置填充颜色

  # 坐标轴设置
  scale_x_continuous(breaks = sort(unique(df_stat$Time)),
                    expand = c(0.05, 0.05),
                    labels = function(x) paste0(x)) +
  labs(
    title = sub("g__", "", target_genus),  # CNS期刊通常不在图内放标题，而是用图注
    x = "Time (day)",
    y = paste("Abundance")
  ) +

  theme_pubr(base_size = 16) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 18, hjust = 0.5, face = 'bold'),
    legend.position = c(0.8, 0.9),
    legend.background = element_rect(fill = "transparent"), # 图例整体背景透明
    # legend.box.background = element_rect(fill = "transparent"), # 图例框背景透明
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )
print(p)


# 保存图形
ggsave(paste0(output_dir, "Fig3e_", sub("g__", "", target_genus), "_TimeSeries.pdf"),
       p, width = 4.5, height = 4.5, dpi = 300)
ggsave(paste0(output_dir, "Fig3e_", sub("g__", "", target_genus), "_TimeSeries.png"),
       p, width = 4.5, height = 4.5, dpi = 300)


# ===== 5. 输出统计数据 =====
# 将处理后的统计数据保存为文件，便于后续分析
write.csv(df_stat, 
          paste0(output_dir, sub("g__", "", target_genus), "_abundance_stats.csv"), 
          row.names = FALSE)

# 输出完成信息
cat("分析完成! 结果保存在:", output_dir, "\n")