rm(list=ls())
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)

# ========== 1. 设置分类学层级选项 ==========
# 可选："Genus" 或 "Phylum"
tax_level <- "Genus" # 修改为 "Phylum" 可切换分析层级

# 1. 读取数据
path.AnnotatedIR <- './data/results_data/Annotated_IROptogenetics.xlsx'
raw.AnnotatedIR <- data.frame(read_excel(path.AnnotatedIR))

# 2. 需要统计的丰度列
abundance_cols <- c("Total", "IR1", "IR2", "IR3", "Dark1", "Dark2", "Dark3")

# 3. 统计
summary.counts <- as.data.frame(table(raw.AnnotatedIR[[tax_level]]))
colnames(summary.counts) <- c(tax_level, 'Gene_counts')
path.save.counts <- './data/ex-fig4/IROptogenetics_Genus_Counts.xlsx'
writexl::write_xlsx(summary.counts, path.save.counts)

summary.abundance <- aggregate(raw.AnnotatedIR[, abundance_cols], 
                              by = setNames(list(raw.AnnotatedIR[[tax_level]]), tax_level), 
                              FUN = sum, na.rm = TRUE)

summary.df <- merge(summary.counts, summary.abundance, by = tax_level)
summary.df$GeneCounts_normed <- summary.df$Gene_counts / summary.df$Total

# 4. 仅可视化gene counts为top20的
summary.df <- summary.df[order(summary.df$Gene_counts, decreasing = TRUE), ]
summary.df <- summary.df[1:min(20, nrow(summary.df)), ]
summary.df[[tax_level]] <- factor(summary.df[[tax_level]], levels = summary.df[[tax_level]])

# 5. CNS正刊风格气泡图美化（可自定义美化颜色）
custom_colors <- colorRampPalette(c("#344893", "#fee090",'#FF4066', "#D53555", '#AA2B44'))(100)

p <- ggplot(summary.df, aes(x = Gene_counts, 
                            y = !!as.name(tax_level), 
                            size = Total, 
                            color = GeneCounts_normed)) +
  geom_point(alpha = 0.6) +
  geom_text(aes(label = !!as.name(tax_level)), 
            hjust = -0.1, 
            vjust = 0.5, 
            fontface = "italic", 
            size = 4.5,
            color = "black") +
  scale_size_continuous(range = c(4, 18), name = "Abundance (RPKM)") +
  scale_color_gradientn(colors = custom_colors, name = "Normalized gene counts") +
  labs(
    x = "IR gene counts",
    y = NULL,
    title = paste0("IR optogenetic switch under ", tax_level)
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent", color = NA)
  ) +
  coord_cartesian(clip = "off")
p

# 6. 输出气泡图
ggsave(filename = paste0("Fig2b_IRswitch_Bubbleplot_", tax_level, ".pdf"), plot = p, width = 8, height = 8, limitsize = FALSE)

# 7. 绘制所选分类学层级gene counts占比的饼图

# 先统计所选层级的gene counts总和（用全数据，不只是summary.df的top20）
pie_df_all <- raw.AnnotatedIR %>%
  group_by(!!as.name(tax_level)) %>%
  summarise(Gene_counts = n()) %>%
  arrange(desc(Gene_counts))

# 计算占比
pie_df_all$Fraction <- pie_df_all$Gene_counts / sum(pie_df_all$Gene_counts)

# 只保留前20，其他合并为Other
top_n <- 20
if (nrow(pie_df_all) > top_n) {
  pie_df_top <- pie_df_all[1:top_n, ]
  other_sum <- sum(pie_df_all$Gene_counts[(top_n+1):nrow(pie_df_all)])
  other_fraction <- other_sum / sum(pie_df_all$Gene_counts)
  pie_df_other <- data.frame(
    label = "Other",
    Gene_counts = other_sum,
    Fraction = other_fraction
  )
  pie_df_top$label <- as.character(pie_df_top[[tax_level]])
  pie_df <- rbind(
    pie_df_top[, c("label", "Gene_counts", "Fraction")],
    pie_df_other
  )
} else {
  pie_df <- pie_df_all
  pie_df$label <- as.character(pie_df[[tax_level]])
  pie_df <- pie_df[, c("label", "Gene_counts", "Fraction")]
}

# 按照占比从高到低排序
pie_df <- pie_df %>% arrange(desc(Fraction))

# 生成带百分比的图例标签
pie_df$Label <- paste0(pie_df$label, " (", percent(pie_df$Fraction, accuracy = 0.1), ")")

# 根据tax_level选择渐变色盘
if (tax_level == "Phylum") {
  # Phylum层面用蓝-绿-黄渐变
  pie_colors <- rev(colorRampPalette(c("#1b9e77", "#a6d854", "#673AB7"))(nrow(pie_df)))
} else if (tax_level == "Genus") {
  # Genus层面用原有蓝-黄-红渐变
  pie_colors <- rev(colorRampPalette(c("#344893", "#fee090", "#D53555"))(nrow(pie_df)))
} else {
  # 其他层级用灰-紫-粉渐变
  pie_colors <- rev(colorRampPalette(c("#bdbdbd", "#8073ac", "#f1b6da"))(nrow(pie_df)))
}

# 设置Label为有序因子，确保图例和饼图顺序一致（从高到低）
pie_df$Label <- factor(pie_df$Label, levels = pie_df$Label)

# 绘制饼图，图例用Label
pie_plot <- ggplot(pie_df, aes(x = "", y = Gene_counts, fill = Label)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = pie_colors) +
  labs(title = "", fill = tax_level) +
  theme_void(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.title = element_text(size = 16,face = "bold"),
    legend.text = element_text(size = 14)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))

pie_plot

# 输出饼图
ggsave(filename = paste0("Fig2b_IRswitch_", tax_level, "_Pieplot.pdf"), plot = pie_plot, width = 14.5, height = 7)
