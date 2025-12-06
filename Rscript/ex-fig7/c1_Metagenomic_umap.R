rm(list=ls())
# 加载必要包
library(readxl)
library(dplyr)
library(vegan)
library(umap)
library(factoextra)
library(NbClust)
library(ggplot2)
library(RColorBrewer)
library(Rtsne)

# 1. 读取数据
path.metageno.kegglevel3 <- "./data/ex-fig7/kegg_level3_profile.xlsx"
rawdata.metageno <- data.frame(read_excel(path.metageno.kegglevel3))
head(rawdata.metageno)

path.metadata <- "./ex-fig7/metadata.xlsx"
raw.metadata <- data.frame(read_excel(path.metadata))
head(raw.metadata)

col.name <- c('Level3', 'Level2', 'Level1', raw.metadata$Sample, 'Initial')
data.metageno <- rawdata.metageno[,col.name]
head(data.metageno)
# Level3                   Level2                               Level1
# 1                           Metabolic pathways Global and overview maps                           Metabolism
# 2        Biosynthesis of secondary metabolites Global and overview maps                           Metabolism
# 3 Microbial metabolism in diverse environments Global and overview maps                           Metabolism
# 4                    Biosynthesis of cofactors Global and overview maps                           Metabolism
# 5                  Biosynthesis of amino acids Global and overview maps                           Metabolism
# 6                         Two-component system      Signal transduction Environmental Information Processing
# Dark1    Dark2    Dark3      IR1      IR2      IR3  Initial
# 1 10820564 11484780 12612242 11307064 11260134 12381520 10737014
# 2  4581514  4871118  5283198  4766600  4780884  5259244  4569202
# 3  3068044  3178536  3451974  3085932  3088010  3330178  3086174
# 4  1830192  2000574  2308630  2102100  1992414  2305562  1747702
# 5  1568272  1657642  1929848  1744198  1632056  1882302  1543600
# 6  1501020  1611196  1824202  1639728  1539104  1745284  1772194

# 2. 构建丰度矩阵
abun_mat <- as.matrix(data.metageno[, c(raw.metadata$Sample,'Initial')])
rownames(abun_mat) <- data.metageno$Level3

# 计算分组平均丰度
# 假设metadata有一列Group，内容为"Dark"或"IR"
# 如果没有Group列，请先在metadata中添加
if(!"Group" %in% colnames(raw.metadata)) {
  stop("metadata中缺少Group列，请添加分组信息（如'Dark'/'IR'）")
}
dark_samples <- raw.metadata$Sample[raw.metadata$Group == "Dark"]
ir_samples <- raw.metadata$Sample[raw.metadata$Group == "IR"]

# 计算平均值
mean_dark <- rowMeans(abun_mat[, dark_samples, drop=FALSE], na.rm = TRUE)
mean_ir <- rowMeans(abun_mat[, ir_samples, drop=FALSE], na.rm = TRUE)
mean_all <- rowMeans(abun_mat, na.rm = TRUE)

# 数据预处理：过滤低丰度pathway
# 计算每个pathway的平均丰度，过滤掉丰度较低的
abundance_threshold <- quantile(mean_all, 0.01)  # 保留前99%的pathway，可调整
keep_pathways <- mean_all >= abundance_threshold
abun_mat_filtered <- abun_mat[keep_pathways, ]

cat("原始pathway数量:", nrow(abun_mat), "\n")
cat("过滤后pathway数量:", nrow(abun_mat_filtered), "\n")
cat("丰度阈值:", abundance_threshold, "\n")


# ==== 1. 先对All分组做降维和聚类，获得Cluster标签 ====
# 可选："All"（全部样品，包括Initial）、"Dark"、"IR"
sample_group_for_umap <- "All"  # 修改为 "Dark" 或 "IR" 可只用对应分组

# 先对All分组做降维和聚类，获得Cluster标签
selected_samples_all <- colnames(abun_mat_filtered)
abun_mat_log_all <- log2(abun_mat_filtered[, selected_samples_all, drop = FALSE] + 1)

dim_reduction_method <- "UMAP"  # 可选: "UMAP" 或 "TSNE"
if(dim_reduction_method == "UMAP") {
  library(umap)
  set.seed(0)
  umap_out_all <- umap(
    abun_mat_log_all,
    n_components = 2,
    n_neighbors = 30,
    min_dist = 11,
    spread = 16.5
  )
  coords_all <- umap_out_all$layout
  colnames(coords_all) <- c('Dim1', 'Dim2')
} else if(dim_reduction_method == "TSNE") {
  library(Rtsne)
  set.seed(0)
  tsne_out_all <- Rtsne(
    abun_mat_log_all,
    dims = 2,
    perplexity = 30,
    theta = 0.5,
    max_iter = 1000
  )
  coords_all <- tsne_out_all$Y
  colnames(coords_all) <- c('Dim1', 'Dim2')
} else {
  stop("dim_reduction_method must be either 'UMAP' or 'TSNE'")
}

# 聚类（只对All分组做）
df.reduced_all <- data.frame(coords_all)
colnames(df.reduced_all) <- c('UMAP1', 'UMAP2')
umap_dist_all <- dist(df.reduced_all)
nbclust_result_all <- NbClust(df.reduced_all, distance = "euclidean", min.nc = 2, max.nc = 8, method = "kmeans")
optimal_clusters_all <- as.numeric(nbclust_result_all$Best.nc[1])
set.seed(0)
kmeans_result_all <- kmeans(df.reduced_all, centers = optimal_clusters_all)
df.reduced_all$Cluster <- as.factor(kmeans_result_all$cluster)

# ========= 基于Cluster label 调整UMAP最佳参数 START（只对All分组） =========
## --- 已调整好后的
best_params <- data.frame(n_neighbors = 30)
best_params$min_dist <- 11
best_params$spread <- 16.5

cat('Best params:\n')
print(best_params)
# cat('Best silhouette score:', best_score, '\n')

# ========= 2. 根据用户选择的分组进行降维（聚类标签始终用All分组的） =========

if (sample_group_for_umap == "All") {
  # 直接使用All分组的降维和聚类结果
  df.reduced <- df.reduced_all
  dim_coords <- as.matrix(df.reduced_all[, c("UMAP1", "UMAP2")])
  cluster_labels <- df.reduced_all$Cluster
} else if (sample_group_for_umap == "Dark") {
  selected_samples <- c('Dark1', 'Dark2', 'Dark3', 'Initial')
  abun_mat_log <- log2(abun_mat_filtered[, selected_samples, drop = FALSE] + 1)
  # 降维（不聚类）
  if(dim_reduction_method == "UMAP") {
    set.seed(0)
    umap_out <- umap(
      abun_mat_log,
      n_components = 2,
      n_neighbors = best_params$n_neighbors,
      min_dist = best_params$min_dist,
      spread = best_params$spread
    )
    coords <- umap_out$layout
    colnames(coords) <- c('UMAP1', 'UMAP2')
  } else if(dim_reduction_method == "TSNE") {
    set.seed(0)
    tsne_out <- Rtsne(
      abun_mat_log,
      dims = 2,
      perplexity = 30,
      theta = 0.5,
      max_iter = 1000
    )
    coords <- tsne_out$Y
    colnames(coords) <- c('UMAP1', 'UMAP2')
  }
  df.reduced <- data.frame(coords)
  # 聚类标签用All分组的
  df.reduced$Cluster <- df.reduced_all$Cluster
  dim_coords <- as.matrix(df.reduced[, c("UMAP1", "UMAP2")])
  cluster_labels <- df.reduced$Cluster
} else if (sample_group_for_umap == "IR") {
  selected_samples <- c('IR1', 'IR2', 'IR3')
  abun_mat_log <- log2(abun_mat_filtered[, selected_samples, drop = FALSE] + 1)
  # 降维（不聚类）
  if(dim_reduction_method == "UMAP") {
    set.seed(0)
    umap_out <- umap(
      abun_mat_log,
      n_components = 2,
      n_neighbors = best_params$n_neighbors,
      min_dist = best_params$min_dist,
      spread = best_params$spread
    )
    coords <- umap_out$layout
    colnames(coords) <- c('UMAP1', 'UMAP2')
  } else if(dim_reduction_method == "TSNE") {
    set.seed(0)
    tsne_out <- Rtsne(
      abun_mat_log,
      dims = 2,
      perplexity = 30,
      theta = 0.5,
      max_iter = 1000
    )
    coords <- tsne_out$Y
    colnames(coords) <- c('UMAP1', 'UMAP2')
  }
  df.reduced <- data.frame(coords)
  # 聚类标签用All分组的
  df.reduced$Cluster <- df.reduced_all$Cluster
  dim_coords <- as.matrix(df.reduced[, c("UMAP1", "UMAP2")])
  cluster_labels <- df.reduced$Cluster
} else {
  stop("sample_group_for_umap 只能为 'All', 'Dark', 'IR'")
}

# ========= 7. 合并标签
filtered_indices <- which(rownames(abun_mat) %in% rownames(abun_mat_filtered))
df.reduced$pathways <- rownames(abun_mat_filtered)
df.reduced$Level1 <- data.metageno$Level1[filtered_indices]
df.reduced$Level2 <- data.metageno$Level2[filtered_indices]
df.reduced$Level3 <- data.metageno$Level3[filtered_indices]

# 新增三列丰度统计
df.reduced$Mean_Dark <- mean_dark[rownames(df.reduced)]
df.reduced$Mean_IR <- mean_ir[rownames(df.reduced)]
df.reduced$Mean_All <- mean_all[rownames(df.reduced)]
df.reduced$Mean_Initial <- data.frame(abun_mat_filtered)$Initial
df.reduced$Cluster <- paste('Cluster',df.reduced$Cluster)
# 8. 可视化（CNS正刊风格美化，仿照Cell/Nature新细胞分型UMAP风格）

library(ggtext)    # 富文本标题和图例
library(showtext)  # 矢量字体
showtext_auto()

# 自定义正刊风格配色（Cluster分组，推荐5-8色，柔和但有区分度）
cluster_colors <- c(
  "#6A8DCD", # 柔和蓝
  "#F17C67", # 柔和橙红
  "#FFD491", # 柔和黄
  "#A98BC5", # 柔和紫
  "#6EC6CA", # 柔和蓝绿
  "#F6A6B2", # 柔和粉
  "#A3A7C2", # 柔和灰蓝
  "#E6A0C4"  # 柔和粉紫
)
n_clusters <- length(levels(df.reduced_all$Cluster))
if (n_clusters > length(cluster_colors)) {
  cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
} else {
  cluster_colors <- cluster_colors[1:n_clusters]
}

# (1) Cluster分组上色（Cell/Nature风格，实心点，柔和配色，极简美化）
p <- ggplot(df.reduced, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 20, face = "plain", hjust = 0.5, family = "Arial"),
    axis.title = element_text(size = 16, face = "plain", family = "Arial"),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 16, family = "Arial", face = "plain"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = paste0(sample_group_for_umap, " (Metagenomics)"),
    x = NULL,
    y = NULL,
    color = "Cluster"
  )

print(p)

# (2) 丰度（Mean_Initial）上色，仿照Cell/Nature UMAP风格（实心点，连续色带，极简）
df.reduced$log_Mean_Dark <- log10(df.reduced$Mean_Dark + 1)
df.reduced$log_Mean_All <- log10(df.reduced$Mean_All + 1)
df.reduced$log_Mean_IR <- log10(df.reduced$Mean_IR + 1)
df.reduced$log_Mean_Initial <- log10(df.reduced$Mean_Initial + 1)

# 自定义连续色带（正刊常用蓝-紫-橙色带，或viridis/cividis）
abun_palette <- colorRampPalette(c('#FFFBEF', '#FFE79E', '#FF9800', "#E54091",  "#93206B" ))(100)

p_abun <- ggplot(df.reduced, aes(x = UMAP1, y = UMAP2, color = log_Mean_Initial)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_color_gradientn(
    colors = abun_palette,
    name = expression("log (RPKM)")
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 20, face = "plain", hjust = 0.5, family = "Arial"),
    axis.title = element_text(size = 16, face = "plain", family = "Arial"),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 16, family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = paste0(sample_group_for_umap, " abundance"),
    x = NULL,
    y = NULL
  )

print(p_abun)

# 9. 保存结果（将行名保存为'pathway'列）
write.csv(df.reduced, "./ex-fig7/umap_kmeans_metageno.csv", row.names = FALSE)
ggsave("./results_data/Fig3/umap_kmeans_metageno.pdf", p, width = 7, height = 6)


# --- Method 2. Merge four figures ---
library(ggpubr)

# 分别为三种丰度绘制UMAP图，并将标题居中
# 使用CNS期刊常用的Viridis配色方案，对比度更强
p_dark <- ggplot(df.reduced, aes(x = UMAP1, y = UMAP2, color = log_Mean_Dark)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "log (RPKM)") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)  # 居中，不加粗
  ) +
  labs(title = "Mean_Dark", color = "log")

p_ir <- ggplot(df.reduced, aes(x = UMAP1, y = UMAP2, color = log_Mean_IR)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "log (RPKM)") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)  # 居中，不加粗
  ) +
  labs(title = "Mean_IR", color = "log")

p_all <- ggplot(df.reduced, aes(x = UMAP1, y = UMAP2, color = log_Mean_All)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "log (RPKM)") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)  # 居中，不加粗
  ) +
  labs(title = "Mean_All", color = "log")

p_initial <- ggplot(df.reduced, aes(x = UMAP1, y = UMAP2, color = log_Mean_Initial)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", direction = -1, name = "log (RPKM)") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)  # 居中，不加粗
  ) +
  labs(title = "Mean_Initial", color = "log")

# 合并四图
p_abun_combined <- ggarrange(
  p_initial, p_dark, p_ir, p_all, 
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "right"
)

print(p_abun_combined)


## ================ Metabolomics ================
### ----------- Intracellular -----------------
# 1. 读取代谢组数据
path.metabolomic.kegglevel3 <- "./data/ex-fig7/pathway_abundance.txt" # pathway_abundance_processed.txt
rawdata.metabolo <- read.csv(path.metabolomic.kegglevel3, sep = "\t", check.names = FALSE)
# pathway, pathway_description, D2.3, D2.2, D2.1, F2.3, I2.2, I2.1, ...

# 2. 读取代谢组metadata
path.metadata.metabolomics <- "./data/ex-fig7/metadata.metabolite.csv"
rawdata.metadata.metabolomics <- read.csv(path.metadata.metabolomics, stringsAsFactors = FALSE)
# samples, group

# 3. 计算三种均值（所有、Dark、IR）
# 提取样本列名
sample_cols <- intersect(colnames(rawdata.metabolo), rawdata.metadata.metabolomics$samples)
# 按metadata分组
group_info <- rawdata.metadata.metabolomics
rownames(group_info) <- group_info$samples

# 计算均值
all_mean <- rowMeans(rawdata.metabolo[, sample_cols, drop = FALSE], na.rm = TRUE)
dark_samples <- group_info$samples[group_info$group == "Dark"]
ir_samples <- group_info$samples[group_info$group == "IR"]
mean_dark <- if(length(dark_samples) > 0) rowMeans(rawdata.metabolo[, intersect(dark_samples, sample_cols), drop = FALSE], na.rm = TRUE) else rep(NA, nrow(rawdata.metabolo))
mean_ir   <- if(length(ir_samples) > 0) rowMeans(rawdata.metabolo[, intersect(ir_samples, sample_cols), drop = FALSE], na.rm = TRUE) else rep(NA, nrow(rawdata.metabolo))

rawdata.metabolo$Mean_All  <- all_mean
rawdata.metabolo$Mean_Dark <- mean_dark
rawdata.metabolo$Mean_IR   <- mean_ir

# 3.1 取对数，避免颜色差异过大
# 加1防止log(0)
rawdata.metabolo$log_Mean_All  <- log10(rawdata.metabolo$Mean_All + 1)
rawdata.metabolo$log_Mean_Dark <- log10(rawdata.metabolo$Mean_Dark + 1)
rawdata.metabolo$log_Mean_IR   <- log10(rawdata.metabolo$Mean_IR + 1)

# 4. 筛选丰度最高的前20个代谢通路
# 基于所有样本的平均丰度进行排序
top10_pathways <- rawdata.metabolo[order(rawdata.metabolo$Mean_All, decreasing = TRUE), ]
top10_pathways <- top10_pathways[1:20, ]

# 5. 合并代谢组和宏基因组数据
# 以pathway_description为key合并
metabolo_plotdata <- rawdata.metabolo[, c("pathway_description", "log_Mean_All", "log_Mean_Dark", "log_Mean_IR")]
colnames(metabolo_plotdata) <- c("pathway_description", "Metab_log_Mean_All", "Metab_log_Mean_Dark", "Metab_log_Mean_IR")

# 宏基因组绘图数据（假设df.reduced有pathway_description, Mean_All, Mean_Dark, Mean_IR, UMAP1, UMAP2等）
# merge，指定宏基因组的Level3列和代谢组的pathway_description列对应
plotdata_merged <- merge(
  df.reduced, 
  metabolo_plotdata, 
  by.x = "Level3", 
  by.y = "pathway_description", 
  all.x = TRUE, all.y = FALSE
)
# 若代谢组没有匹配到的通路，Metab_log_Mean_*为NA
# 若需要NA填0，可加：
plotdata_merged$Metab_log_Mean_All[is.na(plotdata_merged$Metab_log_Mean_All)] <- 0
plotdata_merged$Metab_log_Mean_Dark[is.na(plotdata_merged$Metab_log_Mean_Dark)] <- 0
plotdata_merged$Metab_log_Mean_IR[is.na(plotdata_merged$Metab_log_Mean_IR)] <- 0

# 6. 为前10个通路添加标签标识
plotdata_merged$is_top10 <- plotdata_merged$Level3 %in% top10_pathways$pathway_description
plotdata_merged$label <- ifelse(plotdata_merged$is_top10, plotdata_merged$Level3, "")

# plotdata_merged 现在可用于后续联合绘图

# 计算Metab_log_Mean_Dark和Metab_log_Mean_IR的全局颜色范围
metab_min <- min(plotdata_merged$Metab_log_Mean_Dark, plotdata_merged$Metab_log_Mean_IR, na.rm = TRUE)
metab_max <- max(plotdata_merged$Metab_log_Mean_Dark, plotdata_merged$Metab_log_Mean_IR, na.rm = TRUE)

library(ggpubr)
library(ggrepel)  # 添加ggrepel包用于智能标签定位

# 绘制Metab_log_Mean_Dark
p_abun_metabolite_dark <- ggplot(plotdata_merged, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Metab_log_Mean_Dark), size = 2.5, alpha = 0.8) +
  # 为前10个通路添加标签，使用ggrepel避免重叠
  geom_text_repel(aes(label = label), 
                  data = subset(plotdata_merged, is_top10 == TRUE),
                  size = 3, color = "black", 
                  fontface = "plain",  # 移除加粗
                  max.overlaps = 20,   # 允许更多重叠以避免标签被移除
                  min.segment.length = 0.1,  # 最小线段长度
                  segment.color = "gray50",  # 延长线颜色
                  segment.size = 0.5,        # 延长线粗细
                  box.padding = 0.5,         # 标签框内边距
                  point.padding = 0.2) +     # 点到标签的距离
  scale_color_gradientn(
    colors = c("#F8F4FF",  "#9370DB", "#8A2BE2", "#4B0082"),  # 7个紫色渐变 "#E6E6FA", "#D8BFD8","#C8A2C8",
    name = "log Abundance (a.u.)", 
    limits = c(metab_min, metab_max), 
    oob = scales::squish
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Metabolites of Dark", 
    x = "UMAP1", y = "UMAP2"
  )

# 绘制Metab_log_Mean_IR
p_abun_metabolite_ir <- ggplot(plotdata_merged, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Metab_log_Mean_IR), size = 2.5, alpha = 0.8) +
  # 为前10个通路添加标签，使用ggrepel避免重叠
  geom_text_repel(aes(label = label), 
                  data = subset(plotdata_merged, is_top10 == TRUE),
                  size = 3, color = "black", 
                  fontface = "plain",  # 移除加粗
                  max.overlaps = 20,   # 允许更多重叠以避免标签被移除
                  min.segment.length = 0.1,  # 最小线段长度
                  segment.color = "gray50",  # 延长线颜色
                  segment.size = 0.5,        # 延长线粗细
                  box.padding = 0.5,         # 标签框内边距
                  point.padding = 0.2) +     # 点到标签的距离
  scale_color_gradientn(
    colors = c("#F8F4FF",   "#9370DB", "#8A2BE2", "#4B0082"),  # 7个紫色渐变 "#E6E6FA","#D8BFD8","#C8A2C8",
    name = "log Abundance (a.u.)", 
    limits = c(metab_min, metab_max), 
    oob = scales::squish
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Metabolites of IR", 
    x = "UMAP1", y = "UMAP2"
  )

# 并排展示两张图，方便对比
p_abun_metabolite_combined <- ggarrange(
  p_abun_metabolite_dark, 
  p_abun_metabolite_ir, 
  ncol = 2, nrow = 1, 
  common.legend = TRUE, legend = "right"
)

print(p_abun_metabolite_combined)
