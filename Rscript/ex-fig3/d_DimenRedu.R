rm(list=ls())
# source('./Rscript/functions.R')

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(ggsci)

# load data
path <- './data/rawdata.shared/asv_taxon.xlsx'
sheet.name <- 'asv_taxon'
# path <- "./data/Fig3/16S/asv_taxon.xlsx"
# sheet_name <- "asv_taxon"
rawdata <- data.frame(read_excel(path, sheet = sheet_name))

## --- filtering samples ---
target_cols <- c('ASV_ID','Dark1', 'Dark2', 'Dark3', 'IR1', 'IR2', 'IR3', 'Initial',"Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rawdata <- rawdata[,target_cols]

# 1. 质控：去除丰度极低的ASV（如均值小于0.1）
abundance_cols <- which(
  !(colnames(rawdata) %in% c(
    "ASV_ID", "Domain", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  ))
)
rawdata$mean_abundance <- rowMeans(rawdata[, abundance_cols], na.rm = TRUE)
rawdata_qc <- rawdata %>% filter(mean_abundance > 1) # 0.5
core_mags <- rawdata_qc

# 2. 计算距离矩阵 (Bray-Curtis)
set.seed(0)
core_mags_numeric <- core_mags[, abundance_cols]
core_mags_numeric <- as.data.frame(sapply(core_mags_numeric, as.numeric))
dist_matrix <- vegdist(core_mags_numeric, method = "bray")

# 3. 主坐标分析 (PCoA)
set.seed(0)
pcoa <- cmdscale(dist_matrix, k = 3, eig = TRUE)
scores <- pcoa$points  # 降维坐标

# 4. UMAP进一步降维 (基于PCoA结果)
library(umap)
set.seed(0)
# ## --- ##
# best_score <- -Inf
# best_params <- list()
# for (n_neighbors in seq(5, 50, by = 5)) {
#   for (n_components in 2:5) {
#     set.seed(0)
#     umap_out <- umap(
#       scores,
#       n_neighbors = n_neighbors,
#       n_components = n_components
#     )
#     umap_coords <- umap_out$layout[, 1:2] # 只用前两维聚类和评估
#     # k-means 聚类，遍历不同的聚类数
#     for (k in 2:10) {
#       kmeans_result <- kmeans(data.frame(umap_coords), centers = k)
#       # 计算 silhouette 得分
#       sil <- silhouette(kmeans_result$cluster, dist(umap_coords))
#       avg_sil <- mean(sil[, 3])
#       if (avg_sil > best_score) {
#         best_score <- avg_sil
#         best_params <- list(
#           n_neighbors = n_neighbors,
#           n_components = n_components,
#           centers = k
#         )
#       }
#     }
#   }
# }
# print(best_params)
# n_comp <- best_params$n_components
# n_nei <- best_params$n_neighbors
# n_cen <- best_params$centers
## --- ##
n_comp <- 2
n_nei <- 8 # the higher n_nei indicators how many neighber we need to form a cluster, which directly impact the closeness and distribution of different clusters. 
umap_out <- umap(scores, n_components=n_comp, n_neighbors = n_nei)
# umap_out <- umap(as.matrix(dist_matrix), n_components=8)
umap_coords <- umap_out$layout

# # 5. 聚类 (HDBSCAN处理密度不均)
# library(dbscan)
# set.seed(0)
# clusters <- hdbscan(umap_coords, minPts=25)$cluster

df.reduced <- data.frame(umap_coords)
colnames(df.reduced) <- c('umap1', 'umap2')
# --- Perform clustering on UMAP results ---
library(factoextra)
library(NbClust)
set.seed(0)
umap_dist <- dist(df.reduced[, c('umap1', 'umap2')])  # Calculate distance matrix

# Use NbClust to determine the optimal number of clusters
set.seed(0)
nbclust_result <- NbClust(df.reduced[, c('umap1', 'umap2')], distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
optimal_clusters <- as.numeric(nbclust_result$Best.nc[1])  # Extract the optimal number of clusters

# Debugging: Check optimal_clusters
print(paste("Optimal number of clusters:", optimal_clusters))

# Ensure optimal_clusters is valid
if (is.null(optimal_clusters) || optimal_clusters <= 0) {
  stop("Error: optimal_clusters is NULL or invalid.")
}

# Perform k-means clustering
set.seed(0)
kmeans_result <- kmeans(df.reduced[, c('umap1', 'umap2')], centers = optimal_clusters)
df.reduced$Cluster <- as.factor(kmeans_result$cluster)  # Add cluster information to the dataframe


# 6. 可视化
df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Cluster = as.factor(as.factor(kmeans_result$cluster))
)

# Only keep annotation rows that match the number of rows in df
annotation_cols <- c("ASV_ID", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (nrow(df) != nrow(rawdata_qc)) {
  annotation_df <- rawdata_qc[seq_len(nrow(df)), annotation_cols]
} else {
  annotation_df <- rawdata_qc[, annotation_cols]
}
umap_df <- cbind(df, annotation_df)
# umap_df <- cbind(df, rawdata_qc[, c("ASV_ID", "phylum", "class", "order", "family", "genus", "species")])

# --- 6.2 Annotation/Labeling ---
path.Anno <- "./data/Fig2/labels/Annotated_IROptogenetics.xlsx"
label_cols <- c("Phylum", "Genus", 'Total', 'Dark1', 'Dark2', 'Dark3', 'IR1', 'IR2', 'IR3')
umap_df <- labeling(umap_df, path.Anno, label_cols)

umap_df.mean <- merge(umap_df, rawdata_qc[,c('ASV_ID', 'mean_abundance', 'Dark1', 'Dark2', 'Dark3', 'IR1', 'IR2', 'IR3')], by = 'ASV_ID')

# 判断是否保存 umap_df.mean 至本地 results 文件夹
results_dir <- "./results_data/Fig2"
if (!dir.exists(results_dir)) dir.create(results_dir)
umap_df_mean_path <- file.path(results_dir, "umap_df.mean.csv")
if (!file.exists(umap_df_mean_path)) {
  write.csv(umap_df.mean, umap_df_mean_path, row.names = FALSE)
  message("umap_df.mean 已保存至: ", umap_df_mean_path)
} else {
  message("umap_df.mean 已存在，跳过保存。")
}

p1 <- ggplot(umap_df.mean, aes(x = UMAP1, y = UMAP2, color = Opto)) + # Cluster
  geom_point(size = 2, alpha = 0.8) +
  theme_classic()  + scale_color_manual(values = c('#C7253E',  pal_igv()(80), '#E5E5EA')) + # '#C7253E' ; pal_igv()(80)
  theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.background = element_rect(fill = 'transparent'),
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = 'right',
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 11),
    legend.background = element_blank(),
    legend.key.width = unit(0.75, "cm"),
    legend.key = element_blank()
  )
# + scale_color_brewer(palette = "Set1") 
p1
umap_df.mean$log_mean_abundance <- log1p(umap_df.mean$mean_abundance)
p2 <- ggplot(umap_df.mean, aes(UMAP1, UMAP2, colour = log_mean_abundance)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "", colour = "Log (Mean abundance)") + 
  theme_classic() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.background = element_rect(fill = 'transparent'),
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = 'right',
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    legend.key.width = unit(0.75, "cm"),
    legend.key = element_blank()
  )
p2

p3 <- ggplot(umap_df.mean, aes(x = UMAP1, y = UMAP2, color = Cluster)) + # Cluster
  geom_point(size = 2, alpha = 0.8) +
  theme_classic()  + scale_color_manual(values = pal_igv()(80)) + # '#C7253E' ; pal_igv()(80)
  theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.background = element_rect(fill = 'transparent'),
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = 'right',
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    legend.key.width = unit(0.75, "cm"),
    legend.key = element_blank()
  )
# + scale_color_brewer(palette = "Set1") 
p3
# 统计每个Cluster中Opto为'IR'的ASV占比
library(dplyr)

# 计算每个Cluster中Opto为'IR'的数量和总数
# cluster_ir_stats <- umap_df.mean %>%
#   group_by(Cluster) %>%
#   summarise(
#     IR_count = sum(Opto == "IR", na.rm = TRUE),
#     Total_count = n(),
#     IR_ratio = IR_count / Total_count
#   )
cluster_ir_stats <- umap_df.mean %>%
  group_by(Cluster) %>%
  summarise(
    IR_count = sum(Opto == "IR", na.rm = TRUE),
    Total_count = nrow(.),
    IR_ratio = IR_count / Total_count
  )


# 计算每个Cluster的平均丰度
cluster_ir_stats <- umap_df.mean %>%
  group_by(Cluster) %>%
  summarise(
    IR_count = sum(Opto == "IR", na.rm = TRUE),
    Total_count = nrow(.),
    IR_ratio = IR_count / Total_count,
    mean_abundance = mean(mean_abundance, na.rm = TRUE)
  )
print(cluster_ir_stats)
# IR_count Total_count  IR_ratio mean_abundance
# 1      164         506 0.3241107       68.53275

# --- test ---
# 确保Cluster列是因子类型
umap_df.mean$Cluster <- as.factor(umap_df.mean$Cluster)

# 获取所有唯一的Cluster值
unique_clusters <- unique(umap_df.mean$Cluster)

# 创建一个空的数据框来存储结果
cluster_ir_stats <- data.frame(
  Cluster = unique_clusters,
  IR_count = numeric(length(unique_clusters)),
  Total_count = numeric(length(unique_clusters)),
  IR_ratio = numeric(length(unique_clusters)),
  mean_abundance = numeric(length(unique_clusters))
)

# 为每个Cluster计算统计量
for (i in seq_along(unique_clusters)) {
  cluster_data <- umap_df.mean[umap_df.mean$Cluster == unique_clusters[i], ]
  cluster_ir_stats$IR_count[i] <- sum(cluster_data$Opto == "IR", na.rm = TRUE)
  cluster_ir_stats$Total_count[i] <- nrow(cluster_data)
  cluster_ir_stats$IR_ratio[i] <- cluster_ir_stats$IR_count[i] / cluster_ir_stats$Total_count[i]
  cluster_ir_stats$mean_abundance[i] <- mean(cluster_data$mean_abundance, na.rm = TRUE)
}

# 打印结果
print(cluster_ir_stats)

# --- test ---


# 可视化各Cluster含有IR的占比
library(ggplot2)
p_ir_ratio <- ggplot(cluster_ir_stats, aes(x = Cluster, y = IR_ratio * 100, fill = mean_abundance)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Mean abundance") +
  labs(
    title = "",
    x = "Cluster",
    y = "Proportion of IR (%)"
  ) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.background = element_rect(fill = 'transparent'),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'right',
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.background = element_blank(),
    legend.key.width = unit(0.7, "cm"),
    legend.key = element_blank()
  )
print(p_ir_ratio)
