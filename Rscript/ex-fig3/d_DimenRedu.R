rm(list=ls())

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(ggsci)

# load data
path <- './data/rawdata.shared/asv_taxon.xlsx'
sheet.name <- 'asv_taxon'
rawdata <- data.frame(read_excel(path, sheet = sheet_name))

## --- filtering samples ---
target_cols <- c('ASV_ID','Dark1', 'Dark2', 'Dark3', 'IR1', 'IR2', 'IR3', 'Initial',"Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rawdata <- rawdata[,target_cols]

# 1. quality control
abundance_cols <- which(
  !(colnames(rawdata) %in% c(
    "ASV_ID", "Domain", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  ))
)
rawdata$mean_abundance <- rowMeans(rawdata[, abundance_cols], na.rm = TRUE)
rawdata_qc <- rawdata %>% filter(mean_abundance > 1) # 0.5
core_mags <- rawdata_qc

# 2. distence matrix
set.seed(0)
core_mags_numeric <- core_mags[, abundance_cols]
core_mags_numeric <- as.data.frame(sapply(core_mags_numeric, as.numeric))
dist_matrix <- vegdist(core_mags_numeric, method = "bray")

# 3. PCoA
set.seed(0)
pcoa <- cmdscale(dist_matrix, k = 3, eig = TRUE)
scores <- pcoa$points  # 降维坐标

# 4. UMAP
library(umap)
set.seed(0)


n_comp <- 2
n_nei <- 8 # the higher n_nei indicators how many neighber we need to form a cluster, which directly impact the closeness and distribution of different clusters. 
umap_out <- umap(scores, n_components=n_comp, n_neighbors = n_nei)
# umap_out <- umap(as.matrix(dist_matrix), n_components=8)
umap_coords <- umap_out$layout

# # 5. Clustering
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

# 6. Visualization
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

results_dir <- "./results_data/Fig2"
if (!dir.exists(results_dir)) dir.create(results_dir)
umap_df_mean_path <- file.path(results_dir, "umap_df.mean.csv")
if (!file.exists(umap_df_mean_path)) {
  write.csv(umap_df.mean, umap_df_mean_path, row.names = FALSE)
  message("umap_df.mean saved to: ", umap_df_mean_path)
} else {
  message("umap_df.mean existed，skip。")
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

p3

library(dplyr)

cluster_ir_stats <- umap_df.mean %>%
  group_by(Cluster) %>%
  summarise(
    IR_count = sum(Opto == "IR", na.rm = TRUE),
    Total_count = nrow(.),
    IR_ratio = IR_count / Total_count
  )


cluster_ir_stats <- umap_df.mean %>%
  group_by(Cluster) %>%
  summarise(
    IR_count = sum(Opto == "IR", na.rm = TRUE),
    Total_count = nrow(.),
    IR_ratio = IR_count / Total_count,
    mean_abundance = mean(mean_abundance, na.rm = TRUE)
  )
print(cluster_ir_stats)


umap_df.mean$Cluster <- as.factor(umap_df.mean$Cluster)

unique_clusters <- unique(umap_df.mean$Cluster)

cluster_ir_stats <- data.frame(
  Cluster = unique_clusters,
  IR_count = numeric(length(unique_clusters)),
  Total_count = numeric(length(unique_clusters)),
  IR_ratio = numeric(length(unique_clusters)),
  mean_abundance = numeric(length(unique_clusters))
)

for (i in seq_along(unique_clusters)) {
  cluster_data <- umap_df.mean[umap_df.mean$Cluster == unique_clusters[i], ]
  cluster_ir_stats$IR_count[i] <- sum(cluster_data$Opto == "IR", na.rm = TRUE)
  cluster_ir_stats$Total_count[i] <- nrow(cluster_data)
  cluster_ir_stats$IR_ratio[i] <- cluster_ir_stats$IR_count[i] / cluster_ir_stats$Total_count[i]
  cluster_ir_stats$mean_abundance[i] <- mean(cluster_data$mean_abundance, na.rm = TRUE)
}

print(cluster_ir_stats)

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
