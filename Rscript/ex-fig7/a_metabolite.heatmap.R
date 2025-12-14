rm(list=ls())
# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(grid)

par(family = "Arial")

# load data
path <- './data/ex-fig7/a_membrane_metabolites.xlsx' 
sheet.name <- 'Intra'  # 'Intra' 
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

rawdata <- subset(rawdata, rawdata$Group != 'Y')
data <- as.matrix(rawdata[,c(3:length(rawdata))])
rownames(data) <- rawdata$Label
data <- t(data)
data <- log10(data)
summary(data)

# Perform clustering and save the result
library(pheatmap)
p <- pheatmap(data, cluster_row = FALSE, cluster_col = FALSE, show_rownames= TRUE,
              color=colorRampPalette(rev(c('#A7414A' ,"#FE7875",'#FFD0C9',"#A882C1")))(1000), # 
              # scale = 'column', 
              fontsize = 12,
              display_numbers = matrix(sprintf("%.2e", 10^data), ncol=ncol(data)),   # 注意这里改为ncol
              number_color = "white", 
              cutree_rows = 3,
              angle_col = 0,  # Rotate column names by 45 degrees
              legend = TRUE,
              legend_title = "log10(Norm.A)",  # Add legend title
              border_color = NA,  # Remove border color
              bg = "transparent",  # Set background to transparent
              gaps_col = NULL)  # Add gaps_col to adjust margins

# Extract clustering results
cluster_results <- cutree(p$tree_col, k = 10)  # Extract clusters for 10 groups
clustered_metabolites <- split(colnames(data), cluster_results)  # Group metabolites by cluster

# Add cluster information to the original data
cluster_mapping <- data.frame(Metabolite = colnames(data), Cluster = cluster_results)

# Save the updated data with cluster information to an xlsx file
library(writexl)
write_xlsx(cluster_mapping, "./data/Fig2a.metabolites.intra_clusters.xlsx")

### =========== Change font to Arial ===========
if (!dir.exists("./Figures/Fig4")) {
  dir.create("./Figures/Fig4", recursive = TRUE)
}

png("./results/Fig4/Fig4c.metabolite.heatmap.png", width = 12, height = 5, units = "in", res = 300, type = "cairo")

par(family = "Arial")

print(p)

