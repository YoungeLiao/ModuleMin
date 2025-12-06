# # ------------- 3rd try - Tecent ----------------

library(tidyverse)
library(ggraph)
library(igraph)
library(GGally)

# 读取Python输出的数据
edges <- read.csv("./scr/Visualization/Fig3/network_data/ANM-Cluster6-20250618/network_edges.csv")
nodes <- read.csv("./scr/Visualization/Fig3/network_data/ANM-Cluster6-20250618/network_nodes.csv")
metadata <- read.csv("./scr/Visualization/Fig3/network_data/ANM-Cluster6-20250618/metadata.csv", row.names = "Sample")
# 创建igraph图对象并保证必要属性
net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# ensure weight exists on edges and absolute weight for layout/width
if (!"weight" %in% edge_attr_names(net)) {
  E(net)$weight <- edges$weight
} else {
  E(net)$weight <- as.numeric(E(net)$weight)
}
E(net)$abs_weight <- abs(E(net)$weight)

# edge color: 正相关为蓝（positive），负相关为红（negative）
edges$cor_type <- ifelse(edges$weight >= 0, "positive", "negative")
E(net)$cor_type <- edges$cor_type

# 输出目录（与脚本中其他保存路径一致）
out_dir <- "./scr/Visualization/Fig3/network_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 统一美学设置（CNS 风格：简洁、清晰、高分辨率）
base_text_size <- 12
title_size <- 14
subtitle_size <- 11

# 预处理节点属性，确保列存在且类型正确
V(net)$module <- as.character(V(net)$module)
if (is.null(V(net)$betweenness)) V(net)$betweenness <- igraph::betweenness(net, weights = NA)
V(net)$betweenness <- as.numeric(V(net)$betweenness)
if (is.null(V(net)$is_keystone)) V(net)$is_keystone <- FALSE
V(net)$is_keystone <- as.logical(V(net)$is_keystone)

# 选择布局（使用Fruchterman-Reingold，按边权重绝对值加权以兼容负权）
set.seed(123)
layout_df <- create_layout(net, layout = "fr", weights = E(net)$abs_weight)

# 选择labels：仅标注基石与高中心性节点（>=90%）以保持图面简洁
betw_thresh <- quantile(V(net)$betweenness, 0.90, na.rm = TRUE)
label_nodes <- which(V(net)$is_keystone | V(net)$betweenness >= betw_thresh)
labels_df <- layout_df[label_nodes, , drop = FALSE]

# 颜色方案（对比强烈且对色盲相对友好）
edge_colors <- c("positive" = "#2b83ba", "negative" = "#d7191c") # blue / red
# 模块填充色使用Brewer Set3（若模块过多需自定义）
module_colors <- scales::hue_pal()(max(1, length(unique(na.omit(V(net)$module)))))

# CNS标准配色方案
edge_colors <- c("positive" = "#1f77b4", "negative" = "#d62728")  # 深蓝色和深红色
module_colors <- c("#8785C1", "#A5DBD9", "#FFB17D", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")  # 10色CNS标准配色



# 主图：整网（适度显示边宽/透明度，突出基石节点）
p_net_cns <- ggraph(layout_df) +
  # 边：按绝对权重调整宽度与透明度
  geom_edge_link(aes(color = cor_type, width = abs(weight), alpha = abs(weight)),
                 show.legend = TRUE, lineend = "round") +
  scale_edge_color_manual(values = edge_colors, name = "Correlation") +
  scale_edge_width(range = c(0.2, 1.2), guide = guide_legend(order = 1), name = "Edge strength") +
  scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
  # 节点底层：以模块着色，黑色细边框
  geom_node_point(aes(fill = factor(module), size = betweenness),
                  shape = 21, color = "black", stroke = 0.35) +
  scale_fill_manual(values = module_colors, na.value = "grey80", name = "Module") +
  scale_size_continuous(range = c(2, 8), name = "Betweenness") +
  # 基石节点强调：更粗边框与稍大尺寸
  geom_node_point(data = layout_df %>% dplyr::filter(as.logical(is_keystone)),
                  aes(x = x, y = y),
                  shape = 21, fill = NA, color = "black", stroke = 1, size = 6, inherit.aes = FALSE) +
  # 标签：只标注基石与高中心性节点，避免过度拥挤，使用 repel
  geom_node_text(data = labels_df,
                 aes(x = x, y = y, label = ifelse(!is.na(Genus) & Genus != "", Genus, name)),
                 size = 3.2, repel = TRUE, min.segment.length = 0, seed = 42) +
  # 主题与图注布局
  theme_void(base_size = base_text_size) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = base_text_size),
    plot.title = element_text(size = title_size, face = "bold"),
    plot.subtitle = element_text(size = subtitle_size),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Microbial Interaction Network",
    subtitle = paste0("Nodes: ", vcount(net), "    Edges: ", ecount(net))
  )
p_net_cns

# 保存高分辨率文件（PDF + TIFF）
ggsave(filename = file.path(out_dir, "full_network_CNS.pdf"),
       plot = p_net_cns, width = 12, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = file.path(out_dir, "full_network_CNS.tiff"),
       plot = p_net_cns, width = 12, height = 10, units = "in", dpi = 600)

# 另存一份：更强调所有标签（仅在节点不太多时使用）
max_labels_allowed <- 120
if (vcount(net) <= max_labels_allowed) {
  p_net_labels <- p_net_cns +
    geom_node_text(aes(label = ifelse(!is.na(Genus) & Genus != "", Genus, name)),
                   size = 2.6, repel = TRUE, max.overlaps = 40)
  ggsave(filename = file.path(out_dir, "full_network_all_labels_CNS.pdf"),
         plot = p_net_labels, width = 14, height = 12, units = "in", device = cairo_pdf)
}

# 打印提示
message("CNS-style network figures saved to: ", normalizePath(out_dir))
