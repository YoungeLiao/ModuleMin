rm(list=ls())
library(tidyverse)
library(ggraph)
library(igraph)
library(GGally)

edges <- read.csv("./scr/Visualization/Fig3/network_data/ANM-Cluster6-20250618/network_edges.csv")
nodes <- read.csv("./scr/Visualization/Fig3/network_data/ANM-Cluster6-20250618/network_nodes.csv")
metadata <- read.csv("./scr/Visualization/Fig3/network_data/ANM-Cluster6-20250618/metadata.csv", row.names = "Sample")

net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
# ensure weight exists on edges and absolute weight for layout/width
if (!"weight" %in% edge_attr_names(net)) {
  E(net)$weight <- edges$weight
} else {
  E(net)$weight <- as.numeric(E(net)$weight)
}
E(net)$abs_weight <- abs(E(net)$weight)

edges$cor_type <- ifelse(edges$weight >= 0, "positive", "negative")
E(net)$cor_type <- edges$cor_type

out_dir <- "./scr/Visualization/Fig3/network_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

base_text_size <- 12
title_size <- 14
subtitle_size <- 11

V(net)$module <- as.character(V(net)$module)
if (is.null(V(net)$betweenness)) V(net)$betweenness <- igraph::betweenness(net, weights = NA)
V(net)$betweenness <- as.numeric(V(net)$betweenness)
if (is.null(V(net)$is_keystone)) V(net)$is_keystone <- FALSE
V(net)$is_keystone <- as.logical(V(net)$is_keystone)

set.seed(123)
set.seed(2)
layout_df <- create_layout(net, layout = "fr", weights = E(net)$abs_weight)

betw_thresh <- quantile(V(net)$betweenness, 0.90, na.rm = TRUE)
label_nodes <- which(V(net)$is_keystone | V(net)$betweenness >= betw_thresh)
labels_df <- layout_df[label_nodes, , drop = FALSE]

edge_colors <- c("positive" = "#2b83ba", "negative" = "#d7191c") 
module_colors <- scales::hue_pal()(max(1, length(unique(na.omit(V(net)$module)))))
edge_colors <- c("positive" = "#BCDDF7", "negative" = "#DA7696") 
module_colors <- c( "#1C3885","#9D90CB", "#F4A25C", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")  # 10色CNS标准配色

p_net_cns <- ggraph(layout_df) +
  geom_edge_link(aes(color = cor_type, width = abs(weight), alpha = abs(weight)),
                 show.legend = TRUE, lineend = "round") +
  scale_edge_color_manual(values = edge_colors, name = "Correlation") +
  scale_edge_width(range = c(0.2, 1.2), guide = guide_legend(order = 1), name = "Edge strength") +
  scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
  
  geom_node_point(aes(fill = factor(module), size = betweenness),
                  shape = 21, color = "black", stroke = 0.35) +
  scale_fill_manual(values = module_colors, na.value = "grey80", name = "Module") +
  scale_size_continuous(range = c(2, 8), name = "Betweenness") +
  
  geom_node_point(data = layout_df %>% dplyr::filter(as.logical(is_keystone)),
                  aes(x = x, y = y),
                  shape = 21, fill = NA, color = "black", stroke = 1, size = 6, inherit.aes = FALSE) +
  
  geom_node_text(data = labels_df,
                 aes(x = x, y = y, label = ifelse(!is.na(Genus) & Genus != "", Genus, name)),
                 size = 3.2, repel = TRUE, min.segment.length = 0, seed = 42) +
  
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

ggsave(filename = file.path(out_dir, "full_network_CNS.pdf"),
       plot = p_net_cns, width = 12, height = 10, units = "in", device = cairo_pdf)
ggsave(filename = file.path(out_dir, "full_network_CNS.tiff"),
       plot = p_net_cns, width = 12, height = 10, units = "in", dpi = 600)

max_labels_allowed <- 120
if (vcount(net) <= max_labels_allowed) {
  p_net_labels <- p_net_cns +
    geom_node_text(aes(label = ifelse(!is.na(Genus) & Genus != "", Genus, name)),
                   size = 2.6, repel = TRUE, max.overlaps = 40)
  ggsave(filename = file.path(out_dir, "full_network_all_labels_CNS.pdf"),
         plot = p_net_labels, width = 14, height = 12, units = "in", device = cairo_pdf)
}


