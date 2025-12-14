rm(list=ls())

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(ggsci)

# load data
path <- './data/ex-fig7/b_Carbon_Metabolism.xlsx'
sheet.name <- 'IR'
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

agg_by_function <- rawdata %>%
    group_by(Function) %>%
    summarise(
        Gene_Count = sum(as.numeric(Gene_Count), na.rm = TRUE),
        Total_reads = sum(as.numeric(Total_reads), na.rm = TRUE)
    ) %>%
    arrange(desc(Gene_Count)) %>%
    ungroup() %>%
    as.data.frame()

head(agg_by_function)
# Function Gene_Count Total_reads
# 1   Reductive citrate cycle (Arnon-Buchanan cycle)       7141      804522
# 2              Dicarboxylate-hydroxybutyrate cycle       4102      469424
# 3                     3-Hydroxypropionate bi-cycle       3581      452440
# 4 Reductive pentose phosphate cycle (Calvin cycle)       2316      355704
# 5               Incomplete reductive citrate cycle       2151      267430
# 6          Hydroxypropionate-hydroxybutylate cycle       1820      170530

# prepare plotting data: order Functions by Gene_Count (largest first) and make y-axis show largest at top
agg_by_function <- agg_by_function %>%
    arrange(desc(Gene_Count))
# ggplot draws factor levels from bottom->top in the given order, so reverse to put largest at top
agg_by_function$Function <- factor(agg_by_function$Function, levels = rev(agg_by_function$Function))

# ensure output folder exists
if(!dir.exists("./plots")) dir.create("./plots", recursive = TRUE)

# size legend breaks (choose a few representative values)
size_breaks <- pretty(agg_by_function$Gene_Count, n = 4)

# leave space on the right so largest circles aren't clipped
x_max <- max(agg_by_function$Gene_Count, na.rm = TRUE)
x_limits <- c(0, x_max * 1.12)

# bubble plot: x = Gene_Count, y = Function, size = Gene_Count, fill = Total_reads
p <- ggplot(agg_by_function, aes(x = Gene_Count, y = Function, size = Gene_Count, fill = Total_reads)) +
    geom_point(shape = 21, color = "grey20", stroke = 0.9, alpha = 0.9) +
    # size legend with explicit breaks
    scale_size_area(max_size = 14, breaks = size_breaks,
                    guide = guide_legend(title = "Gene count",
                                         title.position = "top",
                                         label.position = "top",
                                         direction = "vertical",
                                         order = 1,
                                         override.aes = list(alpha = 0.9, shape = 21, color = "grey20", stroke = 0.9))) +
    # continuous color (fill) bar for Total_reads
    scale_fill_gradient(low = "#fde725", high = "#440154", name = "Total reads",
                        guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                               barwidth = unit(0.5, "cm"), barheight = unit(4, "cm"),
                                               direction = "vertical",
                                               order = 2)) +
    # expand x limits to avoid clipping of largest/smallest circles
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    labs(x = "Gene count", y = "", title = "") +
    # publication-ready theme: clean, larger readable sans font (use generic 'sans' to avoid missing system fonts)
    theme_bw(base_size = 14, base_family = "sans") +
    theme(
        panel.grid = element_blank(),
        # thicker panel border
        panel.border = element_rect(color = "black", size = 0.9, fill = NA),
        axis.text.x = element_text(size = 12, family = "sans", color = "black"),
        axis.text.y = element_text(size = 12, family = "sans", color = "black"),
        axis.title = element_text(size = 13, family = "sans"),
        plot.title = element_text(hjust = 0.5, size = 14, family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 12, family = "sans"),
        legend.text = element_text(size = 10, family = "sans"),
        # show legend keys with a thin border
        # legend.key = element_rect(fill = NA, colour = "black", size = 0.6),
        legend.spacing.y = unit(0.2, "cm"),
        plot.margin = unit(c(6, 8, 6, 6), "pt")
    )

print(p)

# save high-resolution figure as PDF suitable for journal submission
# use the base PDF device to avoid requiring Cairo/X11 on macOS where cairo may be unavailable
ggsave(filename = "./plots/Fig4_a_CarbonFixation_bubble.pdf", plot = p,
    width = 8, height = max(4, 0.2 * nrow(agg_by_function)), units = "in",
    device = "pdf", dpi = 300)
