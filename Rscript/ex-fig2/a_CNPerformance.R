rm(list=ls())
# ===== default setting =====
library(readxl)
library(ggplot2)
library(ggpubr)
library(scales)
library(tidyr)
library(dplyr)
library(grid) # for unit()

# ===== 1. load data =====
path <- './data/ex-fig2/a_CNPerformance.xlsx'
sheet.name <- 'Performance'
rawdata <- data.frame(read_excel(path, sheet = sheet.name))
head(rawdata)

df <- rawdata[, c("Time", "Group", "NO3", "NO2", "NH4", "HCO3")]

df_long <- pivot_longer(df,
                        cols = c("NO3", "NO2", "NH4", "HCO3"),
                        names_to = "Indicator",
                        values_to = "Concentration")

head(df_long)
# # A tibble: 6 × 4
# Time Group Indicator Concentration
# <dbl> <chr> <chr>             <dbl>
#   1     0 Dark  NO3               0.543
# 2     0 Dark  NO2             127.   
# 3     0 Dark  NH4              89.8  
# 4     0 Dark  HCO3           1454.   
# 5     0 Dark  NO3               0.372
# 6     0 Dark  NO2             123. 

df_summary <- df_long %>%
  group_by(Time, Group, Indicator) %>%
  mutate(
    mean = mean(Concentration, na.rm = TRUE),
    sd = sd(Concentration, na.rm = TRUE)
  ) %>%
  distinct(Time, Group, Indicator, .keep_all = TRUE) %>%
  ungroup()


head(df_summary)

scale_factor <- max(df$HCO3, na.rm = TRUE) / max(c(df$NO3, df$NO2, df$NH4), na.rm = TRUE)

df_summary <- df_summary %>%
  mutate(
    mean_plot = ifelse(Indicator == "HCO3", mean / scale_factor, mean),
    sd_plot = ifelse(Indicator == "HCO3", sd / scale_factor, sd)
  )

custom_colors <- c(
  "NO3" = "#1B9E77",   
  "NO2" = "#D95F02",  
  "NH4" = "#7570B3", 
  "HCO3" = "#E7298A" 
)

custom_shapes <- c("NO3"=21, "NO2"=22, "NH4"=24)
custom_linetypes <- c("NO3"="solid", "NO2"="dashed", "NH4"="dotdash")

groups <- unique(df_summary$Group)

n_data <- df_summary[df_summary$Indicator %in% c("NO3", "NO2", "NH4"), ]
n_min_all <- min(n_data$mean_plot - 2*n_data$sd_plot, na.rm = TRUE)
n_max_all <- max(n_data$mean_plot + 2*n_data$sd_plot, na.rm = TRUE)
n_range_all <- n_max_all - n_min_all
n_min_plot_all <- max(0, n_min_all - 0.08 * n_range_all)
n_max_plot_all <- n_max_all + 0.15 * n_range_all

c_data <- df_summary[df_summary$Indicator == "HCO3", ]
c_min_all <- min(c_data$mean - 2*c_data$sd, na.rm = TRUE)
c_max_all <- max(c_data$mean + 2*c_data$sd, na.rm = TRUE)
c_range_all <- c_max_all - c_min_all
c_min_plot_all <- max(0, c_min_all - 0.08 * c_range_all)
c_max_plot_all <- c_max_all + 0.15 * c_range_all

plots <- list()
for (grp in groups) {
  df_grp <- df_summary[df_summary$Group == grp, ]
  
  p <- ggplot() +

    geom_col(
      data = subset(df_grp, Indicator == "HCO3"),
      aes(x = Time, y = mean_plot, fill = Indicator),
      width = 0.7,
      alpha = 0.25
    ) +

    geom_line(
      data = subset(df_grp, Indicator %in% c("NO3", "NO2", "NH4")),
      aes(x = Time, y = mean_plot, color = Indicator, group = Indicator, linetype = Indicator),
      size = 1.2
    ) +
    geom_point(
      data = subset(df_grp, Indicator %in% c("NO3", "NO2", "NH4")),
      aes(x = Time, y = mean_plot, color = Indicator, shape = Indicator),
      size = 3.5, stroke = 1.2, fill = "white"
    ) +
    geom_errorbar(
      data = subset(df_grp, Indicator %in% c("NO3", "NO2", "NH4")),
      aes(x = Time, ymin = mean_plot - sd_plot, ymax = mean_plot + sd_plot, color = Indicator),
      width = 0.18,
      size = 0.9
    ) +
    geom_errorbar(
      data = subset(df_grp, Indicator == "HCO3"),
      aes(x = Time, ymin = mean_plot - sd_plot, ymax = mean_plot + sd_plot),
      width = 0.5,
      size = 0.9,
      color = custom_colors["HCO3"]
    ) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_shape_manual(values = custom_shapes) +
    scale_y_continuous(
      name = expression(bold("N (mg·L"^{-1}*")")),
      limits = c(n_min_plot_all, n_max_plot_all),
      expand = expansion(mult = c(0, 0.02)),
      sec.axis = sec_axis(
        ~ . * scale_factor,
        name = expression(bold("C (mg·L"^{-1}*")")),
        breaks = pretty(c(c_min_plot_all, c_max_plot_all)),
        labels = function(x) sprintf("%.0f", x)
      )
    ) +
    scale_x_continuous(
      breaks = sort(unique(df_grp$Time)),
      labels = sort(unique(df_grp$Time))
    ) +
    labs(
      x = expression(bold("Time (day)")),
      color = NULL,
      fill = NULL,
      shape = NULL,
      linetype = NULL,
      title = grp
    ) +
    theme_classic(base_size = 16, base_family = "Arial") +
    theme(
      plot.margin = unit(c(0.7, 0.3, 0.3, 0.7), "cm"),
      plot.background = element_rect(fill = 'transparent', color = NA),
      panel.background = element_rect(fill = 'transparent'),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.line = element_line(size = 1.1, color = "black"),
      axis.text.x = element_text(size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16, color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 1.2)
    )
  
  plots[[as.character(grp)]] <- p
}


if (!requireNamespace("showtext", quietly = TRUE)) {
  install.packages("showtext")
}
library(showtext)
showtext_auto(enable = TRUE)
combined_plot <- NULL
if (all(c("Dark", "IR") %in% names(plots))) {
  combined_plot <- ggarrange(
    plots[["Dark"]], plots[["IR"]],
    ncol = 1, nrow = 2,
    labels = c("Dark", "IR"),
    font.label = list(size = 18, face = "bold", family = "Arial"),
    heights = c(1, 1),
    common.legend = TRUE,
    legend = "right",
    align = "hv"
  )
  print(combined_plot)
  ggsave(
    filename = "./Figures/Fig1/Fig1cd_Performance.pdf",
    plot = combined_plot,
    width = 6, height = 8, units = "in", dpi = 300, bg = "transparent"
  )
} else {
  for (p in plots) print(p)
}

