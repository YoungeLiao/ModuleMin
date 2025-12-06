rm(list=ls())
# source('./Rscript/functions.R')
# load library
library(readxl)
library(ggplot2)
library(ggpubr)

# load data
path <- './data/Fig2/a_CO2.tic.xlsx' # Fig1a.cO2.tic.xlsx # Fig1d.CO2.xlsx
sheet.name <- 'CO2' # MF-1020 #'MD&C'  # 'MAB_MDC' # 'Bottle_5Microbes' # Bottle_MAB12-0505-0525 #'Bottle'# 'Bottle_MAB12' # 'MAB'
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

# Remove rows with NA values
data <- na.omit(rawdata)
head(data)

library(dplyr)
# Calculate statistics by group
group_stats <- data %>%
  group_by(Group) %>%
  summarise(
    mean = mean(CCUS, na.rm = TRUE),
    sd = sd(CCUS, na.rm = TRUE),
    median = median(CCUS, na.rm = TRUE),
    IQR = IQR(CCUS, na.rm = TRUE)
  )
print(group_stats)

# Define pairwise comparisons
# comparisons <- list(c("CO2-Dark", "CO2-IR"), c('NH4-Dark', 'NH4-IR'))
comparisons <- list(c("Dark", "IR"))
result <- compare_means(CCUS ~ Group, data = data, method = "t.test")  
title_name = '' # 'CO2 Capture Rate' # "Substrate Conversion"
manual_fill_colar =  c('Dark'='#9499A6', 'IR'='#D75B66') # c("CO2-Dark" = "#818181", "CO2-IR" = "#E57272", "Yellow" = "#FAC026", "NH4-Dark" = "#EBF2EA", "NH4-IR" = "#C0334D", "Yellow-NC" = "#F3D4A0")

# Create the violin plot
ggplot(data, aes(x = Group, y = CCUS, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +  # Scatter points
  # stat_compare_means(method = "t.test", label = "p.signif", comparisons = comparisons) +  # Add significance values for specific comparisons
  scale_fill_manual(values = manual_fill_colar) +  # Custom fill colors from Set_Visualization
  labs(title = title_name,
       y = "CO2 utilization (mg C L-1)",
       x = '',
       fill = "Group")  + 
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent'),
        # panel.border = element_rect(fill = 'transparent', color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        # axis.ticks = element_line(linewidth = 1),
        legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_blank()
        )

