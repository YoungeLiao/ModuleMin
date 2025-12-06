rm(list=ls())
# source('./Rscript/functions.R')
# load library
library(readxl)
library(ggplot2)
library(ggpubr)

# load data
path <- './data/Fig4/d_Lipid_Energy_Protein_NC.xlsx' # CN_Conversion.xlsx # NH4CO3-MF.xlsx # Optogenetics.xlsx # './data/Hydraulic.xlsx' # './data/MultiSludge.xlsx' # './data/Optogenetics.xlsx' # './data/Sum-Batch.xlsx'
sheet.name <- 'Sheet1' # MF-1020 #'MD&C'  # 'MAB_MDC' # 'Bottle_5Microbes' # Bottle_MAB12-0505-0525 #'Bottle'# 'Bottle_MAB12' # 'MAB'
rawdata <- data.frame(read_excel(path, sheet = sheet.name))
data <- rawdata
# data <- subset(rawdata, rawdata$Batch == 'Opto5')


# Calculate statistics by group
library(dplyr)
group_stats <- data %>%
  group_by(Group) %>%
  summarise(
    mean = mean(Intensity, na.rm = TRUE),
    sd = sd(Intensity, na.rm = TRUE),
    median = median(Intensity, na.rm = TRUE),
    IQR = IQR(Intensity, na.rm = TRUE)
  )
print(group_stats)

# group_stats <- data %>%
#   group_by(Group) %>%
#   summarise(
#     mean = mean(Em598, na.rm = TRUE),
#     sd = sd(Em598, na.rm = TRUE),
#     median = median(Em598, na.rm = TRUE),
#     IQR = IQR(Em598, na.rm = TRUE)
#   )
# print(group_stats)

# # Remove non-finite values
# data <- data[is.finite(data$Intensity), ]

# # Define pairwise comparisons
comparisons <- list(c("DarkML", "IRML"))
result.Intensity <- compare_means(Intensity ~ Group, data = data, method = "t.test")  
result.ML <- compare_means(Intensity ~ Group, data = data, method = "t.test")  

# Create the bar plot with scatter points and error bars
ggplot(data, aes(x = Lipid, y = Intensity, fill = Light)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.7) +  # Bar plot showing mean values
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.7)) +  # Error bars showing standard error
  geom_jitter(size = 2, alpha = 0.4, stroke = 1,
              position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)) +  # Scatter points
  # stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("AnammoxDark", "AnammoxIR"))) +  # Add significance values for specific comparisons
  scale_fill_manual(values = c("Dark" = "#3f7ebd","IR" = "#E57272", 'Dark-2' = "#818181", "Yellow" = "#FAC026", "Dark-NC" = "#EBF2EA", "IR-NC" = "#C0334D", "Yellow-NC" = "#F3D4A0")) +  # Custom fill colors
  labs(title = "",
       x = "",
       y = "Fluorescence intensity (a.u.)",
       fill = "Group") + 
  scale_y_continuous(expand = c(0, 0.1, 0.1, 0)) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = c(0.2,0.92),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.background = element_blank()
  )
