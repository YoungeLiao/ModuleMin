rm(list=ls())

# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

# load data
path <- './data/Fig3/a_BiomassSynthesis.xlsx' # ex.fig2e.BiomassSynthesisRate.xlsx
sheet.name <- 'BS' # BS
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

data <- rawdata[is.finite(rawdata$BiomassSynthesisRate), ] # data <- rawdata[is.finite(rawdata$BiomassSynthesisRate), ]
head(data)

# Calculate statistics by group
## BiomassSynthesisRate
library(dplyr)
# data$Group <- as.factor(data$Group)
group_stats <- data %>%
  group_by(Group) %>%
  summarise(
    mean = mean(BiomassSynthesisRate, na.rm = TRUE),
    sd = sd(BiomassSynthesisRate, na.rm = TRUE),
    median = median(BiomassSynthesisRate, na.rm = TRUE),
    IQR = IQR(BiomassSynthesisRate, na.rm = TRUE)
  )
print(group_stats)


# Define pairwise comparisons
comparisons <- list(c("Dark", "IR"))
result <- compare_means(BiomassSynthesisRate ~ Group, data = data, method = "t.test")  

Cust_palette <- c('#989CC8', '#E59693', '#B387FE', '#FED565') 
# Create the boxplot with scatter points and significance comparison
ggplot(rawdata, aes(x = Group, y = BiomassSynthesisRate, fill = Group)) +
  geom_boxplot(outlier.shape = NA,
               position = position_dodge(width = 0.4)) +  # Boxplot without outliers
  # geom_jitter(width = 0.2, size = 2, alBiomassSynthesisRatea = 0.5, stroke = 1) +  # Scatter points
  # stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c("Dark", "IR"))) +  # Add significance values for specific comparisons
  geom_jitter(data=rawdata,
              mapping=aes(x=Group, y=BiomassSynthesisRate),  # , color = Group
              size = 3, alpha = 0.5,
              BiomassSynthesisRate = 0.5,
              position=position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
  scale_fill_manual(values = Cust_palette) +  # Custom fill colors from Set_Visualization
  labs(title = "",
       x = "",
       y = "Biomass synthesis rate (mg DCW L-1)", # "Biomass Synthesis (mg L-1)"
       fill = "Group") + 
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_blank()
  )

