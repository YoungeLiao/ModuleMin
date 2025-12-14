rm(list=ls())
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

path <- './data/ex-fig1/ac_CO2.MultiSludgeGroup.xlsx'
sheet.name <- 'MultiSludge' # MultiSludge # CN
rawdata <- data.frame(read_excel(path, sheet = sheet.name))
head(rawdata, n=10)

summary_data <- rawdata %>%
  group_by(Group, Time) %>%
  summarise(
    Mean_HCO3 = mean(HCO3),
    SD_HCO3 = sd(HCO3),
    n = n(),
    SE_HCO3 = SD_HCO3 / sqrt(n)
  ) %>%
  ungroup()

theme_set(theme_classic() +
            theme(
              axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 14),
              legend.position = "right",
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              plot.margin = unit(c(1, 1, 1, 1), "lines"),
              plot.background = element_rect(fill = 'transparent', color = NA),
              panel.background = element_rect(fill = 'transparent')
            ))


p <- ggplot() +
  geom_point(data = rawdata, 
             aes(x = Time, y = HCO3, color = Group), 
             size = 2, alpha = 0.5) +

  geom_point(data = summary_data, 
             aes(x = Time, y = Mean_HCO3, color = Group), 
             size = 4) +

  geom_errorbar(data = summary_data,
                aes(x = Time, y = Mean_HCO3, ymin = Mean_HCO3 - SE_HCO3, 
                    ymax = Mean_HCO3 + SE_HCO3, color = Group),
                width = 0.2, size = 0.8) +

  geom_line(data = summary_data, 
            aes(x = Time, y = Mean_HCO3, color = Group), 
            size = 1) +

  scale_color_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  labs(x = "Time (d)", y = expression(HCO[3]^-~""~"(mg·"~L^"-1"~")")) +
  ylim(0, max(rawdata$HCO3) * 1.1) +
  guides(color = guide_legend(override.aes = list(size = 3)))

print(p)

# ===== for sheet 'CN' ======
p <- ggplot() +
  geom_point(data = rawdata, 
             aes(x = Time, y = HCO3, color = Group), 
             size = 2, alpha = 0.5) +
  geom_point(data = summary_data, 
             aes(x = Time, y = Mean_HCO3, color = Group), 
             size = 4) +
  geom_errorbar(data = summary_data,
                aes(x = Time, y = Mean_HCO3, ymin = Mean_HCO3 - SE_HCO3, 
                    ymax = Mean_HCO3 + SE_HCO3, color = Group),
                width = 0.2, size = 0.8) +
  geom_ribbon(data = summary_data,
              aes(x = Time, 
                  ymin = Mean_HCO3 - 1.96 * SE_HCO3, 
                  ymax = Mean_HCO3 + 1.96 * SE_HCO3,
                  fill = Group),
              alpha = 0.2) +
  geom_line(data = summary_data, 
            aes(x = Time, y = Mean_HCO3, color = Group), 
            size = 1) +
  scale_color_manual(values = c("#363062","#B31312", "#67729D", "#E7BCDE", "#BBE2EC",  "#F6EACB")) +
  scale_fill_manual(values = c("#363062","#B31312", "#67729D", "#E7BCDE", "#BBE2EC",  "#F6EACB")) +
  labs(x = "Time (d)", y = expression(HCO[3]^-~""~"(mg·"~L^"-1"~")")) +
  ylim(0, max(rawdata$HCO3) * 1.1) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p
