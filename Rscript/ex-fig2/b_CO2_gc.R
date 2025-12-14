rm(list=ls())
# source('./Rscript/functions.R')
# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot) 

# load data
path <- './data/ex.fig1b.co2.gc.xlsx'
sheet.name <- 'gl_co2_stack' 
rawdata.2 <- data.frame(read_excel(path, sheet = sheet.name))

path <- './data/ex.fig1b.co2.gc.xlsx'
sheet.name <- 'gl_co2' 
rawdata.1 <- data.frame(read_excel(path, sheet = sheet.name))

head(rawdata)
# output:
# Samples Time    Group Batch GC_Area    CO2_ppm       C_m        TIC     TIC_m
# 1      D1    0 Dark-Gas Opto5   808.2 0.09979302 0.2464888  747.66656 22429.997
# 2      D2    0 Dark-Gas Opto6   498.4 0.05915872 0.1461220 1124.53604 33736.081
# 3      D3    0 Dark-Gas Opto6   366.8 0.04189767 0.1034872  954.02520 28620.756
# 4      I1    0   IR-Gas Opto5   514.7 0.06129668 0.1514028  -95.84776 -2875.433
# 5      I2    0   IR-Gas Opto6   484.8 0.05737490 0.1417160  779.77372 23393.212
# 6      I3    0   IR-Gas Opto6   486.8 0.05763723 0.1423640  875.78348 26273.504


plotdata <- rawdata.1 %>%
    select(Group, C_m, TIC_m) %>% 
    pivot_longer(cols = c(C_m, TIC_m), names_to = "Phase", values_to = "Carbon_Mass") %>%  # 转换为长格式
    group_by(Group, Phase) %>% 
    summarise(Carbon_Mass = mean(Carbon_Mass, na.rm = TRUE))  # 计算均值


y.limit <- c(0, 36000)
y.BREAK <- 10000
main_plot <- ggplot(plotdata, aes(x = Group, y = Carbon_Mass, fill = Phase)) +
    stat_summary(fun = sum, geom = "bar", position = "stack", width = 0.7, aes(color = Phase), size = 0.5, fill = "transparent") +  # 边框型柱状图，边框颜色与分组颜色一致
    geom_jitter(data = rawdata.2, aes(x = Group, y = C_m, fill = Phase), 
                            position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), 
                            size = 3, alpha = 0.5, stroke = 1) +  # 添加数据点
    scale_fill_manual(values = c("C_m" = "#E57272", "TIC_m" = "#4489FE")) +  # 自定义填充颜色
    scale_color_manual(values = c("C_m" = "#E57272", "TIC_m" = "#4489FE"), labels = c("C_m" = "Gas", "TIC_m" = "Liquid")) +  # 数据点颜色与边框颜色一致
    scale_y_continuous(breaks = seq(from = y.limit[1], to = y.limit[2], by = y.BREAK)) + 
    labs(title = "",
             x = "",
             y = "CO2 equivalent (mg)",
             fill = "Phase") +
    theme_classic() +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                plot.background = element_rect(fill = 'transparent', color = NA),
                panel.background = element_rect(fill = 'transparent'),
                panel.grid = element_blank(),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 16),
                legend.position = c(0.75, 0.9),
                legend.title = element_blank(),
                legend.text = element_text(size = 16),
                legend.background = element_blank(),
                legend.key.width = unit(1, "cm"),  # Make legend lines longer
                legend.key = element_blank()  # Remove legend background
    )
main_plot


zoom_plot <- ggplot(rawdata.2 %>% filter(Phase == "Gas"), aes(x = Group, y = C_ppm*10000)) + # fill = Group
    stat_summary(fun = mean, geom = "bar", position = "dodge", width = 0.7, color = '#E57272', fill = "transparent") +  # Bar plot showing mean values
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.7)) +  # Error bars showing standard error
    geom_jitter(width = 0.2, size = 3, alpha = 0.5, stroke = 1, color = '#E57272') +  # Scatter points
    scale_fill_manual(values =  c('Gas' = '#E57272')) +  # Custom fill colors
    labs(title = "",
             x = "",
             y = "CO2 (ppm)") +  # , fill = "Group"
    theme_classic() +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                plot.background = element_rect(fill = 'transparent', color = NA),
                panel.background = element_rect(fill = 'transparent'),
                panel.grid = element_blank(),
                axis.title.y.right = element_text(size = 16),  # Adjust right y-axis title size
                axis.text.y.right = element_text(size = 16),  # Adjust right y-axis text size
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 16),
                legend.position = 'none',
                legend.title = element_blank(),
                legend.text = element_text(size = 16),
                legend.background = element_blank()
    ) +
    scale_y_continuous(sec.axis = dup_axis(name = "CO2 (ppm)"))  # Add secondary y-axis on the right
zoom_plot
