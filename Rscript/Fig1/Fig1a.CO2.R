rm(list=ls())
# source('./Rscript/functions.R')
# load library
library(readxl)
library(ggplot2)
library(ggpubr)

# load data
path <- './data/Fig1/Fig1d.CO2.xlsx' # Fig1a.CO2 # Fig1a.co2.tic
sheet.name <- 'CO2' # 'CO2' # Biomass
rawdata <- data.frame(read_excel(path, sheet = sheet.name))


# Set_Visualization <- function(){
#   pale7 <- c('#4489FE','#818181','#00BEA4','#E57272','#ACC5FE','#B387FE', '#FAC026' )
#   pale8 <- c('#B387FE','#818181','#FAC026','#FECDD2','#4489FE','#00BEA4', '#E57272','#B387FE','#90A3AD' )
#   pale9 <- c('#818181','#4489FE','#00BEA4', '#FAC026','#E57272','#FECDD2','#B387FE','#90A3AD', '#DBB4DA')
#   pale8_IR <- c('#818181','#FECDD2','#FAC026','#B387FE','#E57272','#00BEA4','#4489FE','#ACC5FE' )
#   pale7_NaUV <- c('#818181','#4489FE','#00BEA4', '#FAC026','#E57272','#FECDD2','#90A3AD' )
#   
#   top.mar=0.1
#   right.mar=0.1
#   bottom.mar=0.1
#   left.mar=0.1
#   
#   mytheme <- theme_classic() +
#     theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           plot.title = element_text(hjust = 0.5,size = 18),
#           plot.background = element_rect(fill = 'transparent', color = NA),
#           panel.background = element_rect(fill = 'transparent'),
#           panel.grid = element_blank(),
#           axis.title = element_text(size = 16),
#           axis.text = element_text(size = 16),
#           # axis.ticks = element_line(linewidth = 1),
#           legend.position = c(0.21, 0.12),
#           legend.title = element_blank(),
#           legend.text = element_text(size = 16),
#           legend.background = element_blank(),
#           legend.key.width = unit(2, "cm"),  # Make legend lines longer
#           legend.key = element_blank()  # Remove legend background
#     )
#   
#   mytheme1<-theme_test()+
#     theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
#           axis.text = element_text(size = 14),
#           axis.text.x = element_text(size = 14),
#           axis.title = element_text(size = 16, face = 'bold'),
#           legend.text = element_text(size = 14),
#           # legend.title = element_text(size = 16),
#           legend.title = element_blank(),
#           # legend.position = c(1, 0.3), # 调整legend位置, FAB - NO3
#           #legend.position = c(0.14, 0.32), # 调整legend位置, NO2
#           # legend.position = c(0.8, 0.75), # 调整legend位置, NH4
#           legend.position = 'right',
#           # legend.position = "none", # Supporting
#           legend.background = element_blank(),
#           plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
#                            units="inches"))
#   
#   mytheme2<-theme_classic()+
#     theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
#           axis.text = element_text(size = 14),
#           axis.text.x = element_text(size = 14),
#           axis.title = element_text(size = 16, face = 'bold'),
#           legend.text = element_text(size = 14),
#           # legend.title = element_text(size = 16),
#           legend.title = element_blank(),
#           # legend.position = c(0.15, 0.2), # 调整legend位置, FAB - NO3
#           #legend.position = c(0.14, 0.32), # 调整legend位置, NO2
#           # legend.position = c(0.8, 0.75), # 调整legend位置, NH4
#           legend.position = 'right',
#           # legend.position = "none", # Supporting
#           legend.background = element_blank(),
#           plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
#                            units="inches"))
#   
#   visualization_collection <- list(pale7, pale8, pale7_NaUV, pale8_IR, mytheme, mytheme1, mytheme2)
#   return(visualization_collection)
# }
# visua_collect <- Set_Visualization()
# # visualization_collection <- list(pale7, pale9, pale7_NaUV, mytheme1, mytheme2)
# pale7 <- visua_collect[[1]]
# pale8 <- visua_collect[[2]]
# pale7_NaUV <- visua_collect[[3]]
# pale8_IR <- visua_collect[[4]]
# mytheme <- visua_collect[[5]]
# mytheme1 <- visua_collect[[6]]
# mytheme2 <- visua_collect[[7]]

# filtering datasets
# Exp <- 'Opto6'
# data <- subset(rawdata, rawdata$Batch == Exp)
# Remove rows with NA values
# data <- is.na(rawdata)
data <- rawdata[is.finite(rawdata$TIC), ]

# Create the scatter plot with linear trend and confidence interval
ggplot(data, aes(x = Time, y = TIC, color = Group, fill = Group, shape = Group)) +  # Add 'shape' aesthetic
  geom_point(size = 2.5) +  # Adjust point size
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, linetype = "dashed") +  # Polynomial regression (second-degree)
  scale_shape_manual(values = c("Dark" = 16, "IR" = 17, "Yellow" = 18, "Dark-NC" = 15, "IR-NC" = 3, "Yellow-NC" = 4)) +  # Custom shapes for groups
  scale_fill_manual(values = c("Dark" = "#818181", "IR" = "#E57272", "Yellow" = "#FAC026", "Dark-NC" = "#EBF2EA", "IR-NC" = "#C0334D", "Yellow-NC" = "#F3D4A0")) +  # Custom fill colors
  scale_color_manual(values = c("Dark" = "#818181", "IR" = "#E57272", "Yellow" = "#FAC026", "Dark-NC" = "#EBF2EA", "IR-NC" = "#C0334D", "Yellow-NC" = "#F3D4A0")) +  # Custom color for trend lines
  labs(title = "",
       x = "Time (h)",
       y = "mg C L-1",
       color = "Group",
       fill = "Group",
       shape = "Group") +  # Add shape legend
  theme_classic() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        # axis.ticks = element_line(linewidth = 1),
        legend.position = c(0.21, 0.12),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_blank(),
        legend.key.width = unit(2, "cm"),  # Make legend lines longer
        legend.key = element_blank()  # Remove legend background
  )



# plotx.name <- 'Time'
# ploty.name <- 'TIC_All' 
# plotcolor.name <- 'Group'
# axis.y.name <- 'TIC (mg C·L-1)' # 'TIN (mg N L-1)' # 'Ammonia (mg N L-1)' 'Nitrate (mg N L-1)'   #'Nitrite (mg N L-1)'  # 'Sulfate (mg S L-1)'
# axis.x.name <- 'Time (day)'
# fig.name <- 'Batch Opto4'
# x.limit <- c(0, 25)
# y.limit <- c(1, 2000)
# BREAK <- 5
# 
# # filtering datasets
# Exp <- 'Opto4'
# plotdata <- subset(rawdata, rawdata$Batch == Exp)
# 
# l <- ggline(plotdata, x = plotx.name, y = ploty.name, color = plotcolor.name,
#   palette = pale8_IR, 
#   linetype = plotcolor.name,  # 根据组别设置线型
#   width = 0.2,  # 线宽
#   numeric.x.axis = TRUE,  # 将 x 轴设置为数值型
#   point.size = 3,  # 调整点大小
#   add = "mean_se", add.params = list(size = 3, width = 0.16)) + 
#   scale_x_continuous(breaks = seq(from = x.limit[1], to = x.limit[2] + 1, by = BREAK)) + 
#   coord_cartesian(xlim = x.limit) + 
#   labs(y = axis.y.name, x = axis.x.name, title = fig.name) + 
#   # scale_shape_manual(values = c("Dark" = 16, "IR" = 17, "Yellow" = 18, "Dark-NC" = 15, "IR-NC" = 3, "Yellow-NC" = 4)) +  # 自定义点形状
#   scale_linetype_manual(values = c("Dark" = "solid", "IR" = "dashed", "Yellow" = "dotted")) +  # 自定义线型
#   mytheme + 
#   # guides(
#   #   linetype = guide_legend(override.aes = list(size = 1.2)),  # 调整图例中线型的显示
#   #   color = guide_legend(override.aes = list(size = 3))  # 调整图例中颜色的显示
#   # ) 
#   guides(colour = guide_legend(override.aes = list(linetype = c('solid', 'dashed'))))
# 
# l
