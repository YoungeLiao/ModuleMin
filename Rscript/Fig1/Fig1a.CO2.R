rm(list=ls())
# load library
library(readxl)
library(ggplot2)
library(ggpubr)

# load data
path <- './data/Fig1/Fig1d.CO2.xlsx' # Fig1a.CO2 # Fig1a.co2.tic
sheet.name <- 'CO2' # 'CO2' # Biomass
rawdata <- data.frame(read_excel(path, sheet = sheet.name))

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

