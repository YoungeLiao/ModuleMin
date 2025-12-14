rm(list=ls())
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)


# load data
path <- './data/Fig3/d_TaxaGrowth.xlsx'
sheet.name <- 'plot_abundance' 
growth <- data.frame(read_excel(path, sheet = sheet.name))

head(growth)
# Genus  Enrichment Group
# 1 g__Aliidiomarina 0.324732989    IR
# 2 g__Aliidiomarina 0.005310099    IR
# 3 g__Aliidiomarina 0.239738758    IR
# 4 g__Aliidiomarina 0.002881083  Dark
# 5 g__Aliidiomarina 0.008547602  Dark
# 6 g__Aliidiomarina 0.266462529  Dark

# ensure correct types
growth$Enrichment <- as.numeric(as.character(growth$Enrichment))
# growth$Enrichment <- as.numeric(as.character(growth$Enrichment))
growth$Group <- factor(growth$Group, levels = c("IR", "Dark"))

# CNS-like colors
main_colors <- c("IR" = '#E59693', "Dark" = '#989CC8')

library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)


# define genus_order (order genera by overall mean Enrichment, descending; fallback to existing levels)
genus_order <- growth %>%
  group_by(Genus) %>%
  summarise(meanEnrichment = mean(Enrichment, na.rm = TRUE), n = sum(!is.na(Enrichment)), .groups = "drop") %>%
  arrange(desc(-meanEnrichment)) %>%
  pull(Genus) %>%
  as.character()
if (length(genus_order) == 0) genus_order <- unique(as.character(growth$Genus))

# ---------- Growth: ORIGINAL SCALE with error bars (show points; center error bars, only upper part) ----------
summary_growth <- growth %>%
  group_by(Genus, Group) %>%
  summarise(
    meanGrowth = mean(Enrichment, na.rm = TRUE),
    sdGrowth = sd(Enrichment, na.rm = TRUE),
    n = sum(!is.na(Enrichment)),
    seGrowth = ifelse(n > 0, sdGrowth / sqrt(n), NA_real_),
    .groups = "drop"
  )

summary_growth$Genus <- factor(summary_growth$Genus, levels = genus_order)
growth$Genus <- factor(growth$Genus, levels = genus_order)

dodge_width <- 0.75

p_growth <- ggplot() +
  # bars (dodge by Group)
  geom_col(data = summary_growth,
           aes(x = Genus, y = meanGrowth, fill = Group, group = Group),
           position = position_dodge(width = dodge_width),
           width = 0.75,
           color = "grey30") +
  # error bars showing only the upper half (from mean to mean + SE)
  geom_errorbar(data = summary_growth,
                aes(x = Genus, ymin = meanGrowth, ymax = meanGrowth + seGrowth, group = Group),
                position = position_dodge(width = dodge_width),
                width = 0.25,
                color = "black", size = 0.6) +
  # raw data points (jitter within each dodged position)
  geom_point(data = growth,
             aes(x = Genus, y = Enrichment, fill = Group, group = Group),
             shape = 21, color = "grey20", size = 1.6, alpha = 0.9,
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width)) +
  coord_flip() +
  scale_fill_manual(values = main_colors, name = "Group") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = "", y = "Abundance Enrichment", title = "") +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_blank()
  )

print(p_growth)


# ---------- Growth: z-score (within-genus) normalization with error bars ----------
# standardize Enrichment per Genus (z = (x - mean_genus) / sd_genus)
growth <- growth %>%
  group_by(Genus) %>%
  mutate(
    .mean_genus = mean(Enrichment, na.rm = TRUE),
    .sd_genus = sd(Enrichment, na.rm = TRUE),
    Enrichment_z = ifelse(.sd_genus > 0, (Enrichment - .mean_genus) / .sd_genus, NA_real_)
  ) %>%
  ungroup() %>%
  select(-.mean_genus, -.sd_genus)

summary_z <- growth %>%
  group_by(Genus, Group) %>%
  summarise(
    meanEnrichment_z = mean(Enrichment_z, na.rm = TRUE),
    sdEnrichment_z = sd(Enrichment_z, na.rm = TRUE),
    n = sum(!is.na(Enrichment_z)),
    seEnrichment_z = ifelse(n > 0, sdEnrichment_z / sqrt(n), NA_real_),
    .groups = "drop"
  )

summary_z$Genus <- factor(summary_z$Genus, levels = genus_order)
growth$Genus <- factor(growth$Genus, levels = genus_order)

p_z <- ggplot() +
  geom_col(data = summary_z,
           aes(x = Genus, y = meanEnrichment_z, fill = Group),
           position = position_dodge2(width = 0.9, preserve = "single"),
           width = 0.75,
           color = "grey30") +
  geom_errorbar(data = summary_z,
                aes(x = Genus, ymin = meanEnrichment_z - seEnrichment_z, ymax = meanEnrichment_z + seEnrichment_z, group = Group),
                position = position_dodge2(width = 0.9),
                width = 0.25,
                color = "black", size = 0.6) +
  geom_point(data = growth,
             aes(x = Genus, y = Enrichment_z, fill = Group),
             shape = 21, color = "grey20", size = 1.6, alpha = 0.9,
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)) +
  scale_fill_manual(values = main_colors, name = "Group") +
  coord_flip() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = "", y = "Enrichment (z-score within genus)", title = "Genus enrichment by Group (z-score) — means ± SE") +
  theme_pubr(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_blank()
  )

print(p_z)
