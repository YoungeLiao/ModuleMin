srm(list=ls())

# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

# load data
## 1. Edge
path <- './data/Fig3/c_network_edge.xlsx' 
sheet.name <- 'Edge_Clu6' # asv_taxon
rawdata.edge <- data.frame(read_excel(path, sheet = sheet.name))

head(rawdata.edge)
# source  target     weight Source_genus Source_species Target_genus Target_species
# 1 ASV1459 ASV2875  0.9999768         <NA>           <NA>         <NA>           <NA>
#   2 ASV1459   ASV44  0.6942713         <NA>           <NA>         <NA>           <NA>
#   3 ASV1459 ASV5219  0.9999417         <NA>           <NA>         <NA>           <NA>
#   4 ASV1459 ASV1619  0.6390222         <NA>           <NA>         <NA>           <NA>
#   5 ASV1459  ASV378 -0.7759323         <NA>           <NA>         <NA>           <NA>
#   6 ASV1459   ASV17 -0.7777944         <NA>           <NA>         <NA>           <NA>

## 2. asv taxon
path.taxon <- './data/Fig3/c_network_edge.xlsx' 
sheet.name.taxon <- 'asv_taxon' # asv_taxon
rawdata.taxon <- data.frame(read_excel(path.taxon, sheet = sheet.name.taxon))

head(rawdata.edge)
# ASV_ID Dark1 Dark2 Dark3   IR1  IR2  IR3 Initial      Domain               Kingdom            Phylum                  Class
# 1 ASV2457     0     0     0     0    0    0       0 d__Bacteria k__norank_d__Bacteria   p__Bacteroidota         c__Bacteroidia
# 2    ASV1   173   517 11721 14250  241 9914       0 d__Bacteria k__norank_d__Bacteria p__Pseudomonadota c__Gammaproteobacteria
# 3  ASV527     0     0     0     0    0    0       0 d__Bacteria k__norank_d__Bacteria   p__Bacteroidota         c__Bacteroidia
# 4   ASV10  2155  1773   447   601 1466  756    2137 d__Bacteria k__norank_d__Bacteria p__Pseudomonadota c__Gammaproteobacteria
# 5 ASV4984     0     0     0     0    0    0       0 d__Bacteria k__norank_d__Bacteria      p__Bacillota          c__Clostridia
# 6   ASV34     2  1061     2   259 1346 1132       0 d__Bacteria k__norank_d__Bacteria   p__Bacteroidota         c__Bacteroidia
# Order                            Family               Genus                                           Species
# 1                        o__Cytophagales               f__Flammeovirgaceae      g__Xanthovirga                       s__bacterium_g__Xanthovirga
# 2                    o__Enterobacterales                 f__Idiomarinaceae    g__Aliidiomarina                            s__bacterium_N159G.618
# 3                    o__Flavobacteriales              f__Flavobacteriaceae     g__Xanthomarina                   s__unclassified_g__Xanthomarina
# 4                o__Gammaproteobacteria_ f__norank_o__Gammaproteobacteria_ g__Wenzhouxiangella               s__unclassified_g__Wenzhouxiangella
# 5 o__Peptostreptococcales-Tissierellales                f__Caminicellaceae    g__Wukongibacter          s__uncultured_bacterium_g__Wukongibacter
# 6                    o__Flavobacteriales                 f__Cryomorphaceae         g__Vicingus s__uncultured_Bacteroidetes_bacterium_g__Vicingus

# add annotations for source nodes
rawdata.edge <- rawdata.edge %>%
  left_join(rawdata.taxon %>% 
              select(ASV_ID, Genus, Species) %>%
              rename(Source_genus = Genus,
                     Source_species = Species),
            by = c("source" = "ASV_ID"))

# add annotations for target nodes
rawdata.edge <- rawdata.edge %>%
  left_join(rawdata.taxon %>%
              select(ASV_ID, Genus, Species) %>%
              rename(Target_genus = Genus,
                     Target_species = Species),
            by = c("target" = "ASV_ID"))

head(rawdata.edge)
# source  target     weight Source_genus               Source_species
# 1 ASV1459 ASV2875  0.9999768 g__Halomonas s__unclassified_g__Halomonas
# 2 ASV1459   ASV44  0.6942713 g__Halomonas s__unclassified_g__Halomonas
# 3 ASV1459 ASV5219  0.9999417 g__Halomonas s__unclassified_g__Halomonas
# 4 ASV1459 ASV1619  0.6390222 g__Halomonas s__unclassified_g__Halomonas
# 5 ASV1459  ASV378 -0.7759323 g__Halomonas s__unclassified_g__Halomonas
# 6 ASV1459   ASV17 -0.7777944 g__Halomonas s__unclassified_g__Halomonas
# Target_genus                                    Target_species
# 1                            g__Mesonia                        s__unclassified_g__Mesonia
# 2                           g__Vicingus s__uncultured_Bacteroidetes_bacterium_g__Vicingus
# 3                          g__Halomonas                      s__unclassified_g__Halomonas
# 4                          g__Halomonas                      s__unclassified_g__Halomonas
# 5 g__unclassified_k__norank_d__Bacteria             s__unclassified_k__norank_d__Bacteria
# 6                       g__Peredibacter   s__uncultured_Bacteriovorax_sp._g__Peredibacter

# Define interested genus list 
target_genera <- c('g__Aliidiomarina', 'g__Vicingus', 'g__Halomonas', 'g__Kangiella')

# filtering
rawdata.edge.filtered <- rawdata.edge %>%
  filter(Source_genus %in% target_genera)

head(rawdata.edge.filtered)

output.path <- './data/Fig3/c_network_edge_filtered.output.xlsx'

library(writexl)
write_xlsx(rawdata.edge.filtered, 
           path = output.path)
