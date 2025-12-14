rm(list=ls())
# 加载必要包
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# 1. load data
path.metabpathList.info <- "./results_data/Fig3/Fig3b2_metabolomics_umap_plotdata.xlsx"
metabpathList.info <- data.frame(read_excel(path.metabpathList.info))
metabpathList.top <- subset(metabpathList.info, top_pathways == TRUE)

# 2. load pathway
path.pathway <- "./results_data/Fig3/pathway_abundance.xlsx"
pathway <- data.frame(read_excel(path.pathway))

# 3. load metabolites annotation
path.MetabList.anno <- "./results_data/Fig3/metabolite_annotations.xlsx"
MetabList.anno <- data.frame(read_excel(path.MetabList.anno))

# 4. obtain top pathways
top_pathway_names <- metabpathList.top$pathways

# 5. match pathway_description to obtain map id
top_pathway_mapids <- pathway %>%
  filter(pathway_description %in% top_pathway_names) %>%
  select(pathway, pathway_description)

# 6. transform map id into vector
top_mapids <- unique(top_pathway_mapids$pathway)

# 7. MetabList.anno$pathways data processing
MetabList.anno_long <- MetabList.anno %>%
  mutate(pathways = as.character(pathways)) %>%
  filter(!is.na(pathways)) %>%
  separate_rows(pathways, sep = ";")

# 8. extract metabolites that subjected to top pathways
top_metabolites <- MetabList.anno_long %>%
  filter(pathways %in% top_mapids) %>%
  select(metabolite, kegg_id, pathways) %>%
  left_join(top_pathway_mapids, by = c("pathways" = "pathway")) %>%
  select(metabolite, kegg_id, pathway_mapid = pathways, pathway_description) %>%
  arrange(pathway_description, metabolite)

# 9. print results 
print(top_metabolites)
unique(top_metabolites$metabolite)

# 10. save data
write.csv(top_metabolites, file = "./results_data/Fig3/Fig3b3_top_pathways_metabolites.csv", row.names = FALSE)
