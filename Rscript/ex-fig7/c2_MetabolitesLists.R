rm(list=ls())
# 加载必要包
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# 1. 读取数据
path.metabpathList.info <- "./results_data/Fig3/Fig3b2_metabolomics_umap_plotdata.xlsx"
metabpathList.info <- data.frame(read_excel(path.metabpathList.info))
# 只保留top pathways
metabpathList.top <- subset(metabpathList.info, top_pathways == TRUE)

# 2. 读取pathway映射表
path.pathway <- "./results_data/Fig3/pathway_abundance.xlsx"
pathway <- data.frame(read_excel(path.pathway))

# 3. 读取代谢物注释表
path.MetabList.anno <- "./results_data/Fig3/metabolite_annotations.xlsx"
MetabList.anno <- data.frame(read_excel(path.MetabList.anno))

# 4. 获取top pathways的名称
top_pathway_names <- metabpathList.top$pathways

# 5. 根据pathway_description匹配获得对应的map id
top_pathway_mapids <- pathway %>%
  filter(pathway_description %in% top_pathway_names) %>%
  select(pathway, pathway_description)

# 6. 将map id展开为向量
top_mapids <- unique(top_pathway_mapids$pathway)

# 7. 处理MetabList.anno$pathways为长格式，便于匹配
# 有的代谢物对应多个map id，用;分隔
MetabList.anno_long <- MetabList.anno %>%
  mutate(pathways = as.character(pathways)) %>%
  filter(!is.na(pathways)) %>%
  separate_rows(pathways, sep = ";")

# 8. 提取属于top pathways的代谢物
top_metabolites <- MetabList.anno_long %>%
  filter(pathways %in% top_mapids) %>%
  select(metabolite, kegg_id, pathways) %>%
  left_join(top_pathway_mapids, by = c("pathways" = "pathway")) %>%
  select(metabolite, kegg_id, pathway_mapid = pathways, pathway_description) %>%
  arrange(pathway_description, metabolite)

# 9. 输出结果
print(top_metabolites)
unique(top_metabolites$metabolite)
# 10. 可选：保存为csv
write.csv(top_metabolites, file = "./results_data/Fig3/Fig3b3_top_pathways_metabolites.csv", row.names = FALSE)
