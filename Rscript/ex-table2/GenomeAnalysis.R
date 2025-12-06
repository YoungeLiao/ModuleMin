rm(list=ls())

# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

# --- load data ---
# 1. kegg_anno
path.anno.kegg <- './data/ex-table2/kegg_anno.txt'
rawdata.anno.kegg <- read.delim(path.anno.kegg, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(rawdata.anno.kegg)

## 2. MAG taxon
path.MAG.taxon <- './data/rawdata.shared/all.MAG.summary.xlsx'
sheet.name.MAG.taxon <- 'all.MAG.summary' # asv_taxon
rawdata.MAG.taxon <- data.frame(read_excel(path.MAG.taxon, sheet = sheet.name.MAG.taxon))
head(rawdata.MAG.taxon)

# Filtering 
colname <- c('Genome.ID', 'Gene.ID', 'KEGG_gene','KEGG_name','KO','KO_description', 'Pathway', 'Enzyme', 'Module', 'Level3')
anno.kegg <- rawdata.anno.kegg[, colname]
head(anno.kegg)
# Genome.ID        Gene.ID        KEGG_gene KEGG_name     KO
# 1     MAG64 MAG64_gene0002   psn:Pedsa_2486      mdlB K18890
# 2     MAG64 MAG64_gene0003 sinh:LS482_03730      truA K06173
# 3     MAG64 MAG64_gene0005   dmm:dnm_022740         - K14445
# 4     MAG64 MAG64_gene1353    rum:CK1_26570         - K14445
# 5     MAG64 MAG64_gene0006  mod:AS202_15860         - K01969
# 6     MAG64 MAG64_gene0010  pex:IZT61_07865      dapB K00215
# KO_description
# 1                            ATP-binding cassette, subfamily B, multidrug efflux pump
# 2                                     tRNA pseudouridine38-40 synthase [EC:5.4.99.12]
# 3 solute carrier family 13 (sodium-dependent dicarboxylate transporter), member 2/3/5
# 4 solute carrier family 13 (sodium-dependent dicarboxylate transporter), member 2/3/5
# 5                          3-methylcrotonyl-CoA carboxylase beta subunit [EC:6.4.1.4]
# 6                            4-hydroxy-tetrahydrodipicolinate reductase [EC:1.17.1.8]
# Pathway    Enzyme                      Module
# 1                                         ko02010         -                           -
#   2                                               - 5.4.99.12                           -
#   3                                               -         -                           -
#   4                                               -         -                           -
#   5                                 ko01100;ko00280   6.4.1.4                      M00036
# 6 ko00300;ko00261;ko01100;ko01230;ko01120;ko01110  1.17.1.8 M00526;M00016;M00527;M00525
# Level3
# 1                                                                                                                                                              ABC transporters
# 2                                                                                                                                                                             -
#   3                                                                                                                                                                             -
#   4                                                                                                                                                                             -
#   5                                                                                                                 Metabolic pathways;Valine, leucine and isoleucine degradation
# 6 Lysine biosynthesis;Monobactam biosynthesis;Metabolic pathways;Biosynthesis of amino acids;Microbial metabolism in diverse environments;Biosynthesis of secondary metabolites

# ---- KO filtering ----
# 筛选特定的KO号
target_KOs <- c('K00174', 'K00175', 'K03737', 'K01647', 'K01648') # 'K00176'
# ACL <- c('K01676', 'K01677', 'K01679', 'K01681')
# target_KOs <- ACL
filtered.anno.kegg <- anno.kegg %>%
  filter(KO %in% target_KOs)

# 查看筛选结果
head(filtered.anno.kegg)


# ---- Genome filtering ----
## obtain targeted genomes 
# targeted_MAG <- c('MAG90', 'MAG57', 'MAG64', 'MAG42')
targeted_MAG <- c('MAG90', 'MAG57')
Filtered.anno.kegg.genome <- subset(anno.kegg, anno.kegg$Genome.ID == targeted_MAG)

# extract genome of functional taxa
MAG90 <- subset(Filtered.anno.kegg.genome, Filtered.anno.kegg.genome$Genome.ID == 'MAG90')
MAG57 <- subset(Filtered.anno.kegg.genome, Filtered.anno.kegg.genome$Genome.ID == 'MAG57')
head(MAG90)
# Genome.ID        Gene.ID       KEGG_gene KEGG_name     KO                                                         KO_description
# 67847     MAG90 MAG90_gene0001  lcp:LC55x_5582      nqrF K00351 Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1]
# 67849     MAG90 MAG90_gene2397  ili:K734_05255      nqrF K00351 Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1]
# 67851     MAG90 MAG90_gene1187 iab:K5X84_07770      abgT K12942                               aminobenzoyl-glutamate transport protein
# 67853     MAG90 MAG90_gene0005 lal:AT746_10200      glsA K01425                                               glutaminase [EC:3.5.1.2]
# 67855     MAG90 MAG90_gene0007 tcn:H9L16_10080      tcaB K07552             MFS transporter, DHA1 family, multidrug resistance protein
# 67857     MAG90 MAG90_gene0008 iab:K5X84_01640       ssb K03111                                      single-strand DNA-binding protein
# Pathway  Enzyme Module
# 67847                                                                               - 7.2.1.1      -
#   67849                                                                               - 7.2.1.1      -
#   67851                                                                               -       -      -
#   67853 ko00470;ko05206;ko02020;ko04724;ko04727;ko04964;ko00250;ko00220;ko01100;ko05230 3.5.1.2      -
#   67855                                                                               -       -      -
#   67857                                                         ko03430;ko03440;ko03030       -      -
#   Level3
# 67847                                                                                                                                                                                                                                                                         -
#   67849                                                                                                                                                                                                                                                                         -
#   67851                                                                                                                                                                                                                                                                         -
#   67853 D-Amino acid metabolism;MicroRNAs in cancer;Two-component system;Glutamatergic synapse;GABAergic synapse;Proximal tubule bicarbonate reclamation;Alanine, aspartate and glutamate metabolism;Arginine biosynthesis;Metabolic pathways;Central carbon metabolism in cancer
# 67855                                                                                                                                                                                                                                                                         -
#   67857                                                                                                                                                                                                                  Mismatch repair;Homologous recombination;DNA replication

MAG90.keggname <- as.data.frame(table(MAG90$KEGG_name))
MAG57.keggname <- as.data.frame(table(MAG57$KEGG_name))

# 按Freq降序排序并取top5
MAG90.keggname.top5 <- MAG90.keggname %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 5)

MAG57.keggname.top5 <- MAG57.keggname %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 5)

# 查看结果
print("MAG90 top 5 KEGG names:")
print(MAG90.keggname.top5)
# Var1 Freq
# 1    -  150
# 2 hipA    4
# 3 mucR    4
# 4 acrA    3
# 5 narL    3

print("MAG57 top 5 KEGG names:")
print(MAG57.keggname.top5)
# Var1 Freq
# 1      -  150
# 2   bapA   27
# 3   sadA   20
# 4 rsbU_P    7
# 5    pal    5


# 获取top5的基因名称（排除'-'）
MAG90.top5.genes <- MAG90.keggname.top5$Var1[MAG90.keggname.top5$Var1 != '-']
MAG57.top5.genes <- MAG57.keggname.top5$Var1[MAG57.keggname.top5$Var1 != '-']

# 根据top5基因名称过滤原基因组数据
MAG90.abundant <- MAG90[MAG90$KEGG_name %in% MAG90.top5.genes, ]
MAG57.abundant <- MAG57[MAG57$KEGG_name %in% MAG57.top5.genes, ]

# 查看结果
print("MAG90 abundant genes:")
print(head(MAG90.abundant))
print("MAG57 abundant genes:")
print(head(MAG57.abundant))



# 定义输出文件路径
output.path <- './data/ex-table2/TargetedGenome.output.xlsx'

# 将过滤后的数据写入Excel文件
library(writexl)
write_xlsx(Filtered.anno.kegg.genome, 
           path = output.path)

