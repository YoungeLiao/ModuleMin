rm(list=ls())
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)


# load data
path.anno <- './data/rawdata.shared/anno_overview.txt'
anno <- read.delim(path.anno, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(anno)

path.abundance <- './data/rawdata.shared/reads_number_relative.txt'
abundance <- read.delim(path.abundance, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(abundance)
# #Query	Length	Domain	Kingdom	Phylum	Class	Order	Family	Genus	SpeciesNOG	COG_description	COG_Function	COG_Fun_description	COG_Category	KEGG_Gene	KO	KEGG_Name	KO_description	KEGG_Pathway	KEGG_Path_description	KEGG_Enzyme	KEGG_En_description	KEGG_Modules	KEGG_Mo_description	KEGG_Level1	KEGG_Level2	CAZY_Family	CAZY_Class	CAZY_Cl_description	CAZY_Fa_description	ARDB_ARG	ARDB_Type	ARDB_Antibiotic_type	ARDB_Class	ARDB_Cl_description	ARDB_Antibiotic_class	CARD_ARO	CARD_ARO_name	CARD_ARO_description	CARD_AMR_Gene_Family	CARD_Drug_Class	CARD_Resistance_Mechanism	VFDB_ID	VFDB_related_genes	VFDB_product	VFDB_VFs	VFDB_VF_Function	VFDB_Species	VF_Category
# Yellow2_k97_214683_1	399	d__Bacteria	k__unclassified_d__Bacteria	p__Actinomycetota	c__Actinomycetes	o__Geodermatophilales	f__Geodermatophilaceae	g__Blastococcus	s__Blastococcus_sp.	COG0596	2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase MenH and related esterases, alpha/beta hydrolase fold	H;R	Coenzyme transport and metabolism;General function prediction only	METABOLISM;POORLY CHARACTERIZED	-	-	-	-

# Extract Query、Phylum、Genus as new dataframe
anno_sub <- anno %>%
  select(X.Query, Phylum, Genus)
colnames(anno_sub) <- c('GeneID', 'Phylum', 'Genus')
# Check
head(anno_sub)

abun_sub <- abundance %>%
  select(GeneID, IR1, IR2, IR3, Dark1, Dark2, Dark3, Initial)
head(abun_sub)
df <- merge(anno_sub, abun_sub)
# merged_df <- dplyr::left_join(abun_sub, anno_sub, by = "GeneID")
head(df)
# GeneID                        Phylum                                         Genus
# 1      C1_Dark_k97_1_1                             -                                             -
#   2  C1_Dark_k97_10000_1             p__Pseudomonadota                              g__Defluviimonas
# 3 C1_Dark_k97_100000_2    p__Candidatus_Rokubacteria    g__unclassified_p__Candidatus_Rokubacteria
# 4 C1_Dark_k97_100002_1 p__Candidatus_Hydrogenedentes g__unclassified_p__Candidatus_Hydrogenedentes
# 5 C1_Dark_k97_100006_2             p__Pseudomonadota                                g__Marinicella
# 6 C1_Dark_k97_100007_2             p__Pseudomonadota                                 g__Idiomarina
# IR1          IR2          IR3        Dark1        Dark2        Dark3      Initial
# 1 4.271292e-08 4.060595e-08 3.809673e-08 8.552859e-08 0.000000e+00 0.000000e+00 0.000000e+00
# 2 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
# 3 1.708517e-07 1.218178e-07 7.619346e-08 3.848786e-07 1.976041e-07 7.576613e-08 1.258207e-07
# 4 0.000000e+00 1.218178e-07 0.000000e+00 0.000000e+00 3.952083e-08 0.000000e+00 0.000000e+00
# 5 3.417034e-07 2.558175e-06 2.247707e-06 2.266508e-06 1.383229e-06 1.136492e-06 8.388045e-08
# 6 1.067823e-06 1.218178e-07 1.104805e-06 0.000000e+00 4.742499e-07 3.523125e-06 0.000000e+00

# 确保丰度列为数值
ab_cols <- c("IR1","IR2","IR3","Dark1","Dark2","Dark3","Initial")
df[ab_cols] <- lapply(df[ab_cols], function(x) as.numeric(as.character(x)))

# 将缺失或 "-" 的属名统一为 "Unclassified"（可选）
df$Genus[is.na(df$Genus) | df$Genus == ""] <- "Unclassified"
df$Genus[df$Genus == "-"] <- "Unclassified"

# 按属求和得到属水平的相对丰度
library(dplyr)
genus_abundance <- df %>%
  group_by(Genus) %>%
  summarise(across(all_of(ab_cols), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup()

# 保存结果
out_dir <- "./data/Fig2/processed_Metagenomic_abundance"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
readr::write_tsv(genus_abundance, file.path(out_dir, "genus_abundance.tsv"))

# 快览
print(head(genus_abundance, 5))
