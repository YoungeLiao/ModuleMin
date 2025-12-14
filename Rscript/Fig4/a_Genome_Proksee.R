rm(list=ls())

# load library
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)

# --- load data ---
# 1. kegg_anno
MAG <- 'MAG57'
path.MAG <- paste('./data/Fig4/', MAG, '/anno_summary_', MAG, '.xlsx', sep = '')
rawdata.MAG <- data.frame(read_excel(path.MAG))
head(rawdata.MAG)
# Gene.ID              Location Genome Strand Start  End Length..bp.       Swiss_Prot_Hit
# 1 MAG90_gene0001 IR3_k97_191797_ORF001  MAG90      +   367 1020         654                    -
#   2 MAG90_gene0002 IR3_k97_191797_ORF002  MAG90      +  1008 1904         897                    -
#   3 MAG90_gene0003 IR3_k97_191797_ORF003  MAG90      +  2003 3598        1596 sp|P46133|ABGT_ECOLI
# 4 MAG90_gene0004 IR3_k97_191797_ORF004  MAG90      +  3661 4509         849                    -
#   5 MAG90_gene0005 IR3_k97_191797_ORF005  MAG90      -  5427 4486         942 sp|A4XNJ7|GLSA_PSEMY
# 6 MAG90_gene0006 IR3_k97_191797_ORF006  MAG90      -  8342 5508        2835 sp|Q87LA0|UVRA_VIBPA
# Swiss_Prot_Description  COG.ID COG.Type     KO KEGG_name
# 1                                                                                                                   - COG2871        C K00351      nqrF
# 2                                                                                                                   - COG2871        C K00351      nqrF
# 3              p-aminobenzoyl-glutamate transport protein OS=Escherichia coli (strain K12) OX=83333 GN=abgT PE=1 SV=3 COG2978        H K12942      abgT
# 4                                                                                                                   - COG0834      E;T      -         -
#   5                                       Glutaminase OS=Pseudomonas mendocina (strain ymp) OX=399739 GN=glsA PE=3 SV=1 COG2066        E K01425      glsA
# 6 UvrABC system protein A OS=Vibrio parahaemolyticus serotype O3:K6 (strain RIMD 2210633) OX=223926 GN=uvrA PE=3 SV=1 COG0178        L K03701      uvrA
# KO_description
# 1 Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1]
# 2 Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1]
# 3                               aminobenzoyl-glutamate transport protein
# 4                                                                      -
#   5                                               glutaminase [EC:3.5.1.2]
# 6                                             excinuclease ABC subunit A
# Level3
# 1                                                                                                                                                                                                                                                                         -
#   2                                                                                                                                                                                                                                                                         -
#   3                                                                                                                                                                                                                                                                         -
#   4                                                                                                                                                                                                                                                                         -
#   5 D-Amino acid metabolism;MicroRNAs in cancer;Two-component system;Glutamatergic synapse;GABAergic synapse;Proximal tubule bicarbonate reclamation;Alanine, aspartate and glutamate metabolism;Arginine biosynthesis;Metabolic pathways;Central carbon metabolism in cancer
# 6                                                                                                                                                                                                                                                Nucleotide excision repair

# --- extract Location ---
## details: extract location of each role, and split it to remove the last characters '_ORF...'. For example, transform 'IR_k97_191797_ORF001' into 'IR_k97_191797'

# ---- extract targeted colnumn and modify its mae to targeted column name ----- 
col.targeted <- c('KEGG_name', 'Location', 'Genome', 'Start', 'End', 'Strand', 'Gene.ID')
col.name.targeted <- c('name', 'contig', 'type', 'start', 'stop', 'strand', 'legend')

data.MAG <- rawdata.MAG[col.targeted]

data.MAG$Location <- sub("_ORF.*", "", data.MAG$Location)

names(data.MAG) <- col.name.targeted

head(data.MAG)

output.path <- paste('./data/Fig4/', MAG, '/anno_summary_', MAG, '_output.xlsx', sep = '')

library(writexl)
write_xlsx(data.MAG, 
           path = output.path)

output.path <- paste('./data/Fig4/', MAG, '/', MAG, '_processed.tsv', sep = '')
write.table(data.MAG, file = output.path, sep = '\t', quote = FALSE, row.names = FALSE)
