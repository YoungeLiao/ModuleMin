#!/usr/bin/bash
# ===
# This script is used to annotate the filtered optogenetic sequences
# ===

# """
# Input:
#  - filtered_ir_proteins.fasta: the filtered optogenetic sequences
# Output:
# - annotated.fasta: the annotated optogenetic sequences
# - trimmed.fasta
# """
# ===== Seting parameters =====
SAVE_DATA_DIR="../../data/Fig2/database"
INPUT_DATA_DIR="../../data/rawdata.shared"
DATABASE_DIR="../../data/Fig2/database"
OUTPUT_DIR="../../data/Fig2/Fig2bc.IRswitch.pyresults"
# ===== Seting parameters =====

#################### 1. Function Annotation ####################
# #################### 1.1 KASS
# # 下载并安装 KEGG API 的 Python 客户端库：
# pip install bioservices
# # 编写 kaas.py 脚本
# # output kegg annotation 
# # Details can be found in KEGG annotation file
# mv $INPUT_DATA_DIR/anno_overview.xlsx $INPUT_DATA_DIR/anno_overview.txt

# #################### 2. Phynogenetics ####################
# #################### 2.1 Installation 
# brew install mafft # 用于多重序列比对

# conda install -c bioconda trimal #（用于修剪比对结果）

# # brew install raxml # AxML（用于构建进化树）
# conda install -c bioconda raxml # AxML（用于构建进化树）

#################### 2.2 Conudct alignments ####################
# 1. 多重序列比对
mafft $OUTPUT_DIR/filtered_ir_proteins.fasta > $OUTPUT_DIR/aligned.fasta
# 2. 修剪比对结果（使用 TrimAl）：
trimal -in $OUTPUT_DIR/aligned.fasta -out $OUTPUT_DIR/trimmed.fasta -automated1
# 3. 构建进化树（使用 RAxML）： # will take some time
raxmlHPC -s $OUTPUT_DIR/trimmed.fasta -n tree_output -m PROTGAMMAAUTO -p 12345

# # method 2 - clustalo FastTree
# """
# output: 
# - output.tree: tree file
# """

# # installation 
# ## conda install -c bioconda clustalo
# clustalo -i $OUTPUT_DIR/filtered_ir_proteins_all.fasta -o $OUTPUT_DIR/combined.aln # will take some time, but much faster than method 1
# FastTree $OUTPUT_DIR/combined.aln > $OUTPUT_DIR/output.tree
