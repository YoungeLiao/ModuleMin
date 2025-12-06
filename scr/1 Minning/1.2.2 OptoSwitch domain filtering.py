"""
===
This script further filters opto switch candidates based on functional domains.
===
Input:
- bin = '/MAG73' # '' # 'MAG90'
- pfamscan_file = output_dir + "/pfamscan_results_new_3domain.txt". Ensure it exists, which is produced by scripts 1.2.1
Output:
- output_dir + filtered_ir_proteins_all.fasta: obtain the fasta file of 'pfamscan_results' in 1.2.1 
"""
# %%
from Bio import SeqIO


# %%
# Step 1: 过滤包含特定结构域的序列
def filter_by_domain(pfamscan_file, domains, output_file):
    filtered_sequences = []

    with open(pfamscan_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.strip().split()
                if len(cols) > 6 and any(domain in cols[6] for domain in domains):  # Final candidates required detection of at least one photosensory domain (PAS, GAF, or PHY) 
                    filtered_sequences.append(cols[0])

    # 提取过滤后的序列
    candidate_sequences = list(SeqIO.parse( output_dir + "/candidates.faa", "fasta"))
    filtered_sequences = [seq for seq in candidate_sequences if seq.id in filtered_sequences]

    # 保存过滤后的序列
    SeqIO.write(filtered_sequences, output_file, "fasta")

# %%
# 主程序***
if __name__ == "__main__":
    # 输入文件和参数
    bin = '/MAG73' # '' # 'MAG90'
    output_dir = "../../data/Fig2/Fig2bc.IRswitch.pyresults" + bin
    pfamscan_file = output_dir + "/pfamscan_results_new_3domain.txt"
    domains = ["PAS", "GAF", "PHY"]
    output_file = output_dir + "/filtered_ir_proteins_all.fasta"

    # Step 1: 过滤结构域
    filter_by_domain(pfamscan_file, domains, output_file)

#========== output_file explanation ==========
# ../results/filtered_ir_proteins_all.fasta: 未经过结构域过滤的候选序列，但经过Diamond alignment的候选序列
# ../results/filtered_ir_proteins.fasta: 经过结构域过滤后的候选序列，即包含 PAS、GAF 或 PHY 结构域的序列
# %%
