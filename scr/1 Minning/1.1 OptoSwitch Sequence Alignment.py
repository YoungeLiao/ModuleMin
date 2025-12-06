
"""
===
This scripts will build opto switch database based on public database keyworks searching
and then align sequences to the customed opto switch database to obtain candidates sequences
=== 

Input: 
- save_data_dir: the directory to save the database
- keywords_url: the keywords to search for the database
-* input_file: sequences to be aligned to the IR database, one of the following:
  1. gene.uniGeneset_MABAll.faa: all the non-redundant protein sequences of metagenomic data, which could align all the sequences from samples
  2. MAG9_CDS.faa: Only align certain BIN's sequence, i.e. the protein sequences to be aligned of MAG9

Ouput: 
- custom_IR_db.faa: fasta sequence of custom IR optogenetic switch
- custom_IR_db_clustered.faa: clustered fasta sequence of custom IR optogenetic switch
-** custom_IR_db.dmnd: custom IR optogenetic switch database
- diamond_results.tsv: the details of all alignments. Note: this file hasn't filtered  
-** high_confidence_hits.tsv: the details of filtered alignments. 
-** candidates.faa: candidate sequences
-** candidate_ids.txt: the id of candidate sequences

- OUTPUT_DIR="../results/MAG9"

"""
# %%
import os
import subprocess
import pandas as pd

# 配置文件
config = {
    'input_file':'MAG73_CDS.faa', # 'gene.uniGeneset_MABAll.faa' # **  MAG90_CDS.faa
    'input_data_dir': "../../data/rawdata.shared/MAG/MAG73", # "../../data/rawdata.shared' # ../rawdata/MAG/MAG90
    
    'output_dir': "../../data/Fig2/Fig2bc.IRswitch.pyresults/MAG73", # ../../data/Fig2/Fig2bc.IRswitch.pyresults/ # ../results/MAG90
    'output_file':'diamond_results.tsv', # diamond_results_MAG90.tsv

    'keywords_url': 'https://rest.uniprot.org/uniprotkb/stream?query=bacteriophytochrome&format=fasta', # **
    'save_data_dir': "../../data/Fig2/Fig2bc.database", 
    'ref_database_faa': 'bph_refs.faa', 
    'build_db': False  # 是否构建参考数据库
}

def build_custom_reference_db(config):
    save_data_dir = config['save_data_dir']
    ref_faa_path = os.path.join(save_data_dir, config['ref_database_faa'])
    custom_faa_path = os.path.join(save_data_dir, "custom_IR_db.faa")
    clustered_faa_path = os.path.join(save_data_dir, "custom_IR_db_clustered.faa")
    diamond_db_path = os.path.join(save_data_dir, "custom_IR_db")
    keywords_url = config['keywords_url']

    # 1.1 Obtain Targeted Sequence
    print("Downloading bacteriophytochrome sequences from UniProt...")
    # Remove trailing ']' in the URL if present (fix bug)
    if keywords_url.endswith(']'):
        keywords_url = keywords_url[:-1]
    # Ensure save_data_dir exists
    if not os.path.exists(save_data_dir):
        os.makedirs(save_data_dir)
    curl_command = f"curl -f -s -S -o \"{ref_faa_path}\" \"{keywords_url}\""
    result = subprocess.run(curl_command, shell=True)
    if result.returncode != 0 or not os.path.isfile(ref_faa_path) or os.path.getsize(ref_faa_path) == 0:
        raise RuntimeError(f"Failed to download UniProt FASTA to {ref_faa_path}. Please check your internet connection and the URL.")

    # [Optional]合并自定义序列（如已有本地文件）
    print("Merging custom sequences...")
    if not os.path.isfile(ref_faa_path) or os.path.getsize(ref_faa_path) == 0:
        raise FileNotFoundError(f"Reference FASTA file not found: {ref_faa_path}")
    # For now, just copy the downloaded file (no user sequences)
    with open(ref_faa_path, "r") as fin, open(custom_faa_path, "w") as fout:
        fout.write(fin.read())

    # 1.2 去冗余与格式化
    print("Removing redundancy using CD-HIT...")
    cd_hit_command = (
        f"cd-hit -i \"{custom_faa_path}\" "
        f"-o \"{clustered_faa_path}\" "
        f"-c 0.95 -n 5 -M 16000"
    )
    result = subprocess.run(cd_hit_command, shell=True)
    if result.returncode != 0 or not os.path.isfile(clustered_faa_path) or os.path.getsize(clustered_faa_path) == 0:
        raise RuntimeError(f"CD-HIT failed or produced no output: {clustered_faa_path}")

    print("Converting to DIAMOND database format...")
    diamond_makedb_command = (
        f"diamond makedb --in \"{clustered_faa_path}\" -d \"{diamond_db_path}\""
    )
    result = subprocess.run(diamond_makedb_command, shell=True)
    if result.returncode != 0 or not os.path.isfile(diamond_db_path + ".dmnd"):
        raise RuntimeError(f"DIAMOND makedb failed or did not produce {diamond_db_path}.dmnd")

    print("Custom reference database built successfully.")

def diamond_alignment(config):
    input_data_dir = config['input_data_dir']
    input_file = config['input_file']
    output_dir = config['output_dir']
    output_file = config['output_file']
    save_data_dir = config['save_data_dir']

    # 检查output_dir是否存在，若不存在则创建一个新的文件夹
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
        
    # 2. DIAMOND比对参数优化
    print("Running DIAMOND alignment...")
    diamond_blastp_command = (
        f"diamond blastp "
        f"--query {input_data_dir}/{input_file} "
        f"--db {save_data_dir}/custom_IR_db.dmnd "
        f"--out {output_dir}/{output_file} "
        f"--outfmt 6 "
        f"--evalue 1e-5 "
        f"--id 30 "
        f"--query-cover 70 "
        f"--max-target-seqs 100 "
        f"--threads 32"
    )
    subprocess.run(diamond_blastp_command, shell=True, check=True)
    
    # 3. 结果初步筛选：高置信度
    print("Filtering high-confidence alignment results...")
    awk_command = f"awk '$11 < 1e-5 && $4 > 70' {output_dir}/{output_file} > {output_dir}/high_confidence_hits.tsv"
    subprocess.run(awk_command, shell=True, check=True)
    
    print("Extracting unique query IDs...")
    cut_command = f"cut -f1 {output_dir}/high_confidence_hits.tsv | sort | uniq > {output_dir}/candidate_ids.txt"
    subprocess.run(cut_command, shell=True, check=True)
    
    print("Extracting candidate sequences...")
    seqkit_command = f"seqkit grep -f {output_dir}/candidate_ids.txt {input_data_dir}/{input_file} > {output_dir}/candidates.faa"
    subprocess.run(seqkit_command, shell=True, check=True)
    
    print("DIAMOND alignment and filtering completed successfully.")

# %%
def main(config):
    if config['build_db']:
        build_custom_reference_db(config)
    diamond_alignment(config)

if __name__ == "__main__":
    main(config)
# %%
