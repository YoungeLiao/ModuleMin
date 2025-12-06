
"""
===
*Update 20251031
This script further filters opto switch candidates based on functional domains.
===
Input:
- unique_sequences.faa: all the non-redundant protein sequences of metagenomic data, which could align all the sequences from samples
Output:
- unique_sequences.faa: fasta sequence of opto switch candidates
- ** unique_sequence_ids.txt # Key results, will also be analyzed later for genome regulation decription
- pfamscan_results_new_3domain.txt
"""
from Bio import SeqIO
import os
import re
import subprocess

# Configuration
config = {
    'save_data_dir': "../../data/Fig2/database",
    'database_dir': "../../data/Fig2/database", 
    'output_dir': "../../data/Fig2/Fig2bc.IRswitch.pyresults/MAG73",
    'pfam_db_url': "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    'pfam_dat_url': "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz",
    'build_db': False
}

# 运行命令的函数：
def run_command(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}")
        print(e.stderr)
        raise

# 检查并安装 Pfam Scan：
def check_and_install_pfam_scan():
    try:
        subprocess.run("pfam_scan.pl -h", shell=True, check=True, capture_output=True, text=True)
        print("Pfam Scan is already installed.")
    except subprocess.CalledProcessError:
        print("Installing Pfam Scan...")
        run_command("conda install -c bioconda pfam_scan")

# 下载 Pfam 数据库：
def download_pfam_database(config):
    database_dir = config['database_dir']
    pfam_db_path = os.path.join(database_dir, "Pfam-A.hmm")
    pfam_dat_path = os.path.join(database_dir, "Pfam-A.hmm.dat")

    if not os.path.exists(pfam_db_path) or not os.path.exists(pfam_dat_path):
        print("Downloading Pfam database...")
        run_command(f"curl -O {config['pfam_db_url']}")
        run_command("gunzip Pfam-A.hmm.gz")
        run_command(f"curl -O {config['pfam_dat_url']}")
        run_command("gunzip Pfam-A.hmm.dat.gz")
        run_command(f"mv Pfam-A.hmm Pfam-A.hmm.dat {database_dir}")
    else:
        print("Pfam database already exists.")

def build_domain_database(config):
    """
    Build a reduced domain HMM database containing all Pfam models whose NAME starts
    with any of the prefixes (e.g. PAS, GAF, PHY), so subfamilies like PAS_2 are included.
    """
    database_dir = config['database_dir']
    target_db_dir = os.path.join(database_dir, "target_domains_db")
    pfam_hmm = os.path.join(database_dir, "Pfam-A.hmm")
    pfam_dat = os.path.join(database_dir, "Pfam-A.hmm.dat")
    prefixes = ["PAS", "GAF", "PHY"]

    if not os.path.exists(pfam_hmm):
        raise FileNotFoundError(f"{pfam_hmm} not found. Please download Pfam-A.hmm first.")

    # Parse Pfam-A.hmm to collect model names whose NAME starts with any prefix
    model_names = []
    with open(pfam_hmm, 'r') as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    name = parts[1]
                    for p in prefixes:
                        if name.startswith(p):
                            model_names.append(name)
                            break

    if not model_names:
        print("No matching Pfam models found for prefixes:", prefixes)
        return

    os.makedirs(target_db_dir, exist_ok=True)
    temp_hmms = []

    # Fetch each matching model into a temporary .hmm file in database_dir
    for name in model_names:
        tmp_hmm = os.path.join(database_dir, f"{name}.hmm")
        run_command(f"hmmfetch {pfam_hmm} {name} > {tmp_hmm}")
        temp_hmms.append(tmp_hmm)

    # Concatenate all fetched models into single target HMM file
    target_domains_hmm = os.path.join(target_db_dir, "Pfam-A.hmm")
    run_command(f"cat {' '.join(temp_hmms)} > {target_domains_hmm}")

    # Remove temporary individual model files
    for tmp in temp_hmms:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except Exception:
                # fallback to shell remove for consistency with run_command behavior
                run_command(f"rm -f {tmp}")

    # Prepare .dat file if available
    if os.path.exists(pfam_dat):
        run_command(f"mv {pfam_dat} {os.path.join(target_db_dir, 'Pfam-A.hmm.dat')}")

    # Press the concatenated HMM database
    run_command(f"hmmpress {target_domains_hmm}")

    # Move/rename hmmpress output files if they exist (ensure consistent names)
    for suf in [".h3f", ".h3i", ".h3m", ".h3p"]:
        src = target_domains_hmm + suf
        dst = os.path.join(target_db_dir, "Pfam-A.hmm" + suf)
        if os.path.exists(src):
            run_command(f"mv {src} {dst}")

    print(f"Built target domain DB with {len(model_names)} models -> {target_domains_hmm}")

# 运行 Pfam Scan：
def run_pfam_scan(config):
    output_dir = config['output_dir']
    target_db_dir = os.path.join(config['database_dir'], "target_domains_db")
    candidates_faa = os.path.join(output_dir, "candidates.faa")
    pfamscan_results = os.path.join(output_dir, "pfamscan_results_new_3domain.txt")

    print("Running Pfam Scan...")
    run_command(f"pfam_scan.pl -fasta {candidates_faa} -dir {target_db_dir} -outfile {pfamscan_results}")

# 过滤序列：
def filter_sequences(config):
    output_dir = config['output_dir']
    pfamscan_results = os.path.join(output_dir, "pfamscan_results_new_3domain.txt")
    unique_sequence_ids = os.path.join(output_dir, "unique_sequence_ids.txt")
    unique_sequences_faa = os.path.join(output_dir, "unique_sequences.faa")

    print("Filtering sequences...")
    # run_command(f"grep -v '^#' {pfamscan_results} | cut -f1 | sort | uniq | grep -v '^$' > {unique_sequence_ids}")
    run_command(f"grep -v '^#' {pfamscan_results} | cut -d ' ' -f1 | sort | uniq | grep -v '^$' > {unique_sequence_ids}")
    # grep -v "^#" $OUTPUT_DIR/pfamscan_results_new_3domain.txt | awk '{print $1}' | sort | uniq | grep -v '^$' > $OUTPUT_DIR/unique_sequence_ids.txt

    run_command(f"awk 'BEGIN {{ while (getline < \"{unique_sequence_ids}\") ids[$1] = 1 }} "
                f"/^>/ {{ header = substr($1, 2); print_seq = (header in ids) }} "
                f"print_seq' {os.path.join(output_dir, 'candidates.faa')} > {unique_sequences_faa}")

    print("Number of unique genes:")
    run_command(f"grep -c '^>' {unique_sequences_faa}")

def main(config):
    check_and_install_pfam_scan()
    download_pfam_database(config)
    build_domain_database(config)
    run_pfam_scan(config)
    filter_sequences(config)

if __name__ == "__main__":
    main(config)


# # Step 2: Filtering targeted domain
# def filter_by_domain(pfamscan_file, domain_prefixes, candidates_faa, output_file):
#     """
#     Filter candidate sequences whose pfam annotations contain any domain that starts with
#     one of the given prefixes (e.g. PAS, PAS_2, GAF, PHY, PHY_X).
#     - pfamscan_file: path to pfam_scan.pl output (text)
#     - domain_prefixes: iterable of prefixes, e.g. ["PAS","GAF","PHY"]
#     - candidates_faa: path to candidate fasta file
#     - output_file: path to write filtered fasta
#     """
#     if not os.path.exists(pfamscan_file):
#         raise FileNotFoundError(f"pfamscan file not found: {pfamscan_file}")
#     if not os.path.exists(candidates_faa):
#         raise FileNotFoundError(f"candidates fasta not found: {candidates_faa}")

#     # compile regex to match domains that start with any of the prefixes (word boundary + prefix)
#     prefix_pattern = r'\b(' + '|'.join(re.escape(p) for p in domain_prefixes) + r')\w*\b'
#     prog = re.compile(prefix_pattern)

#     matched_ids = set()
#     with open(pfamscan_file, 'r') as fh:
#         for line in fh:
#             if line.startswith("#") or not line.strip():
#                 continue
#             cols = line.strip().split()
#             if len(cols) < 2:
#                 continue
#             seq_id = cols[0]
#             # join remaining columns to search for domain names robustly
#             rest = " ".join(cols[1:])
#             if prog.search(rest):
#                 matched_ids.add(seq_id)

#     # if no matches, write empty file and return
#     if not matched_ids:
#         # ensure output dir exists
#         os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
#         open(output_file, 'w').close()
#         print(f"No domain matches found in {pfamscan_file}; wrote empty {output_file}")
#         return

#     # use SeqIO.index for memory-efficient lookup
#     seq_index = SeqIO.index(candidates_faa, "fasta")
#     filtered_records = []
#     for sid in matched_ids:
#         if sid in seq_index:
#             filtered_records.append(seq_index[sid])
#         else:
#             # try alternative id forms (in case pfam uses full header); match by startswith
#             for key in seq_index:
#                 if key.startswith(sid):
#                     filtered_records.append(seq_index[key])
#                     break

#     # write results
#     SeqIO.write(filtered_records, output_file, "fasta")
#     seq_index.close()
#     print(f"Wrote {len(filtered_records)} filtered sequences to {output_file}")



# # Main***
# if __name__ == "__main__":
#     # 输入文件和参数
#     bin = '/MAG73' # '' # 'MAG90'
#     # normalize output_dir so joining works whether bin has leading slash or not
#     base_dir = "../../data/Fig2/Fig2bc.IRswitch.pyresults"
#     output_dir = os.path.join(base_dir, bin.lstrip('/'))
#     pfamscan_file = os.path.join(output_dir, "pfamscan_results_new_3domain.txt")
#     candidates_faa = os.path.join(output_dir, "candidates.faa")
#     domains = ["PAS", "GAF", "PHY"]
#     output_file = os.path.join(output_dir, "filtered_ir_proteins_all.fasta")

#     # Step 1: 过滤结构域（包括以 PAS/GAF/PHY 开头的子类）
#     filter_by_domain(pfamscan_file, domains, candidates_faa, output_file)
