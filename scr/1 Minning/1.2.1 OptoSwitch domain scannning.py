"""
===
This script further filters opto switch candidates based on functional domains.
===
Input:
- unique_sequences.faa: all the non-redundant protein sequences of metagenomic data, which could align all the sequences from samples
Output:
- unique_sequences.faa: fasta sequence of opto switch candidates
- ** unique_sequence_ids.txt # Key results, will also be analyzed later for genome regulation decription
- pfamscan_results_new_3domain.txt
"""
# %%
import os
import subprocess
import sys

# Configuration
config = {
    'save_data_dir': "../../data/Fig2/database",
    'database_dir': "../../data/Fig2/database", 
    'output_dir': "../../data/Fig2/Fig2bc.IRswitch.pyresults/MAG73",
    'pfam_db_url': "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    'pfam_dat_url': "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz",
    'build_db': False
}

def run_command(command, cwd=None):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True, cwd=cwd)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        raise

def check_and_install_pfam_scan():
    try:
        subprocess.run("pfam_scan.pl -h", shell=True, check=True, capture_output=True, text=True)
        print("Pfam Scan is already installed.")
    except subprocess.CalledProcessError:
        print("Pfam Scan not found. Installing via conda...")
        run_command("conda install -y -c bioconda pfam_scan")

def download_pfam_database(config):
    database_dir = config['database_dir']
    pfam_db_path = os.path.join(database_dir, "Pfam-A.hmm")
    pfam_dat_path = os.path.join(database_dir, "Pfam-A.hmm.dat")

    pfam_db_exists = os.path.exists(pfam_db_path)
    pfam_dat_exists = os.path.exists(pfam_dat_path)

    if pfam_db_exists and pfam_dat_exists:
        print("Pfam database already exists. Skipping download.")
        return
    print("Downloading Pfam database...")
    os.makedirs(database_dir, exist_ok=True)
    if not pfam_db_exists:
        run_command(f"curl -O {config['pfam_db_url']}", cwd=database_dir)
        run_command("gunzip -f Pfam-A.hmm.gz", cwd=database_dir)
    if not pfam_dat_exists:
        run_command(f"curl -O {config['pfam_dat_url']}", cwd=database_dir)
        run_command("gunzip -f Pfam-A.hmm.dat.gz", cwd=database_dir)

def build_domain_database(config):
    database_dir = config['database_dir']
    target_db_dir = os.path.join(database_dir, "target_domains_db")
    os.makedirs(target_db_dir, exist_ok=True)
    target_domains_hmm = os.path.join(target_db_dir, "Pfam-A.hmm")

    if os.path.exists(target_domains_hmm) and all(
        os.path.exists(target_domains_hmm + ext) for ext in [".h3f", ".h3i", ".h3m", ".h3p"]
    ):
        print("Domain database already exists.")
        return

    print("Building domain database for PAS, GAF, and PHY...")
    pas_hmm = os.path.join(database_dir, "PAS.hmm")
    gaf_hmm = os.path.join(database_dir, "GAF.hmm")
    phy_hmm = os.path.join(database_dir, "PHY.hmm")
    run_command(f"hmmfetch {os.path.join(database_dir, 'Pfam-A.hmm')} PAS > {pas_hmm}")
    run_command(f"hmmfetch {os.path.join(database_dir, 'Pfam-A.hmm')} GAF > {gaf_hmm}")
    run_command(f"hmmfetch {os.path.join(database_dir, 'Pfam-A.hmm')} PHY > {phy_hmm}")

    target_hmm_tmp = os.path.join(database_dir, "target_domains.hmm")
    run_command(f"cat {pas_hmm} {gaf_hmm} {phy_hmm} > {target_hmm_tmp}")

    for f in [pas_hmm, gaf_hmm, phy_hmm]:
        if os.path.exists(f):
            os.remove(f)

    run_command(f"mv {target_hmm_tmp} {target_domains_hmm}")
    run_command(f"hmmpress {target_domains_hmm}")

    pfam_dat_src = os.path.join(database_dir, "Pfam-A.hmm.dat")
    pfam_dat_dst = os.path.join(target_db_dir, "Pfam-A.hmm.dat")
    if os.path.exists(pfam_dat_src) and not os.path.exists(pfam_dat_dst):
        run_command(f"mv {pfam_dat_src} {pfam_dat_dst}")

def run_pfam_scan(config):
    output_dir = config['output_dir']
    target_db_dir = os.path.join(config['database_dir'], "target_domains_db")
    candidates_faa = os.path.join(output_dir, "candidates.faa")
    pfamscan_results = os.path.join(output_dir, "pfamscan_results_new_3domain.txt")

    print("Running Pfam Scan...")

    if os.path.exists(pfamscan_results):
        print(f"Output file {pfamscan_results} already exists. Removing it before running Pfam Scan.")
        os.remove(pfamscan_results)

    run_command(
        f"pfam_scan.pl -fasta {candidates_faa} -dir {target_db_dir} -outfile {pfamscan_results}"
    )

def filter_sequences(config):
    output_dir = config['output_dir']
    pfamscan_results = os.path.join(output_dir, "pfamscan_results_new_3domain.txt")
    unique_sequence_ids = os.path.join(output_dir, "unique_sequence_ids.txt")
    unique_sequences_faa = os.path.join(output_dir, "unique_sequences.faa")
    candidates_faa = os.path.join(output_dir, "candidates.faa")

    print("Filtering sequences...")

    run_command(
        f"grep -v '^#' {pfamscan_results} | awk '{{print $1}}' | sort | uniq | grep -v '^$' > {unique_sequence_ids}"
    )

    awk_script = (
        f"awk 'BEGIN {{ while (getline < \"{unique_sequence_ids}\") ids[$1]=1 }} "
        f"/^>/ {{ header=substr($1,2); print_seq=(header in ids) }} "
        f"print_seq' {candidates_faa} > {unique_sequences_faa}"
    )
    run_command(awk_script)

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
# %%
