

### Step 1: Sequence Alignment

---

Input: 
- ***input_file**: sequences to be aligned to the IR database, one of the following:
  1. gene.uniGeneset_MABAll.faa: all the non-redundant protein sequences of metagenomic data, which could align all the sequences from samples
  2. MAG9_CDS.faa: Only align certain BIN's sequence, i.e. the protein sequences to be aligned of MAG9

- **save_data_dir**: the directory to save the database
- **keywords_url:** the keywords to search for the database

Output:

- **\**custom_IR_db.dmnd**: custom IR optogenetic switch database
- custom_IR_db.faa: fasta sequence of custom IR optogenetic switch
- custom_IR_db_clustered.faa: clustered fasta sequence of custom IR optogenetic switch
- 
- diamond_results.tsv: the details of all alignments. Note: this file hasn't filtered  

-** high_confidence_hits.tsv: the details of filtered alignments. 

-** candidates.faa: candidate sequences

-** candidate_ids.txt: the id of candidate sequences

\- OUTPUT_DIR="../results/MAG9"





### Step 2: Domain Scanning 

---

