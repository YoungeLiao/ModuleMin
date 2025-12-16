# ModuleMin

ModuleMin is a toolkit for mining and validating optogenetic elements, identifying keystone taxa, and analyzing modular functions of microbiomes. This project integrates metagenomics, 16S rRNA sequencing, and metabolomics data to provide a comprehensive workflow for discovering and validating functional modules in microbial communities.

## Main Features

1. **Optogenetic Element Mining and Validation**
   - Mine potential optogenetic elements from metagenomic data
   - Validate candidate elements through sequence alignment and domain analysis
   - Construct phylogenetic trees to analyze evolutionary relationships of elements

2. **Keystone Taxa Identification**
   - Analyze the structure-function relationships in microbial communities
   - Identify microbial taxa that play key roles in ecosystems
   - Integrate multi-omics data for correlation analysis

3. **Microbiome Modular Function Analysis**
   - Construct microbial co-occurrence networks
   - Analyze the composition and stability of functional modules
   - Evaluate the impact of environmental factors on module functions

## Project Structure

```
ModuleMin/
├── scr/                          # Main analysis scripts
│   ├── 1 Minning/                # Optogenetic element mining scripts
│   │   ├── 1.1 OptoSwitch Sequence Alignment.py  # Sequence alignment script
│   │   ├── 1.2.1 Opto domain scanning.py         # Domain scanning script
│   │   ├── 1.2.2 OptoSwitch domain filtering.py  # Domain filtering script
│   │   └── 1.5 Visualization_IRDomain.py          # IR domain visualization script
│   ├── 2 Phylogenetics/           # Phylogenetic analysis scripts
│   │   └── 1.3.1 Phylogenetic tree construction.sh  # Phylogenetic tree construction script
│   ├── 3.1 FunTaxMerge.py         # Function and taxa integration analysis script
│   └── Visualization/             # Result visualization scripts
├── data/                          # Data directory
│   ├── Fig1-4/                    # Data for each figure
│   ├── ex-fig1-7/                 # Supplementary figure data
│   ├── rawdata.shared/            # Raw shared data
│   └── results_data/              # Result data
└── Rscript/                       # R language analysis scripts
    ├── Fig1-4/                    # R scripts for each figure
    └── ex-fig1-7/                 # R scripts for supplementary figures
```

## Main Script Descriptions

### Optogenetic Element Mining

1. **Sequence Alignment (1.1 OptoSwitch Sequence Alignment.py)**
   - Download optogenetic element reference sequences from public databases
   - Use DIAMOND tool for sequence alignment
   - Filter high-confidence alignment results

2. **Domain Scanning (1.2.1 Opto domain scanning.py)**
   - Scan candidate sequences for functional domains using Pfam database
   - Identify key photosensory domains such as PAS, GAF, and PHY
   - Build custom functional domain database

3. **Domain Filtering (1.2.2 OptoSwitch domain filtering.py)**
   - Filter candidate sequences based on domain analysis results
   - Retain sequences containing key photosensory domains

4. **IR Domain Visualization (1.5 Visualization_IRDomain.py)**
   - Visualize the distribution of different domains in sequences
   - Generate analysis charts of domain composition and positions

### Phylogenetic Analysis

1. **Phylogenetic Tree Construction (1.3.1 Phylogenetic tree construction.sh)**
   - Perform multiple sequence alignment using MAFFT
   - Trim alignment results using TrimAl
   - Build phylogenetic trees using RAxML

### Keystone Taxa Identification

1. **Function and Taxa Integration Analysis (3.1 FunTaxMerge.py)**
   - Integrate metagenomic, 16S rRNA, and MAG data
   - Analyze abundance changes of different taxa
   - Identify potential keystone taxa

## Usage Instructions

### Requirements

- Python 3.7+
- R 3.6+
- Bioinformatics tools:
  - DIAMOND
  - CD-HIT
  - MAFFT
  - TrimAl
  - RAxML
  - PfamScan
  - seqkit

### Installing Dependencies

```bash
# Python dependencies
pip install biopython pandas matplotlib seaborn numpy

# R dependencies
install.packages(c("ggplot2", "dplyr", "tidyr", "vegan"))

# Bioinformatics tools (using conda)
conda install -c bioconda diamond cd-hit mafft trimal raxml pfam_scan seqkit
```

### Running Examples

1. **Mining Optogenetic Elements**
```bash
# Modify input/output paths in the configuration file
python scr/1\ Minning/1.1\ OptoSwitch\ Sequence\ Alignment.py

# Domain analysis
python scr/1\ Minning/1.2.1\ Opto\ domain\ scanning.py

# Domain filtering
python scr/1\ Minning/1.2.2\ OptoSwitch\ domain\ filtering.py
```

2. **Building Phylogenetic Trees**
```bash
# Modify input/output paths in the script
bash scr/2\ Phylogenetics/1.3.1\ Phylogenetic\ tree\ construction.sh
```

3. **Identifying Keystone Taxa**
```bash
# Modify input/output paths in the configuration file
python scr/3.1\ FunTaxMerge.py
```

## Result Description

- **Candidate Optogenetic Elements**: `results_data/StringentCandidates.fasta`
- **Domain Analysis Results**: `results_data/pfamscan_results.xlsx`
- **Phylogenetic Tree**: `results_data/PhyCandidates.txt`
- **Keystone Taxa Analysis**: `scr/All_IR_Genus_Abundance_BySample.xlsx`

## Citation

If you use ModuleMin in your research, please cite the following literature:

[Add relevant citation information here]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Bug reports and feature requests are welcome. If you want to contribute code, please:

1. Fork this project
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Contact

For questions or suggestions, please contact us through:

liaoy21@mails.tsinghua.edu.cn
