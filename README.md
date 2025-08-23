# METTL5-mediated Translation Control Motif Discovery Pipeline

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pending-orange)](https://doi.org/)
[![GitHub Stars](https://img.shields.io/github/stars/yourusername/mettl5-motif-discovery?style=social)](https://github.com/yourusername/mettl5-motif-discovery)

## Table of Contents
- [Overview](#-overview)
- [Scientific Background](#-scientific-background)
- [Key Features](#-key-features)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Detailed Usage](#-detailed-usage)
- [Input Data Format](#-input-data-format)
- [Output Description](#-output-description)
- [Algorithm Details](#-algorithm-details)
- [Example Analysis](#-example-analysis)
- [Troubleshooting](#-troubleshooting)
- [Citation](#-citation)
- [Contact](#-contact)

## Overview

This repository contains a comprehensive computational framework for discovering and analyzing 5'UTR regulatory motifs that control mRNA translation efficiency. The pipeline was specifically developed to identify sequence elements regulated by METTL5-mediated 18S rRNA methylation in prostate cancer, but can be applied to any translational regulation study.

### Key Discovery
Our pipeline successfully identified the **GCACGN{1,5}CC** motif as a critical regulatory element in 5'UTRs of genes whose translation is controlled by METTL5-mediated rRNA methylation.

## Scientific Background

METTL5 is an RNA methyltransferase that catalyzes N6-methyladenosine (m6A) modification at position 1832 of 18S ribosomal RNA. This modification plays a crucial role in:

- **Cancer Progression**: METTL5 is upregulated in prostate cancer and correlates with poor prognosis
- **Metabolic Regulation**: Controls mitochondrial oxidative phosphorylation
- **Selective Translation**: Regulates specific mRNAs through 5'UTR motifs
- **Molecular Axis**: Functions through the METTL5/IRF7/DNA2 regulatory pathway

## Key Features

### Motif Discovery Capabilities
- **Multiple Algorithm Support**
  - Exact k-mer matching
  - IUPAC wildcard patterns (R, Y, M, K, S, W, H, B, V, D, N)
  - Gapped motif discovery (e.g., GCACGN{2,4}CC)
  - Position Weight Matrix (PWM) analysis
  - Mismatch-tolerant searching

### Statistical Analysis
- **Enrichment Metrics**
  - UP/DOWN regulation ratios
  - Fold change calculations
  - Frequency analysis
- **Statistical Tests**
  - Fisher's exact test
  - Hypergeometric test
  - Binomial test for strand specificity
- **Strand Specificity Analysis**
  - RNA vs DNA-level regulation detection
  - Reverse complement comparison

### Visualization Suite
- Frequency distribution plots
- Enrichment scatter plots
- Strand specificity analysis
- Motif similarity networks
- Heatmaps of significant motifs

### Performance Features
- Parallel processing with multiprocessing
- Optimized for large-scale genomic data
- Memory-efficient sequence handling
- Configurable CPU core usage

## Installation

### System Requirements
- **Operating System**: Linux, macOS, or Windows
- **Python**: 3.7 or higher
- **Memory**: Minimum 8GB RAM (16GB recommended)
- **Storage**: 10GB free space for analysis outputs

### Step 1: Clone Repository
```bash
git clone https://github.com:yangrf2021/5-UTR-Regulatory-Motifs-Discovery-Pipeline.git
```

### Step 2: Create Virtual Environment (Recommended)
```bash
python -m venv motif_env
source motif_env/bin/activate  # On Windows: motif_env\Scripts\activate
```

### Step 3: Install Dependencies
```bash
pip install -r requirements.txt
```

### Dependencies List
```python
pandas>=1.3.0
numpy>=1.20.0
biopython>=1.79
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
networkx>=2.6
tqdm>=4.62.0
multiprocessing
```

## Quick Start

### Basic Analysis
```bash
# Run with default parameters
python src/run_analysis.py
```

### Custom Analysis
```bash
# Run with custom parameters
python src/run_analysis.py \
    --max-utr-length 500 \
    --min-frequency 0.05 \
    --fold-change 1.5 \
    --enable-gapped \
    --enable-pwm
```

## Detailed Usage

### Command Line Interface

```bash
python src/run_analysis.py [OPTIONS]
```

### Parameters

#### Core Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max-utr-length` | 500 | Maximum 5'UTR length to analyze (nt) |
| `--min-frequency` | 0.05 | Minimum motif frequency threshold |
| `--fold-change` | 1.5 | Fold change threshold for enrichment |
| `--min-occur-freq` | 0.1 | Minimum occurrence frequency |

#### Motif Discovery Options
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max-wildcards` | 1 | Maximum IUPAC wildcards allowed |
| `--max-mismatches` | 1 | Maximum mismatches allowed |
| `--enable-gapped` | False | Enable gapped motif discovery |
| `--enable-pwm` | False | Enable PWM-based analysis |
| `--enable-structural` | False | Enable structural motif analysis |
| `--enable-pwm-gapped` | False | Enable PWM with gaps |

#### Analysis Range
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-length` | 4 | Minimum motif length (bp) |
| `--max-length` | 8 | Maximum motif length (bp) |

#### Performance
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n-cores` | Auto | Number of CPU cores (default: all-1) |
| `--seed` | 42 | Random seed for reproducibility |

### Python API Usage

```python
from src.motif_analyzer import MotifAnalyzer

# Initialize analyzer
analyzer = MotifAnalyzer(
    max_utr_length=500,
    min_frequency=0.05,
    fold_change_threshold=1.5,
    enable_gapped=True,
    enable_pwm=True
)

# Run analysis for specific length
results = analyzer.analyze_single_length(
    length=7,
    n_cores=8,
    max_wildcards=1
)

# Access results
significant_motifs = results[results['Significant']]
up_enriched = significant_motifs[significant_motifs['Enrichment'] == 'UP']
```

## Input Data Format

### Required Files
Place these files in the working directory:

1. **te_up_utr.fasta** - 5'UTR sequences of translationally upregulated genes
2. **te_down_utr.fasta** - 5'UTR sequences of translationally downregulated genes

### FASTA Format Specification
```fasta
>GENE_ID | Selection: Transcript_Selection_Method | Additional_Info
ATGGCGCACGNNNNCCTGATCGATCGATCGATCGATCGATCGATCG...
```

#### Header Components
- **GENE_ID**: Unique gene identifier
- **Selection**: Transcript selection method (e.g., "Longest_UTR", "Most_Abundant")
- **Additional_Info**: Optional metadata

### Example Input File
```fasta
>ENSG00000141510 | Selection: Longest_UTR | Gene: TP53
ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAA
>ENSG00000171862 | Selection: Most_Abundant | Gene: PTEN
GCGCACGNNNCCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
```

## Output Description

### Directory Structure
```
motif_analysis_results/
│
├── analysis_summary.csv              # Overall summary statistics
├── analysis_summary.png              # Summary visualization
├── README.md                         # Auto-generated report
│
└── motif_results_length_N/          # Results for each length
    ├── 01_motif_analysis.csv        # All discovered motifs
    ├── 02_significant_motifs.csv    # Filtered significant motifs
    ├── 03_up_enriched_motifs.csv    # UP-enriched motifs
    ├── 03_down_enriched_motifs.csv  # DOWN-enriched motifs
    ├── 04_up_clustered_motifs.csv   # Clustered UP motifs
    ├── 04_down_clustered_motifs.csv # Clustered DOWN motifs
    │
    ├── visualizations/
    │   ├── frequency_distribution.png
    │   ├── enrichment_scatter.png
    │   ├── strand_specificity.png
    │   ├── significant_heatmap.png
    │   └── motif_network.png
    │
    └── motif_details/
        ├── up_motif_GCACGCC_frequencies.csv
        ├── up_motif_GCACGCC_details.txt
        └── up_motif_GCACGCC_genes.txt
```

### Output File Descriptions

#### 1. Summary Files
- **analysis_summary.csv**: Overview of all analyzed lengths
- **README.md**: Detailed analysis report with key findings

#### 2. Motif Analysis Files
- **01_motif_analysis.csv**: Complete list of discovered motifs with statistics
- **02_significant_motifs.csv**: Motifs passing significance thresholds
- **03_*_enriched_motifs.csv**: Direction-specific enriched motifs
- **04_*_clustered_motifs.csv**: Similar motifs grouped together

#### 3. Detailed Analysis Files
- **frequencies.csv**: Position-specific nucleotide frequencies
- **details.txt**: Comprehensive motif statistics
- **genes.txt**: List of genes containing the motif

### Output Data Fields

#### Main Analysis CSV
| Column | Description |
|--------|-------------|
| Pattern | Motif sequence or pattern |
| Type | Motif type (exact/mismatch/gapped/pwm) |
| UP_present | Count in upregulated genes |
| DOWN_present | Count in downregulated genes |
| UP_frequency | Frequency in upregulated set |
| DOWN_frequency | Frequency in downregulated set |
| UP_to_DOWN_ratio | Enrichment ratio |
| Strand_specificity | Strand bias score |
| Strand_pvalue | Statistical significance |
| Significant | Boolean significance flag |
| Enrichment | Direction of enrichment |

## Algorithm Details

### Motif Discovery Methods

#### 1. Exact K-mer Matching
```python
# Traditional k-mer counting
for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i+k]
    count[kmer] += 1
```

#### 2. IUPAC Wildcard Support
```python
IUPAC_MAP = {
    'R': ['A', 'G'],  # puRine
    'Y': ['C', 'T'],  # pYrimidine
    'M': ['A', 'C'],  # aMino
    'K': ['G', 'T'],  # Keto
    'S': ['G', 'C'],  # Strong
    'W': ['A', 'T'],  # Weak
    'N': ['A', 'C', 'G', 'T']  # aNy
}
```

#### 3. Gapped Motif Discovery
```python
# Format: PREFIX{min,max}SUFFIX
# Example: GCACG{2,4}CC matches:
#   GCACGNNCC (gap=2)
#   GCACGNNNCC (gap=3)
#   GCACGNNNNCC (gap=4)
```

#### 4. Position Weight Matrix (PWM)
```python
# PWM scoring for sequence matching
def calculate_pwm_score(pwm, sequence):
    score = 0
    for i, base in enumerate(sequence):
        score += pwm[base][i]
    return score
```

### Statistical Methods

#### Enrichment Calculation
```python
enrichment_ratio = (up_frequency / down_frequency)
significance = fisher_exact_test(contingency_table)
```

#### Strand Specificity
```python
forward_matches = count_matches(motif, sequences)
reverse_matches = count_matches(reverse_complement(motif), sequences)
specificity = forward_matches / reverse_matches
```

## Example Analysis

### Sample Run Output
```
Starting METTL5 motif analysis pipeline
Parameters: {'max_utr_length': 500, 'min_frequency': 0.05, ...}

Analyzing motifs of length 6
Found 127 significant motifs
  UP-enriched: 45
  DOWN-enriched: 82

Analyzing motifs of length 7
Found 156 significant motifs
  UP-enriched: 52
  DOWN-enriched: 104

Top discovered motif: GCACGN{2,4}CC
  - Enrichment: 2.8-fold in DOWN-regulated genes
  - Frequency: 15.3% of target genes
  - Strand specificity: 3.2 (p < 0.001)
  - Genes affected: 147

Analysis complete!
Results saved in motif_analysis_results/
```

### Biological Validation Results
The discovered GCACGN{2,4}CC motif has been validated through:
- Dual-luciferase reporter assays
- Mutation analysis
- Polysome profiling
- Clinical sample validation

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: Memory Error
```
MemoryError: Unable to allocate array
```
**Solution**: Reduce `--max-utr-length` or process in batches

#### Issue 2: No Significant Motifs Found
```
WARNING: No significant motifs identified
```
**Solution**: Adjust parameters:
- Decrease `--min-frequency`
- Decrease `--fold-change`
- Enable `--enable-gapped` and `--enable-pwm`

#### Issue 3: Slow Performance
```
Processing taking too long
```
**Solution**: 
- Increase `--n-cores`
- Reduce motif length range
- Disable structural analysis if not needed

### FAQ

**Q: Can I use this for other organisms?**
A: Yes, the pipeline is organism-agnostic. Just provide appropriate FASTA files.

**Q: What if I don't have UP/DOWN regulated gene sets?**
A: You can compare any two gene sets (e.g., treated vs control).

**Q: How do I interpret strand specificity?**
A: 
- Ratio ≈ 1: No strand preference (DNA-level)
- Ratio > 2: Strong strand preference (RNA-level)

## Citation

This code is associated with a manuscript currently under preparation. Please check back for citation information once published.

## Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/NewFeature`)
3. Commit changes (`git commit -m 'Add NewFeature'`)
4. Push to branch (`git push origin feature/NewFeature`)
5. Open a Pull Request

## Contact

**Principal Investigator**
- Name: Ruifeng Yang
- Email: yangrf1996@163.com

**Technical Support**
- GitHub Issues: [Create Issue](https://github.com:yangrf2021/5-UTR-Regulatory-Motifs-Discovery/issues)
- Email: yangrf1996@163.com

## License

This project is licensed under the MIT License:

```
MIT License

Copyright (c) 2025 Ruifeng Yang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Acknowledgments

- Thanks to all contributors and collaborators
- Computational resources provided by Shanghai Ninth People's Hospital
- Supported by grants: No. 82473192, 82172741
- Special thanks to the Bioinformatics Core Facility

## Disclaimer

This software is provided for research purposes only. Any clinical or diagnostic use requires additional validation and regulatory approval.

## Performance Metrics

| Dataset Size | Runtime | Memory Usage | CPU Cores |
|-------------|---------|--------------|-----------|
| 100 genes | 2 min | 1 GB | 4 |
| 1,000 genes | 15 min | 4 GB | 8 |
| 10,000 genes | 2 hours | 16 GB | 16 |

## Version History

- **v1.0.0** (2025-01): Initial release
  - Core motif discovery algorithms
  - Statistical analysis framework
  - Visualization suite
  
- **v0.9.0** (2023-12): Beta version
  - PWM implementation
  - Gapped motif support

## TODO

- [ ] Add RNA secondary structure prediction
- [ ] Implement machine learning-based motif scoring
- [ ] Add support for paired-end RNA-seq data
- [ ] Create web interface
- [ ] Add Docker container support

---

**Last Updated**: January 2025  
**Version**: 1.0.0  
**Status**: Active Development

[![Made with Python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com:yangrf2021/5-UTR-Regulatory-Motifs-Discovery/graphs/commit-activity)
