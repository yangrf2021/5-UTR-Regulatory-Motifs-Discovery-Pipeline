# 5'UTR Regulatory Motif Discovery Pipeline for METTL5-mediated Translation Control

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/yangrf2021/5-UTR-Regulatory-Motifs-Discovery-Pipeline)

## Table of Contents
- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Module Documentation](#detailed-module-documentation)
- [Algorithm Workflow](#algorithm-workflow)
- [Input Requirements](#input-requirements)
- [Output Description](#output-description)
- [Parameter Configuration](#parameter-configuration)
- [Example Usage](#example-usage)
- [Troubleshooting](#troubleshooting)
- [Performance Considerations](#performance-considerations)
- [Citation](#citation)
- [Contact](#contact)

## Overview

This repository contains a comprehensive bioinformatics pipeline for discovering and analyzing 5'UTR regulatory motifs that control mRNA translation efficiency. The pipeline was specifically developed to identify sequence elements regulated by METTL5-mediated 18S rRNA methylation, ultimately revealing the **GCACGN{2,4}CC** motif as a critical regulatory element in prostate cancer, but can be applied to any translational regulation study.

### Core Discovery
Our pipeline successfully identified that METTL5-mediated 18S rRNA m6A modification at position 1832 selectively regulates the translation of mRNAs containing specific 5'UTR motifs, particularly affecting mitochondrial metabolism genes through the METTL5/IRF7/DNA2 axis.

## Scientific Background

### METTL5 and rRNA Methylation
METTL5 (Methyltransferase Like 5) is an RNA methyltransferase that catalyzes N6-methyladenosine (m6A) modification at position 1832 of 18S ribosomal RNA. This modification:
- Affects ribosome assembly and function
- Selectively regulates mRNA translation efficiency
- Shows upregulation in various cancers
- Correlates with poor prognosis in prostate cancer

### Translation Control Mechanism
The pipeline identifies how rRNA methylation leads to selective translation through:
1. **Ribosome Specialization**: Modified ribosomes show preference for specific mRNA features
2. **5'UTR Recognition**: Certain sequence motifs in 5'UTRs are preferentially translated
3. **Metabolic Reprogramming**: Affects mitochondrial gene expression and OXPHOS
4. **Cancer Progression**: Promotes tumor growth through metabolic adaptation

## Key Features

### 1. Multi-Strategy Motif Discovery
- **Exact k-mer matching**: Traditional sequence counting
- **IUPAC wildcard patterns**: Support for degenerate nucleotides
- **Mismatch tolerance**: Allow specified number of mismatches
- **Gapped motifs**: Discover patterns with variable-length gaps
- **Position Weight Matrix (PWM)**: Probabilistic sequence matching
- **Hybrid approaches**: PWM with gaps for complex patterns

### 2. Statistical Analysis Suite
- Enrichment ratio calculation (UP/DOWN regulation)
- Fisher's exact test for contingency tables
- Hypergeometric test for overrepresentation
- Binomial test for strand specificity
- Multiple testing correction (Bonferroni/FDR)

### 3. Advanced Features
- Parallel processing for large-scale analysis
- Transcript isoform selection tracking (MANE annotations)
- RNA secondary structure prediction (stem-loops, G-quadruplexes)
- Motif clustering based on sequence similarity
- Comprehensive visualization suite

## Installation

### System Requirements
- **Operating System**: Linux, macOS, or Windows
- **Python**: 3.7 or higher
- **Memory**: Minimum 8GB RAM (16GB recommended for large datasets)
- **CPU**: Multi-core processor recommended for parallel processing

### Dependencies
```bash
# Core dependencies
pip install pandas numpy biopython scipy matplotlib seaborn tqdm networkx

# Optional dependencies for RNA structure prediction
# Via pip:
pip install ViennaRNA

# Via conda (recommended):
conda install -c bioconda viennarna
```

### Clone Repository
```bash
git clone https://github.com/yangrf2021/5-UTR-Regulatory-Motifs-Discovery-Pipeline.git
```

## Quick Start

### 1. Prepare Input Files
Place two FASTA files in the working directory:
- `te_up_utr.fasta`: 5'UTR sequences of translationally upregulated genes
- `te_down_utr.fasta`: 5'UTR sequences of translationally downregulated genes

### 2. Run Analysis
```bash
python sequence-analysis.py
```

### 3. Check Results
Results will be generated in `motif_results_length_*` directories.

## Detailed Module Documentation

### Core Classes

#### `MotifAnalyzer`
Main class orchestrating the entire analysis pipeline.

**Key Attributes:**
- `up_sequences`: Dictionary storing upregulated gene sequences
- `down_sequences`: Dictionary storing downregulated gene sequences
- `max_utr_length`: Maximum 5'UTR length to analyze (default: 500nt)
- `min_frequency`: Minimum frequency threshold for initial filtering
- `fold_change_threshold`: Minimum enrichment ratio for significance
- `min_occur_freq`: Minimum occurrence frequency in gene set

**Key Methods:**
- `load_sequences()`: Load and preprocess FASTA files
- `analyze_single_length()`: Analyze motifs of specific length
- `find_motifs()`: Core motif discovery algorithm
- `process_enriched_motifs()`: Cluster and analyze significant motifs

### Sequence Processing Functions

#### `get_matched_sequences()`
Find all sequences matching a given pattern.

**Parameters:**
- `pattern`: Target pattern (supports IUPAC wildcards)
- `sequences`: Dictionary of sequences to search
- `max_mismatches`: Maximum allowed mismatches (default: 0)
- `gap_pattern`: Optional gap pattern format (e.g., "GAG{2,4}GG")

**Returns:**
- `matches`: List of matching sequences
- `match_positions`: List of (gene_id, position) tuples
- `mismatch_details`: Detailed mismatch information

#### `calculate_position_frequencies()`
Calculate nucleotide frequencies at each position of matched sequences.

**Purpose:** Build Position Weight Matrix (PWM) for motif representation

**Returns:** DataFrame with position-specific nucleotide frequencies

#### `scan_sequences_with_pwm()`
Scan sequences using a Position Weight Matrix for flexible matching.

**Parameters:**
- `pwm`: Position Weight Matrix
- `sequences`: Target sequences
- `threshold`: Matching threshold (0-1, default: 0.8)

### Statistical Analysis Functions

#### `calculate_strand_specificity()`
Determine if a motif shows strand preference (RNA vs DNA-level regulation).

**Interpretation:**
- Ratio ≈ 1: No strand preference (likely DNA-level)
- Ratio > 2: Strong forward strand preference (likely RNA-level)
- Ratio < 0.5: Reverse strand preference

**Returns:** (specificity_ratio, p_value)

#### `calculate_enrichment_statistics()`
Comprehensive statistical analysis of motif enrichment.

**Calculates:**
- Frequency in each gene set
- Fold change ratios
- Fisher's exact test p-value
- Odds ratio
- Hypergeometric test p-value

### Motif Discovery Functions

#### `process_sequence_chunk()`
Process a chunk of sequences for parallel motif discovery.

**Workflow:**
1. Extract all exact k-mers
2. Filter by minimum frequency
3. Generate wildcard variants
4. Test mismatch patterns
5. Explore gapped configurations
6. Apply PWM analysis
7. Combine PWM with gaps

#### `find_motifs()`
Main motif discovery orchestrator with parallel processing.

**Features:**
- Automatic CPU core detection
- Progress bar with tqdm
- Memory-efficient chunking
- Result aggregation across chunks

### Clustering and Analysis

#### `cluster_similar_motifs()`
Group similar motifs to reduce redundancy.

**Algorithm:**
- Similarity based on position-wise nucleotide matching
- Handles special formats (gaps, mismatches, PWM)
- Hierarchical clustering approach
- Customizable similarity threshold

#### `calculate_and_save_pattern_info()`
Generate comprehensive analysis report for each significant motif.

**Outputs:**
- Position frequency matrix
- Strand specificity analysis
- Gene list with motif
- Transcript type distribution
- Structural correlation (if enabled)
- Position distribution plots

### Visualization Functions

#### `create_visualizations()`
Generate comprehensive visualization suite.

**Plots Generated:**
1. **Frequency Distribution**: Histogram of motif frequencies
2. **Enrichment Scatter**: Frequency vs fold change
3. **Strand Specificity**: Distribution of strand bias
4. **Significant Heatmap**: Top motifs visualization
5. **Position Distribution**: Motif location in 5'UTRs
6. **Motif Type Analysis**: Distribution by discovery method
7. **Transcript Type Correlation**: MANE selection analysis

### Structural Motif Analysis

#### `find_structural_motifs()`
Identify potential RNA secondary structures.

**Detects:**
- Stem-loop structures
- Minimum stem length: 4nt
- Maximum loop size: 8nt
- Requires complementary base pairing

#### `find_g_quadruplexes()`
Detect potential G-quadruplex forming sequences.

**Pattern:** (G{3,5}[ACGT]{1,7}){3,}G{3,5}

**Biological Relevance:** G-quadruplexes can affect translation efficiency

## Algorithm Workflow

```
1. SEQUENCE LOADING
   ├── Parse FASTA files
   ├── Truncate to max UTR length
   └── Extract transcript metadata

2. MOTIF DISCOVERY (per length)
   ├── K-mer extraction
   ├── Frequency filtering
   ├── Pattern expansion
   │   ├── Wildcard generation
   │   ├── Mismatch variants
   │   ├── Gapped patterns
   │   └── PWM scoring
   └── Parallel processing

3. STATISTICAL ANALYSIS
   ├── Calculate frequencies
   ├── Compute enrichment ratios
   ├── Statistical testing
   └── Significance filtering

4. POST-PROCESSING
   ├── Motif clustering
   ├── Detailed analysis
   ├── Gene list generation
   └── Visualization

5. OUTPUT GENERATION
   ├── CSV files
   ├── Text reports
   ├── Plots
   └── Gene lists
```

## Input Requirements

### FASTA Format Specification
```fasta
>GENE_ID | Selection: TRANSCRIPT_TYPE | Additional_Info
ATGGCGCACGNNNNCCTGATCGATCGATCGATCGATCGATCGATCG...
```

**Header Components:**
- `GENE_ID`: Unique gene identifier (required)
- `Selection`: Transcript selection method (optional, e.g., "MANE_Select", "Longest_UTR")
- `Additional_Info`: Any additional metadata (optional)

### Sequence Requirements
- **Alphabet**: Standard nucleotides (A, C, G, T/U)
- **Special Characters**: N for unknown bases (will be excluded from analysis)
- **Length**: Sequences longer than `max_utr_length` will be truncated
- **Quality**: No quality scores needed (FASTA, not FASTQ)

## Output Description

### Directory Structure
```
motif_results_length_N/
├── 01motif_analysis.csv                    # All discovered motifs
├── 02significant_motifs.csv                # Filtered significant motifs
├── 02significant_motifs_with_details.csv   # Detailed gene information
├── 02significant_motifs_gene_lists.csv     # Gene lists per motif
├── 02up_enriched_motifs.csv               # UP-regulated enriched
├── 02down_enriched_motifs.csv             # DOWN-regulated enriched
├── 03up_clustered_motifs.csv              # Clustered UP motifs
├── 03down_clustered_motifs.csv            # Clustered DOWN motifs
├── 04_[direction]_clustered_motifs.csv    # Final clustered results
├── motif_types_summary.csv                # Statistics by motif type
├── motif_types_distribution.png           # Visualization of types
├── 01frequency_distribution.png           # Frequency histograms
├── 02fold_change_distribution.png         # Enrichment distribution
├── 03enrichment_scatter.png               # Scatter plot analysis
├── 04frequency_heatmap.png                # Heatmap of top motifs
├── 05strand_specificity_by_direction.png  # Strand bias analysis
├── 06motif_type_significance.png          # Significance by type
├── [direction]_motif_[PATTERN]_details.txt     # Detailed report
├── [direction]_motif_[PATTERN]_frequencies.csv # Position frequencies
├── [direction]_motif_[PATTERN]_genes.txt       # Gene list
└── transcript_type_analysis/              # MANE selection analysis
```

### Key Output Files

#### `01motif_analysis.csv`
Complete list of all discovered motifs with full statistics.

**Columns:**
- `Pattern`: Motif sequence or pattern
- `Length`: Motif length
- `Type`: Discovery method (exact/mismatch/gapped/pwm)
- `UP_present`: Count in upregulated genes
- `DOWN_present`: Count in downregulated genes
- `UP_frequency`: Frequency in UP set
- `DOWN_frequency`: Frequency in DOWN set
- `UP_to_DOWN_ratio`: Enrichment ratio
- `P_value`: Fisher's exact test p-value
- `Odds_ratio`: Statistical odds ratio
- `Strand_specificity`: Forward/reverse strand ratio
- `Strand_pvalue`: Strand bias significance
- `Significant`: Boolean significance flag
- `Enrichment`: Direction (UP/DOWN/None)

#### `02significant_motifs.csv`
Filtered motifs meeting significance criteria.

**Significance Criteria:**
- Frequency ≥ `min_occur_freq`
- Fold change ≥ `fold_change_threshold`
- Present in multiple genes

#### Gene Lists (`*_genes.txt`)
Plain text files with one gene ID per line, suitable for:
- Gene Ontology (GO) enrichment analysis
- KEGG pathway analysis
- Integration with other tools

## Parameter Configuration

### Core Parameters
```python
# In main() function of sequence-analysis.py

# Sequence processing
max_utr_length = 500        # Maximum 5'UTR length in nucleotides
                           # Longer sequences will be truncated

# Frequency thresholds
min_frequency = 0.05        # Initial filtering (5% of sequences)
                           # Used during motif discovery phase

min_occur_freq = 0.1       # Final significance threshold (10%)
                           # Used for determining significant motifs

# Enrichment criteria
fold_change_threshold = 1.5 # Minimum enrichment ratio
                           # UP/DOWN or DOWN/UP ratio

# Motif discovery options
max_wildcards = 1          # Maximum IUPAC wildcards per motif
                          # Higher values increase search space

max_mismatches = 1        # Maximum allowed mismatches
                         # For approximate matching

# Advanced features
enable_gapped = True      # Enable gapped motif discovery
                         # Finds patterns like GCACG{2,4}CC

enable_pwm = True        # Enable Position Weight Matrix
                        # For probabilistic matching

enable_structural = False # Enable RNA structure analysis
                        # Requires ViennaRNA package

# Performance
n_cores = None           # Number of CPU cores (None = auto-detect)
```

### Length Range
```python
# Analyze motifs from 4 to 8 nucleotides
for length in range(4, 9):
    analyzer.analyze_single_length(length)
```

## Example Usage

### Basic Usage
```python
# Default parameters
analyzer = MotifAnalyzer()

# Analyze all lengths
for length in range(4, 9):
    results = analyzer.analyze_single_length(length)
```

### Custom Parameters
```python
# More stringent analysis
analyzer = MotifAnalyzer(
    max_utr_length=1000,      # Analyze longer UTRs
    min_frequency=0.02,       # More sensitive
    fold_change_threshold=2.0, # Stricter enrichment
    min_occur_freq=0.15,      # Higher significance threshold
    max_mismatches=2,         # Allow more mismatches
    enable_gapped=True,
    enable_pwm=True,
    enable_structural=True    # Include structure analysis
)
```

### Parallel Processing
```python
# Specify CPU cores for large datasets
results = analyzer.analyze_single_length(
    length=7,
    n_cores=16,  # Use 16 cores
    max_wildcards=2  # More wildcards for deeper search
)
```

### Accessing Results Programmatically
```python
# Load and filter results
import pandas as pd

# Read significant motifs
sig_motifs = pd.read_csv('motif_results_length_7/02significant_motifs.csv')

# Filter for highly enriched motifs
highly_enriched = sig_motifs[sig_motifs['UP_to_DOWN_ratio'] > 3]

# Get DOWN-enriched motifs
down_enriched = sig_motifs[sig_motifs['Enrichment'] == 'DOWN']

# Extract gene lists
up_genes = set()
for _, row in sig_motifs.iterrows():
    if row['Enrichment'] == 'UP':
        genes = row['UP_genes'].split(',')
        up_genes.update(genes)
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Memory Error
**Problem:** `MemoryError: Unable to allocate array`

**Solutions:**
- Reduce `max_utr_length` parameter
- Process lengths sequentially rather than in parallel
- Use fewer CPU cores to reduce memory overhead
- Split input files into smaller batches

#### 2. No Significant Motifs Found
**Problem:** No motifs pass significance thresholds

**Solutions:**
- Decrease `min_frequency` (e.g., 0.02)
- Decrease `fold_change_threshold` (e.g., 1.2)
- Decrease `min_occur_freq` (e.g., 0.05)
- Enable more discovery methods (`enable_gapped`, `enable_pwm`)
- Check input file quality and sequence diversity

#### 3. Slow Performance
**Problem:** Analysis takes too long

**Solutions:**
- Increase `n_cores` parameter
- Reduce motif length range
- Disable structural analysis if not needed
- Use SSD for faster I/O
- Consider cloud computing for large datasets

#### 4. RNA Package Import Error
**Problem:** `ImportError: No module named 'RNA'`

**Solutions:**
```bash
# Install ViennaRNA
conda install -c bioconda viennarna
# or
pip install ViennaRNA

# If still not working, disable structural analysis:
enable_structural = False
```

#### 5. Empty Output Files
**Problem:** Output files are created but empty

**Solutions:**
- Check input FASTA format
- Verify sequence quality (no invalid characters)
- Ensure sequences are DNA (ACGT), not protein
- Check file permissions in output directory

## Performance Considerations

### Computational Complexity
- **Time Complexity:** O(n × m × k) where:
  - n = number of sequences
  - m = average sequence length
  - k = motif length
- **Space Complexity:** O(n × 4^k) for storing all k-mers

### Optimization Tips

#### For Large Datasets (>10,000 sequences)
```python
# Use aggressive filtering
min_frequency = 0.1  # Higher threshold
max_wildcards = 0    # No wildcards initially
enable_structural = False  # Disable structure

# Process in batches
chunk_size = 1000
for i in range(0, len(sequences), chunk_size):
    chunk = sequences[i:i+chunk_size]
    # Process chunk
```

#### For Deep Analysis (High Sensitivity)
```python
# Comprehensive search
min_frequency = 0.01  # Very low threshold
max_wildcards = 2     # More wildcards
max_mismatches = 2    # More mismatches
enable_gapped = True
enable_pwm = True
enable_structural = True
```

### Memory Requirements

| Dataset Size | Memory Usage | Recommended RAM |
|-------------|--------------|-----------------|
| 100 genes | ~500 MB | 4 GB |
| 1,000 genes | ~2 GB | 8 GB |
| 10,000 genes | ~8 GB | 16 GB |
| 50,000 genes | ~32 GB | 64 GB |

### Processing Time Estimates

| Dataset | Length Range | CPU Cores | Estimated Time |
|---------|-------------|-----------|----------------|
| 1,000 genes | 4-8 bp | 4 | 5-10 minutes |
| 10,000 genes | 4-8 bp | 8 | 30-60 minutes |
| 50,000 genes | 4-8 bp | 16 | 2-4 hours |

## Advanced Features

### Custom Motif Patterns

#### IUPAC Wildcards
```python
# Supported IUPAC codes
wildcards = {
    'R': ['A', 'G'],      # puRine
    'Y': ['C', 'T'],      # pYrimidine
    'M': ['A', 'C'],      # aMino
    'K': ['G', 'T'],      # Keto
    'S': ['G', 'C'],      # Strong
    'W': ['A', 'T'],      # Weak
    'H': ['A', 'C', 'T'], # not G
    'B': ['C', 'G', 'T'], # not A
    'V': ['A', 'C', 'G'], # not T
    'D': ['A', 'G', 'T'], # not C
    'N': ['A', 'C', 'G', 'T'] # aNy
}
```

#### Gap Patterns
```python
# Format: PREFIX{min,max}SUFFIX
# Examples:
"GCACG{2,4}CC"    # 2-4 nucleotides between GCACG and CC
"ATG{1,3}TAA"     # 1-3 nucleotides between ATG and TAA
"GGG{3,5}CCC"     # 3-5 nucleotides between GGG and CCC
```

#### PWM Patterns
```python
# Format: PWM:BASE_PATTERN
# Example: "PWM:GCACGCC"
# Creates position weight matrix from GCACGCC matches
# Then uses probabilistic scoring for flexible matching
```

### Integration with Other Tools

#### Export for MEME Suite
```python
# Convert to MEME format
def to_meme_format(motifs_df, output_file):
    with open(output_file, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        # Add motif data...
```

#### Export for HOMER
```python
# Create HOMER motif file
def to_homer_format(motifs_df, output_file):
    with open(output_file, 'w') as f:
        for _, motif in motifs_df.iterrows():
            f.write(f">{motif['Pattern']}\t{motif['P_value']}\n")
            # Add frequency matrix...
```

## Validation Methods

### Experimental Validation Suggestions

1. **Luciferase Reporter Assays**
   - Clone 5'UTRs with/without motifs
   - Measure translation efficiency

2. **Ribosome Profiling**
   - Validate translation changes
   - Confirm ribosome occupancy

3. **CLIP-seq**
   - Identify RNA-binding proteins
   - Validate protein-motif interactions

4. **Mutation Analysis**
   - Mutate key motif positions
   - Measure functional impact

## Future Developments

### Planned Features
- [ ] Machine learning-based motif scoring
- [ ] Integration with RNA-seq data
- [ ] Web interface for easy access
- [ ] Docker containerization
- [ ] Cloud deployment options
- [ ] Real-time analysis updates
- [ ] Interactive visualization dashboard

### Contributing
We welcome contributions! Please submit pull requests or open issues for:
- Bug fixes
- New features
- Documentation improvements
- Performance optimizations

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{yang2025mettl5,
  title={METTL5-mediated 18S rRNA methylation controls mRNA translation through 5'UTR regulatory motifs},
  author={Yang, Ruifeng and others},
  journal={Journal Name},
  year={2025},
  note={Manuscript in preparation}
}
```

## Contact

**Developer:** Ruifeng Yang  
**Email:** yangrf1996@163.com  
**Institution:** Department of Urology, Fudan University Shanghai Cancer Center, Fudan University, Shanghai, China.; The Core Laboratory in Medical Center of Clinical Research, Shanghai Ninth People’s Hospital, State Key Laboratory of Medical Genomics, Shanghai Jiao Tong University School of Medicine, Shanghai, China.  

For questions, bug reports, or collaborations, please:
1. Open an issue on GitHub
2. Email the developer directly
3. Check the FAQ section

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- National Natural Science Foundation of China (Grant No. 82473192, 82172741)
- Fudan University Shanghai Cancer Center
- Shanghai Ninth People's Hospital
- All contributors and collaborators
- Open-source community for tools and libraries

---

**Last Updated:** January 2025  
**Version:** 1.0.0  
**Status:** Active Development
