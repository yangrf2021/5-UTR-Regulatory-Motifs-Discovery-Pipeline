# 5-UTR-Regulatory-Motifs-Discovery-Pipeline

📋 **Overview**
A computational framework for identifying 5'UTR regulatory motifs in translationally regulated genes, particularly those affected by METTL5-mediated rRNA methylation in prostate cancer. This pipeline integrates multiple motif discovery algorithms with advanced statistical analysis to identify sequence elements that control translation efficiency.

🔬 **Background**
METTL5 is an RNA methyltransferase that catalyzes m6A modification at position 1832 of 18S rRNA. Our research demonstrates that METTL5-mediated rRNA methylation selectively regulates the translation of specific mRNAs containing the motif GCACGN{2,4}CC in their 5'UTR regions. This regulatory mechanism plays a crucial role in:
Prostate cancer progression
Mitochondrial oxidative phosphorylation
Translation efficiency control through the METTL5/IRF7/DNA2 axis

🏗️ **Architecture**
src/
├── __init__.py                 # Package initialization
├── motif_analyzer.py           # Main analyzer class
├── sequence_utils.py           # Sequence manipulation utilities
├── motif_discovery.py          # Core motif discovery algorithms
├── statistics.py               # Statistical analysis functions
├── visualization.py            # Visualization manager
└── run_analysis.py            # Main execution script

🚀 **Installation**
Prerequisites：Python 3.7+
CUDA-capable GPU (optional, for faster parallel processing)

Dependencies：pip install -r requirements.txt

📊 **Input Data Format**
Required Files：
1. te_up_utr.fasta - 5'UTR sequences of translationally upregulated genes
2. te_down_utr.fasta - 5'UTR sequences of translationally downregulated genes

FASTA Format Example：>GENE_ID | Selection: Longest_UTR | Additional_Info
ATGGCGCACGNNNCCTGATCG...

🔧 **Usage**
Basic Usage：python src/run_analysis.py

Advanced Usage with Custom Parameters：python src/run_analysis.py \
    --max-utr-length 500 \ # Maximum 5'UTR length to analyze
    --min-frequency 0.05 \ # Minimum frequency threshold for motif selection
    --fold-change 1.5 \ # Fold change threshold for enrichment
    --min-occur-freq 0.1 \ # Minimum occurrence frequency for significance
    --max-wildcards 1 \ # Maximum allowed wildcards in motifs
    --max-mismatches 1 \ # Maximum allowed mismatches
    --enable-gapped \ # Enable gapped motif discovery
    --enable-pwm \ # Enable Position Weight Matrix analysis
    --enable-structural \ # Enable structural motif analysis
    --min-length 4 \ # Minimum motif length
    --max-length 8 \ # Maximum motif length
    --n-cores 8 \ # Number of CPU cores to use
    --seed 42

📁 **Output Structure**
motif_analysis_results/
├── analysis_summary.csv           # Overall summary statistics
├── analysis_summary.png           # Summary visualizations
├── README.md                      # Auto-generated analysis report
└── motif_results_length_N/       # Results for each motif length
    ├── 01_motif_analysis.csv     # All discovered motifs
    ├── 02_significant_motifs.csv # Filtered significant motifs
    ├── 03_up_enriched_motifs.csv # UP-regulated enriched motifs
    ├── 03_down_enriched_motifs.csv # DOWN-regulated enriched motifs
    ├── 04_clustered_motifs.csv   # Clustered similar motifs
    ├── visualizations/            # Analysis plots
    │   ├── frequency_distribution.png
    │   ├── enrichment_scatter.png
    │   ├── strand_specificity.png
    │   ├── significant_heatmap.png
    │   └── motif_network.png
    └── [motif]_details.txt       # Detailed motif analysis

🔍 **Key Features**
1. Motif Discovery Algorithms

Exact matching: Traditional k-mer analysis
IUPAC wildcards: Support for degenerate nucleotides
Gapped motifs: Discovery of motifs with variable gaps (e.g., GCACGN{2,4}CC)
PWM-based: Position Weight Matrix scoring
Mismatch tolerance: Allow specified number of mismatches

2. Statistical Analysis

Enrichment calculation: UP/DOWN regulation ratios
Strand specificity: Detection of RNA vs DNA-level regulation
Fisher's exact test: Statistical significance testing
Hypergeometric test: Alternative statistical validation

3. Advanced Features

Parallel processing: Multi-core support for large datasets
Motif clustering: Group similar motifs together
Structural analysis: G-quadruplex and stem-loop detection
Comprehensive visualization: Multiple plot types for result interpretation

📈 **Example Results**
Discovered Motif: GCACGN{2,4}CC

Enrichment: 2.8-fold in translationally downregulated genes
Frequency: Present in 15% of target genes
Strand specificity: 3.2 (indicates RNA-level regulation)
Clinical relevance: Associated with poor prognosis in prostate cancer

🧬 **Biological Validation**
The computational predictions from this pipeline have been validated through:

Dual-luciferase reporter assays
Polysome profiling
CLIP-seq experiments
Clinical sample analysis

📚 **Citation**
If you use this pipeline in your research, please cite:
bibtex@article{yang2025mettl5,
  title={METTL5-mediated 18S rRNA methylation controls translation of specific mRNAs via 5'UTR motifs in prostate cancer},
  author={Yang, Ruifeng and Others},
  journal={Nature Metabolism},
  year={2025},
  doi={10.1038/xxxxx}
}

🤝 **Contributing**
We welcome contributions! Please see CONTRIBUTING.md for details.

📧 **Contact**
Author: Ruifeng Yang
Email: yangrf1996@163.com
Institution: [Your Institution]

📄 **License**
This project is licensed under the MIT License - see the LICENSE file for details.

🙏 **Acknowledgments
**
Thanks to all contributors and collaborators
Supported by [Grant Numbers]
Computational resources provided by [Institution]

⚠️ **Disclaimer**
This software is for research purposes only. Clinical applications require additional validation.

Last Updated: January 2025
Version: 1.0.0

