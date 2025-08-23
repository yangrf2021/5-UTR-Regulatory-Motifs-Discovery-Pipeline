# 5-UTR-Regulatory-Motifs-Discovery-Pipeline

**Overview**
A computational framework for identifying 5'UTR regulatory motifs in translationally regulated genes, particularly those affected by METTL5-mediated rRNA methylation in prostate cancer. This pipeline integrates multiple motif discovery algorithms with advanced statistical analysis to identify sequence elements that control translation efficiency.

**Background**
METTL5 is an RNA methyltransferase that catalyzes m6A modification at position 1832 of 18S rRNA. Our research demonstrates that METTL5-mediated rRNA methylation selectively regulates the translation of specific mRNAs containing the motif GCACGN{2,4}CC in their 5'UTR regions. This regulatory mechanism plays a crucial role in:
Prostate cancer progression
Mitochondrial oxidative phosphorylation
Translation efficiency control through the METTL5/IRF7/DNA2 axis

**Architecture**
src/
├── __init__.py                 # Package initialization
├── motif_analyzer.py           # Main analyzer class
├── sequence_utils.py           # Sequence manipulation utilities
├── motif_discovery.py          # Core motif discovery algorithms
├── statistics.py               # Statistical analysis functions
├── visualization.py            # Visualization manager
└── run_analysis.py            # Main execution script

**Installation**

