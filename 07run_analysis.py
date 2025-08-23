#!/usr/bin/env python
"""
Main script to run the complete motif analysis pipeline
"""

import os
import sys
import argparse
import logging
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.motif_analyzer import MotifAnalyzer


def set_random_seeds(seed=42):
    """Set random seeds for reproducibility"""
    import random
    import numpy as np

    random.seed(seed)
    np.random.seed(seed)


def main():
    """Main analysis function"""

    parser = argparse.ArgumentParser(
        description='METTL5-mediated Translation Control Motif Discovery'
    )

    parser.add_argument(
        '--max-utr-length', type=int, default=500,
        help='Maximum 5\'UTR length (default: 500)'
    )
    parser.add_argument(
        '--min-frequency', type=float, default=0.05,
        help='Minimum frequency threshold (default: 0.05)'
    )
    parser.add_argument(
        '--fold-change', type=float, default=1.5,
        help='Fold change threshold (default: 1.5)'
    )
    parser.add_argument(
        '--min-occur-freq', type=float, default=0.1,
        help='Minimum occurrence frequency (default: 0.1)'
    )
    parser.add_argument(
        '--max-wildcards', type=int, default=1,
        help='Maximum wildcards in motifs (default: 1)'
    )
    parser.add_argument(
        '--max-mismatches', type=int, default=1,
        help='Maximum mismatches allowed (default: 1)'
    )
    parser.add_argument(
        '--enable-gapped', action='store_true',
        help='Enable gapped motif discovery'
    )
    parser.add_argument(
        '--enable-pwm', action='store_true',
        help='Enable PWM analysis'
    )
    parser.add_argument(
        '--enable-structural', action='store_true',
        help='Enable structural motif analysis'
    )
    parser.add_argument(
        '--min-length', type=int, default=4,
        help='Minimum motif length (default: 4)'
    )
    parser.add_argument(
        '--max-length', type=int, default=8,
        help='Maximum motif length (default: 8)'
    )
    parser.add_argument(
        '--n-cores', type=int, default=None,
        help='Number of CPU cores to use (default: auto)'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducibility (default: 42)'
    )

    args = parser.parse_args()

    # Set random seeds
    set_random_seeds(args.seed)

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)

    logger.info("Starting METTL5 motif analysis pipeline")
    logger.info(f"Parameters: {vars(args)}")

    # Initialize analyzer
    analyzer = MotifAnalyzer(
        max_utr_length=args.max_utr_length,
        min_frequency=args.min_frequency,
        fold_change_threshold=args.fold_change,
        min_occur_freq=args.min_occur_freq,
        max_mismatches=args.max_mismatches,
        enable_gapped=args.enable_gapped,
        enable_pwm=args.enable_pwm,
        enable_structural=args.enable_structural,
        enable_pwm_gapped=args.enable_gapped and args.enable_pwm
    )

    # Create results directory
    os.makedirs("motif_analysis_results", exist_ok=True)

    # Run analysis for each length
    summary_data = []

    for length in range(args.min_length, args.max_length + 1):
        logger.info(f"\nAnalyzing motifs of length {length}")

        results = analyzer.analyze_single_length(
            length,
            n_cores=args.n_cores,
            max_wildcards=args.max_wildcards
        )

        if not results.empty:
            significant = results[results['Significant']]
            up_enriched = significant[significant['Enrichment'] == 'UP']
            down_enriched = significant[significant['Enrichment'] == 'DOWN']

            summary_data.append({
                'Length': length,
                'Total_motifs': len(results),
                'Significant_motifs': len(significant),
                'UP_enriched': len(up_enriched),
                'DOWN_enriched': len(down_enriched),
                'Max_UP_enrichment': results['UP_to_DOWN_ratio'].replace(
                    [float('inf'), float('-inf')], float('nan')
                ).max(),
                'Max_DOWN_enrichment': results['DOWN_to_UP_ratio'].replace(
                    [float('inf'), float('-inf')], float('nan')
                ).max()
            })

            logger.info(f"Found {len(significant)} significant motifs")
            logger.info(f"  UP-enriched: {len(up_enriched)}")
            logger.info(f"  DOWN-enriched: {len(down_enriched)}")

    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(
            "motif_analysis_results/analysis_summary.csv",
            index=False
        )

        # Create summary plots
        create_summary_plots(summary_df)

        # Write README
        write_analysis_readme(args, summary_df)

    logger.info("\nAnalysis complete!")
    logger.info("Results saved in motif_analysis_results/")


def create_summary_plots(summary_df):
    """Create summary visualizations"""

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Total vs Significant motifs
    axes[0, 0].bar(
        summary_df['Length'] - 0.2,
        summary_df['Total_motifs'],
        width=0.4,
        label='Total',
        color='gray',
        alpha=0.7
    )
    axes[0, 0].bar(
        summary_df['Length'] + 0.2,
        summary_df['Significant_motifs'],
        width=0.4,
        label='Significant',
        color='green',
        alpha=0.7
    )
    axes[0, 0].set_xlabel('Motif Length')
    axes[0, 0].set_ylabel('Number of Motifs')
    axes[0, 0].set_title('Motif Discovery Summary')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: UP vs DOWN enriched
    axes[0, 1].plot(
        summary_df['Length'],
        summary_df['UP_enriched'],
        marker='o',
        color='red',
        label='UP-enriched'
    )
    axes[0, 1].plot(
        summary_df['Length'],
        summary_df['DOWN_enriched'],
        marker='s',
        color='blue',
        label='DOWN-enriched'
    )
    axes[0, 1].set_xlabel('Motif Length')
    axes[0, 1].set_ylabel('Number of Enriched Motifs')
    axes[0, 1].set_title('Enrichment Direction by Length')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Maximum enrichment values
    axes[1, 0].plot(
        summary_df['Length'],
        summary_df['Max_UP_enrichment'],
        marker='o',
        color='red',
        label='Max UP enrichment'
    )
    axes[1, 0].plot(
        summary_df['Length'],
        summary_df['Max_DOWN_enrichment'],
        marker='s',
        color='blue',
        label='Max DOWN enrichment'
    )
    axes[1, 0].set_xlabel('Motif Length')
    axes[1, 0].set_ylabel('Maximum Enrichment Ratio')
    axes[1, 0].set_title('Peak Enrichment by Length')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Proportion significant
    axes[1, 1].bar(
        summary_df['Length'],
        summary_df['Significant_motifs'] / summary_df['Total_motifs'],
        color='purple',
        alpha=0.7
    )
    axes[1, 1].set_xlabel('Motif Length')
    axes[1, 1].set_ylabel('Proportion Significant')
    axes[1, 1].set_title('Discovery Efficiency')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('motif_analysis_results/analysis_summary.png', dpi=300)
    plt.close()


def write_analysis_readme(args, summary_df):
    """Write README file with analysis details"""

    with open('motif_analysis_results/README.md', 'w') as f:
        f.write("# Motif Analysis Results\n\n")
        f.write("## Analysis Parameters\n\n")

        f.write(f"- Maximum UTR length: {args.max_utr_length} nt\n")
        f.write(f"- Minimum frequency threshold: {args.min_frequency}\n")
        f.write(f"- Fold change threshold: {args.fold_change}\n")
        f.write(f"- Minimum occurrence frequency: {args.min_occur_freq}\n")
        f.write(f"- Maximum wildcards: {args.max_wildcards}\n")
        f.write(f"- Maximum mismatches: {args.max_mismatches}\n")
        f.write(f"- Gapped motifs: {'Enabled' if args.enable_gapped else 'Disabled'}\n")
        f.write(f"- PWM analysis: {'Enabled' if args.enable_pwm else 'Disabled'}\n")
        f.write(f"- Structural analysis: {'Enabled' if args.enable_structural else 'Disabled'}\n\n")

        f.write("## Summary Statistics\n\n")
        f.write(summary_df.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Key Findings\n\n")

        # Find best length
        best_length = summary_df.loc[
            summary_df['Significant_motifs'].idxmax(), 'Length'
        ]
        f.write(f"- Optimal motif length: {best_length} bp\n")

        total_sig = summary_df['Significant_motifs'].sum()
        total_up = summary_df['UP_enriched'].sum()
        total_down = summary_df['DOWN_enriched'].sum()

        f.write(f"- Total significant motifs discovered: {total_sig}\n")
        f.write(f"  - UP-enriched: {total_up} ({total_up / total_sig * 100:.1f}%)\n")
        f.write(f"  - DOWN-enriched: {total_down} ({total_down / total_sig * 100:.1f}%)\n\n")

        f.write("## Directory Structure\n\n")
        f.write("```\n")
        f.write("motif_results_length_N/\n")
        f.write("├── 01_motif_analysis.csv         # All discovered motifs\n")
        f.write("├── 02_significant_motifs.csv     # Filtered significant motifs\n")
        f.write("├── 03_up_enriched_motifs.csv     # UP-regulated enriched\n")
        f.write("├── 03_down_enriched_motifs.csv   # DOWN-regulated enriched\n")
        f.write("├── 04_clustered_motifs.csv       # Clustered similar motifs\n")
        f.write("├── visualizations/               # Analysis plots\n")
        f.write("└── [motif]_details.txt           # Detailed motif analysis\n")
        f.write("```\n\n")

        f.write("## Next Steps\n\n")
        f.write("1. Review significant motifs in each length directory\n")
        f.write("2. Perform GO/KEGG enrichment on gene lists (*_genes.txt)\n")
        f.write("3. Validate top motifs experimentally\n")
        f.write("4. Investigate motif-containing genes for biological relevance\n")


if __name__ == "__main__":
    main()