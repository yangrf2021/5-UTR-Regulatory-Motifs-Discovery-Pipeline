"""
Main analyzer class for motif discovery and analysis
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

from .sequence_utils import (
    get_matched_sequences,
    get_matched_sequences_with_pwm_gaps,
    calculate_position_frequencies,
    scan_sequences_with_pwm
)
from .motif_discovery import (
    find_structural_motifs,
    find_g_quadruplexes,
    cluster_similar_motifs
)
from .statistics import calculate_strand_specificity
from .visualization import VisualizationManager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MotifAnalyzer:
    """
    Main class for analyzing regulatory motifs in 5'UTR sequences

    This analyzer identifies sequence motifs that are enriched in translationally
    regulated genes, particularly those affected by METTL5-mediated rRNA methylation.

    Attributes:
        max_utr_length (int): Maximum 5'UTR length to analyze
        min_frequency (float): Minimum frequency threshold for motif selection
        fold_change_threshold (float): Fold change threshold for enrichment
        min_occur_freq (float): Minimum occurrence frequency for significance
        max_mismatches (int): Maximum allowed mismatches in motif matching
        enable_gapped (bool): Enable gapped motif discovery
        enable_pwm (bool): Enable Position Weight Matrix analysis
        enable_structural (bool): Enable structural motif analysis
    """

    def __init__(
            self,
            max_utr_length: int = 500,
            min_frequency: float = 0.05,
            fold_change_threshold: float = 1.5,
            min_occur_freq: float = 0.1,
            max_mismatches: int = 0,
            enable_gapped: bool = False,
            enable_pwm: bool = False,
            enable_structural: bool = False,
            enable_pwm_gapped: bool = False
    ):
        """Initialize the motif analyzer with specified parameters"""

        self.max_utr_length = max_utr_length
        self.min_frequency = min_frequency
        self.fold_change_threshold = fold_change_threshold
        self.min_occur_freq = min_occur_freq
        self.max_mismatches = max_mismatches
        self.enable_gapped = enable_gapped
        self.enable_pwm = enable_pwm
        self.enable_structural = enable_structural
        self.enable_pwm_gapped = enable_pwm_gapped

        # Initialize sequence storage
        self.up_sequences = {}
        self.down_sequences = {}
        self.gene_ids = {}

        # Transcript selection information
        self.up_transcript_selection = {}
        self.down_transcript_selection = {}

        # Structural motifs storage
        self.up_structural_motifs = []
        self.down_structural_motifs = []

        # Visualization manager
        self.viz_manager = VisualizationManager()

        # Load sequences
        self._load_sequences()

        # Analyze structural motifs if enabled
        if self.enable_structural:
            self._analyze_structural_motifs()

    def _load_sequences(self):
        """Load sequences from FASTA files"""

        try:
            # Load upregulated sequences
            logger.info("Loading up-regulated sequences...")
            if os.path.exists("te_up_utr.fasta"):
                self._load_fasta("te_up_utr.fasta", "up")
            else:
                logger.error("File te_up_utr.fasta not found")

            # Load downregulated sequences
            logger.info("Loading down-regulated sequences...")
            if os.path.exists("te_down_utr.fasta"):
                self._load_fasta("te_down_utr.fasta", "down")
            else:
                logger.error("File te_down_utr.fasta not found")

            logger.info(f"Loaded {len(self.up_sequences)} up-regulated sequences")
            logger.info(f"Loaded {len(self.down_sequences)} down-regulated sequences")

            # Summarize transcript selection
            self._summarize_transcript_selection()

        except Exception as e:
            logger.error(f"Error loading sequences: {str(e)}")
            raise

    def _load_fasta(self, filename: str, regulation_type: str):
        """Load sequences from a single FASTA file"""

        sequences = self.up_sequences if regulation_type == "up" else self.down_sequences
        transcript_selection = (
            self.up_transcript_selection if regulation_type == "up"
            else self.down_transcript_selection
        )

        for record in SeqIO.parse(filename, "fasta"):
            gene_id = record.id
            seq = str(record.seq)[:self.max_utr_length]
            sequences[gene_id] = seq
            self.gene_ids[gene_id] = gene_id

            # Extract transcript selection logic if present
            if "Selection:" in record.description:
                selection_logic = record.description.split("Selection:")[1].split("|")[0].strip()
                transcript_selection[gene_id] = selection_logic

    def _summarize_transcript_selection(self):
        """Summarize and visualize transcript selection distribution"""

        if not self.up_transcript_selection and not self.down_transcript_selection:
            logger.warning("No transcript selection information found")
            return

        # Count selection types
        up_counts = defaultdict(int)
        for selection in self.up_transcript_selection.values():
            up_counts[selection] += 1

        down_counts = defaultdict(int)
        for selection in self.down_transcript_selection.values():
            down_counts[selection] += 1

        # Log summary
        logger.info("Transcript selection summary for UP-regulated genes:")
        for selection, count in sorted(up_counts.items()):
            percentage = count / len(self.up_transcript_selection) * 100 if self.up_transcript_selection else 0
            logger.info(f"  {selection}: {count} ({percentage:.1f}%)")

        logger.info("Transcript selection summary for DOWN-regulated genes:")
        for selection, count in sorted(down_counts.items()):
            percentage = count / len(self.down_transcript_selection) * 100 if self.down_transcript_selection else 0
            logger.info(f"  {selection}: {count} ({percentage:.1f}%)")

        # Save summary
        os.makedirs("transcript_selection", exist_ok=True)

        selection_data = []
        for gene_id, selection in self.up_transcript_selection.items():
            selection_data.append({
                "Gene": gene_id,
                "Selection": selection,
                "Regulation": "UP"
            })

        for gene_id, selection in self.down_transcript_selection.items():
            selection_data.append({
                "Gene": gene_id,
                "Selection": selection,
                "Regulation": "DOWN"
            })

        if selection_data:
            pd.DataFrame(selection_data).to_csv(
                "transcript_selection/transcript_selection_summary.csv",
                index=False
            )

    def _analyze_structural_motifs(self):
        """Analyze structural motifs in sequences"""

        logger.info("Analyzing structural motifs...")

        # Analyze upregulated sequences
        if self.up_sequences:
            stem_loops, sl_positions = find_structural_motifs(self.up_sequences)
            g_quads, gq_positions = find_g_quadruplexes(self.up_sequences)

            self.up_structural_motifs = {
                'stem_loops': {'motifs': stem_loops, 'positions': sl_positions},
                'g_quadruplexes': {'motifs': g_quads, 'positions': gq_positions}
            }

            logger.info(f"Found {len(stem_loops)} stem-loops in UP-regulated genes")
            logger.info(f"Found {len(g_quads)} G-quadruplexes in UP-regulated genes")

        # Analyze downregulated sequences
        if self.down_sequences:
            stem_loops, sl_positions = find_structural_motifs(self.down_sequences)
            g_quads, gq_positions = find_g_quadruplexes(self.down_sequences)

            self.down_structural_motifs = {
                'stem_loops': {'motifs': stem_loops, 'positions': sl_positions},
                'g_quadruplexes': {'motifs': g_quads, 'positions': gq_positions}
            }

            logger.info(f"Found {len(stem_loops)} stem-loops in DOWN-regulated genes")
            logger.info(f"Found {len(g_quads)} G-quadruplexes in DOWN-regulated genes")

    def analyze_single_length(
            self,
            length: int,
            n_cores: Optional[int] = None,
            max_wildcards: int = 1
    ) -> pd.DataFrame:
        """
        Analyze motifs of a specific length

        Args:
            length: Motif length to analyze
            n_cores: Number of CPU cores to use
            max_wildcards: Maximum number of wildcard positions

        Returns:
            DataFrame containing analysis results
        """

        logger.info(f"Analyzing motifs of length {length}")

        # Create result directory
        result_dir = f'motif_results_length_{length}'
        os.makedirs(result_dir, exist_ok=True)

        # Find motifs in both sequence sets
        from .motif_discovery import find_motifs

        up_patterns, up_counts, up_positions = find_motifs(
            self.up_sequences, length, self.min_frequency,
            n_cores, max_wildcards, self.max_mismatches,
            self.enable_gapped, self.enable_pwm, self.enable_pwm_gapped
        )

        down_patterns, down_counts, down_positions = find_motifs(
            self.down_sequences, length, self.min_frequency,
            n_cores, max_wildcards, self.max_mismatches,
            self.enable_gapped, self.enable_pwm, self.enable_pwm_gapped
        )

        # Analyze all patterns
        results = self._analyze_patterns(
            up_patterns, down_patterns, up_positions, down_positions
        )

        # Save results
        if results:
            results_df = pd.DataFrame(results)
            results_df = self._determine_significance(results_df)

            # Save all results
            results_df.to_csv(f'{result_dir}/01_motif_analysis.csv', index=False)

            # Save significant results
            significant_df = results_df[results_df['Significant']]
            if not significant_df.empty:
                significant_df.to_csv(f'{result_dir}/02_significant_motifs.csv', index=False)

                # Process enriched motifs
                self._process_enriched_motifs(
                    significant_df, up_positions, down_positions, result_dir
                )

            # Create visualizations
            self.viz_manager.create_analysis_visualizations(
                results_df, significant_df, length, result_dir
            )

            return results_df

        return pd.DataFrame()

    def _analyze_patterns(
            self,
            up_patterns: Dict,
            down_patterns: Dict,
            up_positions: Dict,
            down_positions: Dict
    ) -> List[Dict]:
        """Analyze patterns for enrichment and significance"""

        results = []
        all_patterns = set(up_patterns) | set(down_patterns)

        for pattern in all_patterns:
            up_present = len(up_patterns.get(pattern, set()))
            down_present = len(down_patterns.get(pattern, set()))
            up_absent = len(self.up_sequences) - up_present
            down_absent = len(self.down_sequences) - down_present

            # Calculate frequencies
            up_frequency = up_present / len(self.up_sequences) if self.up_sequences else 0
            down_frequency = down_present / len(self.down_sequences) if self.down_sequences else 0

            # Calculate enrichment ratios
            up_to_down_ratio = (
                up_frequency / down_frequency if down_frequency > 0
                else float('inf')
            )
            down_to_up_ratio = (
                down_frequency / up_frequency if up_frequency > 0
                else float('inf')
            )

            # Calculate strand specificity
            combined_sequences = {**self.up_sequences, **self.down_sequences}
            strand_specificity, strand_pval = calculate_strand_specificity(
                pattern, combined_sequences
            )

            # Determine pattern type
            pattern_type = self._determine_pattern_type(pattern)

            results.append({
                'Pattern': pattern,
                'Type': pattern_type,
                'UP_present': up_present,
                'UP_absent': up_absent,
                'DOWN_present': down_present,
                'DOWN_absent': down_absent,
                'UP_frequency': up_frequency,
                'DOWN_frequency': down_frequency,
                'UP_to_DOWN_ratio': up_to_down_ratio,
                'DOWN_to_UP_ratio': down_to_up_ratio,
                'Strand_specificity': strand_specificity,
                'Strand_pvalue': strand_pval
            })

        return results

    def _determine_pattern_type(self, pattern: str) -> str:
        """Determine the type of a pattern"""

        if '(' in pattern:
            return 'mismatch'
        elif '{' in pattern:
            return 'gapped'
        elif pattern.startswith('PWM:'):
            return 'pwm'
        else:
            return 'exact'

    def _determine_significance(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """Determine significant motifs based on enrichment criteria"""

        results_df['Significant'] = False
        results_df['Enrichment'] = 'None'

        # UP-enriched motifs
        up_mask = (
                (results_df['UP_frequency'] >= self.min_occur_freq) &
                (results_df['UP_to_DOWN_ratio'] >= self.fold_change_threshold)
        )
        results_df.loc[up_mask, 'Significant'] = True
        results_df.loc[up_mask, 'Enrichment'] = 'UP'

        # DOWN-enriched motifs
        down_mask = (
                (results_df['DOWN_frequency'] >= self.min_occur_freq) &
                (results_df['DOWN_to_UP_ratio'] >= self.fold_change_threshold)
        )
        results_df.loc[down_mask, 'Significant'] = True
        results_df.loc[down_mask, 'Enrichment'] = 'DOWN'

        return results_df

    def _process_enriched_motifs(
            self,
            significant_df: pd.DataFrame,
            up_positions: Dict,
            down_positions: Dict,
            result_dir: str
    ):
        """Process and cluster enriched motifs"""

        # Separate UP and DOWN enriched motifs
        up_enriched = significant_df[significant_df['Enrichment'] == 'UP']
        down_enriched = significant_df[significant_df['Enrichment'] == 'DOWN']

        # Process UP-enriched motifs
        if not up_enriched.empty:
            up_enriched.to_csv(f'{result_dir}/03_up_enriched_motifs.csv', index=False)
            self._cluster_and_analyze_motifs(
                up_enriched, self.up_sequences, up_positions,
                result_dir, 'up'
            )

        # Process DOWN-enriched motifs
        if not down_enriched.empty:
            down_enriched.to_csv(f'{result_dir}/03_down_enriched_motifs.csv', index=False)
            self._cluster_and_analyze_motifs(
                down_enriched, self.down_sequences, down_positions,
                result_dir, 'down'
            )

    def _cluster_and_analyze_motifs(
            self,
            enriched_df: pd.DataFrame,
            sequences: Dict,
            positions: Dict,
            result_dir: str,
            direction: str
    ):
        """Cluster similar motifs and perform detailed analysis"""

        # Create motif dictionary for clustering
        motifs_dict = dict(zip(
            enriched_df['Pattern'],
            enriched_df[f'{direction.upper()}_to_{("DOWN" if direction == "up" else "UP")}_ratio']
        ))

        # Cluster motifs
        clusters = cluster_similar_motifs(motifs_dict)

        # Save clustered results
        cluster_data = []
        for i, cluster in enumerate(clusters):
            for motif in cluster:
                row = enriched_df[enriched_df['Pattern'] == motif].iloc[0].to_dict()
                row['Cluster'] = i + 1
                row['Cluster_representative'] = cluster[0]
                cluster_data.append(row)

        if cluster_data:
            pd.DataFrame(cluster_data).to_csv(
                f'{result_dir}/04_{direction}_clustered_motifs.csv',
                index=False
            )

            # Analyze representative motifs
            for cluster in clusters:
                self._analyze_motif_details(
                    cluster[0], sequences, positions.get(cluster[0], []),
                    result_dir, direction
                )

    def _analyze_motif_details(
            self,
            pattern: str,
            sequences: Dict,
            pattern_positions: List,
            result_dir: str,
            direction: str
    ):
        """Perform detailed analysis of a specific motif"""

        output_prefix = f'{result_dir}/{direction}_motif_{pattern.replace(":", "_")}'

        # Get matches
        matches, positions, _ = get_matched_sequences(
            pattern, sequences, max_mismatches=self.max_mismatches
        )

        if not matches:
            return

        # Calculate position frequencies
        freq_matrix = calculate_position_frequencies(pattern, matches)
        freq_matrix.to_csv(f'{output_prefix}_frequencies.csv')

        # Calculate strand specificity
        strand_spec, p_value = calculate_strand_specificity(pattern, sequences)

        # Save detailed report
        with open(f'{output_prefix}_details.txt', 'w') as f:
            f.write(f'Pattern: {pattern}\n')
            f.write(f'Number of matching genes: {len(set(p[0] for p in positions))}\n')
            f.write(f'Total matches: {len(matches)}\n')
            f.write(f'Strand specificity: {strand_spec:.2f} (p-value: {p_value:.6f})\n\n')

            f.write('Position frequencies:\n')
            f.write(str(freq_matrix))
            f.write('\n\nMatching genes:\n')
            matching_genes = sorted(set(p[0] for p in positions))
            f.write(', '.join(matching_genes))

        # Save gene list for enrichment analysis
        with open(f'{output_prefix}_genes.txt', 'w') as f:
            for gene in matching_genes:
                f.write(f"{gene}\n")