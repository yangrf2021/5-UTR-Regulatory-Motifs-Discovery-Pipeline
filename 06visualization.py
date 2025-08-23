"""
Visualization utilities for motif analysis
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
from typing import Dict, List, Optional


class VisualizationManager:
    """Manager class for creating various visualizations"""

    def __init__(self):
        """Initialize visualization settings"""
        plt.style.use('seaborn-v0_8-darkgrid')
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.weight'] = 'bold'
        self.colors = {
            'up': '#FF6B6B',
            'down': '#4ECDC4',
            'neutral': '#95A5A6'
        }

    def create_analysis_visualizations(
            self,
            results_df: pd.DataFrame,
            significant_df: pd.DataFrame,
            length: int,
            output_dir: str
    ):
        """Create comprehensive visualization suite"""

        if results_df.empty:
            return

        # Create visualization directory
        viz_dir = os.path.join(output_dir, 'visualizations')
        os.makedirs(viz_dir, exist_ok=True)

        # Generate various plots
        self.plot_frequency_distribution(results_df, length, viz_dir)
        self.plot_enrichment_scatter(results_df, length, viz_dir)
        self.plot_strand_specificity(results_df, length, viz_dir)

        if not significant_df.empty:
            self.plot_significant_heatmap(significant_df, length, viz_dir)
            self.plot_motif_network(significant_df, length, viz_dir)

    def plot_frequency_distribution(
            self,
            results_df: pd.DataFrame,
            length: int,
            output_dir: str
    ):
        """Plot frequency distribution of motifs"""

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # UP frequency distribution
        axes[0].hist(
            results_df['UP_frequency'],
            bins=30,
            color=self.colors['up'],
            alpha=0.7,
            edgecolor='black'
        )
        axes[0].set_xlabel('Frequency in UP-regulated genes')
        axes[0].set_ylabel('Number of motifs')
        axes[0].set_title(f'UP Frequency Distribution (Length {length})')
        axes[0].grid(True, alpha=0.3)

        # DOWN frequency distribution
        axes[1].hist(
            results_df['DOWN_frequency'],
            bins=30,
            color=self.colors['down'],
            alpha=0.7,
            edgecolor='black'
        )
        axes[1].set_xlabel('Frequency in DOWN-regulated genes')
        axes[1].set_ylabel('Number of motifs')
        axes[1].set_title(f'DOWN Frequency Distribution (Length {length})')
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'frequency_distribution.png'), dpi=300)
        plt.close()

    def plot_enrichment_scatter(
            self,
            results_df: pd.DataFrame,
            length: int,
            output_dir: str
    ):
        """Plot enrichment scatter plot"""

        fig, ax = plt.subplots(figsize=(10, 8))

        # Prepare data
        x = results_df['UP_frequency']
        y = results_df['UP_to_DOWN_ratio'].apply(
            lambda val: min(val, 10) if val != float('inf') else 10
        )

        # Color by significance
        colors = []
        for _, row in results_df.iterrows():
            if row.get('Significant', False):
                if row.get('Enrichment') == 'UP':
                    colors.append(self.colors['up'])
                elif row.get('Enrichment') == 'DOWN':
                    colors.append(self.colors['down'])
                else:
                    colors.append(self.colors['neutral'])
            else:
                colors.append(self.colors['neutral'])

        # Create scatter plot
        scatter = ax.scatter(x, y, c=colors, alpha=0.6, s=50, edgecolors='black')

        ax.set_xlabel('Frequency in UP-regulated genes')
        ax.set_ylabel('UP/DOWN Enrichment Ratio (capped at 10)')
        ax.set_title(f'Motif Enrichment Analysis (Length {length})')
        ax.grid(True, alpha=0.3)

        # Add threshold lines
        ax.axhline(y=1.5, color='red', linestyle='--', alpha=0.5, label='Enrichment threshold')
        ax.axvline(x=0.1, color='blue', linestyle='--', alpha=0.5, label='Frequency threshold')

        ax.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'enrichment_scatter.png'), dpi=300)
        plt.close()

    def plot_strand_specificity(
            self,
            results_df: pd.DataFrame,
            length: int,
            output_dir: str
    ):
        """Plot strand specificity distribution"""

        fig, ax = plt.subplots(figsize=(10, 6))

        # Filter out infinite values
        strand_spec = results_df['Strand_specificity'].apply(
            lambda x: min(x, 5) if x != float('inf') else 5
        )

        # Create violin plot
        parts = ax.violinplot(
            [strand_spec],
            positions=[1],
            widths=0.7,
            showmeans=True,
            showmedians=True
        )

        for pc in parts['bodies']:
            pc.set_facecolor(self.colors['neutral'])
            pc.set_alpha(0.7)

        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax.set_ylabel('Strand Specificity (capped at 5)')
        ax.set_title(f'Strand Specificity Distribution (Length {length})')
        ax.set_xticks([1])
        ax.set_xticklabels([f'Length {length}'])
        ax.grid(True, alpha=0.3)

        # Add interpretation text
        ax.text(0.5, 4.5, 'High specificity:\nRNA-level regulation', fontsize=10)
        ax.text(0.5, 0.5, 'Low specificity:\nDNA-level regulation', fontsize=10)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'strand_specificity.png'), dpi=300)
        plt.close()

    def plot_significant_heatmap(
            self,
            significant_df: pd.DataFrame,
            length: int,
            output_dir: str
    ):
        """Plot heatmap of significant motifs"""

        if len(significant_df) > 50:
            # Limit to top 50 motifs
            significant_df = significant_df.nlargest(50, 'UP_to_DOWN_ratio')

        fig, ax = plt.subplots(figsize=(8, max(6, len(significant_df) * 0.3)))

        # Prepare data for heatmap
        heatmap_data = significant_df[['UP_frequency', 'DOWN_frequency']].copy()
        heatmap_data.index = significant_df['Pattern']

        # Create heatmap
        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt='.3f',
            cmap='RdBu_r',
            center=0.5,
            cbar_kws={'label': 'Frequency'},
            ax=ax
        )

        ax.set_title(f'Significant Motif Frequencies (Length {length})')
        ax.set_xlabel('Gene Set')
        ax.set_ylabel('Motif Pattern')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'significant_heatmap.png'), dpi=300)
        plt.close()

    def plot_motif_network(
            self,
            significant_df: pd.DataFrame,
            length: int,
            output_dir: str
    ):
        """Create network visualization of motif relationships"""

        if len(significant_df) < 2:
            return

        fig, ax = plt.subplots(figsize=(12, 10))

        # Create graph
        G = nx.Graph()

        # Add nodes for each significant motif
        for _, row in significant_df.iterrows():
            G.add_node(
                row['Pattern'],
                enrichment=row['Enrichment'],
                frequency=row['UP_frequency'] if row['Enrichment'] == 'UP' else row['DOWN_frequency']
            )

        # Add edges based on sequence similarity
        from .motif_discovery import calculate_motif_similarity

        for i, row1 in significant_df.iterrows():
            for j, row2 in significant_df.iterrows():
                if i < j:
                    similarity = calculate_motif_similarity(
                        row1['Pattern'], row2['Pattern']
                    )
                    if similarity > 0.5:
                        G.add_edge(
                            row1['Pattern'],
                            row2['Pattern'],
                            weight=similarity
                        )

        # Layout
        pos = nx.spring_layout(G, k=2, iterations=50)

        # Draw nodes
        node_colors = [
            self.colors['up'] if G.nodes[n]['enrichment'] == 'UP'
            else self.colors['down']
            for n in G.nodes()
        ]

        node_sizes = [
            G.nodes[n]['frequency'] * 5000
            for n in G.nodes()
        ]

        nx.draw_networkx_nodes(
            G, pos,
            node_color=node_colors,
            node_size=node_sizes,
            alpha=0.7,
            ax=ax
        )

        # Draw edges
        if G.edges():
            edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
            nx.draw_networkx_edges(
                G, pos,
                width=[w * 3 for w in edge_weights],
                alpha=0.5,
                ax=ax
            )

        # Draw labels
        nx.draw_networkx_labels(
            G, pos,
            font_size=8,
            font_weight='bold',
            ax=ax
        )

        ax.set_title(f'Motif Similarity Network (Length {length})')
        ax.axis('off')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'motif_network.png'), dpi=300)
        plt.close()