"""
Statistical analysis utilities
"""

import numpy as np
from scipy import stats
from typing import Dict, Tuple


def calculate_strand_specificity(
        motif: str,
        sequences: Dict[str, str]
) -> Tuple[float, float]:
    """
    Calculate strand specificity of a motif

    Strand specificity indicates whether a motif shows preference
    for one DNA strand over its reverse complement, which can
    indicate RNA-level vs DNA-level regulation.

    Args:
        motif: Motif pattern
        sequences: Dictionary of sequences

    Returns:
        Tuple of (specificity_ratio, p_value)
    """

    # Handle special format motifs
    base_motif = extract_base_motif(motif)

    # Get reverse complement
    rc_motif = reverse_complement(base_motif)

    # Count matches
    from .sequence_utils import get_matched_sequences

    forward_matches, _, _ = get_matched_sequences(base_motif, sequences)
    reverse_matches, _, _ = get_matched_sequences(rc_motif, sequences)

    # Calculate specificity
    if len(reverse_matches) == 0:
        return float('inf'), 0

    specificity = len(forward_matches) / len(reverse_matches)

    # Calculate p-value using binomial test
    total = len(forward_matches) + len(reverse_matches)
    if total > 0:
        p_value = stats.binom_test(
            len(forward_matches), total, p=0.5, alternative='two-sided'
        )
    else:
        p_value = 1.0

    return specificity, p_value


def extract_base_motif(motif: str) -> str:
    """
    Extract base motif from special format patterns

    Args:
        motif: Motif pattern (may include special markers)

    Returns:
        Base motif sequence
    """
    import re

    base_motif = motif

    # Remove PWM markers
    if '(pwm)' in motif:
        base_motif = motif.replace('(pwm)', '')

    # Handle mismatch patterns
    if '(' in motif and ')' in motif:
        base_motif = motif.split('(')[0]

    # Handle gapped patterns
    if '{' in motif:
        parts = re.split(r'\{.*?\}', motif)
        base_motif = ''.join(parts)

    # Remove PWM prefix
    if motif.startswith('PWM:'):
        base_motif = motif[4:]

    return base_motif


def reverse_complement(sequence: str) -> str:
    """
    Get reverse complement of a DNA sequence

    Args:
        sequence: DNA sequence

    Returns:
        Reverse complement sequence
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
        'S': 'S', 'W': 'W', 'H': 'D', 'B': 'V',
        'V': 'B', 'D': 'H', 'N': 'N'
    }

    return ''.join([
        complement_map.get(base, base)
        for base in sequence[::-1]
    ])


def calculate_enrichment_statistics(
        up_count: int,
        up_total: int,
        down_count: int,
        down_total: int
) -> Dict[str, float]:
    """
    Calculate enrichment statistics for a motif

    Args:
        up_count: Count in upregulated genes
        up_total: Total upregulated genes
        down_count: Count in downregulated genes
        down_total: Total downregulated genes

    Returns:
        Dictionary of statistics
    """

    # Calculate frequencies
    up_freq = up_count / up_total if up_total > 0 else 0
    down_freq = down_count / down_total if down_total > 0 else 0

    # Calculate fold changes
    up_to_down = up_freq / down_freq if down_freq > 0 else float('inf')
    down_to_up = down_freq / up_freq if up_freq > 0 else float('inf')

    # Fisher's exact test
    oddsratio, pvalue = stats.fisher_exact([
        [up_count, up_total - up_count],
        [down_count, down_total - down_count]
    ])

    # Hypergeometric test
    total = up_total + down_total
    total_with_motif = up_count + down_count

    hypergeom_pval = stats.hypergeom.sf(
        up_count - 1, total, total_with_motif, up_total
    )

    return {
        'up_frequency': up_freq,
        'down_frequency': down_freq,
        'up_to_down_ratio': up_to_down,
        'down_to_up_ratio': down_to_up,
        'fisher_pvalue': pvalue,
        'odds_ratio': oddsratio,
        'hypergeom_pvalue': hypergeom_pval
    }