"""
Utilities for sequence manipulation and pattern matching
"""

import re
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


def get_iupac_wildcards() -> Dict[str, List[str]]:
    """
    Return IUPAC wildcard mappings

    Returns:
        Dictionary mapping IUPAC codes to nucleotide lists
    """
    return {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'],
        'M': ['A', 'C'], 'K': ['G', 'T'],
        'S': ['G', 'C'], 'W': ['A', 'T'],
        'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
        'V': ['A', 'C', 'G'], 'D': ['A', 'G', 'T'],
        'N': ['A', 'C', 'G', 'T']
    }


def get_matched_sequences(
        pattern: str,
        sequences: Dict[str, str],
        max_mismatches: int = 0,
        gap_pattern: Optional[str] = None
) -> Tuple[List[str], List[Tuple[str, int]], List]:
    """
    Find all sequences matching a pattern

    Args:
        pattern: Pattern to search for
        sequences: Dictionary of sequences
        max_mismatches: Maximum allowed mismatches
        gap_pattern: Optional gap pattern (e.g., "GAG{2,4}GG")

    Returns:
        Tuple of (matches, positions, mismatch_details)
    """
    matches = []
    match_positions = []
    mismatch_details = []
    iupac_map = get_iupac_wildcards()

    # Handle gap patterns
    if gap_pattern:
        gap_match = re.match(r"(.+)\{(\d+),(\d+)\}(.+)", gap_pattern)
        if gap_match:
            prefix, min_gap, max_gap, suffix = gap_match.groups()
            min_gap, max_gap = int(min_gap), int(max_gap)

            # Create regex pattern
            regex_pattern = ''
            for char in prefix:
                if char in iupac_map:
                    regex_pattern += f'[{"".join(iupac_map[char])}]'
                else:
                    regex_pattern += char

            regex_pattern += f'.{{{min_gap},{max_gap}}}'

            for char in suffix:
                if char in iupac_map:
                    regex_pattern += f'[{"".join(iupac_map[char])}]'
                else:
                    regex_pattern += char

            compiled_regex = re.compile(regex_pattern)

            for seq_id, seq in sequences.items():
                for match in compiled_regex.finditer(seq):
                    matches.append(match.group(0))
                    match_positions.append((seq_id, match.start()))

            return matches, match_positions, []

    # Regular pattern matching
    regex_pattern = ''
    for char in pattern:
        if char in iupac_map:
            regex_pattern += f'[{"".join(iupac_map[char])}]'

    if max_mismatches == 0:
        compiled_pattern = re.compile(regex_pattern)
        for seq_id, seq in sequences.items():
            for i in range(len(seq) - len(pattern) + 1):
                subsequence = seq[i:i + len(pattern)]
                if compiled_pattern.match(subsequence):
                    matches.append(subsequence)
                    match_positions.append((seq_id, i))
    else:
        # Handle mismatches
        for seq_id, seq in sequences.items():
            for i in range(len(seq) - len(pattern) + 1):
                subsequence = seq[i:i + len(pattern)]
                mismatches = 0
                mismatch_pos = []

                for j in range(len(pattern)):
                    pattern_char = pattern[j]
                    subseq_char = subsequence[j]

                    if pattern_char in iupac_map:
                        if subseq_char not in iupac_map[pattern_char]:
                            mismatches += 1
                            mismatch_pos.append((j, pattern_char, subseq_char))
                    else:
                        if pattern_char != subseq_char:
                            mismatches += 1
                            mismatch_pos.append((j, pattern_char, subseq_char))

                if mismatches <= max_mismatches:
                    matches.append(subsequence)
                    match_positions.append((seq_id, i))
                    if mismatch_pos:
                        mismatch_details.append((seq_id, i, mismatch_pos))

    return matches, match_positions, mismatch_details


def calculate_position_frequencies(
        pattern: str,
        matches: List[str]
) -> pd.DataFrame:
    """
    Calculate nucleotide frequencies at each position

    Args:
        pattern: Reference pattern
        matches: List of matching sequences

    Returns:
        DataFrame with position frequencies
    """
    length = len(pattern)
    freq_matrix = pd.DataFrame(
        0,
        index=['A', 'C', 'G', 'T'],
        columns=range(length)
    )

    for seq in matches:
        for pos in range(min(length, len(seq))):
            base = seq[pos]
            if base in 'ACGT':
                freq_matrix.loc[base, pos] += 1

    # Convert to frequencies
    return freq_matrix / len(matches) if matches else freq_matrix


def calculate_pwm_score(pwm: pd.DataFrame, sequence: str) -> float:
    """
    Calculate PWM score for a sequence

    Args:
        pwm: Position Weight Matrix
        sequence: Sequence to score

    Returns:
        PWM score
    """
    score = 0
    for i, base in enumerate(sequence):
        if i < len(pwm.columns) and base in pwm.index:
            score += pwm.loc[base, i]
    return score


def scan_sequences_with_pwm(
        pwm: pd.DataFrame,
        sequences: Dict[str, str],
        threshold: float = 0.8
) -> Tuple[List[str], List[Tuple[str, int]], List[float]]:
    """
    Scan sequences using a Position Weight Matrix

    Args:
        pwm: Position Weight Matrix
        sequences: Dictionary of sequences
        threshold: Matching threshold (0-1)

    Returns:
        Tuple of (matches, positions, scores)
    """
    matches = []
    match_positions = []
    match_scores = []

    motif_length = len(pwm.columns)
    max_score = motif_length
    threshold_score = max_score * threshold

    for seq_id, seq in sequences.items():
        for i in range(len(seq) - motif_length + 1):
            subsequence = seq[i:i + motif_length]
            score = calculate_pwm_score(pwm, subsequence)

            if score >= threshold_score:
                matches.append(subsequence)
                match_positions.append((seq_id, i))
                match_scores.append(score)

    return matches, match_positions, match_scores


def get_matched_sequences_with_pwm_gaps(
        prefix: str,
        suffix: str,
        min_gap: int,
        max_gap: int,
        sequences: Dict[str, str],
        prefix_pwm: Optional[pd.DataFrame] = None,
        suffix_pwm: Optional[pd.DataFrame] = None,
        threshold: float = 0.75
) -> Tuple[List[str], List[Tuple[str, int]], List]:
    """
    Find sequences with gapped motifs using PWM

    Args:
        prefix: Prefix pattern
        suffix: Suffix pattern
        min_gap: Minimum gap length
        max_gap: Maximum gap length
        sequences: Dictionary of sequences
        prefix_pwm: Optional prefix PWM
        suffix_pwm: Optional suffix PWM
        threshold: PWM matching threshold

    Returns:
        Tuple of (matches, positions, details)
    """
    matches = []
    match_positions = []

    # Create PWMs if not provided
    if prefix_pwm is None:
        prefix_matches, _, _ = get_matched_sequences(prefix, sequences)
        if len(prefix_matches) >= 10:
            prefix_pwm = calculate_position_frequencies(prefix, prefix_matches)
        else:
            prefix_pwm = pd.DataFrame(
                0, index=['A', 'C', 'G', 'T'],
                columns=range(len(prefix))
            )
            for i, base in enumerate(prefix):
                if base in 'ACGT':
                    prefix_pwm.loc[base, i] = 1.0

    if suffix_pwm is None:
        suffix_matches, _, _ = get_matched_sequences(suffix, sequences)
        if len(suffix_matches) >= 10:
            suffix_pwm = calculate_position_frequencies(suffix, suffix_matches)
        else:
            suffix_pwm = pd.DataFrame(
                0, index=['A', 'C', 'G', 'T'],
                columns=range(len(suffix))
            )
            for i, base in enumerate(suffix):
                if base in 'ACGT':
                    suffix_pwm.loc[base, i] = 1.0

    # Set thresholds
    prefix_threshold = len(prefix) * threshold
    suffix_threshold = len(suffix) * threshold

    # Scan sequences
    for seq_id, seq in sequences.items():
        for i in range(len(seq) - len(prefix) + 1):
            prefix_candidate = seq[i:i + len(prefix)]
            prefix_score = calculate_pwm_score(prefix_pwm, prefix_candidate)

            if prefix_score >= prefix_threshold:
                for gap_len in range(min_gap, max_gap + 1):
                    suffix_start = i + len(prefix) + gap_len
                    if suffix_start + len(suffix) <= len(seq):
                        suffix_candidate = seq[suffix_start:suffix_start + len(suffix)]
                        suffix_score = calculate_pwm_score(suffix_pwm, suffix_candidate)

                        if suffix_score >= suffix_threshold:
                            full_match = seq[i:suffix_start + len(suffix)]
                            matches.append(full_match)
                            match_positions.append((seq_id, i))
                            break

    return matches, match_positions, []