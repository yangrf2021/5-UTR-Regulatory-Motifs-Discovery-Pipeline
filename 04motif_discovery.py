"""
Motif discovery algorithms and utilities
"""

import re
import multiprocessing
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional
from itertools import product
from tqdm import tqdm
import numpy as np

from .sequence_utils import (
    get_matched_sequences,
    get_matched_sequences_with_pwm_gaps,
    calculate_position_frequencies,
    scan_sequences_with_pwm,
    get_iupac_wildcards
)


def chunk_sequences(sequences: Dict[str, str], n_chunks: int) -> List[Dict[str, str]]:
    """Split sequences into chunks for parallel processing"""
    items = list(sequences.items())
    chunk_size = max(1, len(items) // n_chunks)
    for i in range(0, len(items), chunk_size):
        yield dict(items[i:i + chunk_size])


def process_sequence_chunk(args: Tuple) -> Tuple[Dict, Dict, Dict]:
    """Process a chunk of sequences for motif discovery"""

    (chunk_sequences, length, min_frequency, max_wildcards,
     max_mismatches, enable_gapped, enable_pwm, enable_pwm_gapped) = args

    local_pattern_counts = defaultdict(int)
    sequences_with_pattern = defaultdict(set)
    pattern_positions = defaultdict(list)

    # Find exact words
    exact_words = set()
    for seq_id, seq in chunk_sequences.items():
        for i in range(len(seq) - length + 1):
            word = seq[i:i + length]
            if 'N' not in word:
                exact_words.add(word)
                sequences_with_pattern[word].add(seq_id)
                pattern_positions[word].append((seq_id, i))
                local_pattern_counts[word] += 1

    # Process significant patterns
    significant_patterns = []
    iupac_map = get_iupac_wildcards()

    for word in exact_words:
        count = len(sequences_with_pattern[word])
        if count / len(chunk_sequences) >= min_frequency:
            significant_patterns.append(word)

            # Add wildcards
            if max_wildcards > 0:
                for num_wildcards in range(1, min(max_wildcards + 1, length)):
                    for positions in product(range(length), repeat=num_wildcards):
                        if len(set(positions)) == num_wildcards:
                            new_pattern = list(word)
                            for pos in positions:
                                base = word[pos]
                                for wildcard, bases in iupac_map.items():
                                    if base in bases and len(bases) > 1:
                                        new_pattern[pos] = wildcard
                                        break

                            new_pattern = ''.join(new_pattern)
                            if new_pattern != word:
                                sequences_with_pattern[new_pattern].update(
                                    sequences_with_pattern[word]
                                )
                                pattern_positions[new_pattern].extend(
                                    pattern_positions[word]
                                )

            # Add mismatch variants
            if max_mismatches > 0:
                for pos in range(length):
                    orig_base = word[pos]
                    for base in "ACGT":
                        if base != orig_base:
                            mismatch_variant = word[:pos] + base + word[pos + 1:]
                            matches, positions, _ = get_matched_sequences(
                                mismatch_variant, chunk_sequences
                            )
                            if matches:
                                mismatch_key = f"{word}(m1:{pos}>{base})"
                                sequences_with_pattern[mismatch_key].update(
                                    set(p[0] for p in positions)
                                )
                                pattern_positions[mismatch_key].extend(positions)

            # Add gapped patterns
            if enable_gapped and length >= 5:
                for split_point in range(1, length):
                    prefix = word[:split_point]
                    suffix = word[split_point:]

                    for gap_min, gap_max in [(1, 3), (2, 4), (3, 5)]:
                        gap_pattern = f"{prefix}{{{gap_min},{gap_max}}}{suffix}"
                        matches, positions, _ = get_matched_sequences(
                            None, chunk_sequences, gap_pattern=gap_pattern
                        )
                        if matches and len(matches) / len(chunk_sequences) >= min_frequency:
                            gap_key = f"{prefix}N{{{gap_min},{gap_max}}}{suffix}"
                            sequences_with_pattern[gap_key].update(
                                set(p[0] for p in positions)
                            )
                            pattern_positions[gap_key].extend(positions)

            # Add PWM-based patterns
            if enable_pwm:
                matches, _, _ = get_matched_sequences(word, chunk_sequences)
                if len(matches) >= 10:
                    pwm = calculate_position_frequencies(word, matches)
                    pwm_matches, pwm_positions, _ = scan_sequences_with_pwm(
                        pwm, chunk_sequences, threshold=0.75
                    )
                    if pwm_matches:
                        pwm_key = f"PWM:{word}"
                        sequences_with_pattern[pwm_key].update(
                            set(p[0] for p in pwm_positions)
                        )
                        pattern_positions[pwm_key].extend(pwm_positions)

            # Add PWM-gapped patterns
            if enable_gapped and enable_pwm and enable_pwm_gapped and length >= 5:
                for split_point in range(1, length):
                    prefix = word[:split_point]
                    suffix = word[split_point:]

                    for gap_min, gap_max in [(1, 3), (2, 4), (3, 5)]:
                        matches, positions, _ = get_matched_sequences_with_pwm_gaps(
                            prefix, suffix, gap_min, gap_max,
                            chunk_sequences, threshold=0.75
                        )

                        if matches and len(matches) / len(chunk_sequences) >= min_frequency:
                            pwm_gap_key = f"{prefix}N{{{gap_min},{gap_max}}}{suffix}(pwm)"
                            sequences_with_pattern[pwm_gap_key].update(
                                set(p[0] for p in positions)
                            )
                            pattern_positions[pwm_gap_key].extend(positions)

    return (
        dict(sequences_with_pattern),
        dict(local_pattern_counts),
        dict(pattern_positions)
    )


def find_motifs(
        sequences: Dict[str, str],
        length: int,
        min_frequency: float = 0.05,
        n_cores: Optional[int] = None,
        max_wildcards: int = 1,
        max_mismatches: int = 0,
        enable_gapped: bool = False,
        enable_pwm: bool = False,
        enable_pwm_gapped: bool = False
) -> Tuple[Dict, Dict, Dict]:
    """
    Find motifs in sequences using parallel processing

    Args:
        sequences: Dictionary of sequences
        length: Motif length
        min_frequency: Minimum frequency threshold
        n_cores: Number of CPU cores
        max_wildcards: Maximum wildcards allowed
        max_mismatches: Maximum mismatches allowed
        enable_gapped: Enable gapped motif discovery
        enable_pwm: Enable PWM analysis
        enable_pwm_gapped: Enable PWM-gapped analysis

    Returns:
        Tuple of (sequences_with_pattern, pattern_counts, pattern_positions)
    """

    if n_cores is None:
        n_cores = max(1, multiprocessing.cpu_count() - 1)

    sequence_chunks = list(chunk_sequences(sequences, n_cores))
    args_list = [
        (chunk, length, min_frequency, max_wildcards, max_mismatches,
         enable_gapped, enable_pwm, enable_pwm_gapped)
        for chunk in sequence_chunks
    ]

    with multiprocessing.Pool(processes=n_cores) as pool:
        chunk_results = list(tqdm(
            pool.imap(process_sequence_chunk, args_list),
            total=len(args_list),
            desc=f"Finding motifs of length {length}"
        ))

    # Merge results
    all_sequences_with_pattern = defaultdict(set)
    pattern_counts = defaultdict(int)
    all_pattern_positions = defaultdict(list)

    for sequences_dict, counts_dict, positions_dict in chunk_results:
        for pattern, seqs in sequences_dict.items():
            all_sequences_with_pattern[pattern].update(seqs)

        for pattern, count in counts_dict.items():
            pattern_counts[pattern] += count

        for pattern, positions in positions_dict.items():
            all_pattern_positions[pattern].extend(positions)

    return (
        dict(all_sequences_with_pattern),
        dict(pattern_counts),
        dict(all_pattern_positions)
    )


def find_structural_motifs(
        sequences: Dict[str, str],
        min_stem_length: int = 4,
        max_loop_size: int = 8
) -> Tuple[List[Dict], List[Tuple]]:
    """
    Find structural motifs (stem-loops) in sequences

    Note: This is a simplified version. For actual RNA structure prediction,
    use packages like ViennaRNA.

    Args:
        sequences: Dictionary of sequences
        min_stem_length: Minimum stem length
        max_loop_size: Maximum loop size

    Returns:
        Tuple of (motifs, positions)
    """
    structural_motifs = []
    positions = []

    # Simplified stem-loop detection
    for seq_id, seq in sequences.items():
        for i in range(len(seq) - 2 * min_stem_length - max_loop_size):
            for stem_len in range(min_stem_length, min(10, (len(seq) - i) // 2)):
                for loop_len in range(3, max_loop_size + 1):
                    j = i + stem_len + loop_len
                    if j + stem_len <= len(seq):
                        stem1 = seq[i:i + stem_len]
                        stem2 = seq[j:j + stem_len]

                        # Check for complementarity (simplified)
                        complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                                      'T': 'U'}  # T treated as U for RNA

                        is_complement = True
                        for k in range(stem_len):
                            if stem1[k] not in complement:
                                is_complement = False
                                break
                            if complement[stem1[k]] != stem2[-(k + 1)]:
                                is_complement = False
                                break

                        if is_complement:
                            structural_motifs.append({
                                'type': 'stem_loop',
                                'sequence': seq[i:j + stem_len],
                                'stem_length': stem_len,
                                'loop_length': loop_len
                            })
                            positions.append((seq_id, i))

    return structural_motifs, positions


def find_g_quadruplexes(sequences: Dict[str, str]) -> Tuple[List[Dict], List[Tuple]]:
    """
    Find potential G-quadruplex structures in sequences

    Args:
        sequences: Dictionary of sequences

    Returns:
        Tuple of (quadruplexes, positions)
    """
    # G-quadruplex pattern: (G{3,5}[ACGT]{1,7}){3,}G{3,5}
    g_quad_pattern = r'(G{3,5}[ACGT]{1,7}){3,}G{3,5}'
    compiled_pattern = re.compile(g_quad_pattern)

    quadruplexes = []
    positions = []

    for seq_id, seq in sequences.items():
        for match in compiled_pattern.finditer(seq):
            quadruplexes.append({
                'type': 'g_quadruplex',
                'sequence': match.group(0),
                'start': match.start(),
                'end': match.end()
            })
            positions.append((seq_id, match.start()))

    return quadruplexes, positions


def cluster_similar_motifs(
        motifs: Dict[str, float],
        similarity_threshold: float = 0.5
) -> List[List[str]]:
    """
    Cluster similar motifs together

    Args:
        motifs: Dictionary of motifs with their scores
        similarity_threshold: Minimum similarity for clustering

    Returns:
        List of motif clusters
    """
    clusters = []
    assigned = set()

    sorted_motifs = sorted(motifs.keys(), key=lambda x: motifs[x], reverse=True)

    for motif in sorted_motifs:
        if motif in assigned or motif.startswith("PWM:"):
            continue

        cluster = [motif]
        assigned.add(motif)

        for other in sorted_motifs:
            if other != motif and other not in assigned and not other.startswith("PWM:"):
                similarity = calculate_motif_similarity(motif, other)

                if similarity >= similarity_threshold:
                    cluster.append(other)
                    assigned.add(other)

        clusters.append(cluster)

    # Add PWM motifs as individual clusters
    for motif in sorted_motifs:
        if motif.startswith("PWM:"):
            clusters.append([motif])

    return clusters


def calculate_motif_similarity(motif1: str, motif2: str) -> float:
    """
    Calculate similarity between two motifs

    Args:
        motif1: First motif
        motif2: Second motif

    Returns:
        Similarity score (0-1)
    """
    # Handle special formats
    if '{' in motif1 or '{' in motif2:
        # Extract prefix and suffix for gapped motifs
        if '{' in motif1:
            parts1 = re.split(r'\{.*?\}', motif1)
        else:
            parts1 = [motif1]

        if '{' in motif2:
            parts2 = re.split(r'\{.*?\}', motif2)
        else:
            parts2 = [motif2]

        # Compare prefixes
        prefix_sim = 0
        if parts1[0] and parts2[0]:
            min_len = min(len(parts1[0]), len(parts2[0]))
            prefix_sim = sum(parts1[0][i] == parts2[0][i] for i in range(min_len)) / min_len

        # Compare suffixes
        suffix_sim = 0
        if len(parts1) > 1 and len(parts2) > 1:
            if parts1[-1] and parts2[-1]:
                min_len = min(len(parts1[-1]), len(parts2[-1]))
                suffix_sim = sum(parts1[-1][i] == parts2[-1][i] for i in range(min_len)) / min_len

        return (prefix_sim + suffix_sim) / 2

    elif '(' in motif1 or '(' in motif2:
        # Extract base patterns for mismatch motifs
        base1 = motif1.split('(')[0] if '(' in motif1 else motif1
        base2 = motif2.split('(')[0] if '(' in motif2 else motif2

        min_len = min(len(base1), len(base2))
        return sum(base1[i] == base2[i] for i in range(min_len)) / min_len

    else:
        # Regular similarity calculation
        min_len = min(len(motif1), len(motif2))
        similarity = sum(motif1[i] == motif2[i] for i in range(min_len)) / min_len
        return similarity