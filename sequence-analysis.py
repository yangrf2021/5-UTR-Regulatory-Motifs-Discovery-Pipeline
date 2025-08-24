import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
from tqdm import tqdm
import multiprocessing
from itertools import product, combinations
from functools import partial
import re
import logging
import matplotlib.patches as mpatches
import networkx as nx  # For structural analysis
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('sequence_analysis')


# After logging.basicConfig setup and logger definition
def set_random_seeds(seed=42):
    """Set random seeds to ensure reproducible results"""
    import random
    random.seed(seed)

    # NumPy random seed
    np.random.seed(seed)

    # Set seed for matplotlib to ensure consistent figure generation
    plt.rcParams['svg.hashsalt'] = str(seed)

    logger.info(f"Random seed set to {seed} to ensure reproducible results")


# Set font
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'

# Set plotting style
plt.style.use('seaborn')


def get_iupac_wildcards():
    """
    Return IUPAC wildcard mappings

    Purpose: Support wildcard matching in motif searching
    """
    return {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'M': ['A', 'C'],
        'K': ['G', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'H': ['A', 'C', 'T'],
        'B': ['C', 'G', 'T'],
        'V': ['A', 'C', 'G'],
        'D': ['A', 'G', 'T'],
        'N': ['A', 'C', 'G', 'T']
    }


def chunk_sequences(sequences, n_chunks):
    """
    Split sequence dictionary into n chunks

    Purpose: Support parallel processing to improve performance
    """
    items = list(sequences.items())
    chunk_size = max(1, len(items) // n_chunks)  # Ensure at least 1
    for i in range(0, len(items), chunk_size):
        yield dict(items[i:i + chunk_size])


def get_matched_sequences(pattern, sequences, max_mismatches=0, gap_pattern=None):
    """
    Get all sequences matching specified pattern with positions, supporting mismatches and gap patterns

    Parameters:
        pattern: Basic pattern string
        sequences: Sequence dictionary
        max_mismatches: Maximum allowed mismatches
        gap_pattern: Gap pattern (format like: "GAG{2,4}GG" means 2-4 arbitrary bases between GAG and GG)

    Returns:
        matches: List of matching sequence fragments
        match_positions: List of match position information
        mismatch_details: Mismatch details (only meaningful when max_mismatches>0)
    """
    matches = []
    match_positions = []  # Record match positions
    mismatch_details = []  # Record mismatch details
    iupac_map = get_iupac_wildcards()

    # Handle gap patterns
    if gap_pattern:
        # Parse gap pattern, e.g.: "GAG{2,4}GG"
        gap_match = re.match(r"(.+)\{(\d+),(\d+)\}(.+)", gap_pattern)
        if gap_match:
            prefix, min_gap, max_gap, suffix = gap_match.groups()
            min_gap, max_gap = int(min_gap), int(max_gap)

            # Create regex for gap pattern
            regex_pattern = ''
            for char in prefix:
                if char in iupac_map:
                    regex_pattern += f'[{"".join(iupac_map[char])}]'
                else:
                    regex_pattern += char

            # Add gap part
            regex_pattern += f'.{{{min_gap},{max_gap}}}'

            # Add suffix
            for char in suffix:
                if char in iupac_map:
                    regex_pattern += f'[{"".join(iupac_map[char])}]'
                else:
                    regex_pattern += char

            # Compile regex
            compiled_regex = re.compile(regex_pattern)

            # Find all matches
            for seq_id, seq in sequences.items():
                for match in compiled_regex.finditer(seq):
                    matches.append(match.group(0))
                    match_positions.append((seq_id, match.start()))

            return matches, match_positions, []

    # Non-gap pattern processing
    # Pre-compile regex for basic pattern
    regex_pattern = ''
    for char in pattern:
        if char in iupac_map:
            regex_pattern += f'[{"".join(iupac_map[char])}]'
    compiled_pattern = re.compile(regex_pattern)

    # If no mismatches allowed, use regex for direct matching
    if max_mismatches == 0:
        for seq_id, seq in sequences.items():
            for i in range(len(seq) - len(pattern) + 1):
                subsequence = seq[i:i + len(pattern)]
                if compiled_pattern.match(subsequence):
                    matches.append(subsequence)
                    match_positions.append((seq_id, i))
    else:
        # If mismatches allowed, need position-by-position comparison
        for seq_id, seq in sequences.items():
            for i in range(len(seq) - len(pattern) + 1):
                subsequence = seq[i:i + len(pattern)]
                mismatches = 0
                mismatch_pos = []

                for j in range(len(pattern)):
                    pattern_char = pattern[j]
                    subseq_char = subsequence[j]

                    # Check if matches (considering IUPAC wildcards)
                    if pattern_char in iupac_map:
                        if subseq_char not in iupac_map[pattern_char]:
                            mismatches += 1
                            mismatch_pos.append((j, pattern_char, subseq_char))
                    else:
                        if pattern_char != subseq_char:
                            mismatches += 1
                            mismatch_pos.append((j, pattern_char, subseq_char))

                # If mismatches within allowed range, record match
                if mismatches <= max_mismatches:
                    matches.append(subsequence)
                    match_positions.append((seq_id, i))
                    if mismatch_pos:
                        mismatch_details.append((seq_id, i, mismatch_pos))

    return matches, match_positions, mismatch_details


def calculate_position_frequencies(pattern, matches):
    """
    Calculate nucleotide frequency at each position

    Purpose: Determine nucleotide preferences at each position in the motif
    """
    length = len(pattern)
    freq_matrix = pd.DataFrame(0,
                               index=['A', 'C', 'G', 'T'],
                               columns=range(length))

    for seq in matches:
        for pos in range(length):
            if pos < len(seq):  # Prevent errors from unequal length sequences
                base = seq[pos]
                if base in 'ACGT':  # Ensure only standard bases are counted
                    freq_matrix.loc[base, pos] += 1

    # Convert to frequencies
    return freq_matrix / len(matches) if matches else freq_matrix


def calculate_pwm_score(pwm, sequence):
    """
    Calculate sequence score using Position Weight Matrix (PWM)

    Parameters:
        pwm: Position Weight Matrix
        sequence: Sequence to score

    Returns:
        score: PWM score
    """
    score = 0
    for i, base in enumerate(sequence):
        if i < len(pwm.columns) and base in pwm.index:
            score += pwm.loc[base, i]
    return score


def scan_sequences_with_pwm(pwm, sequences, threshold=0.8):
    """
    Scan sequences using PWM to find potential matches

    Parameters:
        pwm: Position Weight Matrix
        sequences: Sequence dictionary
        threshold: Matching threshold (0-1)

    Returns:
        matches: Matching sequences and their scores
        match_positions: Match positions
    """
    matches = []
    match_positions = []
    match_scores = []

    motif_length = len(pwm.columns)
    max_score = motif_length  # Maximum score under ideal conditions
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


def get_matched_sequences_with_pwm_gaps(prefix, suffix, min_gap, max_gap, sequences, prefix_pwm=None, suffix_pwm=None,
                                        threshold=0.75):
    """
    Flexible gap pattern matching using PWM method

    Parameters:
        prefix: Prefix pattern
        suffix: Suffix pattern
        min_gap: Minimum gap length
        max_gap: Maximum gap length
        sequences: Sequences to search
        prefix_pwm: Prefix PWM matrix (optional)
        suffix_pwm: Suffix PWM matrix (optional)
        threshold: PWM matching threshold (0-1)

    Returns:
        matches: List of matching sequence fragments
        match_positions: List of match position information
    """
    matches = []
    match_positions = []

    # If PWM not provided, calculate PWM based on exact matches
    if prefix_pwm is None:
        # First find exact matches of prefix
        prefix_matches, _, _ = get_matched_sequences(prefix, sequences)
        if len(prefix_matches) >= 10:  # Ensure sufficient samples
            prefix_pwm = calculate_position_frequencies(prefix, prefix_matches)
        else:
            # Insufficient samples, use simple diagonal matrix
            prefix_pwm = pd.DataFrame(0, index=['A', 'C', 'G', 'T'], columns=range(len(prefix)))
            for i, base in enumerate(prefix):
                prefix_pwm.loc[base, i] = 1.0

    # Similarly process suffix PWM
    if suffix_pwm is None:
        suffix_matches, _, _ = get_matched_sequences(suffix, sequences)
        if len(suffix_matches) >= 10:
            suffix_pwm = calculate_position_frequencies(suffix, suffix_matches)
        else:
            suffix_pwm = pd.DataFrame(0, index=['A', 'C', 'G', 'T'], columns=range(len(suffix)))
            for i, base in enumerate(suffix):
                suffix_pwm.loc[base, i] = 1.0

    # Set score thresholds
    prefix_max_score = len(prefix)
    suffix_max_score = len(suffix)
    prefix_threshold = prefix_max_score * threshold
    suffix_threshold = suffix_max_score * threshold

    # Scan sequences
    for seq_id, seq in sequences.items():
        for i in range(len(seq) - len(prefix) + 1):
            # Calculate prefix PWM score
            prefix_candidate = seq[i:i + len(prefix)]
            prefix_score = calculate_pwm_score(prefix_pwm, prefix_candidate)

            if prefix_score >= prefix_threshold:
                # Prefix match successful, check various possible gap lengths
                for gap_len in range(min_gap, max_gap + 1):
                    suffix_start = i + len(prefix) + gap_len
                    if suffix_start + len(suffix) <= len(seq):
                        # Calculate suffix PWM score
                        suffix_candidate = seq[suffix_start:suffix_start + len(suffix)]
                        suffix_score = calculate_pwm_score(suffix_pwm, suffix_candidate)

                        if suffix_score >= suffix_threshold:
                            # Complete match successful
                            full_match = seq[i:suffix_start + len(suffix)]
                            matches.append(full_match)
                            match_positions.append((seq_id, i))
                            # Break after finding one match
                            break

    return matches, match_positions, []


def predict_rna_structure(sequence):
    """
    Predict secondary structure of RNA sequence

    Using Vienna RNA package for prediction

    Parameters:
        sequence: RNA sequence

    Returns:
        structure: Dot-bracket notation of structure
        mfe: Minimum free energy
    """
    structure, mfe = RNA.fold(sequence)
    return structure, mfe


def find_structural_motifs(sequences, min_stem_length=4, max_loop_size=8):
    """
    Find structural motifs (like stem-loop structures) in sequences

    Parameters:
        sequences: Sequence dictionary
        min_stem_length: Minimum stem length
        max_loop_size: Maximum loop size

    Returns:
        motifs: Discovered structural motifs
        positions: Motif positions
    """
    structural_motifs = []
    positions = []

    for seq_id, seq in sequences.items():
        # Predict structure for entire sequence
        structure, _ = predict_rna_structure(seq)

        # Parse structure to find stem-loops
        i = 0
        while i < len(structure):
            # Find opening bracket, representing stem start
            if structure[i] == '(':
                # Find corresponding closing bracket
                stack = 1
                j = i + 1
                while j < len(structure) and stack > 0:
                    if structure[j] == '(':
                        stack += 1
                    elif structure[j] == ')':
                        stack -= 1
                    j += 1

                # Check if complete stem-loop found
                if j - i - 1 <= max_loop_size + 2 * min_stem_length:  # Limit total length
                    stem_loop = structure[i:j]
                    # Calculate stem length
                    stem_length = sum(1 for c in stem_loop if c == '(')
                    if stem_length >= min_stem_length:
                        structural_motifs.append({
                            'type': 'stem_loop',
                            'structure': stem_loop,
                            'sequence': seq[i:j],
                            'stem_length': stem_length
                        })
                        positions.append((seq_id, i))
            i += 1

    return structural_motifs, positions


def find_g_quadruplexes(sequences):
    """
    Find G-quadruplex structures in sequences

    Parameters:
        sequences: Sequence dictionary

    Returns:
        quadruplexes: Discovered G-quadruplexes
        positions: G-quadruplex positions
    """
    # G-quadruplex pattern: (G{3,5}[ACGT]{1,7}){3,}G{3,5}
    # This matches at least 3 groups of consecutive G (3-5) plus 1-7 arbitrary bases, plus a final group of G
    g_quad_pattern = r'(G{3,5}[ACGT]{1,7}){3,}G{3,5}'
    compiled_pattern = re.compile(g_quad_pattern)

    quadruplexes = []
    positions = []

    for seq_id, seq in sequences.items():
        for match in compiled_pattern.finditer(seq):
            quadruplexes.append({
                'type': 'g_quadruplex',
                'sequence': match.group(0)
            })
            positions.append((seq_id, match.start()))

    return quadruplexes, positions


def calculate_strand_specificity(motif, sequences):
    """
    Calculate strand specificity of pattern (forward vs reverse complement)

    Add handling for special characters
    """
    # Handle specially formatted patterns
    # If pattern contains special format markers (like brackets, braces), extract basic pattern part
    base_motif = motif
    if '(' in motif:  # Mismatch pattern or PWM gap pattern
        if '(pwm)' in motif:
            base_motif = motif.replace('(pwm)', '')  # Remove pwm marker
        else:
            base_motif = motif.split('(')[0]  # Mismatch pattern
    elif '{' in motif:  # Gap pattern
        parts = re.split(r'\{.*?\}', motif)
        base_motif = ''.join(parts)
    elif motif.startswith('PWM:'):  # PWM pattern
        base_motif = motif[4:]  # Remove "PWM:" prefix

    # Get reverse complement sequence
    rc_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R',
              'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W', 'H': 'D', 'B': 'V',
              'V': 'B', 'D': 'H', 'N': 'N'}

    # Safely get reverse complement, return themselves for non-standard characters
    rc_motif = ''.join([rc_map.get(c, c) for c in base_motif[::-1]])

    # Calculate forward and reverse match counts
    forward_matches, _, _ = get_matched_sequences(base_motif, sequences)
    reverse_matches, _, _ = get_matched_sequences(rc_motif, sequences)

    # Calculate strand specificity
    if len(reverse_matches) == 0:
        return float('inf'), 0  # Avoid division by zero

    specificity = len(forward_matches) / len(reverse_matches)

    # Calculate p-value (binomial distribution approximation)
    total = len(forward_matches) + len(reverse_matches)
    p_value = stats.binom_test(len(forward_matches), total, p=0.5, alternative='two-sided')

    return specificity, p_value


def process_sequence_chunk(args):
    """
    Process a sequence chunk, search for motifs

    Purpose: Process a group of sequences in parallel, find possible motifs and allow generalization using wildcards

    Extension: Add support for mismatches, gaps and PWM
    """
    chunk_sequences, length, min_frequency, max_wildcards, max_mismatches, enable_gapped, enable_pwm, enable_pwm_gapped = args  # Modified this line
    local_pattern_counts = defaultdict(int)
    sequences_with_pattern = defaultdict(set)
    pattern_positions = defaultdict(list)

    # Step 1: Find all exact words (no wildcards, no mismatches)
    exact_words = set()
    for seq_id, seq in chunk_sequences.items():
        for i in range(len(seq) - length + 1):
            word = seq[i:i + length]
            if 'N' not in word:  # Exclude sequences containing N
                exact_words.add(word)
                sequences_with_pattern[word].add(seq_id)
                pattern_positions[word].append((seq_id, i))
                local_pattern_counts[word] += 1

    # Step 2: Wildcard and mismatch generalization for significant words
    significant_patterns = []
    iupac_map = get_iupac_wildcards()

    for word in exact_words:
        count = len(sequences_with_pattern[word])
        if count / len(chunk_sequences) >= min_frequency:
            significant_patterns.append(word)

            # Try adding wildcards (up to max_wildcards)
            if max_wildcards > 0:
                # Find all possible position combinations
                if max_wildcards >= length:
                    max_wildcards = length - 1  # Avoid all wildcards

                for num_wildcards in range(1, max_wildcards + 1):
                    for positions in product(range(length), repeat=num_wildcards):
                        if len(set(positions)) == num_wildcards:  # Ensure positions don't repeat
                            new_pattern = list(word)
                            for pos in positions:
                                # Find suitable wildcard for this position
                                base = word[pos]
                                for wildcard, bases in iupac_map.items():
                                    if base in bases and len(bases) > 1:
                                        new_pattern[pos] = wildcard
                                        break

                            new_pattern = ''.join(new_pattern)
                            if new_pattern != word:  # Ensure pattern has changed
                                sequences_with_pattern[new_pattern].update(
                                    sequences_with_pattern[word]
                                )
                                # Copy position information
                                pattern_positions[new_pattern].extend(pattern_positions[word])

            # Try mismatch generalization (if enabled)
            if max_mismatches > 0:
                # Generate variants for single mismatch
                for pos in range(length):
                    orig_base = word[pos]
                    for base in "ACGT":
                        if base != orig_base:
                            mismatch_variant = word[:pos] + base + word[pos + 1:]
                            # Check occurrence of this variant in sequences
                            matches, positions, _ = get_matched_sequences(mismatch_variant, chunk_sequences)
                            if matches:
                                mismatch_key = f"{word}(m1:{pos}>{base})"
                                sequences_with_pattern[mismatch_key].update(set(p[0] for p in positions))
                                pattern_positions[mismatch_key].extend(positions)

            # Try generating gap patterns (if enabled)
            if enable_gapped and length >= 5:  # Try gaps for longer patterns
                # Try multiple split points, not just midpoint split
                for split_point in range(1, length):  # Keep at least 2 bases on each side
                    prefix = word[:split_point]
                    suffix = word[split_point:]

                    # Try different gap lengths
                    for gap_min, gap_max in [(1, 3), (2, 4), (3, 5)]:
                        gap_pattern = f"{prefix}{{{gap_min},{gap_max}}}{suffix}"
                        matches, positions, _ = get_matched_sequences(None, chunk_sequences, gap_pattern=gap_pattern)
                        if matches and len(matches) / len(chunk_sequences) >= min_frequency:
                            gap_key = f"{prefix}N{{{gap_min},{gap_max}}}{suffix}"
                            sequences_with_pattern[gap_key].update(set(p[0] for p in positions))
                            pattern_positions[gap_key].extend(positions)

            # Try generating gap patterns using PWM method (if enabled)
            if enable_gapped and enable_pwm and enable_pwm_gapped and length >= 5:
                for split_point in range(1, length):  # Keep at least 2 bases on each side
                    prefix = word[:split_point]
                    suffix = word[split_point:]

                    # Try different gap lengths
                    for gap_min, gap_max in [(1, 3), (2, 4), (3, 5)]:
                        # Match gap pattern using PWM method
                        matches, positions, _ = get_matched_sequences_with_pwm_gaps(
                            prefix, suffix, gap_min, gap_max, chunk_sequences,
                            threshold=0.75  # PWM matching threshold
                        )

                        if matches and len(matches) / len(chunk_sequences) >= min_frequency:
                            # Use special marker to indicate this is PWM gap pattern
                            pwm_gap_key = f"{prefix}N{{{gap_min},{gap_max}}}{suffix}(pwm)"
                            sequences_with_pattern[pwm_gap_key].update(set(p[0] for p in positions))
                            pattern_positions[pwm_gap_key].extend(positions)

    # Step 3: Try PWM method (if enabled)
    if enable_pwm:
        # Create PWM for each significant pattern
        for pattern in significant_patterns:
            # Get matching sequences
            matches, _, _ = get_matched_sequences(pattern, chunk_sequences)
            if len(matches) >= 10:  # Ensure sufficient data to create PWM
                # Calculate PWM
                pwm = calculate_position_frequencies(pattern, matches)
                # Scan using PWM
                pwm_matches, pwm_positions, _ = scan_sequences_with_pwm(pwm, chunk_sequences, threshold=0.75)
                if pwm_matches:
                    pwm_key = f"PWM:{pattern}"
                    sequences_with_pattern[pwm_key].update(set(p[0] for p in pwm_positions))
                    pattern_positions[pwm_key].extend(pwm_positions)

    return dict(sequences_with_pattern), dict(local_pattern_counts), dict(pattern_positions)


class MotifAnalyzer:
    def __init__(self, max_utr_length=500, min_frequency=0.05, fold_change_threshold=1.5, min_occur_freq=0.1,
                 max_mismatches=0, enable_gapped=False, enable_pwm=False, enable_structural=False,
                 enable_pwm_gapped=False):
        """
        Initialize motif analyzer

        Purpose: Set motif analysis parameters and load sequence data
        Extension: Add advanced motif analysis options
        """
        self.up_sequences = {}  # 5'UTR sequences of upregulated genes
        self.down_sequences = {}  # 5'UTR sequences of downregulated genes
        self.gene_ids = {}  # Store gene IDs
        self.max_utr_length = max_utr_length  # Limit 5'UTR length
        self.min_frequency = min_frequency  # Minimum frequency threshold
        self.fold_change_threshold = fold_change_threshold  # Enrichment fold change threshold
        self.min_occur_freq = min_occur_freq  # Minimum occurrence frequency threshold

        # New advanced motif analysis parameters
        self.max_mismatches = max_mismatches  # Maximum allowed mismatches
        self.enable_gapped = enable_gapped  # Enable gap motif analysis
        self.enable_pwm = enable_pwm  # Enable PWM analysis
        self.enable_structural = enable_structural  # Enable structural motif analysis
        self.enable_pwm_gapped = enable_pwm_gapped

        # Add support for MANE information
        self.up_transcript_selection = {}  # Store transcript selection logic for upregulated genes
        self.down_transcript_selection = {}  # Store transcript selection logic for downregulated genes

        # Initialize structural motif storage
        self.up_structural_motifs = []
        self.down_structural_motifs = []

        # Load sequences
        self.load_sequences()

        # If structural analysis enabled, pre-calculate structural motifs
        if self.enable_structural:
            self.analyze_structural_motifs()

    def load_sequences(self):
        """
        Load sequences from FASTA files

        Purpose: Read 5'UTR sequences generated by the first script, based on main transcripts
        Improvement: Added extraction and analysis of transcript selection logic
        """
        try:
            logger.info("Loading up-regulated sequences...")
            if os.path.exists("te_up_utr.fasta"):
                for record in SeqIO.parse("te_up_utr.fasta", "fasta"):
                    # Get gene ID and sequence
                    gene_id = record.id
                    seq = str(record.seq)[:self.max_utr_length]  # Limit length to 500nt
                    self.up_sequences[gene_id] = seq
                    self.gene_ids[gene_id] = gene_id  # Save gene ID mapping

                    # Extract transcript selection logic
                    if "Selection:" in record.description:
                        selection_logic = record.description.split("Selection:")[1].split("|")[0].strip()
                        self.up_transcript_selection[gene_id] = selection_logic
            else:
                logger.error("File te_up_utr.fasta not found")

            logger.info("Loading down-regulated sequences...")
            if os.path.exists("te_down_utr.fasta"):
                for record in SeqIO.parse("te_down_utr.fasta", "fasta"):
                    gene_id = record.id
                    seq = str(record.seq)[:self.max_utr_length]
                    self.down_sequences[gene_id] = seq
                    self.gene_ids[gene_id] = gene_id

                    # Extract transcript selection logic
                    if "Selection:" in record.description:
                        selection_logic = record.description.split("Selection:")[1].split("|")[0].strip()
                        self.down_transcript_selection[gene_id] = selection_logic
            else:
                logger.error("File te_down_utr.fasta not found")

            logger.info(f"Loaded {len(self.up_sequences)} up-regulated gene sequences")
            logger.info(f"Loaded {len(self.down_sequences)} down-regulated gene sequences")

            # Summarize and visualize transcript selection
            self.summarize_transcript_selection()

        except Exception as e:
            logger.error(f"Error loading sequences: {str(e)}")
            raise

    def analyze_structural_motifs(self):
        """
        Analyze structural motifs in sequences

        Purpose: Identify possible secondary structure-related regulatory elements
        """
        logger.info("Analyzing structural motifs in sequences...")

        # Analyze structural motifs in upregulated genes
        if self.up_sequences:
            logger.info("Analyzing stem-loops in UP-regulated genes...")
            stem_loops, positions = find_structural_motifs(self.up_sequences)
            logger.info(f"Found {len(stem_loops)} potential stem-loops in UP-regulated genes")

            logger.info("Analyzing G-quadruplexes in UP-regulated genes...")
            g_quads, g_positions = find_g_quadruplexes(self.up_sequences)
            logger.info(f"Found {len(g_quads)} potential G-quadruplexes in UP-regulated genes")

            # Store results
            self.up_structural_motifs = {
                'stem_loops': {'motifs': stem_loops, 'positions': positions},
                'g_quadruplexes': {'motifs': g_quads, 'positions': g_positions}
            }

        # Analyze structural motifs in downregulated genes
        if self.down_sequences:
            logger.info("Analyzing stem-loops in DOWN-regulated genes...")
            stem_loops, positions = find_structural_motifs(self.down_sequences)
            logger.info(f"Found {len(stem_loops)} potential stem-loops in DOWN-regulated genes")

            logger.info("Analyzing G-quadruplexes in DOWN-regulated genes...")
            g_quads, g_positions = find_g_quadruplexes(self.down_sequences)
            logger.info(f"Found {len(g_quads)} potential G-quadruplexes in DOWN-regulated genes")

            # Store results
            self.down_structural_motifs = {
                'stem_loops': {'motifs': stem_loops, 'positions': positions},
                'g_quadruplexes': {'motifs': g_quads, 'positions': g_positions}
            }

        # Compare and analyze structural motifs
        self.compare_structural_motifs()

    def compare_structural_motifs(self):
        """
        Compare structural motifs between upregulated and downregulated genes

        Purpose: Identify possible differential structural regulatory elements
        """
        # Create results directory
        os.makedirs("structural_motifs", exist_ok=True)

        # Prepare comparison data
        up_counts = {
            'stem_loops': len(self.up_structural_motifs.get('stem_loops', {}).get('motifs', [])),
            'g_quadruplexes': len(self.up_structural_motifs.get('g_quadruplexes', {}).get('motifs', []))
        }

        down_counts = {
            'stem_loops': len(self.down_structural_motifs.get('stem_loops', {}).get('motifs', [])),
            'g_quadruplexes': len(self.down_structural_motifs.get('g_quadruplexes', {}).get('motifs', []))
        }

        # Calculate frequency of each structure in each gene group
        up_freq = {
            'stem_loops': up_counts['stem_loops'] / len(self.up_sequences) if self.up_sequences else 0,
            'g_quadruplexes': up_counts['g_quadruplexes'] / len(self.up_sequences) if self.up_sequences else 0
        }

        down_freq = {
            'stem_loops': down_counts['stem_loops'] / len(self.down_sequences) if self.down_sequences else 0,
            'g_quadruplexes': down_counts['g_quadruplexes'] / len(self.down_sequences) if self.down_sequences else 0
        }

        # Calculate enrichment fold change
        enrichment = {
            'stem_loops': up_freq['stem_loops'] / down_freq['stem_loops'] if down_freq['stem_loops'] > 0 else float(
                'inf'),
            'g_quadruplexes': up_freq['g_quadruplexes'] / down_freq['g_quadruplexes'] if down_freq[
                                                                                             'g_quadruplexes'] > 0 else float(
                'inf')
        }

        # Save results
        results = pd.DataFrame({
            'Structure_Type': ['Stem-Loops', 'G-Quadruplexes'],
            'UP_Count': [up_counts['stem_loops'], up_counts['g_quadruplexes']],
            'DOWN_Count': [down_counts['stem_loops'], down_counts['g_quadruplexes']],
            'UP_Frequency': [up_freq['stem_loops'], up_freq['g_quadruplexes']],
            'DOWN_Frequency': [down_freq['stem_loops'], down_freq['g_quadruplexes']],
            'UP/DOWN_Ratio': [enrichment['stem_loops'], enrichment['g_quadruplexes']]
        })

        results.to_csv("structural_motifs/structural_enrichment.csv", index=False)

        # Create visualization
        plt.figure(figsize=(10, 6))
        x = range(len(results))
        width = 0.35

        plt.bar([i - width / 2 for i in x], results['UP_Frequency'], width, label='UP regulated', color='red')
        plt.bar([i + width / 2 for i in x], results['DOWN_Frequency'], width, label='DOWN regulated', color='blue')

        plt.xlabel('Structure Type')
        plt.ylabel('Frequency per Gene')
        plt.title('Structural Motif Frequencies')
        plt.xticks(x, results['Structure_Type'])
        plt.legend()

        plt.tight_layout()
        plt.savefig("structural_motifs/structural_frequencies.png")
        plt.close()

        # Save gene lists containing structural motifs
        for motif_type in ['stem_loops', 'g_quadruplexes']:
            if self.up_structural_motifs.get(motif_type, {}).get('positions'):
                genes = set(pos[0] for pos in self.up_structural_motifs[motif_type]['positions'])
                with open(f"structural_motifs/up_{motif_type}_genes.txt", 'w') as f:
                    for gene in sorted(genes):
                        f.write(f"{gene}\n")

            if self.down_structural_motifs.get(motif_type, {}).get('positions'):
                genes = set(pos[0] for pos in self.down_structural_motifs[motif_type]['positions'])
                with open(f"structural_motifs/down_{motif_type}_genes.txt", 'w') as f:
                    for gene in sorted(genes):
                        f.write(f"{gene}\n")

    def summarize_transcript_selection(self):
        """
        Summarize and visualize the distribution of transcript selection logic

        Purpose: Understand the distribution of different transcript types (MANE Select, MANE Plus Clinical, etc.) in the dataset
        """
        # If no transcript selection information, return
        if not self.up_transcript_selection and not self.down_transcript_selection:
            logger.warning("No transcript selection information found in FASTA descriptions.")
            return

        up_selection_counts = {}
        down_selection_counts = {}

        # Count upregulated gene transcript selection logic
        for selection in self.up_transcript_selection.values():
            up_selection_counts[selection] = up_selection_counts.get(selection, 0) + 1

        # Count downregulated gene transcript selection logic
        for selection in self.down_transcript_selection.values():
            down_selection_counts[selection] = down_selection_counts.get(selection, 0) + 1

        logger.info("Transcript selection summary for UP-regulated genes:")
        for selection, count in sorted(up_selection_counts.items()):
            percentage = count / len(self.up_transcript_selection) * 100 if self.up_transcript_selection else 0
            logger.info(f"  {selection}: {count} ({percentage:.1f}%)")

        logger.info("Transcript selection summary for DOWN-regulated genes:")
        for selection, count in sorted(down_selection_counts.items()):
            percentage = count / len(self.down_transcript_selection) * 100 if self.down_transcript_selection else 0
            logger.info(f"  {selection}: {count} ({percentage:.1f}%)")

        # Create folder to save results
        os.makedirs("transcript_selection", exist_ok=True)

        # Save as CSV
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
            selection_df = pd.DataFrame(selection_data)
            selection_df.to_csv("transcript_selection/transcript_selection_summary.csv", index=False)

            # Create visualization
            plt.figure(figsize=(12, 6))

            # Prepare data
            all_selections = set(self.up_transcript_selection.values()) | set(self.down_transcript_selection.values())

            # Sort (ensure MANE-related ones are in front)
            selection_order = sorted(all_selections)
            if "MANE_Select" in selection_order:
                selection_order.remove("MANE_Select")
                selection_order.insert(0, "MANE_Select")
            if "MANE_Plus_Clinical" in selection_order:
                selection_order.remove("MANE_Plus_Clinical")
                selection_order.insert(1, "MANE_Plus_Clinical")

            up_counts = [up_selection_counts.get(s, 0) for s in selection_order]
            down_counts = [down_selection_counts.get(s, 0) for s in selection_order]

            x = range(len(selection_order))
            width = 0.35

            # Plot
            plt.bar([i - width / 2 for i in x], up_counts, width, label='UP regulated', color='red')
            plt.bar([i + width / 2 for i in x], down_counts, width, label='DOWN regulated', color='blue')

            plt.xlabel('Transcript Selection Method')
            plt.ylabel('Number of Genes')
            plt.title('Transcript Selection Methods Distribution')
            plt.xticks(x, selection_order, rotation=45, ha='right')
            plt.legend()

            plt.tight_layout()
            plt.savefig("transcript_selection/transcript_selection_distribution.png")
            plt.close()

    def find_motifs(self, sequences, length, min_frequency=None, n_cores=None, max_wildcards=1):
        """
        DREME-style motif search

        Purpose: Find potential regulatory motifs in sequence collections
        Extension: Add support for mismatches, gaps and PWM
        """
        if min_frequency is None:
            min_frequency = self.min_frequency

        if n_cores is None:
            n_cores = max(1, multiprocessing.cpu_count() - 1)  # Reserve one core for system

        sequence_chunks = list(chunk_sequences(sequences, n_cores))
        args_list = [(chunk, length, min_frequency, max_wildcards, self.max_mismatches,
                      self.enable_gapped, self.enable_pwm, self.enable_pwm_gapped) for chunk in
                     sequence_chunks]  # Modified this line

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

        return dict(all_sequences_with_pattern), dict(pattern_counts), dict(all_pattern_positions)

    def calculate_and_save_pattern_info(self, pattern, sequences, output_prefix, pattern_positions=None,
                                        regulation_type=None, mismatch_details=None):
        """
        Analyze and save pattern sequence information

        Purpose: Generate detailed reports for each significant motif, including position frequency, strand specificity and transcript type distribution
        Improvements:
        1. Added analysis of motif distribution in different transcript types
        2. Added support for mismatch analysis
        3. Added analysis of structural correlation
        """
        matches, positions, _ = get_matched_sequences(pattern, sequences, max_mismatches=self.max_mismatches)
        if not matches:
            return None

        # Calculate position frequencies
        freq_matrix = calculate_position_frequencies(pattern, matches)

        # Save frequency matrix
        freq_matrix.to_csv(f'{output_prefix}_frequencies.csv')

        # Calculate strand specificity
        strand_specificity, p_value = calculate_strand_specificity(pattern, sequences)

        # Analyze mismatch situation (if provided)
        mismatch_analysis = None
        if mismatch_details:
            # Count mismatches at each position
            pos_mismatch_counts = defaultdict(Counter)
            for _, _, mismatches in mismatch_details:
                for pos, orig, new in mismatches:
                    pos_mismatch_counts[pos][new] += 1

            mismatch_analysis = dict(pos_mismatch_counts)

        # Analyze motif distribution in different transcript types
        transcript_selection_dist = None
        if regulation_type and pattern_positions:
            matching_genes = set(p[0] for p in positions)
            selection_dict = self.up_transcript_selection if regulation_type == 'UP' else self.down_transcript_selection

            # Only analyze when transcript selection information is available
            if selection_dict:
                # Count number of genes containing this motif in different transcript types
                selection_counts = {}
                for gene in matching_genes:
                    if gene in selection_dict:
                        selection = selection_dict[gene]
                        selection_counts[selection] = selection_counts.get(selection, 0) + 1

                # Create transcript type distribution results
                transcript_selection_dist = selection_counts

                # Generate transcript type distribution plot
                if selection_counts:
                    plt.figure(figsize=(10, 6))

                    # Prepare data
                    labels = list(selection_counts.keys())
                    values = list(selection_counts.values())

                    # Sort
                    sorted_indices = sorted(range(len(labels)), key=lambda i: values[i], reverse=True)
                    labels = [labels[i] for i in sorted_indices]
                    values = [values[i] for i in sorted_indices]

                    # Create pie chart
                    plt.pie(values, labels=labels, autopct='%1.1f%%', startangle=90)
                    plt.axis('equal')  # Ensure pie chart is circular
                    plt.title(
                        f'Distribution of {pattern} in Different Transcript Types\n({regulation_type} regulated genes)')
                    plt.savefig(f'{output_prefix}_transcript_dist.png')
                    plt.close()

        # Analyze relationship with structural motifs (if enabled)
        structure_correlation = None
        if self.enable_structural and pattern_positions:
            matching_genes = set(p[0] for p in positions)

            # Get corresponding structural motifs
            structural_motifs = self.up_structural_motifs if regulation_type == 'UP' else self.down_structural_motifs

            # Calculate co-occurrence with different structural motifs
            structure_correlation = {}

            for motif_type in ['stem_loops', 'g_quadruplexes']:
                if motif_type in structural_motifs:
                    structural_genes = set(pos[0] for pos in structural_motifs[motif_type]['positions'])

                    # Calculate co-occurring genes
                    overlap_genes = matching_genes.intersection(structural_genes)

                    # Calculate hypergeometric test p-value
                    total_genes = len(sequences)
                    motif_genes = len(matching_genes)
                    struct_genes = len(structural_genes)
                    overlap_count = len(overlap_genes)

                    # Calculate p-value using hypergeometric distribution
                    p_value = stats.hypergeom.sf(overlap_count - 1, total_genes, struct_genes, motif_genes)

                    structure_correlation[motif_type] = {
                        'overlap_count': overlap_count,
                        'total_motif_genes': motif_genes,
                        'total_struct_genes': struct_genes,
                        'p_value': p_value,
                        'overlap_genes': list(overlap_genes)
                    }

        # Generate detailed text report
        with open(f'{output_prefix}_details.txt', 'w') as f:
            f.write(f'Pattern: {pattern}\n')
            f.write(f'Number of matching genes: {len(set(p[0] for p in positions))}\n')
            f.write(f'Total matches: {len(matches)}\n')
            f.write(f'Strand specificity: {strand_specificity:.2f} (p-value: {p_value:.6f})\n\n')

            # Explain meaning of strand specificity
            f.write("Strand Specificity Interpretation:\n")
            if strand_specificity > 3:
                f.write("- Very high strand specificity (>3): Strongly suggests RNA-level regulation\n")
            elif strand_specificity > 1.5:
                f.write("- High strand specificity (1.5-3): Likely RNA-level regulatory element\n")
            elif strand_specificity >= 0.8 and strand_specificity <= 1.2:
                f.write("- Near equal strand distribution (0.8-1.2): Possible DNA-level regulatory element\n")
            elif strand_specificity < 0.7:
                f.write("- Low strand specificity (<0.7): Complementary strand bias, unusual for 5'UTR elements\n")

            # Add mismatch analysis information
            if mismatch_analysis:
                f.write('\nMismatch Analysis:\n')
                for pos, base_counts in sorted(mismatch_analysis.items()):
                    f.write(f'Position {pos}: Original base {pattern[pos]}, Mismatches: ')
                    for base, count in base_counts.most_common():
                        f.write(f'{base}({count}) ')
                    f.write('\n')

            # Add transcript type distribution information
            if transcript_selection_dist:
                f.write('\nDistribution across transcript selection methods:\n')
                for selection, count in sorted(transcript_selection_dist.items(), key=lambda x: x[1], reverse=True):
                    percentage = count / len(set(p[0] for p in positions)) * 100
                    f.write(f'- {selection}: {count} genes ({percentage:.1f}%)\n')

            # Add structural correlation information
            if structure_correlation:
                f.write('\nStructural Correlation Analysis:\n')
                for motif_type, data in structure_correlation.items():
                    f.write(f'- {motif_type.replace("_", " ").title()}:\n')
                    f.write(f'  Overlap: {data["overlap_count"]} genes\n')
                    f.write(f'  P-value: {data["p_value"]:.6f}\n')
                    if data["p_value"] < 0.05:
                        f.write(f'  Significant correlation detected (p<0.05)\n')

                    # List top 10 co-occurring genes
                    top_genes = sorted(data["overlap_genes"])[:10]
                    if top_genes:
                        f.write(f'  Example overlap genes: {", ".join(top_genes)}\n')

            f.write('\nPosition frequencies:\n')
            f.write(str(freq_matrix))

            # Add position specificity information
            if pattern_positions:
                f.write('\n\nSequence match positions:\n')
                pos_counter = defaultdict(int)
                for seq_id, pos in pattern_positions[:100]:  # Only show first 100 to prevent file being too large
                    pos_counter[pos] += 1
                    f.write(f'{seq_id}: position {pos}\n')

                # Add position distribution statistics
                f.write('\nPosition distribution:\n')
                for pos, count in sorted(pos_counter.items()):
                    f.write(f'Position {pos}: {count} matches\n')

            # List of matching genes - important for gene function enrichment analysis
            f.write('\n\nMatching genes:\n')
            matching_genes = sorted(set(p[0] for p in positions))
            f.write(', '.join(matching_genes))

            # Write matching gene list to separate file for subsequent GO/KEGG enrichment analysis
            with open(f'{output_prefix}_matching_genes.txt', 'w') as gene_file:
                for gene in matching_genes:
                    gene_file.write(f"{gene}\n")

            f.write('\n\nExample matching sequences:\n')
            for seq in matches[:10]:  # Show first 10 matching sequences
                f.write(f'{seq}\n')

        # Create position distribution plot
        if pattern_positions:
            pos_data = [pos for _, pos in pattern_positions]
            plt.figure(figsize=(10, 5))
            plt.hist(pos_data, bins=30, alpha=0.7)
            plt.title(f'Position Distribution of Pattern {pattern}')
            plt.xlabel('Position in 5\'UTR')
            plt.ylabel('Frequency')
            plt.grid(True, alpha=0.3)
            plt.savefig(f'{output_prefix}_position_dist.png')
            plt.close()

        # Plot position frequency graph
        plt.figure(figsize=(10, 5))
        sns.heatmap(freq_matrix, cmap='YlGnBu', annot=True, fmt='.2f')
        plt.title(f'Position Frequencies for Pattern {pattern}')
        plt.savefig(f'{output_prefix}_position_freq.png')
        plt.close()

        # If there's mismatch analysis, plot mismatch heatmap
        if mismatch_analysis:
            # Prepare data
            rows = ['A', 'C', 'G', 'T']
            cols = list(range(len(pattern)))
            mismatch_data = np.zeros((len(rows), len(cols)))

            for pos, base_counts in mismatch_analysis.items():
                for base, count in base_counts.items():
                    if base in rows:
                        row_idx = rows.index(base)
                        mismatch_data[row_idx, pos] = count

            # Plot heatmap
            plt.figure(figsize=(10, 5))
            sns.heatmap(mismatch_data, cmap='YlOrRd', annot=True, fmt='g',
                        xticklabels=cols, yticklabels=rows)
            plt.title(f'Mismatch Distribution for Pattern {pattern}')
            plt.xlabel('Position')
            plt.ylabel('Mismatched Base')
            plt.savefig(f'{output_prefix}_mismatch_dist.png')
            plt.close()

        # If there's structural correlation, create network diagram
        if structure_correlation:
            significant_correlations = {k: v for k, v in structure_correlation.items() if v['p_value'] < 0.05}
            if significant_correlations:
                plt.figure(figsize=(8, 6))
                G = nx.Graph()

                # Add nodes
                G.add_node(pattern, type='sequence_motif')
                for motif_type in significant_correlations.keys():
                    G.add_node(motif_type, type='structural_motif')

                # Add edges
                for motif_type, data in significant_correlations.items():
                    weight = -np.log10(data['p_value'])  # Set edge thickness based on p-value
                    G.add_edge(pattern, motif_type, weight=weight,
                               overlap=data['overlap_count'], pvalue=data['p_value'])

                # Set node colors
                node_colors = ['red' if G.nodes[n]['type'] == 'sequence_motif' else 'blue'
                               for n in G.nodes()]

                # Set edge widths
                edge_widths = [G[u][v]['weight'] * 2 for u, v in G.edges()]

                # Draw network
                pos = nx.spring_layout(G)
                nx.draw_networkx(G, pos, node_color=node_colors, width=edge_widths,
                                 font_size=10, node_size=500, font_weight='bold')

                # Add edge labels
                edge_labels = {(u, v): f"Overlap: {G[u][v]['overlap']}\np:{G[u][v]['pvalue']:.1e}"
                               for u, v in G.edges()}
                nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

                plt.title(f'Structural Correlations for {pattern}')
                plt.axis('off')
                plt.savefig(f'{output_prefix}_structural_network.png')
                plt.close()

        return freq_matrix, strand_specificity, p_value

    def cluster_similar_motifs(self, motifs, similarity_threshold=0.5):
        """
        Cluster similar patterns

        Purpose: Group similar motifs into families to reduce redundant results
        Improvement: Consider mismatches and gap patterns
        """
        clusters = []
        assigned = set()

        # Sort motifs by frequency or other metrics
        sorted_motifs = sorted(motifs.keys(), key=lambda x: motifs[x], reverse=True)

        for motif in sorted_motifs:
            # Skip assigned and specially formatted patterns
            if motif in assigned or motif.startswith("PWM:"):
                continue

            # Create new cluster
            cluster = [motif]
            assigned.add(motif)

            # Find similar patterns
            for other in sorted_motifs:
                if other != motif and other not in assigned and not other.startswith("PWM:"):
                    # Handle gap patterns
                    if '{' in other or '{' in motif:
                        # Extract prefix and suffix of gap patterns for comparison
                        motif_parts = re.split(r'\{.*?\}', motif) if '{' in motif else [motif]
                        other_parts = re.split(r'\{.*?\}', other) if '{' in other else [other]

                        # Remove PWM gap pattern markers
                        if motif_parts and len(motif_parts) > 1 and '(pwm)' in motif_parts[-1]:
                            motif_parts[-1] = motif_parts[-1].replace('(pwm)', '')
                        if other_parts and len(other_parts) > 1 and '(pwm)' in other_parts[-1]:
                            other_parts[-1] = other_parts[-1].replace('(pwm)', '')

                        # Rest of similarity calculation code remains unchanged...
                        # If prefix and suffix are similar, consider patterns similar
                        prefix_similarity = 0
                        if motif_parts[0] and other_parts[0]:
                            prefix_len = min(len(motif_parts[0]), len(other_parts[0]))
                            prefix_similarity = sum(
                                motif_parts[0][i] == other_parts[0][i] for i in range(prefix_len)) / prefix_len

                        suffix_similarity = 0
                        if len(motif_parts) > 1 and len(other_parts) > 1 and motif_parts[-1] and other_parts[-1]:
                            suffix_len = min(len(motif_parts[-1]), len(other_parts[-1]))
                            suffix_similarity = sum(
                                motif_parts[-1][i] == other_parts[-1][i] for i in range(suffix_len)) / suffix_len

                        combined_similarity = (prefix_similarity + suffix_similarity) / 2

                        if combined_similarity >= similarity_threshold:
                            cluster.append(other)
                            assigned.add(other)
                    # Handle mismatch patterns
                    elif '(' in other or '(' in motif:
                        # Extract basic pattern
                        base_motif = motif.split('(')[0] if '(' in motif else motif
                        base_other = other.split('(')[0] if '(' in other else other

                        # Calculate similarity of basic patterns
                        min_len = min(len(base_motif), len(base_other))
                        similarity = sum(base_motif[i] == base_other[i] for i in range(min_len)) / min_len

                        if similarity >= similarity_threshold:
                            cluster.append(other)
                            assigned.add(other)
                    # Regular pattern similarity calculation
                    else:
                        # Simple similarity calculation: proportion of bases at same positions
                        similarity = 0
                        for i in range(min(len(motif), len(other))):
                            if motif[i] == other[i] or (motif[i] in "RYSWKMBDHVN" and other[i] in "RYSWKMBDHVN"):
                                similarity += 1

                        similarity /= min(len(motif), len(other))
                        if similarity >= similarity_threshold:
                            cluster.append(other)
                            assigned.add(other)

            clusters.append(cluster)

        return clusters

    def analyze_single_length(self, length, n_cores=None, max_wildcards=1):
        """
        Analyze motifs based on specific length, using frequency and enrichment fold change to determine significance

        Purpose: Analyze motifs of specific length, identify significantly different motifs based on frequency and enrichment fold change
        Improvements:
        1. Integrate mismatch, gap and PWM analysis
        2. Add structural correlation analysis
        """
        # Create results directory
        result_dir = f'motif_results_length_{length}'
        os.makedirs(result_dir, exist_ok=True)

        # Get motifs
        up_patterns, up_counts, up_positions = self.find_motifs(
            self.up_sequences, length, n_cores=n_cores, max_wildcards=max_wildcards
        )
        down_patterns, down_counts, down_positions = self.find_motifs(
            self.down_sequences, length, n_cores=n_cores, max_wildcards=max_wildcards
        )

        # Analyze all patterns
        results = []
        all_patterns = set(up_patterns) | set(down_patterns)

        logger.info(f"Analyzing {len(all_patterns)} patterns")
        for pattern in all_patterns:
            up_present = len(up_patterns.get(pattern, set()))
            down_present = len(down_patterns.get(pattern, set()))
            up_absent = len(self.up_sequences) - up_present
            down_absent = len(self.down_sequences) - down_present

            # Calculate frequency
            up_frequency = up_present / len(self.up_sequences) if len(self.up_sequences) > 0 else 0
            down_frequency = down_present / len(self.down_sequences) if len(self.down_sequences) > 0 else 0

            # Calculate enrichment fold change
            if down_frequency > 0:
                up_to_down_ratio = up_frequency / down_frequency
            else:
                up_to_down_ratio = float('inf')

            if up_frequency > 0:
                down_to_up_ratio = down_frequency / up_frequency
            else:
                down_to_up_ratio = float('inf')

            # Calculate strand specificity - separately for upregulated and downregulated genes
            up_strand_spec, up_strand_pval = (0, 1)
            down_strand_spec, down_strand_pval = (0, 1)

            if up_present > 0:
                up_strand_spec, up_strand_pval = calculate_strand_specificity(pattern, self.up_sequences)
            if down_present > 0:
                down_strand_spec, down_strand_pval = calculate_strand_specificity(pattern, self.down_sequences)

            # Calculate overall strand specificity
            combined_sequences = {**self.up_sequences, **self.down_sequences}
            strand_specificity, strand_pval = calculate_strand_specificity(pattern, combined_sequences)

            # Analyze characteristics of this pattern
            is_mismatch = '(' in pattern  # Mismatch pattern
            is_gapped = '{' in pattern  # Gap pattern
            is_pwm = pattern.startswith('PWM:')  # PWM pattern

            # Determine pattern category
            if is_mismatch:
                pattern_type = 'mismatch'
            elif is_gapped:
                pattern_type = 'gapped'
            elif is_pwm:
                pattern_type = 'pwm'
            else:
                pattern_type = 'exact'

            # Keep Fisher test calculation, though not used for significance determination
            oddsratio, pvalue = stats.fisher_exact([
                [up_present, up_absent],
                [down_present, down_absent]
            ])

            results.append({
                'Pattern': pattern,
                'Length': length,
                'Type': pattern_type,
                'UP_present': up_present,
                'UP_absent': up_absent,
                'DOWN_present': down_present,
                'DOWN_absent': down_absent,
                'UP_frequency': up_frequency,
                'DOWN_frequency': down_frequency,
                'UP_to_DOWN_ratio': up_to_down_ratio,
                'DOWN_to_UP_ratio': down_to_up_ratio,
                'P_value': pvalue,  # Keep but not used for significance determination
                'Odds_ratio': oddsratio,  # Keep but not used for significance determination
                'UP_strand_specificity': up_strand_spec,
                'UP_strand_pvalue': up_strand_pval,
                'DOWN_strand_specificity': down_strand_spec,
                'DOWN_strand_pvalue': down_strand_pval,
                'Overall_strand_specificity': strand_specificity,
                'Overall_strand_pvalue': strand_pval
            })

        # Save results
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            # Determine significance based on new criteria
            results_df['Significant'] = False
            results_df['Enrichment'] = 'None'

            # Upregulation enriched patterns: UP frequency>min_occur_freq AND UP/DOWN>fold_change_threshold
            up_enriched_mask = (results_df['UP_frequency'] >= self.min_occur_freq) & \
                               (results_df['UP_to_DOWN_ratio'] >= self.fold_change_threshold)
            results_df.loc[up_enriched_mask, 'Significant'] = True
            results_df.loc[up_enriched_mask, 'Enrichment'] = 'UP'

            # Downregulation enriched patterns: DOWN frequency>min_occur_freq AND DOWN/UP>fold_change_threshold
            down_enriched_mask = (results_df['DOWN_frequency'] >= self.min_occur_freq) & \
                                 (results_df['DOWN_to_UP_ratio'] >= self.fold_change_threshold)
            results_df.loc[down_enriched_mask, 'Significant'] = True
            results_df.loc[down_enriched_mask, 'Enrichment'] = 'DOWN'

            # Save all results
            results_df.to_csv(f'{result_dir}/01motif_analysis.csv', index=False)

            # Save significant patterns
            significant_df = results_df[results_df['Significant']]
            significant_df.to_csv(f'{result_dir}/02significant_motifs.csv', index=False)

            # Add: Create CSV file containing detailed gene information
            if not significant_df.empty:
                # Need to import json module
                import json
                self.create_motif_details_file(significant_df, up_positions, down_positions, result_dir)

            # Statistics by type classification
            type_counts = results_df.groupby(['Type', 'Significant']).size().unstack(fill_value=0)
            if not type_counts.empty:
                type_counts.to_csv(f'{result_dir}/motif_types_summary.csv')

                # Create type distribution plot
                plt.figure(figsize=(10, 6))
                type_counts.plot(kind='bar', stacked=True)
                plt.title(f'Motif Types Distribution (Length {length})')
                plt.xlabel('Motif Type')
                plt.ylabel('Count')
                plt.savefig(f'{result_dir}/motif_types_distribution.png')
                plt.close()

            # Save upregulated and downregulated enriched motifs separately
            up_enriched = significant_df[significant_df['Enrichment'] == 'UP']
            down_enriched = significant_df[significant_df['Enrichment'] == 'DOWN']

            if not up_enriched.empty:
                up_enriched = up_enriched.sort_values('UP_to_DOWN_ratio', ascending=False)
                up_enriched.to_csv(f'{result_dir}/02up_enriched_motifs.csv', index=False)
                logger.info(f"Found {len(up_enriched)} UP-enriched motifs")

            if not down_enriched.empty:
                down_enriched = down_enriched.sort_values('DOWN_to_UP_ratio', ascending=False)
                down_enriched.to_csv(f'{result_dir}/02down_enriched_motifs.csv', index=False)
                logger.info(f"Found {len(down_enriched)} DOWN-enriched motifs")

            # Cluster similar motifs
            if not significant_df.empty:
                # Cluster separately for upregulated and downregulated
                if not up_enriched.empty:
                    # Create dictionary, using UP_to_DOWN_ratio as value
                    up_motifs_dict = dict(zip(up_enriched['Pattern'], up_enriched['UP_to_DOWN_ratio']))
                    self.process_enriched_motifs(up_enriched, self.up_sequences, up_positions, result_dir, "up",
                                                 up_motifs_dict)

                if not down_enriched.empty:
                    # Create dictionary, using DOWN_to_UP_ratio as value
                    down_motifs_dict = dict(zip(down_enriched['Pattern'], down_enriched['DOWN_to_UP_ratio']))
                    self.process_enriched_motifs(down_enriched, self.down_sequences, down_positions, result_dir, "down",
                                                 down_motifs_dict)

            # Create data visualization
            self.create_visualizations(results_df, significant_df, length, result_dir)

        return results_df

    def process_enriched_motifs(self, enriched_df, sequences, positions, result_dir, direction, motifs_dict):
        """
        Process enriched motifs, generate clustering and detailed information

        Purpose: Cluster enriched motifs and generate detailed reports
        Improvement: Support clustering of mismatch and gap patterns
        """
        # Separate different types of motifs
        exact_motifs = {p: v for p, v in motifs_dict.items()
                        if '(' not in p and '{' not in p and not p.startswith('PWM:')}
        mismatch_motifs = {p: v for p, v in motifs_dict.items() if '(' in p}
        gapped_motifs = {p: v for p, v in motifs_dict.items() if '{' in p}
        pwm_motifs = {p: v for p, v in motifs_dict.items() if p.startswith('PWM:')}

        # Cluster each type
        exact_clusters = self.cluster_similar_motifs(exact_motifs)
        mismatch_clusters = self.cluster_similar_motifs(mismatch_motifs)
        gapped_clusters = self.cluster_similar_motifs(gapped_motifs)

        # Merge all clusters
        all_clusters = exact_clusters + mismatch_clusters + gapped_clusters
        # Handle PWM separately - each is independent
        for pwm_motif in pwm_motifs:
            all_clusters.append([pwm_motif])

        # Save clustering results
        cluster_data = []
        for i, cluster in enumerate(all_clusters):
            for motif in cluster:
                row = enriched_df[enriched_df['Pattern'] == motif].iloc[0].to_dict()
                row['Cluster'] = i + 1
                row['Cluster_representative'] = cluster[0]  # First one as representative
                if '(' in motif:
                    row['Motif_type'] = 'Mismatch'
                elif '{' in motif:
                    row['Motif_type'] = 'Gapped'
                elif motif.startswith('PWM:'):
                    row['Motif_type'] = 'PWM'
                else:
                    row['Motif_type'] = 'Exact'
                cluster_data.append(row)

        if cluster_data:
            cluster_df = pd.DataFrame(cluster_data)
            cluster_df.to_csv(f'{result_dir}/03{direction}_clustered_motifs.csv', index=False)

            # Create detailed information for representative motif of each cluster
            representatives = [cluster[0] for cluster in all_clusters]
            for rep in representatives:
                output_prefix = f'{result_dir}/{direction}_motif_{rep}'
                # Pass detailed information to calculate_and_save_pattern_info
                mismatch_details = None
                if '(' in rep:
                    # Extract mismatch information
                    mismatch_details = []
                    # Here need to extract detailed information based on actual mismatch format

                self.calculate_and_save_pattern_info(
                    rep, sequences, output_prefix,
                    pattern_positions=positions.get(rep, []),
                    regulation_type=direction.upper(),
                    mismatch_details=mismatch_details
                )

    def create_visualizations(self, results_df, significant_df, length, result_dir):
        """
        Create data visualization

        Purpose: Generate graphical representation of results to make analysis results more intuitive
        Improvement: Add visualization support for different motif types
        """
        if not significant_df.empty:
            # 1. Frequency distribution plot
            plt.figure(figsize=(12, 6))
            plt.subplot(1, 2, 1)
            sns.histplot(data=results_df, x='UP_frequency', bins=30, color='red', alpha=0.6)
            plt.axvline(x=self.min_occur_freq, color='r', linestyle='--', label=f'Min Frequency={self.min_occur_freq}')
            plt.title(f'UP Frequency Distribution (Length {length})')
            plt.xlabel('UP Frequency')
            plt.ylabel('Count')
            plt.legend()

            plt.subplot(1, 2, 2)
            sns.histplot(data=results_df, x='DOWN_frequency', bins=30, color='blue', alpha=0.6)
            plt.axvline(x=self.min_occur_freq, color='b', linestyle='--', label=f'Min Frequency={self.min_occur_freq}')
            plt.title(f'DOWN Frequency Distribution (Length {length})')
            plt.xlabel('DOWN Frequency')
            plt.ylabel('Count')
            plt.legend()

            plt.tight_layout()
            plt.savefig(f'{result_dir}/01frequency_distribution.png')
            plt.close()

            # 2. Enrichment fold change distribution plot
            plt.figure(figsize=(12, 6))
            plt.subplot(1, 2, 1)
            # Limit maximum value for plotting
            up_to_down = results_df['UP_to_DOWN_ratio'].apply(lambda x: min(x, 10) if x != float('inf') else 10)
            sns.histplot(up_to_down, bins=30, color='red', alpha=0.6)
            plt.axvline(x=self.fold_change_threshold, color='r', linestyle='--',
                        label=f'Threshold={self.fold_change_threshold}')
            plt.title(f'UP/DOWN Fold Change Distribution (Length {length})')
            plt.xlabel('UP/DOWN Ratio (capped at 10)')
            plt.ylabel('Count')
            plt.legend()

            plt.subplot(1, 2, 2)
            # Limit maximum value for plotting
            down_to_up = results_df['DOWN_to_UP_ratio'].apply(lambda x: min(x, 10) if x != float('inf') else 10)
            sns.histplot(down_to_up, bins=30, color='blue', alpha=0.6)
            plt.axvline(x=self.fold_change_threshold, color='b', linestyle='--',
                        label=f'Threshold={self.fold_change_threshold}')
            plt.title(f'DOWN/UP Fold Change Distribution (Length {length})')
            plt.xlabel('DOWN/UP Ratio (capped at 10)')
            plt.ylabel('Count')
            plt.legend()

            plt.tight_layout()
            plt.savefig(f'{result_dir}/02fold_change_distribution.png')
            plt.close()

            # 3. Scatter plot: frequency vs enrichment fold change
            plt.figure(figsize=(12, 10))

            # UP enrichment scatter plot
            plt.subplot(2, 1, 1)

            # Use different colors for different types of motifs
            type_colors = {'exact': 'red', 'mismatch': 'orange', 'gapped': 'green', 'pwm': 'purple'}
            default_color = 'gray'

            # Create mapping of type to color
            color_map = results_df['Type'].map(lambda t: type_colors.get(t, default_color))

            # Plot scatter plot and color by type
            plt.scatter(
                results_df['UP_frequency'],
                results_df['UP_to_DOWN_ratio'].apply(lambda x: min(x, 10) if x != float('inf') else 10),
                c=color_map,
                alpha=0.7
            )

            # Add legend
            for t, c in type_colors.items():
                if t in results_df['Type'].values:
                    plt.scatter([], [], c=c, label=f'{t.capitalize()} motif')

            # Add threshold lines
            plt.axhline(y=self.fold_change_threshold, color='r', linestyle='--',
                        label=f'Fold Change={self.fold_change_threshold}')
            plt.axvline(x=self.min_occur_freq, color='r', linestyle='--', label=f'Min Frequency={self.min_occur_freq}')

            # Add background color in significant region
            plt.axvspan(self.min_occur_freq, 1,
                        ymin=0, ymax=self.fold_change_threshold / 10,
                        alpha=0.1, color='red', label='Significant UP Region')

            plt.title(f'UP Frequency vs UP/DOWN Fold Change (Length {length})')
            plt.xlabel('UP Frequency')
            plt.ylabel('UP/DOWN Ratio (capped at 10)')
            plt.ylim(0, 10)
            plt.legend()

            # DOWN enrichment scatter plot
            plt.subplot(2, 1, 2)
            # Create mapping of type to color
            color_map = results_df['Type'].map(lambda t: type_colors.get(t, default_color))

            plt.scatter(
                results_df['DOWN_frequency'],
                results_df['DOWN_to_UP_ratio'].apply(lambda x: min(x, 10) if x != float('inf') else 10),
                c=color_map,
                alpha=0.7
            )

            # Add legend
            for t, c in type_colors.items():
                if t in results_df['Type'].values:
                    plt.scatter([], [], c=c, label=f'{t.capitalize()} motif')

            # Add threshold lines
            plt.axhline(y=self.fold_change_threshold, color='b', linestyle='--',
                        label=f'Fold Change={self.fold_change_threshold}')
            plt.axvline(x=self.min_occur_freq, color='b', linestyle='--', label=f'Min Frequency={self.min_occur_freq}')

            # Add background color in significant region
            plt.axvspan(self.min_occur_freq, 1,
                        ymin=0, ymax=self.fold_change_threshold / 10,
                        alpha=0.1, color='blue', label='Significant DOWN Region')

            plt.title(f'DOWN Frequency vs DOWN/UP Fold Change (Length {length})')
            plt.xlabel('DOWN Frequency')
            plt.ylabel('DOWN/UP Ratio (capped at 10)')
            plt.ylim(0, 10)
            plt.legend()

            plt.tight_layout()
            plt.savefig(f'{result_dir}/03enrichment_scatter.png')
            plt.close()

            # 4. Heatmap: frequency of enriched motifs
            if len(significant_df) <= 30:  # Limit heatmap size
                plt.figure(figsize=(12, max(8, len(significant_df) * 0.3)))

                # Prepare heatmap data
                heatmap_data = significant_df.copy()
                heatmap_data['Pattern'] = heatmap_data.apply(
                    lambda row: f"{row['Pattern']} ({row['Enrichment']}|{row['Type']})", axis=1
                )

                # Sort by enrichment direction and fold change
                heatmap_data = heatmap_data.sort_values(['Enrichment', 'UP_to_DOWN_ratio'],
                                                        ascending=[True, False])

                # Plot heatmap
                sns.heatmap(
                    heatmap_data.set_index('Pattern')[['UP_frequency', 'DOWN_frequency']],
                    annot=True,
                    cmap='YlGnBu',
                    fmt=".3f"
                )

                plt.title(f'Motif Frequencies in UP vs DOWN Regulated Genes (Length {length})')
                plt.tight_layout()
                plt.savefig(f'{result_dir}/04frequency_heatmap.png')
                plt.close()

            # 5. Strand specificity distribution plot - by enrichment direction
            plt.figure(figsize=(10, 6))

            # Separate data
            for direction, color in {'UP': 'red', 'DOWN': 'blue'}.items():
                subset = significant_df[significant_df['Enrichment'] == direction]
                if not subset.empty:
                    sns.kdeplot(
                        data=subset['Overall_strand_specificity'].apply(
                            lambda x: min(x, 5) if x != float('inf') else 5
                        ),
                        color=color,
                        label=f"{direction} enriched motifs"
                    )

            plt.axvline(x=1, color='black', linestyle='--', label='Equal strand distribution')
            plt.title(f'Distribution of Strand Specificity by Enrichment Direction (Length {length})')
            plt.xlabel('Strand Specificity (capped at 5)')
            plt.ylabel('Density')
            plt.legend()
            plt.savefig(f'{result_dir}/05strand_specificity_by_direction.png')
            plt.close()

            # 6. Analyze significance proportion by motif type
            if not significant_df.empty and 'Type' in significant_df.columns:
                # Calculate proportion of significant motifs for each type
                type_sig_ratio = {}
                for motif_type in results_df['Type'].unique():
                    type_total = sum(results_df['Type'] == motif_type)
                    if type_total > 0:
                        type_sig = sum((significant_df['Type'] == motif_type))
                        type_sig_ratio[motif_type] = type_sig / type_total

                if type_sig_ratio:
                    plt.figure(figsize=(8, 6))
                    plt.bar(type_sig_ratio.keys(), type_sig_ratio.values(), color='teal')
                    plt.title(f'Proportion of Significant Motifs by Type (Length {length})')
                    plt.xlabel('Motif Type')
                    plt.ylabel('Significant Proportion')
                    plt.ylim(0, 1)

                    # Add value labels
                    for x, y in type_sig_ratio.items():
                        plt.text(x, y + 0.02, f'{y:.2f}', ha='center')

                    plt.tight_layout()
                    plt.savefig(f'{result_dir}/06motif_type_significance.png')
                    plt.close()

            # 7. [New] Analyze motif distribution by transcript type
            self.analyze_motifs_by_transcript_type(significant_df, length, result_dir)

    def analyze_motifs_by_transcript_type(self, significant_df, length, result_dir):
        """
        Analyze motif distribution in different transcript types

        Purpose: Reveal potential transcript type-specific motif patterns
        """
        # Only execute when transcript selection information exists
        if not self.up_transcript_selection and not self.down_transcript_selection:
            return

        # Collect all transcript types
        all_transcript_types = set(self.up_transcript_selection.values()) | set(self.down_transcript_selection.values())
        if not all_transcript_types:
            return

        # Create transcript type analysis directory
        trans_dir = f"{result_dir}/transcript_type_analysis"
        os.makedirs(trans_dir, exist_ok=True)

        # Analyze upregulated genes
        up_enriched = significant_df[significant_df['Enrichment'] == 'UP']
        if not up_enriched.empty and self.up_transcript_selection:
            # Transcript type dictionary: key is type, value is list of genes containing this type
            transcript_type_genes = defaultdict(list)
            for gene, trans_type in self.up_transcript_selection.items():
                transcript_type_genes[trans_type].append(gene)

            # For each transcript type, calculate frequency of enriched motifs in it
            type_motif_freq = defaultdict(dict)

            for trans_type, genes in transcript_type_genes.items():
                if len(genes) < 5:  # Skip types with too few samples
                    continue

                # Calculate frequency of each motif in this transcript type
                for _, row in up_enriched.iterrows():
                    pattern = row['Pattern']
                    # Get genes containing this motif
                    pattern_genes = self.get_genes_with_motif(pattern, self.up_sequences)
                    # Calculate proportion of genes in this transcript type containing this motif
                    common_genes = set(genes) & set(pattern_genes)
                    freq = len(common_genes) / len(genes) if genes else 0
                    type_motif_freq[trans_type][pattern] = freq

            # Create heatmap
            if type_motif_freq:
                # Limit image size to prevent exceeding matplotlib limits
                width = min(max(8, len(up_enriched) * 0.4), 50)  # Maximum width 50
                height = min(max(6, len(type_motif_freq) * 0.4), 30)  # Maximum height 30
                plt.figure(figsize=(width, height))

                # If data is too large, split into multiple plots
                if len(up_enriched) > 150 or len(type_motif_freq) > 100:
                    logger.warning(f"Data too large, will split images")

                    # Split large plot into multiple small plots
                    max_motifs_per_fig = 100
                    motif_keys = list(next(iter(type_motif_freq.values())).keys())

                    for i in range(0, len(motif_keys), max_motifs_per_fig):
                        end_idx = min(i + max_motifs_per_fig, len(motif_keys))
                        subset_motifs = motif_keys[i:end_idx]

                        # Create subset data
                        subset_data = {
                            trans_type: {motif: freq for motif, freq in freqs.items() if motif in subset_motifs}
                            for trans_type, freqs in type_motif_freq.items()}

                        # Continue normal plotting process, but use subset_data
                        subset_heatmap_data = pd.DataFrame(subset_data).T

                        plt.figure(figsize=(min(max(8, len(subset_motifs) * 0.4), 50), height))
                        sns.heatmap(subset_heatmap_data, cmap='YlOrRd', annot=True, fmt='.2f')
                        plt.title(
                            f'Motif Frequencies (Part {i // max_motifs_per_fig + 1})\n(UP-regulated genes, Length {length})')
                        plt.xlabel('Motif Patterns')
                        plt.ylabel('Transcript Type')
                        plt.tight_layout()
                        plt.savefig(f'{trans_dir}/up_motif_transcript_heatmap_part{i // max_motifs_per_fig + 1}.png')
                        plt.close()

                    return  # Skip creation of original large plot

                # Convert to DataFrame
                heatmap_data = pd.DataFrame(type_motif_freq).T  # Transpose to make transcript types as rows

                # Create heatmap
                sns.heatmap(heatmap_data, cmap='YlOrRd', annot=True, fmt='.2f')
                plt.title(f'Motif Frequencies Across Transcript Types\n(UP-regulated genes, Length {length})')
                plt.xlabel('Motif Patterns')
                plt.ylabel('Transcript Type')
                plt.tight_layout()
                plt.savefig(f'{trans_dir}/up_motif_transcript_heatmap.png')
                plt.close()

                # Save data
                heatmap_data.to_csv(f'{trans_dir}/up_motif_transcript_frequencies.csv')

        # Analyze downregulated genes (similar logic)
        down_enriched = significant_df[significant_df['Enrichment'] == 'DOWN']
        if not down_enriched.empty and self.down_transcript_selection:
            # Similar to above logic, but for downregulated genes
            transcript_type_genes = defaultdict(list)
            for gene, trans_type in self.down_transcript_selection.items():
                transcript_type_genes[trans_type].append(gene)

            type_motif_freq = defaultdict(dict)

            for trans_type, genes in transcript_type_genes.items():
                if len(genes) < 5:
                    continue

                for _, row in down_enriched.iterrows():
                    pattern = row['Pattern']
                    pattern_genes = self.get_genes_with_motif(pattern, self.down_sequences)
                    common_genes = set(genes) & set(pattern_genes)
                    freq = len(common_genes) / len(genes) if genes else 0
                    type_motif_freq[trans_type][pattern] = freq

            if type_motif_freq:
                # Limit image size to prevent exceeding matplotlib limits
                width = min(max(8, len(up_enriched) * 0.4), 50)  # Maximum width 50
                height = min(max(6, len(type_motif_freq) * 0.4), 30)  # Maximum height 30
                plt.figure(figsize=(width, height))

                # If data is too large, split into multiple plots
                if len(up_enriched) > 150 or len(type_motif_freq) > 100:
                    logger.warning(f"Data too large, will split images")

                    # Split large plot into multiple small plots
                    max_motifs_per_fig = 100
                    motif_keys = list(next(iter(type_motif_freq.values())).keys())

                    for i in range(0, len(motif_keys), max_motifs_per_fig):
                        end_idx = min(i + max_motifs_per_fig, len(motif_keys))
                        subset_motifs = motif_keys[i:end_idx]

                        # Create subset data
                        subset_data = {
                            trans_type: {motif: freq for motif, freq in freqs.items() if motif in subset_motifs}
                            for trans_type, freqs in type_motif_freq.items()}

                        # Continue normal plotting process, but use subset_data
                        subset_heatmap_data = pd.DataFrame(subset_data).T

                        plt.figure(figsize=(min(max(8, len(subset_motifs) * 0.4), 50), height))
                        sns.heatmap(subset_heatmap_data, cmap='YlOrRd', annot=True, fmt='.2f')
                        plt.title(
                            f'Motif Frequencies (Part {i // max_motifs_per_fig + 1})\n(UP-regulated genes, Length {length})')
                        plt.xlabel('Motif Patterns')
                        plt.ylabel('Transcript Type')
                        plt.tight_layout()
                        plt.savefig(f'{trans_dir}/up_motif_transcript_heatmap_part{i // max_motifs_per_fig + 1}.png')
                        plt.close()

                    return  # Skip creation of original large plot

                heatmap_data = pd.DataFrame(type_motif_freq).T

                sns.heatmap(heatmap_data, cmap='YlGnBu', annot=True, fmt='.2f')
                plt.title(f'Motif Frequencies Across Transcript Types\n(DOWN-regulated genes, Length {length})')
                plt.xlabel('Motif Patterns')
                plt.ylabel('Transcript Type')
                plt.tight_layout()
                plt.savefig(f'{trans_dir}/down_motif_transcript_heatmap.png')
                plt.close()

                heatmap_data.to_csv(f'{trans_dir}/down_motif_transcript_frequencies.csv')

    def get_genes_with_motif(self, pattern, sequences):
        """
        Get list of genes containing specific motif

        Purpose: Auxiliary function to determine which genes contain specific motif
        Improvement: Support patterns with mismatches and gaps
        """
        # Check if there are special formats
        if '{' in pattern:  # Gap pattern
            # Check if it's PWM gap pattern
            if '(pwm)' in pattern:
                # Extract basic pattern part
                base_pattern = pattern.replace('(pwm)', '')
                pattern_parts = re.split(r'N\{(\d+),(\d+)\}', base_pattern)
                prefix = pattern_parts[0]
                suffix = pattern_parts[-1]
                gap_info = re.search(r'N\{(\d+),(\d+)\}', base_pattern)
                if gap_info:
                    gap_min, gap_max = int(gap_info.group(1)), int(gap_info.group(2))

                    # Use PWM matching
                    matches, positions, _ = get_matched_sequences_with_pwm_gaps(
                        prefix, suffix, gap_min, gap_max, sequences,
                        threshold=0.75
                    )
                    return [pos[0] for pos in positions] if positions else []
            else:
                # Parse standard gap pattern
                matches, positions, _ = get_matched_sequences(None, sequences, gap_pattern=pattern)
                return [pos[0] for pos in positions] if positions else []

        elif '(' in pattern:  # Mismatch pattern
            # Extract basic pattern
            base_pattern = pattern.split('(')[0]
            # Get mismatch information
            mismatch_info = pattern.split('(')[1].split(')')[0] if '(' in pattern else None

            # Currently simple processing: use basic pattern, further analyze in result processing
            matches, positions, _ = get_matched_sequences(base_pattern, sequences,
                                                          max_mismatches=1 if mismatch_info else 0)
            return [pos[0] for pos in positions] if positions else []
        elif pattern.startswith('PWM:'):  # PWM pattern
            # Extract basic pattern for creating PWM
            base_pattern = pattern[4:]  # Remove "PWM:" prefix
            # Get matching sequences
            matches, _, _ = get_matched_sequences(base_pattern, sequences)
            if matches:
                # Create PWM
                pwm = calculate_position_frequencies(base_pattern, matches)
                # Scan using PWM
                _, pwm_positions, _ = scan_sequences_with_pwm(pwm, sequences, threshold=0.75)
                return [pos[0] for pos in pwm_positions] if pwm_positions else []
            return []
        else:  # Regular pattern
            matches, positions, _ = get_matched_sequences(pattern, sequences)
            return [pos[0] for pos in positions] if positions else []

    def get_motif_gene_details(self, pattern, sequences, positions):
        """
        Get detailed gene information for motif, including position and count of motif occurrence in each gene

        Parameters:
            pattern: motif pattern
            sequences: sequence dictionary
            positions: position dictionary

        Returns:
            gene_details: dictionary containing detailed gene information
        """
        gene_details = {}

        # Get position information for this motif
        pattern_positions = positions.get(pattern, [])

        # Count position and occurrences of motif in each gene
        for gene_id, pos in pattern_positions:
            if gene_id not in gene_details:
                gene_details[gene_id] = {
                    'count': 0,
                    'positions': []
                }
            gene_details[gene_id]['count'] += 1
            gene_details[gene_id]['positions'].append(pos)

        return gene_details

    def create_motif_details_file(self, significant_df, up_positions, down_positions, result_dir):
        """
        Create detailed information CSV file for significant motifs

        Parameters:
            significant_df: significant motif DataFrame
            up_positions: upregulated gene position information
            down_positions: downregulated gene position information
            result_dir: results directory
        """
        # Prepare detailed information data
        details_data = []

        for _, row in significant_df.iterrows():
            pattern = row['Pattern']
            enrichment = row['Enrichment']

            # Get UP group gene details
            up_gene_details = self.get_motif_gene_details(pattern, self.up_sequences, up_positions)

            # Get DOWN group gene details
            down_gene_details = self.get_motif_gene_details(pattern, self.down_sequences, down_positions)

            # Basic motif information
            motif_info = {
                'Pattern': pattern,
                'Enrichment': enrichment,
                'Length': row['Length'],
                'Type': row['Type'],
                'UP_frequency': row['UP_frequency'],
                'DOWN_frequency': row['DOWN_frequency'],
                'UP_to_DOWN_ratio': row['UP_to_DOWN_ratio'],
                'DOWN_to_UP_ratio': row['DOWN_to_UP_ratio'],
                'UP_genes_count': len(up_gene_details),
                'DOWN_genes_count': len(down_gene_details),
                'UP_genes': ','.join(sorted(up_gene_details.keys())),
                'DOWN_genes': ','.join(sorted(down_gene_details.keys())),
                'UP_genes_detail': json.dumps(up_gene_details),
                'DOWN_genes_detail': json.dumps(down_gene_details),
            }

            details_data.append(motif_info)

        # Create DataFrame and save
        if details_data:
            details_df = pd.DataFrame(details_data)
            details_df.to_csv(f'{result_dir}/02significant_motifs_with_details.csv', index=False)

            # Also save a simplified version containing only gene lists
            simple_df = details_df[['Pattern', 'Enrichment', 'UP_genes', 'DOWN_genes']]
            simple_df.to_csv(f'{result_dir}/02significant_motifs_gene_lists.csv', index=False)

            # Add to original significant motif file
            significant_df['UP_genes_count'] = details_df['UP_genes_count']
            significant_df['DOWN_genes_count'] = details_df['DOWN_genes_count']
            significant_df.to_csv(f'{result_dir}/02significant_motifs.csv', index=False)

        # Create separate gene detailed information file for each motif
        for pattern in set(details_df['Pattern']):
            pattern_df = details_df[details_df['Pattern'] == pattern]
            if len(pattern_df) > 0:
                row = pattern_df.iloc[0]

                # Create UP gene detailed information
                if row['UP_genes_count'] > 0:
                    up_details = json.loads(row['UP_genes_detail'])
                    up_data = []
                    for gene, details in up_details.items():
                        up_data.append({
                            'Gene': gene,
                            'Count': details['count'],
                            'Positions': ';'.join(map(str, details['positions']))
                        })
                    pd.DataFrame(up_data).to_csv(f'{result_dir}/{pattern}_UP_genes_details.csv', index=False)

                # Create DOWN gene detailed information
                if row['DOWN_genes_count'] > 0:
                    down_details = json.loads(row['DOWN_genes_detail'])
                    down_data = []
                    for gene, details in down_details.items():
                        down_data.append({
                            'Gene': gene,
                            'Count': details['count'],
                            'Positions': ';'.join(map(str, details['positions']))
                        })
                    pd.DataFrame(down_data).to_csv(f'{result_dir}/{pattern}_DOWN_genes_details.csv', index=False)


def main():
    # Set random seed to ensure reproducible results
    set_random_seeds(42)  # You can use any integer as seed

    # Configure parameters
    max_utr_length = 500  # According to paper, only take up to 500 nucleotides downstream of TSS
    min_frequency = 0.05  # Minimum occurrence frequency threshold for individual motifs (for initial screening)
    fold_change_threshold = 1.5  # Enrichment fold change threshold
    min_occur_freq = 0.1  # Minimum occurrence frequency threshold (for determining significance)
    n_cores = None  # Auto-detect
    max_wildcards = 1  # Maximum number of wildcards allowed in motifs

    # New advanced motif analysis parameters
    max_mismatches = 1  # Maximum number of mismatches allowed
    enable_gapped = True  # Enable gap motif analysis
    enable_pwm = True  # Enable PWM analysis
    enable_structural = False  # Enable structural motif analysis
    enable_pwm_gapped = True  # Add this line

    # Initialize analyzer
    analyzer = MotifAnalyzer(
        max_utr_length=max_utr_length,
        min_frequency=min_frequency,
        fold_change_threshold=fold_change_threshold,
        min_occur_freq=min_occur_freq,
        max_mismatches=max_mismatches,
        enable_gapped=enable_gapped,
        enable_pwm=enable_pwm,
        enable_structural=enable_structural,
        enable_pwm_gapped=enable_pwm_gapped
    )

    # Create results root directory
    os.makedirs("motif_analysis_results", exist_ok=True)

    # Overall results summary
    summary_data = []

    # Analyze 4-8bp sequences
    for length in range(4, 9):
        logger.info(f"\nAnalyzing motifs of length {length}")
        results = analyzer.analyze_single_length(
            length,
            n_cores=n_cores,
            max_wildcards=max_wildcards
        )

        # Add to summary
        if not results.empty:
            significant = results[results['Significant']]
            up_enriched = significant[significant['Enrichment'] == 'UP']
            down_enriched = significant[significant['Enrichment'] == 'DOWN']

            # Statistics by type
            motif_types = results['Type'].value_counts().to_dict()
            significant_types = significant['Type'].value_counts().to_dict() if not significant.empty else {}

            summary_data.append({
                'Length': length,
                'Total_motifs': len(results),
                'Significant_motifs': len(significant),
                'UP_enriched': len(up_enriched),
                'DOWN_enriched': len(down_enriched),
                'Max_UP_to_DOWN_ratio': results['UP_to_DOWN_ratio'].replace([np.inf, -np.inf], np.nan).max(),
                'Max_DOWN_to_UP_ratio': results['DOWN_to_UP_ratio'].replace([np.inf, -np.inf], np.nan).max(),
                'Avg_UP_frequency': up_enriched['UP_frequency'].mean() if not up_enriched.empty else 0,
                'Avg_DOWN_frequency': down_enriched['DOWN_frequency'].mean() if not down_enriched.empty else 0,
                'Avg_strand_specificity': results['Overall_strand_specificity'].replace([np.inf, -np.inf],
                                                                                        np.nan).mean(),
                'Exact_motifs': motif_types.get('exact', 0),
                'Mismatch_motifs': motif_types.get('mismatch', 0),
                'Gapped_motifs': motif_types.get('gapped', 0),
                'PWM_motifs': motif_types.get('pwm', 0),
                'Significant_exact': significant_types.get('exact', 0),
                'Significant_mismatch': significant_types.get('mismatch', 0),
                'Significant_gapped': significant_types.get('gapped', 0),
                'Significant_pwm': significant_types.get('pwm', 0)
            })

    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv("motif_analysis_results/analysis_summary.csv", index=False)

        # Plot summary figures
        plt.figure(figsize=(15, 15))

        # 1. Motif count statistics
        plt.subplot(2, 2, 1)
        width = 0.3
        x = np.arange(len(summary_df))

        plt.bar(x, summary_df['Total_motifs'], width, label='Total')
        plt.bar(x + width, summary_df['Significant_motifs'], width, label='Significant')

        # Add up and down regulation counts
        plt.bar(x + 2 * width, summary_df['UP_enriched'], width, label='UP enriched', color='red')
        plt.bar(x + 3 * width, summary_df['DOWN_enriched'], width, label='DOWN enriched', color='blue')

        plt.xlabel('Motif Length')
        plt.ylabel('Number of Motifs')
        plt.title('Motif Counts by Length')
        plt.xticks(x + 1.5 * width, summary_df['Length'])
        plt.legend()

        # 2. Average frequency statistics
        plt.subplot(2, 2, 2)
        plt.plot(summary_df['Length'], summary_df['Avg_UP_frequency'], marker='o', color='red', label='UP frequency')
        plt.plot(summary_df['Length'], summary_df['Avg_DOWN_frequency'], marker='s', color='blue',
                 label='DOWN frequency')
        plt.axhline(y=min_occur_freq, color='r', linestyle='--', label=f'Min Frequency={min_occur_freq}')
        plt.xlabel('Motif Length')
        plt.ylabel('Average Frequency')
        plt.title('Average Frequency by Length')
        plt.legend()

        # 3. Enrichment fold change statistics
        plt.subplot(2, 2, 3)
        plt.plot(summary_df['Length'], summary_df['Max_UP_to_DOWN_ratio'], marker='o', label='UP enrichment',
                 color='red')
        plt.plot(summary_df['Length'], summary_df['Max_DOWN_to_UP_ratio'], marker='s', label='DOWN enrichment',
                 color='blue')
        plt.axhline(y=fold_change_threshold, color='r', linestyle='--', label=f'Fold Change={fold_change_threshold}')
        plt.xlabel('Motif Length')
        plt.ylabel('Maximum Fold Enrichment')
        plt.title('Maximum Fold Enrichment by Length')
        plt.legend()

        # 4. Motif type distribution
        plt.subplot(2, 2, 4)
        motif_types = ['Exact_motifs', 'Mismatch_motifs', 'Gapped_motifs', 'PWM_motifs']
        sig_types = ['Significant_exact', 'Significant_mismatch', 'Significant_gapped', 'Significant_pwm']

        # Calculate proportion of each type
        motif_props = summary_df[motif_types].div(summary_df['Total_motifs'], axis=0)
        sig_props = summary_df[sig_types].div(summary_df[motif_types].values, axis=1)

        # Stacked bar chart
        bottom = np.zeros(len(summary_df))
        for i, col in enumerate(motif_types):
            plt.bar(summary_df['Length'], summary_df[col], bottom=bottom, label=col.split('_')[0])
            bottom += summary_df[col].values

        plt.xlabel('Motif Length')
        plt.ylabel('Number of Motifs')
        plt.title('Motif Types Distribution')
        plt.legend()

        plt.tight_layout()
        plt.savefig("motif_analysis_results/analysis_summary.png")
        plt.close()

        # Add advanced analysis summary plot
        plt.figure(figsize=(15, 10))

        # 1. Mismatch statistics
        plt.subplot(2, 2, 1)
        plt.plot(summary_df['Length'], summary_df['Mismatch_motifs'], marker='o', color='orange',
                 label='Total Mismatch')
        plt.plot(summary_df['Length'], summary_df['Significant_mismatch'], marker='s', color='red',
                 label='Significant Mismatch')
        plt.xlabel('Motif Length')
        plt.ylabel('Count')
        plt.title('Mismatch Motifs by Length')
        plt.legend()

        # 2. Gap statistics
        plt.subplot(2, 2, 2)
        plt.plot(summary_df['Length'], summary_df['Gapped_motifs'], marker='o', color='green',
                 label='Total Gapped')
        plt.plot(summary_df['Length'], summary_df['Significant_gapped'], marker='s', color='lime',
                 label='Significant Gapped')
        plt.xlabel('Motif Length')
        plt.ylabel('Count')
        plt.title('Gapped Motifs by Length')
        plt.legend()

        # 3. PWM statistics
        plt.subplot(2, 2, 3)
        plt.plot(summary_df['Length'], summary_df['PWM_motifs'], marker='o', color='purple',
                 label='Total PWM')
        plt.plot(summary_df['Length'], summary_df['Significant_pwm'], marker='s', color='magenta',
                 label='Significant PWM')
        plt.xlabel('Motif Length')
        plt.ylabel('Count')
        plt.title('PWM Motifs by Length')
        plt.legend()

        # 4. Type significance proportion
        plt.subplot(2, 2, 4)
        for i, (col, label) in enumerate(zip(sig_types, ['Exact', 'Mismatch', 'Gapped', 'PWM'])):
            plt.plot(summary_df['Length'], sig_props[col].fillna(0), marker='o', label=f'{label}')

        plt.xlabel('Motif Length')
        plt.ylabel('Proportion Significant')
        plt.title('Proportion of Significant Motifs by Type')
        plt.ylim(0, 1)
        plt.legend()

        plt.tight_layout()
        plt.savefig("motif_analysis_results/advanced_analysis_summary.png")
        plt.close()

        # Generate instructions for GO enrichment analysis
        with open("motif_analysis_results/README.txt", "w") as f:
            f.write("Motif Analysis Results Guide\n")
            f.write("==========================\n\n")
            f.write("Each motif directory contains a '*_matching_genes.txt' file with genes containing that motif.\n")
            f.write("These gene lists can be used for GO/KEGG enrichment analysis using tools like:\n")
            f.write("1. DAVID (https://david.ncifcrf.gov/)\n")
            f.write("2. Metascape (https://metascape.org/)\n")
            f.write("3. g:Profiler (https://biit.cs.ut.ee/gprofiler/)\n")
            f.write("4. Enrichr (https://maayanlab.cloud/Enrichr/)\n\n")
            f.write("Significant motifs were determined using the following criteria:\n")
            f.write(f"- UP enriched: UP frequency >= {min_occur_freq} and UP/DOWN ratio >= {fold_change_threshold}\n")
            f.write(
                f"- DOWN enriched: DOWN frequency >= {min_occur_freq} and DOWN/UP ratio >= {fold_change_threshold}\n\n")
            f.write(
                "For each significant motif, analyze the matching gene list to determine potential biological functions.\n")

            # Add MANE related information
            f.write("\nMANE Transcript Selection Analysis\n")
            f.write("===============================\n\n")
            f.write("This analysis includes additional insights based on MANE transcript selections:\n")
            f.write(
                "- Transcript selection distribution shows how different transcript types were used in the analysis\n")
            f.write(
                "- Motif frequency across transcript types reveals potential transcript-specific regulatory elements\n")
            f.write(
                "- Individual motif reports include information about their distribution across transcript types\n\n")
            f.write(
                "These analyses may help identify motifs that function specifically in certain transcript isoforms.\n")

            # Add advanced motif analysis information
            f.write("\nAdvanced Motif Analysis\n")
            f.write("======================\n\n")
            f.write("This analysis includes several advanced motif types:\n\n")

            f.write("1. Mismatch Motifs\n")
            f.write("----------------\n")
            f.write("Patterns with format 'PATTERN(mX:Y>Z)' allow for X mismatches where position Y\n")
            f.write("can be base Z instead of the original base. This helps identify less stringent\n")
            f.write("regulatory elements that tolerate some variation.\n\n")

            f.write("2. Gapped Motifs\n")
            f.write("--------------\n")
            f.write("Patterns with format 'PREFIX{min,max}SUFFIX' represent motifs with variable-length\n")
            f.write("gaps between conserved regions. These are common in RNA-binding protein recognition\n")
            f.write("sites and structural motifs.\n\n")

            f.write("3. PWM-based Motifs\n")
            f.write("-----------------\n")
            f.write("Patterns with format 'PWM:PATTERN' use position weight matrices to identify less\n")
            f.write("stringent matches based on the positional nucleotide frequencies of the original pattern.\n\n")

            f.write("4. Structural Motifs\n")
            f.write("------------------\n")
            f.write("The analysis also identifies potential RNA secondary structures like stem-loops and\n")
            f.write("G-quadruplexes that may play roles in post-transcriptional regulation.\n\n")

            f.write("For each motif, correlation with structural elements is analyzed to identify\n")
            f.write("potential structure-sequence relationships in regulatory mechanisms.\n")

    # If structural analysis is enabled, create overall structural report
    if analyzer.enable_structural:
        # Create structural analysis summary directory
        struct_dir = "structural_analysis_summary"
        os.makedirs(struct_dir, exist_ok=True)

        # Merge correlation analysis of motifs and structures from all lengths
        all_struct_correlations = {}

        # Find structural correlation analysis results for all significant motifs
        for length in range(4, 9):
            result_dir = f'motif_results_length_{length}'

            # Scan UP and DOWN results
            for direction in ['up', 'down']:
                for filename in os.listdir(result_dir):
                    if filename.startswith(f"{direction}_motif_") and filename.endswith("_details.txt"):
                        motif_name = filename.replace(f"{direction}_motif_", "").replace("_details.txt", "")

                        # Read file to find structural correlation
                        with open(f"{result_dir}/{filename}", 'r') as f:
                            content = f.read()

                            if "Structural Correlation Analysis" in content:
                                # Extract correlation data
                                struct_section = content.split("Structural Correlation Analysis:")[1].split("\n\n")[0]

                                # Parse correlation data
                                stem_loop_data = {}
                                g_quad_data = {}

                                if "Stem Loops" in struct_section:
                                    stem_loop_text = struct_section.split("- Stem Loops:")[1].split("-")[0]

                                    # Extract p-value
                                    p_match = re.search(r"P-value: ([0-9.e-]+)", stem_loop_text)
                                    if p_match:
                                        p_value = float(p_match.group(1))
                                        stem_loop_data['p_value'] = p_value
                                        stem_loop_data['significant'] = p_value < 0.05

                                if "G Quadruplexes" in struct_section:
                                    g_quad_text = struct_section.split("- G Quadruplexes:")[1].split("\n\n")[0]

                                    # Extract p-value
                                    p_match = re.search(r"P-value: ([0-9.e-]+)", g_quad_text)
                                    if p_match:
                                        p_value = float(p_match.group(1))
                                        g_quad_data['p_value'] = p_value
                                        g_quad_data['significant'] = p_value < 0.05

                                # Save data
                                all_struct_correlations[f"{motif_name}_{direction}_{length}"] = {
                                    'motif': motif_name,
                                    'direction': direction,
                                    'length': length,
                                    'stem_loops': stem_loop_data,
                                    'g_quadruplexes': g_quad_data
                                }

        # Save structural correlation summary
        if all_struct_correlations:
            # Create DataFrame
            struct_corr_data = []

            for key, data in all_struct_correlations.items():
                row = {
                    'Motif': data['motif'],
                    'Direction': data['direction'].upper(),
                    'Length': data['length'],
                    'StemLoop_pvalue': data['stem_loops'].get('p_value', np.nan),
                    'StemLoop_significant': data['stem_loops'].get('significant', False),
                    'GQuad_pvalue': data['g_quadruplexes'].get('p_value', np.nan),
                    'GQuad_significant': data['g_quadruplexes'].get('significant', False)
                }
                struct_corr_data.append(row)

            # Create DataFrame and save
            struct_corr_df = pd.DataFrame(struct_corr_data)
            struct_corr_df.to_csv(f"{struct_dir}/structural_correlations.csv", index=False)

            # Create heatmap
            if len(struct_corr_df) > 0:
                # Prepare heatmap data
                heatmap_data = struct_corr_df.pivot_table(
                    index='Motif',
                    columns='Direction',
                    values=['StemLoop_pvalue', 'GQuad_pvalue'],
                    aggfunc='min'  # Use minimum p-value
                )

                # Plot heatmap
                plt.figure(figsize=(12, max(6, len(heatmap_data) * 0.3)))

                # Log transform p-values
                log_pvals = -np.log10(heatmap_data)

                # Set heatmap mask
                mask = np.isnan(log_pvals)

                # Plot heatmap
                sns.heatmap(log_pvals, cmap='YlOrRd', annot=True, fmt='.2f', mask=mask)
                plt.title('Structure-Sequence Correlations (-log10 p-value)')
                plt.tight_layout()
                plt.savefig(f"{struct_dir}/structural_correlation_heatmap.png")
                plt.close()

                # Create significant correlation network diagram
                # Filter significant correlations
                sig_correlations = struct_corr_df[(struct_corr_df['StemLoop_significant']) |
                                                  (struct_corr_df['GQuad_significant'])]

                if len(sig_correlations) > 0:
                    plt.figure(figsize=(12, 10))
                    G = nx.Graph()

                    # Add nodes
                    motifs = sig_correlations['Motif'].unique()
                    structures = ['Stem-Loop', 'G-Quadruplex']

                    for motif in motifs:
                        G.add_node(motif, type='motif')

                    for struct in structures:
                        G.add_node(struct, type='structure')

                    # Add edges
                    for _, row in sig_correlations.iterrows():
                        motif = row['Motif']

                        if row['StemLoop_significant']:
                            weight = -np.log10(row['StemLoop_pvalue'])
                            G.add_edge(motif, 'Stem-Loop', weight=weight, pvalue=row['StemLoop_pvalue'])

                        if row['GQuad_significant']:
                            weight = -np.log10(row['GQuad_pvalue'])
                            G.add_edge(motif, 'G-Quadruplex', weight=weight, pvalue=row['GQuad_pvalue'])

                    # Set node colors
                    node_colors = ['red' if G.nodes[n]['type'] == 'motif' else 'blue'
                                   for n in G.nodes()]

                    # Set edge widths
                    edge_widths = [G[u][v]['weight'] for u, v in G.edges()]

                    # Draw network
                    pos = nx.spring_layout(G, k=0.8)
                    nx.draw_networkx(G, pos, node_color=node_colors, width=edge_widths,
                                     font_size=9, node_size=600, font_weight='bold',
                                     alpha=0.8)

                    # Add edge labels
                    edge_labels = {(u, v): f"p:{G[u][v]['pvalue']:.1e}"
                                   for u, v in G.edges()}
                    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=7)

                    plt.title('Significant Structure-Sequence Correlations Network')
                    plt.axis('off')
                    plt.savefig(f"{struct_dir}/structural_correlation_network.png")
                    plt.close()

    logger.info("Motif analysis complete!")


if __name__ == "__main__":
    main()