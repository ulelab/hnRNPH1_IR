#!/usr/bin/env python3
"""
Process MAF files: extract FASTA with msa_view, run SpliceAI on 50-nt windows around pos1/pos2
Version 3.0: Uses positional information from v3.1 TSV to extract local 50-nt windows
             and recalculate max1 and max2 scores for each species
"""

import os
import sys
import argparse
import subprocess
import tempfile
import gc
import re
import numpy as np
from collections import defaultdict
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

# Silence TensorFlow output
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
tf.get_logger().setLevel('ERROR')

def load_species_list(species_map_file):
    """Load species IDs from tree_to_clade_mapping.tsv"""
    species_list = []
    with open(species_map_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                species_id = fields[2]  # column 3: species_id
                species_list.append(species_id)
    return species_list

def load_strand_map(decoy_bed_file):
    """Load strand information from decoyIDs.bed: name (col4) -> strand (col6)
    Handles both tab and space delimited BED files"""
    strand_map = {}
    with open(decoy_bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Try tab first, then space (BED files can use either)
            if '\t' in line:
                fields = line.split('\t')
            else:
                fields = line.split()
            if len(fields) >= 6:
                name = fields[3]  # column 4: name
                strand = fields[5]  # column 6: strand
                strand_map[name] = strand
    return strand_map

def load_spliceai_models():
    """Load SpliceAI models"""
    context = 10000
    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
    models = [load_model(resource_filename('spliceai', x), compile=False) for x in paths]
    predict_fn = [model.predict for model in models]
    return predict_fn, context

def load_v31_tsv(v31_tsv_file):
    """
    Load v3.1 TSV file and extract pos1 and pos2 for each decoy/species combination
    Returns: dict[decoy_id][species_id] = (max1, pos1, max2, pos2, d)
    """
    positions = defaultdict(dict)
    
    with open(v31_tsv_file, 'r') as f:
        header = next(f).strip().split('\t')
        species_columns = header[1:]  # Skip 'decoyID' column
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            
            decoy_id = fields[0]
            
            # Parse each species column: [max1][pos1][max2][pos2][d]
            for i, species_id in enumerate(species_columns):
                if i + 1 >= len(fields):
                    continue
                
                score_str = fields[i + 1]
                # Parse format: [max1][pos1][max2][pos2][d]
                match = re.match(r'\[([-\d.]+)\]\[([-\d]+)\]\[([-\d.]+)\]\[([-\d]+)\]\[([-\d.]+)\]', score_str)
                if match:
                    max1 = float(match.group(1))
                    pos1 = int(match.group(2))
                    max2 = float(match.group(3))
                    pos2 = int(match.group(4))
                    d = float(match.group(5))
                    
                    positions[decoy_id][species_id] = (max1, pos1, max2, pos2, d)
                else:
                    # Invalid format, use -1 for all
                    positions[decoy_id][species_id] = (-1.0, -1, -1.0, -1, -1.0)
    
    return positions

def extract_fasta_from_maf(maf_file, msa_view_bin, strand, temp_dir, debug=False):
    """Extract FASTA from MAF using msa_view with strand-aware logic"""
    fasta_file = os.path.join(temp_dir, f"{os.path.basename(maf_file)}.fa")
    
    # Build msa_view command
    cmd = [msa_view_bin, '--in-format', 'MAF', '--out-format', 'FASTA']
    
    # Add -V flag for reverse complement if negative strand
    if strand == '-':
        cmd.append('-V')
    
    cmd.append(maf_file)
    
    # Run msa_view and write to fasta file
    try:
        with open(fasta_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        return fasta_file
    except subprocess.CalledProcessError as e:
        # Invalid/empty MAF files are expected for some regions - only print in debug mode
        if debug:
            error_msg = e.stderr.decode().strip()
            print(f"WARNING: msa_view failed for {os.path.basename(maf_file)}: {error_msg}", file=sys.stderr)
        return None

def get_sequence_from_fasta(fasta_file, species_id):
    """
    Extract the full sequence (with gaps) for a specific species from FASTA file
    Returns: sequence string (with gaps, -, *, spaces preserved)
    """
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip()  # Keep trailing spaces if any
            if not line:
                continue
            if line.startswith('>'):
                # Save previous sequence if it matches target species
                if current_header is not None and current_sequence:
                    # Extract species_id from header
                    if '.' in current_header:
                        header_species_id = current_header.split('.')[0]
                    elif ' ' in current_header:
                        header_species_id = current_header.split()[0]
                    else:
                        header_species_id = current_header
                    
                    if header_species_id == species_id:
                        return ''.join(current_sequence)
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                if current_header is not None:
                    current_sequence.append(line)
        
        # Process last sequence
        if current_header is not None and current_sequence:
            if '.' in current_header:
                header_species_id = current_header.split('.')[0]
            elif ' ' in current_header:
                header_species_id = current_header.split()[0]
            else:
                header_species_id = current_header
            
            if header_species_id == species_id:
                return ''.join(current_sequence)
    
    return None  # Species not found

def ungapped_position_to_gapped(sequence, ungapped_pos):
    """
    Convert an ungapped position (0-based) to a gapped position in the sequence
    Ignores gaps (-, *, spaces) when counting
    Returns: gapped position (0-based) or -1 if position is beyond sequence
    """
    if ungapped_pos < 0:
        return -1
    
    ungapped_count = 0
    for i, char in enumerate(sequence):
        if char not in ['-', '*', ' ', '\n', '\r']:
            if ungapped_count == ungapped_pos:
                return i
            ungapped_count += 1
    
    # Position is beyond the sequence
    return -1

def extract_window_around_position(sequence, ungapped_pos, window_size=50):
    """
    Extract a window of size window_size (default 50) centered on ungapped_pos
    The window will contain at least window_size ungapped nucleotides
    Returns: (window_sequence, start_gapped_pos, end_gapped_pos)
             Returns (None, -1, -1) if position is invalid or window cannot be extracted
    """
    if ungapped_pos < 0:
        return (None, -1, -1)
    
    # Find the gapped position corresponding to ungapped_pos
    center_gapped = ungapped_position_to_gapped(sequence, ungapped_pos)
    if center_gapped < 0:
        return (None, -1, -1)
    
    # Calculate half window size (25 nt on each side)
    half_window = window_size // 2
    
    # Find start position: go back half_window ungapped positions
    start_ungapped = max(0, ungapped_pos - half_window)
    start_gapped = ungapped_position_to_gapped(sequence, start_ungapped)
    if start_gapped < 0:
        start_gapped = 0
    
    # Find end position: go forward half_window ungapped positions
    # First, find the total ungapped length
    total_ungapped = sum(1 for c in sequence if c not in ['-', '*', ' ', '\n', '\r'])
    end_ungapped = min(total_ungapped - 1, ungapped_pos + half_window)
    end_gapped = ungapped_position_to_gapped(sequence, end_ungapped)
    if end_gapped < 0:
        # If we can't find exact end, use sequence end
        end_gapped = len(sequence)
    else:
        # Include the character at end_gapped
        end_gapped += 1
    
    # Extract window
    window_sequence = sequence[start_gapped:end_gapped]
    
    # Verify we have at least window_size ungapped nucleotides
    # If not, try to extend the window
    ungapped_in_window = sum(1 for c in window_sequence if c not in ['-', '*', ' ', '\n', '\r'])
    if ungapped_in_window < window_size:
        # Try to extend the window to get more ungapped nucleotides
        # Extend backward first
        while start_gapped > 0 and ungapped_in_window < window_size:
            start_gapped -= 1
            window_sequence = sequence[start_gapped:end_gapped]
            ungapped_in_window = sum(1 for c in window_sequence if c not in ['-', '*', ' ', '\n', '\r'])
        
        # Then extend forward if needed
        while end_gapped < len(sequence) and ungapped_in_window < window_size:
            end_gapped += 1
            window_sequence = sequence[start_gapped:end_gapped]
            ungapped_in_window = sum(1 for c in window_sequence if c not in ['-', '*', ' ', '\n', '\r'])
    
    return (window_sequence, start_gapped, end_gapped)

def process_window_with_spliceai(window_sequence, predict_fn, context, debug=False):
    """
    Process a 50-nt window sequence with SpliceAI and return the max donor score
    Returns: max_donor_score (float) or -1.0 if processing fails
    """
    if window_sequence is None:
        return -1.0
    
    # Remove all gap characters: hyphens, spaces, and asterisks
    sequence_clean = window_sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
    
    # Skip if sequence is too short or contains only gaps/N's
    if len(sequence_clean) == 0:
        if debug:
            print(f"      Empty window after removing gaps", file=sys.stderr)
        return -1.0
    
    # Check if sequence has any actual nucleotides (not just N's)
    sequence_no_n = sequence_clean.replace('N', '').replace('n', '')
    if len(sequence_no_n) == 0:
        if debug:
            print(f"      Window contains only N's", file=sys.stderr)
        return -1.0
    
    try:
        # Pad with context
        padded_seq = 'N' * (context // 2) + sequence_clean + 'N' * (context // 2)
        x = one_hot_encode(padded_seq)[None, :]
        
        # Run predictions (ensemble of 5 models)
        y = np.mean([fn(x, verbose=0) for fn in predict_fn], axis=0)
        scores = y[0]  # shape (len + context, 4)
        
        # Extract donor probabilities for the actual sequence (trimming context)
        if len(sequence_clean) >= context:
            start = context // 2
            end = start + len(sequence_clean)
            donor_probs = scores[start:end, 2]
        else:
            donor_probs = scores[:len(sequence_clean), 2]
        
        # Find maximum donor score in the window
        max_donor = float(np.max(donor_probs))
        
        if debug:
            print(f"      Window length: {len(sequence_clean)}, max donor score: {max_donor:.6f}", file=sys.stderr)
        
        return max_donor
        
    except Exception as e:
        if debug:
            print(f"      WARNING: SpliceAI failed for window: {e}", file=sys.stderr)
        return -1.0

def process_maf_file(maf_file, msa_view_bin, strand_map, predict_fn, context, all_species, 
                     v31_positions, temp_dir, debug=False):
    """
    Process a single MAF file: extract FASTA, extract 50-nt windows around pos1/pos2,
    run SpliceAI on windows, return updated scores
    """
    # Get decoy ID from MAF filename (without .maf extension)
    decoy_id = os.path.splitext(os.path.basename(maf_file))[0]
    
    # Look up strand from decoyIDs.bed
    strand = strand_map.get(decoy_id, '+')  # Default to '+' if not found
    
    if debug:
        print(f"  Processing {decoy_id}, strand: {strand}", file=sys.stderr)
    
    # Check if we have position data for this decoy
    if decoy_id not in v31_positions:
        if debug:
            print(f"  WARNING: No position data found for {decoy_id} in v3.1 TSV", file=sys.stderr)
        # Return -1 for all species
        return {species: (-1.0, -1, -1.0, -1, -1.0) for species in all_species}
    
    # Extract FASTA from MAF
    fasta_file = extract_fasta_from_maf(maf_file, msa_view_bin, strand, temp_dir, debug=debug)
    if fasta_file is None:
        # Return -1 for all species if extraction failed
        return {species: (-1.0, -1, -1.0, -1, -1.0) for species in all_species}
    
    # Check if FASTA file has content
    if os.path.getsize(fasta_file) == 0:
        if debug:
            print(f"  WARNING: Empty FASTA file for {decoy_id}", file=sys.stderr)
        return {species: (-1.0, -1, -1.0, -1, -1.0) for species in all_species}
    
    # Process each species
    result = {}
    decoy_positions = v31_positions[decoy_id]
    
    for species in all_species:
        if species not in decoy_positions:
            # No position data for this species, return -1
            result[species] = (-1.0, -1, -1.0, -1, -1.0)
            continue
        
        max1_old, pos1, max2_old, pos2, d = decoy_positions[species]
        
        # If positions are invalid, return -1
        if pos1 < 0 or pos2 < 0:
            result[species] = (-1.0, pos1, -1.0, pos2, d)
            continue
        
        # Get sequence for this species
        sequence = get_sequence_from_fasta(fasta_file, species)
        if sequence is None:
            if debug:
                print(f"    {species}: sequence not found in FASTA", file=sys.stderr)
            result[species] = (-1.0, pos1, -1.0, pos2, d)
            continue
        
        if debug:
            print(f"    {species}: pos1={pos1}, pos2={pos2}, d={d:.6f}", file=sys.stderr)
        
        # Extract 50-nt window around pos1
        window1_seq, start1, end1 = extract_window_around_position(sequence, pos1, window_size=50)
        if window1_seq is None:
            if debug:
                print(f"      Could not extract window around pos1={pos1}", file=sys.stderr)
            max1_new = -1.0
        else:
            max1_new = process_window_with_spliceai(window1_seq, predict_fn, context, debug=debug)
        
        # Extract 50-nt window around pos2
        window2_seq, start2, end2 = extract_window_around_position(sequence, pos2, window_size=50)
        if window2_seq is None:
            if debug:
                print(f"      Could not extract window around pos2={pos2}", file=sys.stderr)
            max2_new = -1.0
        else:
            max2_new = process_window_with_spliceai(window2_seq, predict_fn, context, debug=debug)
        
        # Store result with same pos1, pos2, d but updated max1 and max2
        result[species] = (max1_new, pos1, max2_new, pos2, d)
        
        if debug:
            print(f"    {species}: max1={max1_old:.6f}->{max1_new:.6f}, max2={max2_old:.6f}->{max2_new:.6f}", file=sys.stderr)
    
    return result

def format_score_output(max1, pos1, max2, pos2, d):
    """Format output as [max1][pos1][max2][pos2][d]"""
    return f"[{max1:.6f}][{pos1}][{max2:.6f}][{pos2}][{d:.6f}]"

def main():
    parser = argparse.ArgumentParser(description='Process MAF files with msa_view and SpliceAI (v3.0: local 50-nt windows around pos1/pos2 from v3.1 TSV)')
    parser.add_argument('--maf_dir', required=True, help='Directory containing MAF files (or single MAF file for testing)')
    parser.add_argument('--decoy_bed', required=True, help='BED file with decoy IDs and strand info')
    parser.add_argument('--species_map', required=True, help='TSV file with species mapping')
    parser.add_argument('--v31_tsv', required=True, help='v3.1 TSV file with position information')
    parser.add_argument('--msa_view', required=True, help='Path to msa_view binary')
    parser.add_argument('--output_tsv', required=True, help='Output TSV file')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    print(f"Loading v3.1 TSV positions from: {args.v31_tsv}")
    v31_positions = load_v31_tsv(args.v31_tsv)
    print(f"Loaded position data for {len(v31_positions)} decoy IDs")
    
    print(f"Loading species list from: {args.species_map}")
    all_species = load_species_list(args.species_map)
    print(f"Loaded {len(all_species)} species")
    
    print(f"Loading strand map from: {args.decoy_bed}")
    strand_map = load_strand_map(args.decoy_bed)
    print(f"Loaded strand info for {len(strand_map)} decoy IDs")
    
    print("Loading SpliceAI models...")
    predict_fn, context = load_spliceai_models()
    print("SpliceAI models loaded")
    
    # Find all MAF files (handle both directory and single file)
    maf_files = []
    if os.path.isfile(args.maf_dir) and args.maf_dir.endswith('.maf'):
        # Single MAF file provided
        maf_files = [args.maf_dir]
    else:
        # Directory provided
        for root, dirs, files in os.walk(args.maf_dir):
            for file in files:
                if file.endswith('.maf'):
                    maf_files.append(os.path.join(root, file))
    
    print(f"Found {len(maf_files)} MAF files to process")
    
    # Create temporary directory for FASTA files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Open output file for writing (write incrementally to save memory)
        print(f"Writing results to: {args.output_tsv}")
        with open(args.output_tsv, 'w') as f:
            # Write header: decoyID, then each species with 5 fields
            header_parts = ["decoyID"]
            for species in all_species:
                header_parts.append(species)
            f.write("\t".join(header_parts) + "\n")
            
            # Process each MAF file and write immediately
            failed_files = 0
            processed_count = 0
            
            # Sort MAF files for consistent output order
            maf_files_sorted = sorted(maf_files)
            
            for i, maf_file in enumerate(maf_files_sorted, 1):
                if i % 100 == 0:
                    print(f"Processing {i}/{len(maf_files_sorted)}...")
                
                decoy_id = os.path.splitext(os.path.basename(maf_file))[0]
                
                # Skip if we don't have position data for this decoy
                if decoy_id not in v31_positions:
                    if args.debug:
                        print(f"Skipping {decoy_id}: no position data in v3.1 TSV", file=sys.stderr)
                    continue
                
                scores = process_maf_file(maf_file, args.msa_view, strand_map, predict_fn, 
                                         context, all_species, v31_positions, temp_dir, debug=args.debug)
                
                # Write immediately to save memory
                score_values = [format_score_output(scores[species][0], scores[species][1], 
                                                   scores[species][2], scores[species][3], scores[species][4]) 
                               for species in all_species]
                f.write(f"{decoy_id}\t" + "\t".join(score_values) + "\n")
                f.flush()  # Ensure data is written to disk
                
                processed_count += 1
                
                # Track files that failed (all scores are -1.0)
                if all(score[0] == -1.0 for score in scores.values()):
                    failed_files += 1
                
                # Periodic garbage collection to free memory (every 50 files)
                if i % 50 == 0:
                    gc.collect()
        
        print(f"✅ Done! Processed {processed_count} decoy IDs")
        if failed_files > 0:
            print(f"Note: {failed_files} MAF files were invalid/empty (assigned -1.0 for all species)")
        print(f"Output written to: {args.output_tsv}")

if __name__ == '__main__':
    main()
