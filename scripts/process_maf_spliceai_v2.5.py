#!/usr/bin/env python3
"""
Process MAF files: extract FASTA with msa_view, run SpliceAI, output TSV with top 2 scores and positions per species
Version 2.5: Length-adjusted canonical splice site detection for sequences with insertions/deletions
             (max1 search extended to first 100 nt instead of 60 nt)
"""

import os
import sys
import argparse
import subprocess
import tempfile
import gc
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

def parse_fasta_headers(fasta_file):
    """Parse FASTA headers to extract species IDs"""
    species_in_fasta = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # FASTA header format from msa_view: >species_id (no dot) or >species_id.chrom:start-end(strand)
                header = line[1:]  # Remove '>'
                # Extract species_id (everything before first '.' or first space, or whole header if no dot/space)
                if '.' in header:
                    species_id = header.split('.')[0]
                elif ' ' in header:
                    species_id = header.split()[0]
                else:
                    species_id = header  # Header is just the species ID
                # Only add non-empty species IDs
                if species_id:
                    species_in_fasta[species_id] = header
    return species_in_fasta

def find_top2_scores_length_adjusted(donor_probs, hg38_length, current_seq_length, debug=False):
    """
    Length-adjusted algorithm to find max1 (canonical splice site) and max2 (decoy splice site):
    
    Variables:
    - v = hg38_length (reference length)
    - n = current_seq_length (input sequence length)
    - |n| = |n - v| (absolute value of sequence length change)
    - d = |n| / v (normalized length change)
    
    max1 (Canonical Splice Site):
    - If n <= v: search positions [0, 99] (first 100 nt)
    - If n > v: search positions [0, adjusted_interval - 1] where adjusted_interval = 100 * (n / v)
    
    max2 (Decoy Splice Site):
    - Search from pos1 + 20 to end of sequence
    - If (pos1 + 20) >= sequence length, return -1 for both max2 and pos2
    
    Returns: (max1, pos1, max2, pos2, d) where positions are 0-based relative to sequence start
    """
    if len(donor_probs) == 0:
        return (-1.0, -1, -1.0, -1, -1.0)
    
    # Convert to numpy array if not already
    donor_probs = np.array(donor_probs)
    
    # Calculate normalized length change d = |n| / v
    abs_length_change = abs(current_seq_length - hg38_length)
    d = abs_length_change / hg38_length if hg38_length > 0 else 0.0
    
    # Step 1: Find max1 (canonical splice site) with length-adjusted interval
    if current_seq_length <= hg38_length:
        # Case 1: n <= v, search first 100 nt
        search_end = min(100, len(donor_probs))
        if debug:
            print(f"      max1 search: n={current_seq_length} <= v={hg38_length}, using interval [0, {search_end-1}]", file=sys.stderr)
    else:
        # Case 2: n > v, calculate adjusted interval
        adjusted_interval = int(100 * (current_seq_length / hg38_length))
        search_end = min(adjusted_interval, len(donor_probs))
        if debug:
            print(f"      max1 search: n={current_seq_length} > v={hg38_length}, adjusted_interval={adjusted_interval}, using [0, {search_end-1}]", file=sys.stderr)
    
    if search_end == 0:
        # Empty sequence
        return (-1.0, -1, -1.0, -1, d)
    
    canonical_region_scores = donor_probs[:search_end]
    if len(canonical_region_scores) == 0:
        max1 = -1.0
        pos1 = -1
    else:
        pos1 = int(np.argmax(canonical_region_scores))
        max1 = float(canonical_region_scores[pos1])
    
    # Step 2: Find max2 (decoy splice site) starting from pos1 + 20
    search_start = pos1 + 20
    if search_start >= len(donor_probs):
        # pos1 + 20 exceeds sequence length
        max2 = -1.0
        pos2 = -1
        if debug:
            print(f"      max2 search: pos1={pos1}, search_start={search_start} >= seq_len={len(donor_probs)}, returning -1", file=sys.stderr)
    else:
        # Search remaining sequence from pos1 + 20 to end
        decoy_region_scores = donor_probs[search_start:]
        if len(decoy_region_scores) == 0:
            max2 = -1.0
            pos2 = -1
        else:
            pos2_relative = int(np.argmax(decoy_region_scores))
            pos2 = search_start + pos2_relative  # Convert to absolute position
            max2 = float(decoy_region_scores[pos2_relative])
        if debug:
            print(f"      max2 search: pos1={pos1}, search_start={search_start}, found max2={max2:.6f}@pos{pos2}", file=sys.stderr)
    
    if debug:
        print(f"      Length-adjusted scores: max1={max1:.6f}@pos{pos1}, max2={max2:.6f}@pos{pos2}, d={d:.6f}", file=sys.stderr)
    
    return (max1, pos1, max2, pos2, d)

def process_sequence(sequence, species_id, predict_fn, context, fasta_file, strand, hg38_length, debug=False):
    """Process a single sequence with SpliceAI, return top 2 scores, positions, and length change"""
    # Remove all gap characters: hyphens, spaces, and asterisks (which can appear in alignments)
    sequence_clean = sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
    
    if debug:
        print(f"    Processing {species_id}: original_len={len(sequence)}, clean_len={len(sequence_clean)}, strand={strand}, hg38_len={hg38_length}", file=sys.stderr)
        if len(sequence_clean) > 0:
            # Show preview of first 50 bases (for display only - full sequence is used)
            preview_len = min(50, len(sequence_clean))
            print(f"      Preview (first {preview_len} of {len(sequence_clean)} bases): {sequence_clean[:preview_len]}", file=sys.stderr)
    
    # Skip if sequence is too short or contains only gaps/N's
    if len(sequence_clean) == 0:
        if debug:
            print(f"    Skipping {species_id}: empty sequence after removing gaps", file=sys.stderr)
        return (-1.0, -1, -1.0, -1, -1.0)
    
    # Check if sequence has any actual nucleotides (not just N's)
    sequence_no_n = sequence_clean.replace('N', '').replace('n', '')
    if len(sequence_no_n) == 0:
        if debug:
            print(f"    Skipping {species_id}: sequence contains only N's", file=sys.stderr)
        return (-1.0, -1, -1.0, -1, -1.0)
    
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
        
        # Find top 2 scores with length-adjusted logic
        max1, pos1, max2, pos2, d = find_top2_scores_length_adjusted(donor_probs, hg38_length, len(sequence_clean), debug=debug)
        
        if debug:
            print(f"    {species_id}: max1={max1:.6f}@pos{pos1}, max2={max2:.6f}@pos{pos2}, d={d:.6f}", file=sys.stderr)
        
        return (max1, pos1, max2, pos2, d)
        
    except Exception as e:
        print(f"WARNING: SpliceAI failed for {species_id} in {fasta_file}: {e}", file=sys.stderr)
        return (-1.0, -1, -1.0, -1, -1.0)

def get_sequence_length(fasta_file, species_id):
    """Get the length of a sequence for a specific species from FASTA file"""
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
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
                        sequence = ''.join(current_sequence)
                        sequence_clean = sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
                        return len(sequence_clean)
                
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
                sequence = ''.join(current_sequence)
                sequence_clean = sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
                return len(sequence_clean)
    
    return None  # Species not found

def run_spliceai_on_fasta(fasta_file, predict_fn, context, strand, hg38_length, debug=False):
    """Run SpliceAI on FASTA file and return top 2 scores, positions, and length change per species"""
    max_scores = {}
    
    if debug:
        print(f"  Reading FASTA file: {fasta_file}", file=sys.stderr)
        print(f"  Using hg38_length={hg38_length} as reference", file=sys.stderr)
    
    # Parse FASTA file properly (handles multi-line sequences)
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines
            if not line:
                continue
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None and current_sequence:
                    sequence = ''.join(current_sequence)
                    # Extract species_id from header (handle both formats)
                    if '.' in current_header:
                        species_id = current_header.split('.')[0]
                    elif ' ' in current_header:
                        species_id = current_header.split()[0]
                    else:
                        species_id = current_header  # Header is just the species ID
                    # Only process if we have a valid species ID
                    if species_id:
                        max_scores[species_id] = process_sequence(sequence, species_id, predict_fn, context, fasta_file, strand, hg38_length, debug=debug)
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                # Add to current sequence (don't uppercase here, let process_sequence handle it)
                if current_header is not None:
                    current_sequence.append(line)
        
        # Process last sequence
        if current_header is not None and current_sequence:
            sequence = ''.join(current_sequence)
            # Extract species_id from header (handle both formats)
            if '.' in current_header:
                species_id = current_header.split('.')[0]
            elif ' ' in current_header:
                species_id = current_header.split()[0]
            else:
                species_id = current_header  # Header is just the species ID
            # Only process if we have a valid species ID
            if species_id:
                max_scores[species_id] = process_sequence(sequence, species_id, predict_fn, context, fasta_file, strand, hg38_length, debug=debug)
    
    if debug:
        print(f"  Processed {len(max_scores)} species from FASTA", file=sys.stderr)
    
    return max_scores

def process_maf_file(maf_file, msa_view_bin, strand_map, predict_fn, context, all_species, temp_dir, debug=False):
    """Process a single MAF file: extract FASTA, run SpliceAI, return scores"""
    # Get decoy ID from MAF filename (without .maf extension)
    decoy_id = os.path.splitext(os.path.basename(maf_file))[0]
    
    # Look up strand from decoyIDs.bed
    strand = strand_map.get(decoy_id, '+')  # Default to '+' if not found
    
    if debug:
        print(f"  Processing {decoy_id}, strand: {strand}", file=sys.stderr)
    
    # Extract FASTA from MAF
    fasta_file = extract_fasta_from_maf(maf_file, msa_view_bin, strand, temp_dir, debug=debug)
    if fasta_file is None:
        # Return -1 for all species if extraction failed (invalid/empty MAF file)
        return {species: (-1.0, -1, -1.0, -1, -1.0) for species in all_species}
    
    # Check if FASTA file has content
    if os.path.getsize(fasta_file) == 0:
        if debug:
            print(f"  WARNING: Empty FASTA file for {decoy_id}", file=sys.stderr)
        return {species: (-1.0, -1, -1.0, -1, -1.0) for species in all_species}
    
    if debug:
        print(f"  FASTA file size: {os.path.getsize(fasta_file)} bytes", file=sys.stderr)
        # Show first few lines of FASTA
        with open(fasta_file, 'r') as f:
            lines = f.readlines()[:6]
            print(f"  First few FASTA lines:", file=sys.stderr)
            for line in lines:
                print(f"    {line.strip()}", file=sys.stderr)
    
    # Parse FASTA to get species present
    species_in_fasta = parse_fasta_headers(fasta_file)
    if debug:
        print(f"  Species found in FASTA: {len(species_in_fasta)}", file=sys.stderr)
        if len(species_in_fasta) > 0:
            print(f"  Example species: {list(species_in_fasta.keys())[:5]}", file=sys.stderr)
    
    # First, get hg38 sequence length (v) - process hg38 first
    hg38_length = None
    if 'hg38' in species_in_fasta:
        hg38_length = get_sequence_length(fasta_file, 'hg38')
        if hg38_length is None or hg38_length == 0:
            if debug:
                print(f"  WARNING: Could not determine hg38 length, using default 200", file=sys.stderr)
            hg38_length = 200  # Default fallback
        else:
            if debug:
                print(f"  hg38 sequence length (v): {hg38_length}", file=sys.stderr)
    else:
        if debug:
            print(f"  WARNING: hg38 not found in FASTA, using default length 200", file=sys.stderr)
        hg38_length = 200  # Default fallback if hg38 not present
    
    # Run SpliceAI with hg38_length as reference
    max_scores = run_spliceai_on_fasta(fasta_file, predict_fn, context, strand, hg38_length, debug=debug)
    if debug:
        print(f"  SpliceAI scores computed: {len(max_scores)}", file=sys.stderr)
        if len(max_scores) > 0:
            example_scores = {k: v for k, v in list(max_scores.items())[:3]}
            print(f"  Example scores: {example_scores}", file=sys.stderr)
        else:
            print(f"  WARNING: No scores computed! Check FASTA file content.", file=sys.stderr)
    
    # Create result dictionary with all species
    result = {}
    for species in all_species:
        if species in max_scores:
            result[species] = max_scores[species]
        else:
            result[species] = (-1.0, -1, -1.0, -1, -1.0)  # Species not present in FASTA
    
    return result

def format_score_output(max1, pos1, max2, pos2, d):
    """Format output as [max1][pos1][max2][pos2][d]"""
    return f"[{max1:.6f}][{pos1}][{max2:.6f}][{pos2}][{d:.6f}]"

def main():
    parser = argparse.ArgumentParser(description='Process MAF files with msa_view and SpliceAI (v2.5: length-adjusted canonical splice site detection, max1 in first 100 nt)')
    parser.add_argument('--maf_dir', required=True, help='Directory containing MAF files (or single MAF file for testing)')
    parser.add_argument('--decoy_bed', required=True, help='BED file with decoy IDs and strand info')
    parser.add_argument('--species_map', required=True, help='TSV file with species mapping')
    parser.add_argument('--msa_view', required=True, help='Path to msa_view binary')
    parser.add_argument('--output_tsv', required=True, help='Output TSV file')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
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
                scores = process_maf_file(maf_file, args.msa_view, strand_map, predict_fn, 
                                         context, all_species, temp_dir, debug=args.debug)
                
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
