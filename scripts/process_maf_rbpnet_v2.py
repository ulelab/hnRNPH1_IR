#!/usr/bin/env python3
"""
Process MAF files: extract FASTA with msa_view, run RBPnet, output TSV with top scores and positions per species
Output format:
- Top 3 scores from first 100 nt (5' splice site): [css1][cps1][css2][cps2][css3][cps3]
- Top 3 scores from last 100 nt (3' splice site): [css1_3p][cps1_3p][css2_3p][cps2_3p][css3_3p][cps3_3p]
- All scores > 0.02 between position 100 and (length - 100): [int1][ips1][int2][ips2][int3][ips3][int*][ips*] (variable length)
"""

import os
import sys
import argparse
import subprocess
import tempfile
import gc
import numpy as np
from collections import defaultdict

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

def parse_rbpnet_tsv(rbpnet_tsv_file, debug=False):
    """Parse RBPnet TSV output and extract profile_target values
    Returns: list of (position, score) tuples, or None if parsing fails
    """
    try:
        with open(rbpnet_tsv_file, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        
        if len(lines) < 5:
            if debug:
                print(f"  WARNING: RBPnet TSV has fewer than 5 lines", file=sys.stderr)
            return None
        
        # Line 5 (index 4) contains profile_target values (space-separated)
        profile_target_line = lines[4]
        profile_target_values = [float(x) for x in profile_target_line.split()]
        
        # Create list of (position, score) tuples (0-based positions)
        profile_data = [(pos, score) for pos, score in enumerate(profile_target_values)]
        
        if debug:
            print(f"  Parsed {len(profile_data)} profile_target values", file=sys.stderr)
        
        return profile_data
        
    except Exception as e:
        if debug:
            print(f"  ERROR parsing RBPnet TSV: {e}", file=sys.stderr)
        return None

def find_top3_canonical_scores(profile_data, debug=False):
    """
    Find top 3 scores in first 100 nucleotides (5' splice site)
    Returns: list of (score, position) tuples, sorted by score descending
    """
    if profile_data is None or len(profile_data) == 0:
        return [(-1.0, -1), (-1.0, -1), (-1.0, -1)]
    
    # Filter to first 100 positions
    first_100 = [(pos, score) for pos, score in profile_data if pos < 100]
    
    if len(first_100) == 0:
        return [(-1.0, -1), (-1.0, -1), (-1.0, -1)]
    
    # Sort by score descending and take top 3
    first_100_sorted = sorted(first_100, key=lambda x: x[1], reverse=True)
    top3 = first_100_sorted[:3]
    
    # Pad to 3 if needed
    while len(top3) < 3:
        top3.append((-1.0, -1))
    
    # Return as (score, position) tuples
    result = [(score, pos) for pos, score in top3]
    
    if debug:
        print(f"    Top 3 canonical scores (5' splice site, first 100nt): {result}", file=sys.stderr)
    
    return result

def find_top3_canonical_scores_3p(profile_data, debug=False):
    """
    Find top 3 scores in last 100 nucleotides (3' splice site)
    Returns: list of (score, position) tuples, sorted by score descending
    """
    if profile_data is None or len(profile_data) == 0:
        return [(-1.0, -1), (-1.0, -1), (-1.0, -1)]
    
    # Get sequence length
    seq_length = len(profile_data)
    
    if seq_length <= 100:
        # If sequence is <= 100 nt, return -1 for all
        return [(-1.0, -1), (-1.0, -1), (-1.0, -1)]
    
    # Filter to last 100 positions
    last_100_start = seq_length - 100
    last_100 = [(pos, score) for pos, score in profile_data if pos >= last_100_start]
    
    if len(last_100) == 0:
        return [(-1.0, -1), (-1.0, -1), (-1.0, -1)]
    
    # Sort by score descending and take top 3
    last_100_sorted = sorted(last_100, key=lambda x: x[1], reverse=True)
    top3 = last_100_sorted[:3]
    
    # Pad to 3 if needed
    while len(top3) < 3:
        top3.append((-1.0, -1))
    
    # Return as (score, position) tuples
    result = [(score, pos) for pos, score in top3]
    
    if debug:
        print(f"    Top 3 canonical scores (3' splice site, last 100nt): {result}", file=sys.stderr)
    
    return result

def find_intronic_scores(profile_data, threshold=0.02, debug=False):
    """
    Find all scores > threshold between position 100 and (length - 100)
    Returns: list of (score, position) tuples, sorted by position ascending
    """
    if profile_data is None or len(profile_data) == 0:
        return []
    
    # Get sequence length
    seq_length = len(profile_data)
    
    # Calculate end position for intronic region (exclude last 100 nt)
    intronic_end = seq_length - 100
    
    # If sequence is too short (<= 200 nt), no intronic region
    if seq_length <= 200:
        if debug:
            print(f"    Sequence too short ({seq_length} nt) for intronic region (need > 200 nt)", file=sys.stderr)
        return []
    
    # Filter to positions >= 100 and < (length - 100) with score > threshold
    intronic = [(pos, score) for pos, score in profile_data 
                if pos >= 100 and pos < intronic_end and score > threshold]
    
    # Sort by position ascending
    intronic_sorted = sorted(intronic, key=lambda x: x[0])
    
    # Return as (score, position) tuples
    result = [(score, pos) for pos, score in intronic_sorted]
    
    if debug:
        if len(result) > 0:
            print(f"    Found {len(result)} intronic scores > {threshold} (100-{intronic_end}nt): {result[:5]}...", file=sys.stderr)
        else:
            print(f"    No intronic scores > {threshold} found in region 100-{intronic_end}nt", file=sys.stderr)
    
    return result

def format_rbpnet_output(canonical_scores_5p, canonical_scores_3p, intronic_scores):
    """
    Format output as:
    [css1][cps1][css2][cps2][css3][cps3][css1_3p][cps1_3p][css2_3p][cps2_3p][css3_3p][cps3_3p][int1][ips1][int2][ips2][int3][ips3][int*][ips*]...
    """
    parts = []
    
    # Add top 3 canonical scores from 5' splice site (first 100 nt)
    for score, pos in canonical_scores_5p:
        parts.append(f"[{score:.6f}][{pos}]")
    
    # Add top 3 canonical scores from 3' splice site (last 100 nt)
    for score, pos in canonical_scores_3p:
        parts.append(f"[{score:.6f}][{pos}]")
    
    # Add intronic scores (variable length)
    for score, pos in intronic_scores:
        parts.append(f"[{score:.6f}][{pos}]")
    
    # If no data, return all -1 (6 values for 5' and 3' canonical scores)
    if len(parts) == 0 or (all(score == -1.0 for score, _ in canonical_scores_5p) and 
                           all(score == -1.0 for score, _ in canonical_scores_3p)):
        return "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"
    
    return "".join(parts)

def run_rbpnet_on_fasta(fasta_file, rbpnet_bin, rbpnet_model, temp_dir, debug=False):
    """Run RBPnet on FASTA file and return scores per species"""
    species_scores = {}
    
    if debug:
        print(f"  Reading FASTA file: {fasta_file}", file=sys.stderr)
    
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
                    # Extract species_id from header
                    if '.' in current_header:
                        species_id = current_header.split('.')[0]
                    elif ' ' in current_header:
                        species_id = current_header.split()[0]
                    else:
                        species_id = current_header
                    
                    # Only process if we have a valid species ID
                    if species_id:
                        species_scores[species_id] = process_sequence_with_rbpnet(
                            sequence, species_id, rbpnet_bin, rbpnet_model, 
                            temp_dir, fasta_file, debug=debug)
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                # Add to current sequence
                if current_header is not None:
                    current_sequence.append(line)
        
        # Process last sequence
        if current_header is not None and current_sequence:
            sequence = ''.join(current_sequence)
            # Extract species_id from header
            if '.' in current_header:
                species_id = current_header.split('.')[0]
            elif ' ' in current_header:
                species_id = current_header.split()[0]
            else:
                species_id = current_header
            # Only process if we have a valid species ID
            if species_id:
                species_scores[species_id] = process_sequence_with_rbpnet(
                    sequence, species_id, rbpnet_bin, rbpnet_model, 
                    temp_dir, fasta_file, debug=debug)
    
    if debug:
        print(f"  Processed {len(species_scores)} species from FASTA", file=sys.stderr)
    
    return species_scores

def process_sequence_with_rbpnet(sequence, species_id, rbpnet_bin, rbpnet_model, 
                                  temp_dir, fasta_file, debug=False):
    """Process a single sequence with RBPnet, return formatted scores"""
    # Remove all gap characters
    sequence_clean = sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
    
    if debug:
        print(f"    Processing {species_id}: clean_len={len(sequence_clean)}", file=sys.stderr)
    
    # Skip if sequence is too short or contains only gaps/N's
    if len(sequence_clean) == 0:
        if debug:
            print(f"    Skipping {species_id}: empty sequence after removing gaps", file=sys.stderr)
        return "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"
    
    # Check if sequence has any actual nucleotides (not just N's)
    sequence_no_n = sequence_clean.replace('N', '').replace('n', '')
    if len(sequence_no_n) == 0:
        if debug:
            print(f"    Skipping {species_id}: sequence contains only N's", file=sys.stderr)
        return "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"
    
    # Create temporary FASTA file for this sequence
    temp_fasta = os.path.join(temp_dir, f"{species_id}_temp.fa")
    temp_tsv = os.path.join(temp_dir, f"{species_id}_temp.tsv")
    
    try:
        # Write sequence to temporary FASTA
        with open(temp_fasta, 'w') as f:
            f.write(f">{species_id}\n{sequence_clean}\n")
        
        # Run RBPnet
        cmd = [rbpnet_bin, 'predict', 
               '-m', rbpnet_model,
               '-o', temp_tsv,
               '--format', 'fasta',
               temp_fasta]
        
        if debug:
            print(f"    Running RBPnet: {' '.join(cmd)}", file=sys.stderr)
        
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                              check=True, text=True)
        
        # Parse RBPnet output
        profile_data = parse_rbpnet_tsv(temp_tsv, debug=debug)
        
        if profile_data is None:
            return "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"
        
        # Find top 3 canonical scores (5' splice site, first 100 nt)
        canonical_scores_5p = find_top3_canonical_scores(profile_data, debug=debug)
        
        # Find top 3 canonical scores (3' splice site, last 100 nt)
        canonical_scores_3p = find_top3_canonical_scores_3p(profile_data, debug=debug)
        
        # Find intronic scores (> 0.02 between position 100 and (length - 100))
        intronic_scores = find_intronic_scores(profile_data, threshold=0.02, debug=debug)
        
        # Format output
        formatted = format_rbpnet_output(canonical_scores_5p, canonical_scores_3p, intronic_scores)
        
        if debug:
            print(f"    {species_id}: formatted output length={len(formatted)}", file=sys.stderr)
        
        return formatted
        
    except subprocess.CalledProcessError as e:
        if debug:
            print(f"    WARNING: RBPnet failed for {species_id}: {e.stderr}", file=sys.stderr)
        return "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"
    except Exception as e:
        if debug:
            print(f"    WARNING: Error processing {species_id}: {e}", file=sys.stderr)
        return "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"
    finally:
        # Clean up temporary files
        for temp_file in [temp_fasta, temp_tsv]:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except:
                    pass

def process_maf_file(maf_file, msa_view_bin, strand_map, rbpnet_bin, rbpnet_model, 
                     all_species, temp_dir, debug=False):
    """Process a single MAF file: extract FASTA, run RBPnet, return scores"""
    # Get decoy ID from MAF filename (without .maf extension)
    decoy_id = os.path.splitext(os.path.basename(maf_file))[0]
    
    # Look up strand from decoyIDs.bed
    strand = strand_map.get(decoy_id, '+')  # Default to '+' if not found
    
    if debug:
        print(f"  Processing {decoy_id}, strand: {strand}", file=sys.stderr)
    
    # Extract FASTA from MAF
    fasta_file = extract_fasta_from_maf(maf_file, msa_view_bin, strand, temp_dir, debug=debug)
    if fasta_file is None:
        # Return -1 for all species if extraction failed
        return {species: "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]" for species in all_species}
    
    # Check if FASTA file has content
    if os.path.getsize(fasta_file) == 0:
        if debug:
            print(f"  WARNING: Empty FASTA file for {decoy_id}", file=sys.stderr)
        return {species: "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]" for species in all_species}
    
    # Run RBPnet
    species_scores = run_rbpnet_on_fasta(fasta_file, rbpnet_bin, rbpnet_model, temp_dir, debug=debug)
    if debug:
        print(f"  RBPnet scores computed: {len(species_scores)}", file=sys.stderr)
    
    # Create result dictionary with all species
    result = {}
    for species in all_species:
        if species in species_scores:
            result[species] = species_scores[species]
        else:
            result[species] = "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]"  # Species not present in FASTA
    
    return result

def main():
    parser = argparse.ArgumentParser(description='Process MAF files with msa_view and RBPnet (v2)')
    parser.add_argument('--maf_dir', required=True, help='Directory containing MAF files (or single MAF file for testing)')
    parser.add_argument('--decoy_bed', required=True, help='BED file with decoy IDs and strand info')
    parser.add_argument('--species_map', required=True, help='TSV file with species mapping')
    parser.add_argument('--msa_view', required=True, help='Path to msa_view binary')
    parser.add_argument('--rbpnet', required=True, help='Path to rbpnet binary')
    parser.add_argument('--rbpnet_model', required=True, help='Path to RBPnet model file')
    parser.add_argument('--output_tsv', required=True, help='Output TSV file')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    print(f"Loading species list from: {args.species_map}")
    all_species = load_species_list(args.species_map)
    print(f"Loaded {len(all_species)} species")
    
    print(f"Loading strand map from: {args.decoy_bed}")
    strand_map = load_strand_map(args.decoy_bed)
    print(f"Loaded strand info for {len(strand_map)} decoy IDs")
    
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
            # Write header: decoyID, then each species
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
                scores = process_maf_file(maf_file, args.msa_view, strand_map, 
                                         args.rbpnet, args.rbpnet_model,
                                         all_species, temp_dir, debug=args.debug)
                
                # Write immediately to save memory
                score_values = [scores[species] for species in all_species]
                f.write(f"{decoy_id}\t" + "\t".join(score_values) + "\n")
                f.flush()  # Ensure data is written to disk
                
                processed_count += 1
                
                # Track files that failed (all scores are -1)
                if all(score == "[-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1][-1.0][-1]" for score in scores.values()):
                    failed_files += 1
                
                # Periodic garbage collection to free memory (every 50 files)
                if i % 50 == 0:
                    gc.collect()
        
        print(f"✅ Done! Processed {processed_count} decoy IDs")
        if failed_files > 0:
            print(f"Note: {failed_files} MAF files were invalid/empty (assigned -1 for all species)")
        print(f"Output written to: {args.output_tsv}")

if __name__ == '__main__':
    main()
