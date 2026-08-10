#!/usr/bin/env python3
"""
Process MAF files: extract FASTA with msa_view, run SpliceAI, output TSV with max scores per species
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

def run_spliceai_on_fasta(fasta_file, predict_fn, context, debug=False):
    """Run SpliceAI on FASTA file and return max donor score per species"""
    max_scores = {}
    
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
                    # Extract species_id from header (handle both formats)
                    if '.' in current_header:
                        species_id = current_header.split('.')[0]
                    elif ' ' in current_header:
                        species_id = current_header.split()[0]
                    else:
                        species_id = current_header  # Header is just the species ID
                    # Only process if we have a valid species ID
                    if species_id:
                        max_scores[species_id] = process_sequence(sequence, species_id, predict_fn, context, fasta_file, debug=debug)
                
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
                max_scores[species_id] = process_sequence(sequence, species_id, predict_fn, context, fasta_file, debug=debug)
    
    if debug:
        print(f"  Processed {len(max_scores)} species from FASTA", file=sys.stderr)
    
    return max_scores

def process_sequence(sequence, species_id, predict_fn, context, fasta_file, debug=False):
    """Process a single sequence with SpliceAI"""
    # Remove all gap characters: hyphens, spaces, and asterisks (which can appear in alignments)
    sequence_clean = sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
    
    if debug:
        print(f"    Processing {species_id}: original_len={len(sequence)}, clean_len={len(sequence_clean)}", file=sys.stderr)
        if len(sequence_clean) > 0:
            # Show preview of first 50 bases (for display only - full sequence is used)
            preview_len = min(50, len(sequence_clean))
            print(f"      Preview (first {preview_len} of {len(sequence_clean)} bases): {sequence_clean[:preview_len]}", file=sys.stderr)
    
    # Skip if sequence is too short or contains only gaps/N's
    if len(sequence_clean) == 0:
        if debug:
            print(f"    Skipping {species_id}: empty sequence after removing gaps", file=sys.stderr)
        return -1.0
    
    # Check if sequence has any actual nucleotides (not just N's)
    sequence_no_n = sequence_clean.replace('N', '').replace('n', '')
    if len(sequence_no_n) == 0:
        if debug:
            print(f"    Skipping {species_id}: sequence contains only N's", file=sys.stderr)
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
        
        # Get max donor probability
        max_donor = float(np.max(donor_probs))
        if debug:
            print(f"    {species_id}: max_donor={max_donor:.6f}", file=sys.stderr)
        return max_donor
        
    except Exception as e:
        print(f"WARNING: SpliceAI failed for {species_id} in {fasta_file}: {e}", file=sys.stderr)
        return -1.0

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
        return {species: -1.0 for species in all_species}
    
    # Check if FASTA file has content
    if os.path.getsize(fasta_file) == 0:
        if debug:
            print(f"  WARNING: Empty FASTA file for {decoy_id}", file=sys.stderr)
        return {species: -1.0 for species in all_species}
    
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
    
    # Run SpliceAI
    max_scores = run_spliceai_on_fasta(fasta_file, predict_fn, context, debug=debug)
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
            result[species] = -1.0  # Species not present in FASTA
    
    return result

def main():
    parser = argparse.ArgumentParser(description='Process MAF files with msa_view and SpliceAI')
    parser.add_argument('--maf_dir', required=True, help='Directory containing MAF files')
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
    
    # Find all MAF files
    maf_files = []
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
            # Write header
            f.write("decoyID\t" + "\t".join(all_species) + "\n")
            
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
                score_values = [f"{scores[species]:.6f}" for species in all_species]
                f.write(f"{decoy_id}\t" + "\t".join(score_values) + "\n")
                f.flush()  # Ensure data is written to disk
                
                processed_count += 1
                
                # Track files that failed (all scores are -1.0)
                if all(score == -1.0 for score in scores.values()):
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

