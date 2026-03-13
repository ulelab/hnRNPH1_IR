#!/usr/bin/env python3
"""
Process BED file: extract FASTA with bedtools getfasta, run SpliceAI, output TSV with max1, pos1, max2, pos2
Version: Simplified version that processes hg38 sequences directly from BED file (no MAF alignment step)
"""

import os
import sys
import argparse
import subprocess
import tempfile
import gc
import numpy as np
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

# Silence TensorFlow output
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
tf.get_logger().setLevel('ERROR')

def load_spliceai_models():
    """Load SpliceAI models"""
    context = 10000
    paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
    models = [load_model(resource_filename('spliceai', x), compile=False) for x in paths]
    predict_fn = [model.predict for model in models]
    return predict_fn, context

def find_top2_scores(donor_probs, debug=False):
    """
    Find max1 (canonical splice site) and max2 (decoy splice site):
    
    max1 (Canonical Splice Site):
    - Search positions [0, 99] (first 100 nt)
    
    max2 (Decoy Splice Site):
    - Search from pos1 + 20 to end of sequence
    - If (pos1 + 20) >= sequence length, return -1 for both max2 and pos2
    
    Returns: (max1, pos1, max2, pos2) where positions are 0-based relative to sequence start
    """
    if len(donor_probs) == 0:
        return (-1.0, -1, -1.0, -1)
    
    # Convert to numpy array if not already
    donor_probs = np.array(donor_probs)
    
    # Step 1: Find max1 (canonical splice site) in first 100 nt
    search_end = min(100, len(donor_probs))
    
    if search_end == 0:
        # Empty sequence
        return (-1.0, -1, -1.0, -1)
    
    canonical_region_scores = donor_probs[:search_end]
    if len(canonical_region_scores) == 0:
        max1 = -1.0
        pos1 = -1
    else:
        pos1 = int(np.argmax(canonical_region_scores))
        max1 = float(canonical_region_scores[pos1])
    
    if debug:
        print(f"      max1 search: using interval [0, {search_end-1}], found max1={max1:.6f}@pos{pos1}", file=sys.stderr)
    
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
        print(f"      Scores: max1={max1:.6f}@pos{pos1}, max2={max2:.6f}@pos{pos2}", file=sys.stderr)
    
    return (max1, pos1, max2, pos2)

def process_sequence(sequence, decoy_id, predict_fn, context, debug=False):
    """Process a single sequence with SpliceAI, return max1, pos1, max2, pos2"""
    # Remove all gap characters and convert to uppercase
    sequence_clean = sequence.replace('-', '').replace(' ', '').replace('*', '').upper()
    
    if debug:
        print(f"    Processing {decoy_id}: clean_len={len(sequence_clean)}", file=sys.stderr)
        if len(sequence_clean) > 0:
            preview_len = min(50, len(sequence_clean))
            print(f"      Preview (first {preview_len} of {len(sequence_clean)} bases): {sequence_clean[:preview_len]}", file=sys.stderr)
    
    # Skip if sequence is too short or contains only gaps/N's
    if len(sequence_clean) == 0:
        if debug:
            print(f"    Skipping {decoy_id}: empty sequence after removing gaps", file=sys.stderr)
        return (-1.0, -1, -1.0, -1)
    
    # Check if sequence has any actual nucleotides (not just N's)
    sequence_no_n = sequence_clean.replace('N', '').replace('n', '')
    if len(sequence_no_n) == 0:
        if debug:
            print(f"    Skipping {decoy_id}: sequence contains only N's", file=sys.stderr)
        return (-1.0, -1, -1.0, -1)
    
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
        
        # Find top 2 scores
        max1, pos1, max2, pos2 = find_top2_scores(donor_probs, debug=debug)
        
        if debug:
            print(f"    {decoy_id}: max1={max1:.6f}@pos{pos1}, max2={max2:.6f}@pos{pos2}", file=sys.stderr)
        
        return (max1, pos1, max2, pos2)
        
    except Exception as e:
        print(f"WARNING: SpliceAI failed for {decoy_id}: {e}", file=sys.stderr)
        return (-1.0, -1, -1.0, -1)

def extract_fasta_from_bed(bed_file, genome_fasta, temp_dir, debug=False):
    """Extract FASTA from BED file using bedtools getfasta with strand-aware logic"""
    fasta_file = os.path.join(temp_dir, "extracted_sequences.fa")
    
    # Build bedtools getfasta command
    # -s flag: strand-aware (reverse complement for negative strand)
    # -name flag: use BED name column (col4) as FASTA header
    cmd = ['bedtools', 'getfasta',
           '-fi', genome_fasta,
           '-bed', bed_file,
           '-s',  # strand-aware
           '-name',  # use name column as header
           '-fo', fasta_file]
    
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                              check=True, text=True)
        
        # Verify the file was created and has content
        if not os.path.exists(fasta_file):
            if debug:
                print(f"WARNING: FASTA file was not created: {fasta_file}", file=sys.stderr)
            return None
        
        file_size = os.path.getsize(fasta_file)
        if file_size == 0:
            if debug:
                print(f"WARNING: Empty FASTA file created (size: {file_size} bytes)", file=sys.stderr)
            return None
        
        if debug:
            print(f"  Successfully extracted FASTA: {fasta_file} (size: {file_size} bytes)", file=sys.stderr)
        
        return fasta_file
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.strip() if e.stderr else "No error message"
        print(f"WARNING: bedtools getfasta failed: {error_msg}", file=sys.stderr)
        return None
    except Exception as e:
        if debug:
            print(f"WARNING: Unexpected error extracting FASTA: {e}", file=sys.stderr)
        return None

def parse_fasta(fasta_file, debug=False):
    """Parse FASTA file and return dictionary of decoy_id -> sequence"""
    sequences = {}
    current_header = None
    current_sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None and current_sequence:
                    sequences[current_header] = ''.join(current_sequence)
                
                # Start new sequence
                # Header format from bedtools: >name (from BED col4)
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                # Add to current sequence
                if current_header is not None:
                    current_sequence.append(line)
        
        # Process last sequence
        if current_header is not None and current_sequence:
            sequences[current_header] = ''.join(current_sequence)
    
    if debug:
        print(f"  Parsed {len(sequences)} sequences from FASTA", file=sys.stderr)
    
    return sequences

def process_bed_file(bed_file, genome_fasta, predict_fn, context, temp_dir, debug=False):
    """Process BED file: extract FASTA, run SpliceAI, return scores per decoy"""
    # Extract FASTA from BED
    fasta_file = extract_fasta_from_bed(bed_file, genome_fasta, temp_dir, debug=debug)
    if fasta_file is None:
        return {}
    
    # Parse FASTA
    sequences = parse_fasta(fasta_file, debug=debug)
    
    # Process each sequence with SpliceAI
    results = {}
    for decoy_id, sequence in sequences.items():
        if debug:
            print(f"  Processing {decoy_id}...", file=sys.stderr)
        max1, pos1, max2, pos2 = process_sequence(sequence, decoy_id, predict_fn, context, debug=debug)
        results[decoy_id] = (max1, pos1, max2, pos2)
    
    return results

def format_score_output(max1, pos1, max2, pos2):
    """Format output as [max1][pos1][max2][pos2]"""
    return f"[{max1:.6f}][{pos1}][{max2:.6f}][{pos2}]"

def main():
    parser = argparse.ArgumentParser(description='Process BED file with bedtools getfasta and SpliceAI')
    parser.add_argument('--bed_file', required=True, help='BED file with decoy regions')
    parser.add_argument('--genome_fasta', required=True, help='Path to hg38 genome FASTA file')
    parser.add_argument('--output_tsv', required=True, help='Output TSV file')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    # Check if bedtools is available
    try:
        subprocess.run(['bedtools', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: bedtools not found. Please install bedtools.", file=sys.stderr)
        sys.exit(1)
    
    # Check if genome FASTA exists
    if not os.path.exists(args.genome_fasta):
        print(f"ERROR: Genome FASTA file not found: {args.genome_fasta}", file=sys.stderr)
        sys.exit(1)
    
    # Check if BED file exists
    if not os.path.exists(args.bed_file):
        print(f"ERROR: BED file not found: {args.bed_file}", file=sys.stderr)
        sys.exit(1)
    
    print("Loading SpliceAI models...")
    predict_fn, context = load_spliceai_models()
    print("SpliceAI models loaded")
    
    # Create temporary directory for FASTA files
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Processing BED file: {args.bed_file}")
        print(f"Using genome FASTA: {args.genome_fasta}")
        
        # Process BED file
        results = process_bed_file(args.bed_file, args.genome_fasta, predict_fn, context, temp_dir, debug=args.debug)
        
        print(f"Processed {len(results)} decoy IDs")
        
        # Write output TSV
        print(f"Writing results to: {args.output_tsv}")
        with open(args.output_tsv, 'w') as f:
            # Write header: decoyID, hg38
            f.write("decoyID\thg38\n")
            
            # Write results
            for decoy_id in sorted(results.keys()):
                max1, pos1, max2, pos2 = results[decoy_id]
                score_str = format_score_output(max1, pos1, max2, pos2)
                f.write(f"{decoy_id}\t{score_str}\n")
        
        print(f"✅ Done! Wrote {len(results)} entries to {args.output_tsv}")

if __name__ == '__main__':
    main()
