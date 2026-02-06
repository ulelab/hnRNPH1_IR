#!/usr/bin/env python3
"""
Parse phyloP SPH output and calculate summary statistics.

This script processes phyloP SPH output files and calculates:
- Mean phyloP score
- Max phyloP score
- Fraction of constrained sites (phyloP > 1.3)
- Sum of positive phyloP values
- Number of nucleotides with phyloP >= 300
- Fraction of nucleotides with phyloP >= 300 (normalized by sequence length)
"""

import sys
import os
import glob
import pandas as pd
import numpy as np
import argparse

def parse_phylop_sph_file(phylop_file):
    """
    Parse a phyloP SPH output file and extract scores.
    
    phyloP SPH output with --wig-scores format:
    - Header line: "fixedStep chrom=... start=... step=..."
    - One score per line (can be negative for acceleration, positive for conservation)
    - Sentinel value 999999999.000 indicates missing data/gaps
    
    Returns a numpy array of phyloP scores (excluding sentinel values).
    """
    scores = []
    
    try:
        with open(phylop_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines, comments, and warnings
                if not line or line.startswith('#') or line.startswith('WARNING'):
                    continue
                
                # Skip wig format header lines
                if line.startswith('fixedStep') or line.startswith('variableStep'):
                    continue
                
                # Try to extract numeric value
                try:
                    score = float(line)
                    # Filter out NaN, Inf, and sentinel values (999999999.000)
                    if not np.isnan(score) and not np.isinf(score) and abs(score) < 1e8:
                        scores.append(score)
                except ValueError:
                    # Not a number, skip
                    continue
        
        return np.array(scores)
    
    except Exception as e:
        print(f"ERROR reading {phylop_file}: {e}", file=sys.stderr)
        return np.array([])

def calculate_metrics(scores):
    """
    Calculate summary statistics from phyloP scores.
    
    Returns a dictionary with:
    - mean_phylop: mean phyloP score
    - max_phylop: maximum phyloP score
    - fraction_constrained: fraction of sites with phyloP > 1.3
    - fraction_constrained_p01: fraction of sites with phyloP > 2 (-log10(0.01))
    - sum_positive_phylop: sum of all positive phyloP values
    - n_sites_ge_300: number of nucleotides with phyloP >= 300
    - fraction_sites_ge_300: fraction of nucleotides with phyloP >= 300 (normalized by sequence length)
    """
    if len(scores) == 0:
        return {
            'mean_phylop': np.nan,
            'max_phylop': np.nan,
            'fraction_constrained': np.nan,
            'fraction_constrained_p01': np.nan,
            'sum_positive_phylop': np.nan,
            'n_sites_ge_300': np.nan,
            'fraction_sites_ge_300': np.nan
        }
    
    mean_phylop = np.mean(scores)
    max_phylop = np.max(scores)
    
    # Fraction of constrained sites (phyloP > 1.3)
    n_constrained = np.sum(scores > 1.3)
    fraction_constrained = n_constrained / len(scores)
    
    # Fraction of constrained sites (phyloP > 2, which is -log10(0.01))
    n_constrained_p01 = np.sum(scores > 2.0)
    fraction_constrained_p01 = n_constrained_p01 / len(scores)
    
    # Sum of positive phyloP values
    sum_positive = np.sum(scores[scores > 0])
    
    # Number of nucleotides with phyloP >= 300
    n_sites_ge_300 = np.sum(scores >= 300)
    
    # Fraction of nucleotides with phyloP >= 300 (normalized by sequence length)
    fraction_sites_ge_300 = n_sites_ge_300 / len(scores)
    
    return {
        'mean_phylop': mean_phylop,
        'max_phylop': max_phylop,
        'fraction_constrained': fraction_constrained,
        'fraction_constrained_p01': fraction_constrained_p01,
        'sum_positive_phylop': sum_positive,
        'n_sites_ge_300': n_sites_ge_300,
        'fraction_sites_ge_300': fraction_sites_ge_300
    }

def main():
    parser = argparse.ArgumentParser(
        description='Parse phyloP SPH output files and calculate summary statistics'
    )
    parser.add_argument('input_dir', 
                       help='Directory containing phyloP output files (*.phylop)')
    parser.add_argument('output_tsv',
                       help='Output TSV file with summary statistics')
    
    args = parser.parse_args()
    
    # Find all phyloP output files
    phylop_files = glob.glob(os.path.join(args.input_dir, "*.phylop"))
    
    if len(phylop_files) == 0:
        print(f"WARNING: No .phylop files found in {args.input_dir}", file=sys.stderr)
        return 1
    
    print(f"Processing {len(phylop_files)} phyloP output files...")
    
    results = []
    
    for phylop_file in sorted(phylop_files):
        maf_name = os.path.basename(phylop_file).replace('.phylop', '')
        
        # Parse phyloP scores
        scores = parse_phylop_sph_file(phylop_file)
        
        if len(scores) == 0:
            print(f"WARNING: No scores found in {maf_name}", file=sys.stderr)
            continue
        
        # Calculate metrics
        metrics = calculate_metrics(scores)
        
        results.append({
            'MAF_file': maf_name,
            'mean_phylop': metrics['mean_phylop'],
            'max_phylop': metrics['max_phylop'],
            'fraction_constrained': metrics['fraction_constrained'],
            'fraction_constrained_p01': metrics['fraction_constrained_p01'],
            'sum_positive_phylop': metrics['sum_positive_phylop'],
            'n_sites_ge_300': metrics['n_sites_ge_300'],
            'fraction_sites_ge_300': metrics['fraction_sites_ge_300']
        })
    
    # Create DataFrame and save
    if results:
        df = pd.DataFrame(results)
        df.to_csv(args.output_tsv, sep='\t', index=False)
        print(f"Created summary TSV with {len(results)} entries: {args.output_tsv}")
        return 0
    else:
        print("ERROR: No results to write", file=sys.stderr)
        return 1

if __name__ == '__main__':
    sys.exit(main())

