#!/usr/bin/env python3
"""
Combine multiple TSV files from SpliceAI processing into one final file
"""

import os
import sys
import argparse
import pandas as pd
from glob import glob

def main():
    parser = argparse.ArgumentParser(description='Combine SpliceAI TSV files')
    parser.add_argument('--input_dir', required=True, 
                       help='Directory containing *_spliceai_scores.tsv files')
    parser.add_argument('--output_tsv', required=True,
                       help='Output combined TSV file')
    parser.add_argument('--pattern', default='*_spliceai_scores.tsv',
                       help='File pattern to match (default: *_spliceai_scores.tsv)')
    
    args = parser.parse_args()
    
    # Find all TSV files
    pattern = os.path.join(args.input_dir, args.pattern)
    tsv_files = glob(pattern)
    
    if not tsv_files:
        print(f"ERROR: No files found matching pattern: {pattern}")
        sys.exit(1)
    
    print(f"Found {len(tsv_files)} TSV files to combine")
    
    # Read and combine all TSV files
    all_dataframes = []
    
    for tsv_file in sorted(tsv_files):
        print(f"Reading: {os.path.basename(tsv_file)}")
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            all_dataframes.append(df)
            print(f"  Loaded {len(df)} rows, {len(df.columns)-1} species columns")
        except Exception as e:
            print(f"  ERROR reading {tsv_file}: {e}", file=sys.stderr)
            continue
    
    if not all_dataframes:
        print("ERROR: No dataframes loaded")
        sys.exit(1)
    
    # Combine all dataframes
    print("\nCombining dataframes...")
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    
    # Remove duplicates (in case same decoyID appears in multiple files)
    print(f"Before deduplication: {len(combined_df)} rows")
    combined_df = combined_df.drop_duplicates(subset=['decoyID'], keep='first')
    print(f"After deduplication: {len(combined_df)} rows")
    
    # Sort by decoyID
    combined_df = combined_df.sort_values('decoyID')
    
    # Write output with formatting to preserve decimal notation (not scientific)
    # Use float_format to ensure small numbers like 0.000005 are written as decimals
    # instead of scientific notation (5e-06)
    print(f"\nWriting combined TSV to: {args.output_tsv}")
    
    # Format all numeric columns to 6 decimal places, preserving decimal notation
    # We'll write manually to have full control over formatting
    with open(args.output_tsv, 'w') as f:
        # Write header
        f.write('\t'.join(combined_df.columns) + '\n')
        
        # Write data rows with formatted numbers
        for _, row in combined_df.iterrows():
            values = []
            for col in combined_df.columns:
                val = row[col]
                if col == 'decoyID':
                    values.append(str(val))
                else:
                    # Format numeric values to 6 decimal places, preserving decimal notation
                    if pd.isna(val):
                        values.append("-1.000000")
                    else:
                        # Use format to ensure decimal notation (not scientific)
                        values.append(f"{val:.6f}")
            f.write('\t'.join(values) + '\n')
    
    print(f"✅ Done!")
    print(f"   Total decoy IDs: {len(combined_df)}")
    print(f"   Total species columns: {len(combined_df.columns) - 1}")
    print(f"   Output file: {args.output_tsv}")

if __name__ == '__main__':
    main()

