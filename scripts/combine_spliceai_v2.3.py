#!/usr/bin/env python3
"""
Combine all SpliceAI score TSV files from splice_out_v2.3 into a single file.
"""

import os
import glob
from pathlib import Path

# Set up paths
input_dir = Path("splice_out_v2.3")
output_file = Path("all_spliceai_scores_v2.3.tsv")

# Get all TSV files, sorted
tsv_files = sorted(glob.glob(str(input_dir / "*_spliceai_scores.tsv")))

print(f"Found {len(tsv_files)} TSV files to combine")

if not tsv_files:
    print("No TSV files found!")
    exit(1)

# Read first file to get header
with open(tsv_files[0], 'r') as f:
    header = f.readline()

# Write combined file
total_rows = 0
with open(output_file, 'w') as outfile:
    # Write header
    outfile.write(header)
    
    # Process each file
    for i, tsv_file in enumerate(tsv_files, 1):
        print(f"Processing {i}/{len(tsv_files)}: {os.path.basename(tsv_file)}")
        
        with open(tsv_file, 'r') as infile:
            # Skip header line
            next(infile)
            
            # Copy all data rows
            file_rows = 0
            for line in infile:
                outfile.write(line)
                file_rows += 1
                total_rows += 1
        
        print(f"  Added {file_rows} rows")

print(f"\nCombined file created: {output_file}")
print(f"Total rows: {total_rows}")
print(f"Total files processed: {len(tsv_files)}")
