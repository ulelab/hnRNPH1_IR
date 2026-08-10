#!/bin/bash
# Script to split chromosomes into parts with ~50 rows per part
# Usage: ./split_chromosomes.sh <input_bed_file> <output_dir>

INPUT_BED="$1"
OUTPUT_DIR="$2"
ROWS_PER_PART=50

if [ -z "$INPUT_BED" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_bed_file> <output_dir>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Splitting chromosomes from $INPUT_BED"
echo "Output directory: $OUTPUT_DIR"
echo "Target: ~$ROWS_PER_PART rows per part"
echo ""

# Function to split a chromosome based on ~50 rows per part
split_chr() {
    local CHR=$1
    
    # Extract chromosome lines (match first column exactly to avoid matching chr10-19 when searching for chr1)
    awk -v chr="$CHR" '$1 == chr' "$INPUT_BED" > "${OUTPUT_DIR}/${CHR}_all.bed"
    local TOTAL_LINES=$(wc -l < "${OUTPUT_DIR}/${CHR}_all.bed")
    
    # Skip if no lines found
    if [ "$TOTAL_LINES" -eq 0 ]; then
        rm "${OUTPUT_DIR}/${CHR}_all.bed"
        return
    fi
    
    # Calculate number of parts needed (ceiling division: (total + rows_per_part - 1) / rows_per_part)
    local NUM_PARTS=$(( (TOTAL_LINES + ROWS_PER_PART - 1) / ROWS_PER_PART ))
    
    # Ensure at least 1 part
    if [ "$NUM_PARTS" -eq 0 ]; then
        NUM_PARTS=1
    fi
    
    # Calculate lines per part (ceiling division)
    local LINES_PER_PART=$(( (TOTAL_LINES + NUM_PARTS - 1) / NUM_PARTS ))
    
    echo "  Splitting $CHR ($TOTAL_LINES sites) into $NUM_PARTS parts (~$LINES_PER_PART sites each)"
    
    # Split into parts
    for ((i=1; i<=NUM_PARTS; i++)); do
        local START=$(( (i-1) * LINES_PER_PART + 1 ))
        local END=$(( i * LINES_PER_PART ))
        
        if [ $i -eq $NUM_PARTS ]; then
            # Last part gets remaining lines
            sed -n "${START},\$p" "${OUTPUT_DIR}/${CHR}_all.bed" > "${OUTPUT_DIR}/${CHR}_part${i}.bed"
        else
            sed -n "${START},${END}p" "${OUTPUT_DIR}/${CHR}_all.bed" > "${OUTPUT_DIR}/${CHR}_part${i}.bed"
        fi
        
        local PART_LINES=$(wc -l < "${OUTPUT_DIR}/${CHR}_part${i}.bed")
        echo "    Part $i: $PART_LINES sites -> ${OUTPUT_DIR}/${CHR}_part${i}.bed"
    done
    
    # Clean up temporary file
    rm "${OUTPUT_DIR}/${CHR}_all.bed"
}

# Get all unique chromosomes from the BED file
CHROMOSOMES=$(awk '{print $1}' "$INPUT_BED" | sort -u)

# Process each chromosome
for CHR in $CHROMOSOMES; do
    split_chr "$CHR"
    echo ""
done

echo "Done! Split files are in: $OUTPUT_DIR"
