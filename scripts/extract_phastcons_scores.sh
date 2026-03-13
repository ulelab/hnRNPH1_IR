#!/bin/bash
# extract_phastcons_scores.sh
# Extract phastCons 470-way conservation scores from bigWig for BED regions

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Go up from scripts/ to hnRNPH1_IR_MAF/
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIGWIG_FILE="${1:-hg38.phastCons470way.bw}"
BED_FILE="${2:-$PROJECT_DIR/data/DecoysControlsUniqueIDsforRecount_filtered_intron.bed}"
OUTPUT_DIR="$PROJECT_DIR/results/phastCons"
OUTPUT_TSV="${OUTPUT_DIR}/phastCons_scores.tsv"
OUTPUT_MERGED="${OUTPUT_DIR}/bed_with_phastCons.tsv"

# Binary location
ARCH=$(uname -m)
BINARY_PATH="$HOME/bin/$ARCH/bigWigAverageOverBed"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if binary exists
if [ ! -f "$BINARY_PATH" ]; then
    echo "Error: bigWigAverageOverBed not found at $BINARY_PATH"
    echo "Please download it first or ensure it's in your PATH"
    exit 1
fi

# Check if bigWig file exists
if [ ! -f "$BIGWIG_FILE" ]; then
    echo "Error: BigWig file not found: $BIGWIG_FILE"
    echo "Usage: $0 [bigwig_file] [bed_file]"
    echo "  bigwig_file: Path to hg38.phastCons470way.bw (default: hg38.phastCons470way.bw)"
    echo "  bed_file: Path to BED file (default: DecoysControlsUniqueIDsforRecount_filtered_intron.bed)"
    exit 1
fi

# Check if BED file exists
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found: $BED_FILE"
    exit 1
fi

echo "=========================================="
echo "Extracting phastCons scores"
echo "=========================================="
echo "BigWig: $BIGWIG_FILE"
echo "BED file: $BED_FILE"
echo "Output TSV: $OUTPUT_TSV"
echo "Merged output: $OUTPUT_MERGED"
echo "Binary: $BINARY_PATH"
echo ""

# Create a temporary BED file with only required columns (chr, start, end, name)
# bigWigAverageOverBed expects standard BED format and may have issues with float scores
# Remove duplicates (keep first occurrence) since bigWigAverageOverBed requires unique names
TEMP_BED=$(mktemp)
awk -F'\t' '
BEGIN {OFS="\t"}
{
    name = $4
    if (!seen[name]) {
        print $1, $2, $3, $4
        seen[name] = 1
    } else {
        duplicates++
        dup_names[name] = 1
    }
}
END {
    if (duplicates > 0) {
        print "Warning: Removed " duplicates " duplicate entries (keeping first occurrence)" > "/dev/stderr"
        print "Duplicate names:" > "/dev/stderr"
        for (name in dup_names) {
            print "  " name > "/dev/stderr"
        }
    }
}' "$BED_FILE" > "$TEMP_BED"

# Run bigWigAverageOverBed
# Output format: name, size, covered, sum, mean0, mean
echo "Running bigWigAverageOverBed..."
"$BINARY_PATH" "$BIGWIG_FILE" "$TEMP_BED" "$OUTPUT_TSV"

if [ $? -ne 0 ]; then
    echo "Error: bigWigAverageOverBed failed"
    rm -f "$TEMP_BED"
    exit 1
fi

# Clean up temporary file
rm -f "$TEMP_BED"

echo "Successfully extracted scores to: $OUTPUT_TSV"
echo ""

# Merge results back with original BED file
echo "Merging results with original BED file..."
# bigWigAverageOverBed output: name, size, covered, sum, mean0, mean
# We want to add mean (column 6) as a new column to the BED file

awk -F'\t' '
BEGIN {OFS="\t"}
NR==FNR {
    # Read the bigWigAverageOverBed output
    # Output format: name, size, covered, sum, mean0, mean
    # Skip header line if it exists (check if first field is "name")
    if (NR == 1 && $1 == "name") {
        next
    }
    # Store mean_score by name (column 6 is mean - average over covered bases)
    scores[$1] = $6
    next
}
{
    # Process BED file
    # BED format: chr, start, end, name, score, strand, diff
    name = $4
    # Only include entries that have scores (handles duplicates)
    if (name in scores) {
        score = scores[name]
        # Print all BED columns plus mean_phastCons
        print $1, $2, $3, $4, $5, $6, $7, score
    }
}' "$OUTPUT_TSV" "$BED_FILE" > "$OUTPUT_MERGED"

echo "Merged results saved to: $OUTPUT_MERGED"
echo ""

# Summary
TOTAL_REGIONS=$(tail -n +1 "$OUTPUT_TSV" | wc -l)
if [ "$TOTAL_REGIONS" -gt 0 ] && head -1 "$OUTPUT_TSV" | grep -q "name"; then
    TOTAL_REGIONS=$((TOTAL_REGIONS - 1))  # Subtract header
fi

echo "=========================================="
echo "Summary"
echo "=========================================="
echo "  Total regions processed: $TOTAL_REGIONS"
echo "  Output columns: $(head -1 "$OUTPUT_MERGED" | awk -F'\t' '{print NF}')"
echo "  Output file: $OUTPUT_MERGED"
echo ""
echo "Done!"
