#!/bin/bash
# Test script to extract 470-way MAF for HNRNPH1 region using stream_maf_region.py
# Region: chr5:179619358-179620941

set -e

# Configuration
CHROM="chr5"
START=179619358
END=179620941
OUTPUT_MAF="/home/mikej10/advbfx/intronretention/hnRNPH1_IR_MAF/data/HNRNPH1_179620582_470way.maf"
SCRIPT_DIR="/home/mikej10/advbfx/intronretention/hnRNPH1_IR_MAF/scripts"
STREAM_SCRIPT="${SCRIPT_DIR}/stream_maf_region.py"

# Check if script exists
if [ ! -f "$STREAM_SCRIPT" ]; then
    echo "ERROR: stream_maf_region.py not found at $STREAM_SCRIPT"
    exit 1
fi

echo "=========================================="
echo "Testing 470-way MAF extraction for HNRNPH1"
echo "=========================================="
echo "Region: ${CHROM}:${START}-${END}"
echo "Output file: ${OUTPUT_MAF}"
echo "Script: ${STREAM_SCRIPT}"
echo ""
echo "Starting extraction..."
echo ""

# Record start time
START_TIME=$(date +%s)

# Run the stream script
python3 "$STREAM_SCRIPT" "$CHROM" "$START" "$END" "$OUTPUT_MAF" 2>&1 | tee /tmp/stream_maf_test.log

# Record end time
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "Extraction completed!"
echo "=========================================="
echo "Duration: ${DURATION} seconds ($(echo "scale=2; ${DURATION}/60" | bc) minutes)"
echo ""

# Check if output file was created
if [ -f "$OUTPUT_MAF" ]; then
    FILE_SIZE=$(du -h "$OUTPUT_MAF" | cut -f1)
    LINE_COUNT=$(wc -l < "$OUTPUT_MAF")
    
    echo "Output file statistics:"
    echo "  File: ${OUTPUT_MAF}"
    echo "  Size: ${FILE_SIZE}"
    echo "  Lines: ${LINE_COUNT}"
    echo ""
    
    # Count MAF blocks (lines starting with 'a ')
    BLOCK_COUNT=$(grep -c "^a " "$OUTPUT_MAF" 2>/dev/null || echo "0")
    echo "  MAF blocks: ${BLOCK_COUNT}"
    
    # Count species/sequences (lines starting with 's ')
    SEQ_COUNT=$(grep -c "^s " "$OUTPUT_MAF" 2>/dev/null || echo "0")
    echo "  Sequence lines: ${SEQ_COUNT}"
    
    # Show first few lines
    echo ""
    echo "First 20 lines of output:"
    echo "----------------------------------------"
    head -20 "$OUTPUT_MAF"
    echo "----------------------------------------"
    
    # Check for hg38 reference
    if grep -q "^s hg38\.chr5" "$OUTPUT_MAF"; then
        echo ""
        echo "✓ Found hg38.chr5 reference sequence"
        HG38_LINES=$(grep "^s hg38\.chr5" "$OUTPUT_MAF" | head -1)
        echo "  Example: ${HG38_LINES:0:100}..."
    else
        echo ""
        echo "⚠ WARNING: No hg38.chr5 reference found in output"
    fi
    
else
    echo "ERROR: Output file was not created!"
    echo "Check the log above for errors."
    exit 1
fi

echo ""
echo "Test completed successfully!"
