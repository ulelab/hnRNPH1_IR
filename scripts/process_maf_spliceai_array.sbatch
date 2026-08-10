#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --array=1-109
#SBATCH --job-name=msa_spliceai
#SBATCH --output=msa_spliceai_%A_%a.out
#SBATCH --error=msa_spliceai_%A_%a.err

# Activate conda environment
CONDA_BASE="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/mambaforge"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate spliceai-env

# Verify conda environment is activated
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "ERROR: Failed to activate conda environment spliceai-env"
    exit 1
fi
echo "Activated conda environment: $CONDA_DEFAULT_ENV"

# Paths - adjust these for your HPC
MSA_VIEW_BIN="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/kent/src/msa_view"
DECOY_BED="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/consGenesDecoys/DecoysInConservedGenes50thpercentile.bed"
MAF_OUTPUT_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/consGenesDecoys/maf_output"
OUTPUT_DIR="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/consGenesDecoys/splice_out"
SPECIES_MAP="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/tree_to_clade_mapping.tsv"
CHRPARTS_FILE="/scratch/prj/ppn_rnp_networks/users/mike.jones/data/maf/consGenesDecoys/chrparts.txt"
PYTHON_SCRIPT="/scratch/prj/ppn_rnp_networks/users/mike.jones/software/process_maf_spliceai.py"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get the subdirectory name for this array task
BED_PART=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRPARTS_FILE")
SUBDIR_NAME="${BED_PART%.bed}"

if [ -z "$SUBDIR_NAME" ]; then
    echo "ERROR: No subdirectory found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

MAF_SUBDIR="${MAF_OUTPUT_DIR}/${SUBDIR_NAME}"

# Check if MAF subdirectory exists
if [ ! -d "$MAF_SUBDIR" ]; then
    echo "ERROR: MAF subdirectory not found: $MAF_SUBDIR"
    exit 1
fi

# Count MAF files
NUM_MAF_FILES=$(find "$MAF_SUBDIR" -name "*.maf" | wc -l)

echo "=========================================="
echo "Task ${SLURM_ARRAY_TASK_ID}: Processing ${SUBDIR_NAME}"
echo "MAF subdirectory: $MAF_SUBDIR"
echo "Number of MAF files: $NUM_MAF_FILES"
echo "Output directory: $OUTPUT_DIR"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Run Python script to process all MAF files in this subdirectory
# Use python from conda environment
python "$PYTHON_SCRIPT" \
    --maf_dir "$MAF_SUBDIR" \
    --decoy_bed "$DECOY_BED" \
    --species_map "$SPECIES_MAP" \
    --msa_view "$MSA_VIEW_BIN" \
    --output_tsv "${OUTPUT_DIR}/${SUBDIR_NAME}_spliceai_scores.tsv"

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Task ${SLURM_ARRAY_TASK_ID} completed"
    echo "End time: $(date)"
    
    if [ -f "${OUTPUT_DIR}/${SUBDIR_NAME}_spliceai_scores.tsv" ]; then
        OUTPUT_SIZE=$(du -h "${OUTPUT_DIR}/${SUBDIR_NAME}_spliceai_scores.tsv" | cut -f1)
        OUTPUT_LINES=$(wc -l < "${OUTPUT_DIR}/${SUBDIR_NAME}_spliceai_scores.tsv")
        echo "Output TSV file size: $OUTPUT_SIZE"
        echo "Output TSV lines: $OUTPUT_LINES"
    fi
else
    echo "ERROR: Task ${SLURM_ARRAY_TASK_ID} failed with exit code $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE

