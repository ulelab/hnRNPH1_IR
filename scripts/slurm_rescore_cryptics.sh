#!/bin/bash
# Local SpliceAI re-score of cryptics_supported.bed on SLURM, as N parallel slices.
#
# Run on the login node (not with sbatch):
#   bash slurm_rescore_cryptics.sh
#
# It splits the input into N line-based slices, submits one array task per slice
# running SpliceAI_Inference.py, and submits a merge job (afterok) that
# concatenates the slices back in input order. SpliceAI_Inference.py writes rows
# in input order, so the merged file matches an unsliced run row for row.
#
# Any setting below can be overridden from the environment, e.g.
#   N=64 TIME=01:30:00 bash slurm_rescore_cryptics.sh
# PARTITION is unset by default, so jobs go to the cluster's default partition.
# A failed or timed-out task can be resubmitted alone; --resume picks up its
# checkpoint:  sbatch --array=<task id> <WORK>/logs/array.sbatch

set -euo pipefail

WORK="${WORK:-/camp/lab/ulej/home/users/jonesm6/cryptics}"
IN="${IN:-$WORK/cryptics_supported.bed}"
SCRIPT="${SCRIPT:-$WORK/SpliceAI_Inference.py}"
FASTA="${FASTA:-/camp/lab/ulej/home/shared/genomes/hg38/GRCh38.primary_assembly.genome.fa}"
OUT="${OUT:-$WORK/cryptics_supported_splicescores.bed}"
N="${N:-32}"
CONDA_BASE="${CONDA_BASE:-/camp/home/jonesm6/home/users/jonesm6/software/miniconda3}"
CONDA_ENV="${CONDA_ENV:-spliceai-env}"
PARTITION="${PARTITION:-}"
CPUS="${CPUS:-16}"
MEM="${MEM:-16G}"
TIME="${TIME:-03:00:00}"

for f in "$IN" "$SCRIPT" "$FASTA" "$FASTA.fai" "$CONDA_BASE/etc/profile.d/conda.sh"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

SLICES="$WORK/slices"
SCORED="$WORK/scored"
LOGS="$WORK/logs"
mkdir -p "$SLICES" "$SCORED" "$LOGS"

# Chromosome sizes for bedtools slop, taken from the FASTA index so the contig
# names are guaranteed to match the FASTA getfasta reads from.
GENOME="$WORK/genome.sizes"
cut -f1,2 "$FASTA.fai" > "$GENOME"

# A contig-name mismatch (chr1 vs 1) makes bedtools drop rows without failing.
MISSING=$(cut -f1 "$IN" | sort -u | grep -vxFf <(cut -f1 "$GENOME") || true)
if [ -n "$MISSING" ]; then
  echo "ERROR: contigs in $IN not in $FASTA.fai: $(echo $MISSING | head -c 200)" >&2
  exit 1
fi

rm -f "$SLICES"/part_*.bed
split -n "l/$N" -d -a 2 --additional-suffix=.bed "$IN" "$SLICES/part_"
N_SLICES=$(ls "$SLICES"/part_*.bed | wc -l)
N_IN=$(wc -l < "$IN")
echo "input: $N_IN sites -> $N_SLICES slices of ~$((N_IN / N_SLICES)) in $SLICES"

PART_DIRECTIVE=""
PART_ARGS=()
if [ -n "$PARTITION" ]; then
  PART_DIRECTIVE="#SBATCH --partition=$PARTITION"
  PART_ARGS=(--partition="$PARTITION")
fi

cat > "$LOGS/array.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=spliceai_cryptics
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --array=0-$((N_SLICES - 1))
#SBATCH --output=$LOGS/part_%a.log
$PART_DIRECTIVE

set -euo pipefail
source "$CONDA_BASE/etc/profile.d/conda.sh" || { echo "ERROR: could not source conda.sh from $CONDA_BASE"; exit 1; }
conda activate $CONDA_ENV || { echo "ERROR: could not activate conda env $CONDA_ENV"; exit 1; }

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
export TF_NUM_INTRAOP_THREADS=\$SLURM_CPUS_PER_TASK
export TF_NUM_INTEROP_THREADS=2

PART=\$(printf "%02d" "\$SLURM_ARRAY_TASK_ID")
python "$SCRIPT" \\
  --bed "$SLICES/part_\$PART.bed" \\
  --fasta "$FASTA" \\
  --genome "$GENOME" \\
  --batch-size 64 \\
  --resume \\
  --out "$SCORED/part_\$PART.splicescores.bed"
EOF

ARRAY_ID=$(sbatch --parsable "$LOGS/array.sbatch")
echo "submitted array job $ARRAY_ID ($N_SLICES tasks)"

MERGE_ID=$(sbatch --parsable --dependency="afterok:$ARRAY_ID" \
  --job-name=spliceai_merge ${PART_ARGS[@]+"${PART_ARGS[@]}"} --cpus-per-task=1 --mem=2G \
  --time=00:10:00 --output="$LOGS/merge.log" <<EOF
#!/bin/bash
set -euo pipefail
for i in \$(seq -f "%02g" 0 $((N_SLICES - 1))); do
  cat "$SCORED/part_\$i.splicescores.bed"
done > "$OUT"
N_OUT=\$(wc -l < "$OUT")
if [ "\$N_OUT" -ne "$N_IN" ]; then
  echo "ERROR: merged \$N_OUT rows but input had $N_IN" >&2
  exit 1
fi
echo "merged $N_SLICES slices -> $OUT (\$N_OUT rows)"
echo "unscorable (-1): \$(awk -F'\t' '\$5 == -1' "$OUT" | wc -l)"
echo "local SpliceAI >= 0.1: \$(awk -F'\t' '\$5 >= 0.1' "$OUT" | wc -l)"
EOF
)
echo "submitted merge job $MERGE_ID (runs after $ARRAY_ID succeeds) -> $OUT"
echo "watch: squeue -u \$USER    progress: tail -1 $LOGS/part_*.log"
