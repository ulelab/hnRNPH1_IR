#!/bin/bash
# Query the three hg38 Snaptron/Recount junction compilations (TCGA, SRA, GTEx) for the
# final decoy and supported-cryptic sites, one SLURM array task per compilation.
#
# Run on the login node (not with sbatch):
#   bash slurm_query_recount.sh
#
# Needs, in $WORK: query_recount_junctions.py, decoys_final.bed, cryptics_supported_final.bed.
# The junctions.sqlite files are looked for in $DBDIR/<compilation>/junctions.sqlite and
# downloaded from https://snaptron.cs.jhu.edu/data/<compilation>/junctions.sqlite when
# missing (tcgav2 is ~27 GB; set DOWNLOAD=0 to fail instead of downloading).
#
# Outputs: $WORK/recount_<compilation>_junctions.tsv and _summary.tsv, to copy back to
# results/ in the repository.
#
# Override any setting from the environment, e.g. COMPILATIONS="tcgav2 gtexv2" FLANK=3 bash slurm_query_recount.sh

set -euo pipefail

WORK="${WORK:-/camp/lab/ulej/home/users/jonesm6/recount}"
DBDIR="${DBDIR:-$WORK/snaptron}"
SCRIPT="${SCRIPT:-$WORK/query_recount_junctions.py}"
DECOYS="${DECOYS:-$WORK/decoys_final.bed}"
CRYPTICS="${CRYPTICS:-$WORK/cryptics_supported_final.bed}"
COMPILATIONS="${COMPILATIONS:-tcgav2 srav3h gtexv2}"
FLANK="${FLANK:-5}"
DOWNLOAD="${DOWNLOAD:-1}"
MEM="${MEM:-16G}"
TIME="${TIME:-06:00:00}"
PARTITION="${PARTITION:-}"

for f in "$SCRIPT" "$DECOYS" "$CRYPTICS"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
mkdir -p "$DBDIR" "$WORK/logs"

read -r -a COMPS <<< "$COMPILATIONS"
printf "%s\n" "${COMPS[@]}" > "$WORK/logs/compilations.txt"

PART_DIRECTIVE=""
[ -n "$PARTITION" ] && PART_DIRECTIVE="#SBATCH --partition=$PARTITION"

cat > "$WORK/logs/array.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=recount_junctions
#SBATCH --cpus-per-task=2
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --array=0-$(( ${#COMPS[@]} - 1 ))
#SBATCH --output=$WORK/logs/%a.log
$PART_DIRECTIVE

set -euo pipefail
COMP=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" "$WORK/logs/compilations.txt")
DB="$DBDIR/\$COMP/junctions.sqlite"
if [ ! -s "\$DB" ]; then
  if [ "$DOWNLOAD" = "1" ]; then
    mkdir -p "$DBDIR/\$COMP"
    echo "downloading https://snaptron.cs.jhu.edu/data/\$COMP/junctions.sqlite ..."
    curl -fL --retry 5 -C - -o "\$DB.part" "https://snaptron.cs.jhu.edu/data/\$COMP/junctions.sqlite"
    mv "\$DB.part" "\$DB"
  else
    echo "ERROR: \$DB not found and DOWNLOAD=0" >&2; exit 1
  fi
fi
python3 "$SCRIPT" --db "\$DB" --compilation "\$COMP" --flank-bp $FLANK \\
  --bed decoy="$DECOYS" --bed cryptic_supported="$CRYPTICS" \\
  --out-prefix "$WORK/recount_\$COMP"
EOF

JOB=$(sbatch --parsable "$WORK/logs/array.sbatch")
echo "submitted array $JOB for: ${COMPS[*]}"
echo "logs: $WORK/logs/<task>.log   outputs: $WORK/recount_<compilation>_{junctions,summary}.tsv"
