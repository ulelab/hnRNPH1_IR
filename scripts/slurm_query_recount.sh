#!/bin/bash
# Query the three hg38 Snaptron/Recount junction compilations (TCGA, SRA, GTEx) for a set of site BEDs,
# one SLURM array task per compilation.
#
# Run on the login node (not with sbatch):
#   bash slurm_query_recount.sh                         # final decoys and supported cryptics
#   BEDS="altd_reference=$PWD/altd_donors_reference.bed altd_alternative=$PWD/altd_donors_alternative.bed" \
#     OUT_TAG=recount_altd MIN_READS=5 bash slurm_query_recount.sh   # all Vast-DB alternative 5'SS donors
#
# Needs, in $WORK: query_recount_junctions.py and the BEDs named in $BEDS.
# The junctions.sqlite files are looked for in $DBDIR/<compilation>/junctions.sqlite and
# downloaded from https://snaptron.cs.jhu.edu/data/<compilation>/junctions.sqlite when
# missing (tcgav2 is ~27 GB; set DOWNLOAD=0 to fail instead of downloading).
#
# Outputs: $WORK/<OUT_TAG>_<compilation>_junctions.tsv and _summary.tsv, to copy back to results/ in the
# repository. Logs go to $WORK/logs/<OUT_TAG>/.
#
# Override any setting from the environment, e.g. COMPILATIONS="tcgav2 gtexv2" FLANK=3 bash slurm_query_recount.sh

set -euo pipefail

WORK="${WORK:-/camp/lab/ulej/home/users/jonesm6/recount}"
DBDIR="${DBDIR:-$WORK/snaptron}"
SCRIPT="${SCRIPT:-$WORK/query_recount_junctions.py}"
BEDS="${BEDS:-decoy=$WORK/decoys_final.bed cryptic_supported=$WORK/cryptics_supported_final.bed}"
OUT_TAG="${OUT_TAG:-recount}"
COMPILATIONS="${COMPILATIONS:-tcgav2 srav3h gtexv2}"
FLANK="${FLANK:-5}"
MIN_READS="${MIN_READS:-1}"
DOWNLOAD="${DOWNLOAD:-1}"
MEM="${MEM:-16G}"
TIME="${TIME:-06:00:00}"
PARTITION="${PARTITION:-}"

[ -f "$SCRIPT" ] || { echo "ERROR: missing $SCRIPT" >&2; exit 1; }
BED_ARGS=""
for kv in $BEDS; do
  f="${kv#*=}"
  [ "$f" != "$kv" ] || { echo "ERROR: BEDS entries must be LABEL=PATH, got $kv" >&2; exit 1; }
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
  BED_ARGS="$BED_ARGS --bed $kv"
done
LOGS="$WORK/logs/$OUT_TAG"
mkdir -p "$DBDIR" "$LOGS"

read -r -a COMPS <<< "$COMPILATIONS"
printf "%s\n" "${COMPS[@]}" > "$LOGS/compilations.txt"

PART_DIRECTIVE=""
[ -n "$PARTITION" ] && PART_DIRECTIVE="#SBATCH --partition=$PARTITION"

cat > "$LOGS/array.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=$OUT_TAG
#SBATCH --cpus-per-task=2
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --array=0-$(( ${#COMPS[@]} - 1 ))
#SBATCH --output=$LOGS/%a.log
$PART_DIRECTIVE

set -euo pipefail
COMP=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" "$LOGS/compilations.txt")
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
python3 "$SCRIPT" --db "\$DB" --compilation "\$COMP" --flank-bp $FLANK --min-reads $MIN_READS \\
  $BED_ARGS \\
  --out-prefix "$WORK/${OUT_TAG}_\$COMP"
EOF

JOB=$(sbatch --parsable "$LOGS/array.sbatch")
echo "submitted array $JOB for: ${COMPS[*]}"
echo "logs: $LOGS/<task>.log   outputs: $WORK/${OUT_TAG}_<compilation>_{junctions,summary}.tsv"
