#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
WORKSPACE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Usage:
#   bash run_maxentscan_canonical.sh <input_bed> [output_bed] [genome_fasta]
INPUT_BED="${1:?Usage: bash run_maxentscan_canonical.sh <input_bed> [output_bed] [genome_fasta]}"
OUT_BED="${2:-${INPUT_BED%.bed}_maxent_canonical.bed}"
GENOME_FA="${3:-$WORKSPACE_DIR/reference/genomes/Gencode49/GRCh38.primary_assembly.genome.fa}"

# Resolve MaxEntScan perl script without hardcoded absolute paths.
if command -v maxentscan_score5.pl >/dev/null 2>&1; then
  MAXENT_5="$(command -v maxentscan_score5.pl)"
elif [ -n "${CONDA_PREFIX:-}" ] && [ -f "${CONDA_PREFIX}/share/maxentscan-0_2004.04.21-4/maxentscan_score5.pl" ]; then
  MAXENT_5="${CONDA_PREFIX}/share/maxentscan-0_2004.04.21-4/maxentscan_score5.pl"
elif [ -n "${MAMBA_ROOT_PREFIX:-}" ] && [ -f "${MAMBA_ROOT_PREFIX}/share/maxentscan-0_2004.04.21-4/maxentscan_score5.pl" ]; then
  MAXENT_5="${MAMBA_ROOT_PREFIX}/share/maxentscan-0_2004.04.21-4/maxentscan_score5.pl"
else
  echo "Error: maxentscan_score5.pl not found in PATH or active conda/mamba prefix." >&2
  exit 1
fi

if [ ! -f "${INPUT_BED}" ]; then
  echo "Error: input BED not found: ${INPUT_BED}" >&2
  exit 1
fi

if [ ! -f "${GENOME_FA}" ]; then
  echo "Error: genome FASTA not found: ${GENOME_FA}" >&2
  exit 1
fi

# Temp files
TMP_BED="$(mktemp)"
TMP_FA="$(mktemp)"
TMP_MAXENT="$(mktemp)"

echo "Using BED:        ${INPUT_BED}"
echo "Genome FASTA:     ${GENOME_FA}"
echo "MaxEntScan 5'SS:  ${MAXENT_5}"
echo "Output BED:       ${OUT_BED}"
echo

# 1) Build 9-mer 5'SS window BED (3 exonic + 6 intronic bases, strand-aware).
#    The input interval is the FIRST INTRONIC base of the canonical donor
#    (1-based intron_start on +, intron_end on -), so:
#    + strand: 0-based [D-4, D+5)  with D = 1-based first intronic base = $2+1
#    - strand: 0-based [D-6, D+3)  with D = $3
awk 'BEGIN{OFS="\t"} {
  if ($6=="+")      {ss=$2; s=ss-3; e=ss+6}
  else if ($6=="-"){ss=$3; s=ss-6; e=ss+3}
  else             {s=$2;   e=$3}
  if (s<0) s=0;
  $2=s; $3=e;
  print
}' "${INPUT_BED}" > "${TMP_BED}"

# 2) Get 9-mer sequences with bedtools (strand-aware)
bedtools getfasta \
  -fi "${GENOME_FA}" \
  -bed "${TMP_BED}" \
  -s -name \
  -fo "${TMP_FA}"

# 3) Run MaxEntScan 5'SS on the 9-mers
perl "${MAXENT_5}" "${TMP_FA}" > "${TMP_MAXENT}"

# 4) Append MaxEnt score as an extra column to the original BED
#    MaxEnt output: <sequence>\t<score>
#    We map by line order: line N of BED gets line N score.
awk 'NR==FNR {scores[NR]=$2; next} {print $0 "\t" scores[FNR]}' \
  "${TMP_MAXENT}" "${INPUT_BED}" > "${OUT_BED}"

echo "Wrote BED with MaxEnt scores to: ${OUT_BED}"
echo "Columns: original 7 + column 8 = MaxEnt_canonical_score"

# Cleanup
rm -f "${TMP_BED}" "${TMP_FA}" "${TMP_MAXENT}"