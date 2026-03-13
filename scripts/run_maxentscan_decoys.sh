#!/usr/bin/env bash
set -euo pipefail

# Config
INPUT_BED="../data/DecoysControlsUniqueIDsforRecount_filtered.bed"
GENOME_FA="/home/mikej10/advbfx/reference/genomes/GRCh38.primary_assembly.genome.fa"
MAXENT_DIR="/home/mikej10/miniconda3/share/maxentscan-0_2004.04.21-4"
MAXENT_5="${MAXENT_DIR}/maxentscan_score5.pl"

# Output
OUT_BED="${INPUT_BED%.bed}_maxent.bed"

# Temp files
TMP_BED="$(mktemp)"
TMP_FA="$(mktemp)"
TMP_MAXENT="$(mktemp)"

echo "Using BED:        ${INPUT_BED}"
echo "Genome FASTA:     ${GENOME_FA}"
echo "MaxEntScan 5'SS:  ${MAXENT_5}"
echo "Output BED:       ${OUT_BED}"
echo

# 1) Build 9-mer 5'SS window BED (strand-aware; keep original BED untouched)
#    + strand: start-2, end+6
#    - strand: start-6, end+2
awk 'BEGIN{OFS="\t"} {
  if ($6 == "+")      {s=$2-2; e=$3+6}
  else if ($6 == "-") {s=$2-6; e=$3+2}
  else                {s=$2;   e=$3}
  if (s < 0) s = 0;
  $2 = s; $3 = e;
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
echo "Columns: original 6 + column 7 = MaxEnt_5SS_score"

# Cleanup
rm -f "${TMP_BED}" "${TMP_FA}" "${TMP_MAXENT}"