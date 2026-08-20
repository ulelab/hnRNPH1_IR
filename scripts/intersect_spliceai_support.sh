#!/bin/bash
# Intersect SpliceAI splice-site inferences against CLIP/prediction support tracks
# over a range of window sizes.
#
# Logic (the AND/OR union from the flowchart, generalised to N categories and a
# distance window): a SpliceAI site is retained if it falls within W nt of a
# feature in ANY of the supplied category files, on the SAME strand.
#
# Per-category counts are also emitted so the contribution of each track is
# visible, not just the union.
#
# Usage:
#   scripts/intersect_spliceai_support.sh \
#       -a data/Decoys/Splice_All.filtered.10min.bed \
#       -g ~/advbfx/reference/genomes/Gencode49/genome.sizes \
#       -o results/support_sweep \
#       -w "50,25,10" \
#       RBPNET=data/Decoys/rbpnet_...Peaks.bed \
#       PRPF8=data/Decoys/PRPF8_...Peaks.bed \
#       SMB=data/Decoys/SmB_...Peaks.bed
#
# Each trailing argument is NAME=path. Any number of categories is accepted.

set -euo pipefail

WINDOWS="50,25,10"
OUTPREFIX="results/support_sweep"
GENOME=""
AFILE=""

while getopts "a:g:o:w:" opt; do
  case $opt in
    a) AFILE="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    o) OUTPREFIX="$OPTARG" ;;
    w) WINDOWS="$OPTARG" ;;
    *) echo "unknown option" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

[ -n "$AFILE" ] || { echo "ERROR: -a (SpliceAI BED) is required" >&2; exit 1; }
[ -f "$AFILE" ] || { echo "ERROR: -a file not found: $AFILE" >&2; exit 1; }
[ $# -gt 0 ] || { echo "ERROR: supply at least one CATEGORY=path argument" >&2; exit 1; }

NAMES=(); PATHS=()
for kv in "$@"; do
  name="${kv%%=*}"; path="${kv#*=}"
  [ -f "$path" ] || { echo "ERROR: category '$name' file not found: $path" >&2; exit 1; }
  NAMES+=("$name"); PATHS+=("$path")
done

mkdir -p "$(dirname "$OUTPREFIX")"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

echo "SpliceAI input : $AFILE ($(wc -l < "$AFILE") sites)"
for i in "${!NAMES[@]}"; do
  echo "  category ${NAMES[$i]}: ${PATHS[$i]} ($(wc -l < "${PATHS[$i]}") features)"
done
echo

SUMMARY="${OUTPREFIX}_summary.tsv"
{
  printf "window_nt\tspliceai_in\tsupported\tpct"
  for n in "${NAMES[@]}"; do printf "\t%s" "$n"; done
  printf "\tany_two_plus\n"
} > "$SUMMARY"

IFS=',' read -ra WLIST <<< "$WINDOWS"
for W in "${WLIST[@]}"; do
  W=$(echo "$W" | tr -d ' ')
  echo "=== window +/- ${W} nt ==="

  # one count column per category, in A order
  PASTE_ARGS=()
  for i in "${!NAMES[@]}"; do
    out="$TMP/c${i}.txt"
    bedtools window -w "$W" -sm -c -a "$AFILE" -b "${PATHS[$i]}" | cut -f7 > "$out"
    PASTE_ARGS+=("$out")
    echo "    ${NAMES[$i]}: $(awk '$1>0' "$out" | wc -l) sites supported"
  done

  paste "$AFILE" "${PASTE_ARGS[@]}" > "$TMP/joined.bed"

  NCAT=${#NAMES[@]}
  OUTBED="${OUTPREFIX}_w${W}.bed"
  awk -v ncat="$NCAT" 'BEGIN{FS=OFS="\t"}
    {
      hits=0; multi=0
      for (i = NF-ncat+1; i <= NF; i++) { if ($i > 0) { hits++ } }
      if (hits > 0) print $0, hits
    }' "$TMP/joined.bed" > "$OUTBED"

  TOTAL=$(wc -l < "$AFILE")
  SUP=$(wc -l < "$OUTBED")
  TWO=$(awk -v ncat="$NCAT" 'BEGIN{FS="\t"}{if ($NF >= 2) n++} END{print n+0}' "$OUTBED")
  PCT=$(awk -v s="$SUP" -v t="$TOTAL" 'BEGIN{printf "%.1f", 100*s/t}')

  {
    printf "%s\t%s\t%s\t%s" "$W" "$TOTAL" "$SUP" "$PCT"
    for i in "${!NAMES[@]}"; do
      printf "\t%s" "$(awk '$1>0' "$TMP/c${i}.txt" | wc -l)"
    done
    printf "\t%s\n" "$TWO"
  } >> "$SUMMARY"

  echo "    UNION supported: $SUP / $TOTAL (${PCT}%)  -> $OUTBED"
  echo "    supported by >=2 categories: $TWO"
  echo
done

echo "summary -> $SUMMARY"
column -t "$SUMMARY"
