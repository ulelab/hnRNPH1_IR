#!/bin/bash
# Convert a merged crosslink BED into stranded bedGraph tracks for IGV,
# and (if bedGraphToBigWig is available) matching bigWigs.
#
# Stranded, not combined: CLIP crosslinks are strand-specific and this data is
# close to evenly split (PRPF8: 8.18M + / 7.88M -), so a single merged track
# would overlay antisense signal on sense signal and misrepresent both.
# Minus-strand values are written NEGATIVE, the usual convention so IGV draws
# them below the axis and the two strands are legible in one panel.
#
# bigWig is strongly preferred for viewing: the bedGraphs here are ~8M intervals
# and several hundred MB each, which IGV loads slowly or refuses. The bigWigs
# carry identical values and are indexed, so use them unless you specifically
# need the text form.
#
# Usage:
#   scripts/xl_to_bedgraph.sh data/CLIP/PRPF8_merged.xl.bed data/CLIP/bedgraph/PRPF8
#   scripts/xl_to_bedgraph.sh data/CLIP/SmB_merged.xl.bed   data/CLIP/bedgraph/SmB

set -euo pipefail

XL="${1:?usage: xl_to_bedgraph.sh <merged.xl.bed> <out_prefix>}"
PREFIX="${2:?usage: xl_to_bedgraph.sh <merged.xl.bed> <out_prefix>}"
FAI="${3:-reference/GRCh38.primary_assembly.genome.fa.fai}"

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"
[ -f "$XL" ] || { echo "ERROR: missing $XL" >&2; exit 1; }
mkdir -p "$(dirname "$PREFIX")"

NAME=$(basename "$PREFIX")
echo "$XL -> ${PREFIX}.{plus,minus}.bedgraph"

for S in plus minus; do
  if [ "$S" = plus ]; then SYM='+'; SIGN=1; else SYM='-'; SIGN=-1; fi
  OUT="${PREFIX}.${S}.bedgraph"
  {
    echo "track type=bedGraph name=\"${NAME}_${S}\" description=\"${NAME} crosslinks (${SYM})\" visibility=full color=$([ "$S" = plus ] && echo '0,0,180' || echo '180,0,0')"
    awk -v s="$SYM" -v sign="$SIGN" 'BEGIN{FS=OFS="\t"} $6==s {print $1,$2,$3,sign*$5}' "$XL" \
      | sort -k1,1 -k2,2n -S1G
  } > "$OUT"
  N=$(( $(wc -l < "$OUT") - 1 ))
  printf "  %-6s %9d intervals  %s\n" "$S" "$N" "$OUT"
done

# bigWig conversion (optional but recommended for IGV performance)
BGTBW=""
for c in bedGraphToBigWig "$HOME/miniconda3/envs/ucsc-bw/bin/bedGraphToBigWig"; do
  command -v "$c" >/dev/null 2>&1 && { BGTBW="$c"; break; }
  [ -x "$c" ] && { BGTBW="$c"; break; }
done

if [ -z "$BGTBW" ]; then
  echo "  (bedGraphToBigWig not found - skipping bigWig conversion)"
  exit 0
fi
if [ ! -f "$FAI" ]; then
  echo "  (no $FAI - skipping bigWig conversion)"
  exit 0
fi

cut -f1,2 "$FAI" > "${PREFIX}.chrom.sizes"
for S in plus minus; do
  # bedGraphToBigWig rejects the track line, so feed it a headerless copy
  tail -n +2 "${PREFIX}.${S}.bedgraph" > "${PREFIX}.${S}.nohdr.bedgraph"
  if "$BGTBW" "${PREFIX}.${S}.nohdr.bedgraph" "${PREFIX}.chrom.sizes" "${PREFIX}.${S}.bw" 2>/dev/null; then
    printf "  %-6s %9s  %s\n" "$S" "$(du -h "${PREFIX}.${S}.bw" | cut -f1)" "${PREFIX}.${S}.bw"
  else
    echo "  bigWig conversion failed for $S strand" >&2
  fi
  rm -f "${PREFIX}.${S}.nohdr.bedgraph"
done
rm -f "${PREFIX}.chrom.sizes"
