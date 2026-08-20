#!/bin/bash
# Merge per-sample crosslink BEDs from Flow into one strand-aware crosslink track
# per purification target, converting Ensembl contig names to UCSC on the way.
#
# Two things this handles that a plain `cat` does not:
#
#  1. CONTIG NAMING. Flow's .genome.xl.bed files are Ensembl-style ("1", "MT").
#     Every other track in this repo - the SpliceAI inferences, the RBPnet peaks,
#     the Gencode annotation - is UCSC-style ("chr1", "chrM"). Intersecting the
#     two without converting silently returns ZERO overlaps rather than an error.
#     Non-primary scaffolds (GL*, KI*) are dropped: they carry no SpliceAI sites.
#
#  2. SCORE SUMMING. Crosslinks at the same position in different replicates must
#     have their cDNA counts summed, not merely concatenated, or the peak caller
#     sees duplicate single-count entries instead of one high-count position.
#
# Usage:
#   scripts/merge_crosslinks.sh data/CLIP/PRPF8  data/CLIP/PRPF8_merged.xl.bed
#   scripts/merge_crosslinks.sh data/CLIP/SmB    data/CLIP/SmB_merged.xl.bed

set -euo pipefail

INDIR="${1:?usage: merge_crosslinks.sh <indir> <out.bed>}"
OUT="${2:?usage: merge_crosslinks.sh <indir> <out.bed>}"

shopt -s nullglob
FILES=("$INDIR"/*.bed)
[ ${#FILES[@]} -gt 0 ] || { echo "ERROR: no .bed files in $INDIR" >&2; exit 1; }

echo "merging ${#FILES[@]} files from $INDIR"
for f in "${FILES[@]}"; do
  printf "   %10d  %s\n" "$(wc -l < "$f")" "$(basename "$f")"
done

mkdir -p "$(dirname "$OUT")"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# Ensembl -> UCSC, primary contigs only, then sort by position+strand
cat "${FILES[@]}" \
  | awk 'BEGIN{FS=OFS="\t"}
      {
        c = $1
        if (c == "MT") c = "chrM"
        else if (c ~ /^([0-9]+|X|Y)$/) c = "chr" c
        else next                      # drop GL*/KI* scaffolds
        print c, $2, $3, $4, $5, $6
      }' \
  | sort -k1,1 -k2,2n -k6,6 -S2G --parallel=4 > "$TMP/all.bed"

# sum scores across replicates at identical (chrom,start,end,strand)
awk 'BEGIN{FS=OFS="\t"}
  {
    key = $1 FS $2 FS $3 FS $6
    if (key == prev) { sum += $5 }
    else {
      if (NR > 1) { split(prev, p, FS); print p[1], p[2], p[3], ".", sum, p[4] }
      prev = key; sum = $5
    }
  }
  END { if (NR > 0) { split(prev, p, FS); print p[1], p[2], p[3], ".", sum, p[4] } }' \
  "$TMP/all.bed" > "$OUT"

IN=$(cat "${FILES[@]}" | wc -l)
KEPT=$(wc -l < "$TMP/all.bed")
UNIQ=$(wc -l < "$OUT")
echo
echo "  input crosslink rows      : $IN"
echo "  after contig conversion   : $KEPT  (dropped $((IN-KEPT)) on non-primary scaffolds)"
echo "  unique crosslink positions: $UNIQ"
echo "  total cDNA counts         : $(awk -F'\t' '{s+=$5} END{printf "%d", s}' "$OUT")"
echo "  -> $OUT"
