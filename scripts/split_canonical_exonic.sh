#!/bin/bash
# STEP 2 of the reordered dataset generation.
#
# Split the protein-coding SpliceAI set into:
#   deep intronic - sites with NO overlap of a canonical splice-site window
#   exonic        - sites that DO overlap one (canonical + alternative 5'SS)
#
# Wide_canonical_splice_sites.bed is BED6 and every interval is exactly 100 bp
# (+/-50 nt around a canonical splice site), so "exonic" here means "within
# 50 nt of an annotated splice site", not "anywhere in an exon".
#
# Three outputs:
#   <out>/deep_intronic_splice_sites.bed   bedtools intersect -v -s   (BED6)
#   <out>/proteincoding_splice_sites.bed   bedtools intersect -c -s   (BED6 + count)
#   <out>/exonic_splice_sites.bed          awk '$NF > 0' on the above (BED6 + count)
#
# The -c output is kept because it annotates the FULL protein-coding set with
# its canonical-overlap count, which is what the exonic subset is derived from.
#
# Usage:
#   scripts/split_canonical_exonic.sh \
#       [-a data/Decoys/spliceai_05min_proteincoding.bed] \
#       [-b data/Decoys/Wide_canonical_splice_sites.bed] \
#       [-o results]

set -euo pipefail

AFILE="data/Decoys/spliceai_05min_proteincoding.bed"
BFILE="data/Decoys/Wide_canonical_splice_sites.bed"
OUTDIR="results"

while getopts "a:b:o:" opt; do
  case $opt in
    a) AFILE="$OPTARG" ;;
    b) BFILE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    *) echo "usage: $0 [-a proteincoding.bed] [-b wide_canonical.bed] [-o outdir]" >&2; exit 1 ;;
  esac
done

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"

for f in "$AFILE" "$BFILE"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

INTRONIC="$OUTDIR/deep_intronic_splice_sites.bed"
ANNOTATED="$OUTDIR/proteincoding_splice_sites.bed"
EXONIC="$OUTDIR/exonic_splice_sites.bed"

N_IN=$(wc -l < "$AFILE")
echo "protein-coding SpliceAI sites : $N_IN   ($AFILE)"
echo "canonical splice-site windows : $(wc -l < "$BFILE")   ($BFILE)"
echo

# Deep intronic: no canonical-window overlap, same strand.
bedtools intersect -a "$AFILE" -b "$BFILE" -v -s > "$INTRONIC"

# Full set annotated with the number of canonical windows each site overlaps.
bedtools intersect -a "$AFILE" -b "$BFILE" -c -s > "$ANNOTATED"

# Exonic borders: the annotated rows with at least one overlap.
awk -F'\t' '$NF > 0' "$ANNOTATED" > "$EXONIC"

N_INTRONIC=$(wc -l < "$INTRONIC")
N_ANNOT=$(wc -l < "$ANNOTATED")
N_EXONIC=$(wc -l < "$EXONIC")

printf "  %-34s %8d  %s\n" "deep intronic (-v)"       "$N_INTRONIC" "$INTRONIC"
printf "  %-34s %8d  %s\n" "exonic borders (count>0)" "$N_EXONIC"   "$EXONIC"
printf "  %-34s %8d  %s\n" "annotated full set (-c)"  "$N_ANNOT"    "$ANNOTATED"
echo

# The two branches must partition the input exactly. A mismatch means either a
# strand mismatch between the inputs or a malformed BED - fail loudly rather
# than letting a silently truncated dataset flow downstream.
if [ "$N_ANNOT" -ne "$N_IN" ]; then
  echo "ERROR: -c output has $N_ANNOT rows but input had $N_IN" >&2
  exit 1
fi
if [ $((N_INTRONIC + N_EXONIC)) -ne "$N_IN" ]; then
  echo "ERROR: intronic ($N_INTRONIC) + exonic ($N_EXONIC) != input ($N_IN)" >&2
  exit 1
fi
echo "  partition OK: $N_INTRONIC + $N_EXONIC = $N_IN"
echo
echo "Next: knit scripts/decoy_exon_overlaps.Rmd again - Part 2 removes exon overlaps from $INTRONIC"
