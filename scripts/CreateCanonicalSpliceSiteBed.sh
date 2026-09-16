#!/bin/bash
# Canonical 5' splice-site windows: +/-5 nt around each annotated canonical donor.
#
# Built by trimming Wide_canonical_splice_sites.bed (+/-50 nt, 100 bp windows)
# to its central 10 bp, so both files mark exactly the same donors.
#
# The previous version read an introns.bed that still carried the 400 nt exonic
# flanks added for SpliceAI, which put every window 400 nt into the flanking
# exon (+ strand -400, - strand +400) and missed every real donor.
#
# Usage:
#   scripts/CreateCanonicalSpliceSiteBed.sh \
#       [data/Decoys/Wide_canonical_splice_sites.bed] \
#       [data/Decoys/Canonical_splice_sites.bed]

set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"

IN="${1:-data/Decoys/Wide_canonical_splice_sites.bed}"
OUT="${2:-data/Decoys/Canonical_splice_sites.bed}"

[ -f "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }

# The +45/-45 trim only centres the window if every input interval is 100 bp.
BAD=$(awk -F'\t' '$3 - $2 != 100' "$IN" | wc -l)
if [ "$BAD" -ne 0 ]; then
  echo "ERROR: $BAD intervals in $IN are not 100 bp" >&2
  exit 1
fi

awk 'BEGIN{FS=OFS="\t"} {print $1, $2 + 45, $3 - 45, $4, $5, $6}' "$IN" > "$OUT"

N_IN=$(wc -l < "$IN")
N_OUT=$(awk -F'\t' '$3 - $2 == 10' "$OUT" | wc -l)
if [ "$N_OUT" -ne "$N_IN" ]; then
  echo "ERROR: $N_OUT of $N_IN output windows are 10 bp" >&2
  exit 1
fi
echo "wrote $OUT: $N_OUT windows of 10 bp (from $IN)"
