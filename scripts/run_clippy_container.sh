#!/bin/bash
# Call CLIP peaks with Clippy, using the pinned biocontainer.
#
# Use the container, not a local Clippy checkout. The checkout at
# ~/advbfx/software/clippy fails on this machine in every available environment:
# system python hits a NumPy 2.x / numexpr ABI error, the `clippy` conda env has
# no numpy, and the mamba `clippy-dev` env writes a Summits file but then dies in
# the broad-peak merge ("bedtools merge ... -c 11,6 ... only has fields 1 - 0")
# at every threshold including -mg 1 -mb 1. The 1.5.0 container runs the same
# input cleanly. Image is the one pinned by the peak-benchmarking pipeline
# (local-nf-modules/clippy/main.nf).
#
# Parameters match the existing peak files in data/Decoys/:
#   rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5
#
# Usage:
#   scripts/run_clippy_container.sh <crosslinks.bed> <output_prefix> [gtf] [fai]

set -euo pipefail

IMAGE="quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0"
XL="${1:?usage: run_clippy_container.sh <crosslinks.bed> <out_prefix> [gtf] [fai]}"
PREFIX="${2:?usage: run_clippy_container.sh <crosslinks.bed> <out_prefix> [gtf] [fai]}"
GTF="${3:-reference/gencode.v49.annotation.gtf}"
FAI="${4:-reference/GRCh38.primary_assembly.genome.fa.fai}"

WINDOW=${WINDOW:-10}
ADJUST=${ADJUST:-1.0}
MIN_GENE=${MIN_GENE:-5}
THREADS=${THREADS:-4}

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"

for f in "$XL" "$GTF" "$FAI"; do
  [ -f "$f" ] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done

echo "Clippy $IMAGE"
echo "  crosslinks : $XL ($(wc -l < "$XL") positions)"
echo "  annotation : $GTF"
echo "  genome     : $FAI"
echo "  params     : -n $WINDOW -x $ADJUST -mg $MIN_GENE -t $THREADS"
echo

mkdir -p "$(dirname "$PREFIX")"

docker run --rm --user "$(id -u):$(id -g)" \
  -v "$REPO":/data -w /data "$IMAGE" \
  clippy -i "$XL" -o "$PREFIX" -a "$GTF" -g "$FAI" \
         -n "$WINDOW" -x "$ADJUST" -mg "$MIN_GENE" -t "$THREADS"

echo
for f in "${PREFIX}"*_Peaks.bed "${PREFIX}"*_Summits.bed; do
  [ -f "$f" ] && printf "  %10d  %s\n" "$(wc -l < "$f")" "$f"
done
