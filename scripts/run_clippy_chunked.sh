#!/bin/bash
# Call CLIP peaks with Clippy, one chromosome at a time, then concatenate.
#
# Why chunked: a genome-wide run on the merged crosslinks is OOM-killed on this
# machine (exit 137; kernel oom-killer at ~5.6 GB RSS with ~13 GB free of 23 GB).
# The peak-benchmarking pipeline labels Clippy "high_mem" and runs it on HPC
# nodes with far more RAM.
#
# Why this is equivalent, not an approximation: Clippy calls peaks per gene, and
# no gene spans two chromosomes, so partitioning the crosslinks and the
# annotation by chromosome and concatenating the per-chromosome peaks gives the
# same peak set as one genome-wide invocation. (This would NOT hold for the
# intergenic-peak mode, which is off by default: intergenic_peak_threshold=0.)
#
# Usage:
#   scripts/run_clippy_chunked.sh <crosslinks.bed> <out_prefix> [threads]

set -euo pipefail

IMAGE="quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0"
XL="${1:?usage: run_clippy_chunked.sh <crosslinks.bed> <out_prefix> [threads]}"
PREFIX="${2:?usage: run_clippy_chunked.sh <crosslinks.bed> <out_prefix> [threads]}"
THREADS="${3:-2}"

WINDOW=${WINDOW:-10}
ADJUST=${ADJUST:-1.0}
MIN_GENE=${MIN_GENE:-5}
WIDTH=${WIDTH:-0.4}
MIN_HEIGHT=${MIN_HEIGHT:-1.0}
# Clippy stringifies these as Python floats when building its output filename
# ("3" -> "3.0"), so normalise here or SUFFIX will not match the file Clippy
# actually writes and every chunk "fails" with rc=0.
ADJUST=$(awk -v v="$ADJUST" 'BEGIN{printf "%.1f", v}')
MIN_HEIGHT=$(awk -v v="$MIN_HEIGHT" 'BEGIN{printf "%.1f", v}')
# NB: Clippy does NOT encode --width in its output filename (see parse_arguments
# in clip/__init__.py), so two runs differing only in width collide silently.
# Put the width in your output prefix if you vary it.
SUFFIX="rollmean${WINDOW}_minHeightAdjust${MIN_HEIGHT}_minPromAdjust${ADJUST}_minGeneCount${MIN_GENE}"

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"

GTF_FULL="reference/gencode.v49.annotation.gtf"
FAI="reference/GRCh38.primary_assembly.genome.fa.fai"
GTFDIR="reference/gtf_by_chrom"
WORK="data/CLIP/_clippy_chunks/$(basename "$PREFIX")"

[ -f "$XL" ] || { echo "ERROR: missing $XL" >&2; exit 1; }
# The plain GTF is a large derived file and is gitignored, so it may have been
# cleaned up; regenerate it from the shipped .gz rather than failing.
if [ ! -f "$GTF_FULL" ]; then
  if [ -f "${GTF_FULL}.gz" ]; then
    echo "decompressing ${GTF_FULL}.gz ..."
    zcat "${GTF_FULL}.gz" > "$GTF_FULL"
  else
    echo "ERROR: missing $GTF_FULL and ${GTF_FULL}.gz" >&2; exit 1
  fi
fi
mkdir -p "$WORK"

# Split the GTF once (shared across runs) - a 3.3 GB rescan per chromosome would
# otherwise dominate the runtime.
if [ ! -f "$GTFDIR/.complete" ]; then
  echo "splitting GTF by chromosome (one pass)..."
  mkdir -p "$GTFDIR"
  awk -v d="$GTFDIR" 'BEGIN{FS=OFS="\t"} /^#/ {next} {print > (d "/" $1 ".gtf")}' "$GTF_FULL"
  touch "$GTFDIR/.complete"
fi

CHROMS=$(cut -f1 "$XL" | uniq | awk '!seen[$0]++')
echo "Clippy (chunked) on $XL"
echo "  chromosomes: $(echo "$CHROMS" | wc -w)   params: -n $WINDOW -x $ADJUST -mx $MIN_HEIGHT -w $WIDTH -mg $MIN_GENE -t $THREADS"
echo

: > "$WORK/all_peaks.bed"
: > "$WORK/all_summits.bed"

for C in $CHROMS; do
  GTF_C="$GTFDIR/${C}.gtf"
  if [ ! -s "$GTF_C" ]; then
    echo "  $C: no annotation, skipping"
    continue
  fi
  XL_C="$WORK/${C}.xl.bed"
  awk -v c="$C" 'BEGIN{FS=OFS="\t"} $1==c' "$XL" > "$XL_C"
  N=$(wc -l < "$XL_C")
  printf "  %-6s %9d crosslinks ... " "$C" "$N"

  set +e
  docker run --rm --user "$(id -u):$(id -g)" -v "$REPO":/data -w /data "$IMAGE" \
    clippy -i "$XL_C" -o "$WORK/${C}" -a "$GTF_C" -g "$FAI" \
           -n "$WINDOW" -x "$ADJUST" -mx "$MIN_HEIGHT" -mg "$MIN_GENE" -w "$WIDTH" -t "$THREADS" \
    > "$WORK/${C}.log" 2>&1
  RC=$?
  set -e

  PK="$WORK/${C}_${SUFFIX}_Peaks.bed"
  SM="$WORK/${C}_${SUFFIX}_Summits.bed"
  if [ $RC -ne 0 ] || [ ! -f "$PK" ]; then
    echo "FAILED (rc=$RC) - see $WORK/${C}.log"
    [ $RC -eq 137 ] && echo "      (exit 137 = OOM; lower threads or split further)"
    continue
  fi
  cat "$PK" >> "$WORK/all_peaks.bed"
  [ -f "$SM" ] && cat "$SM" >> "$WORK/all_summits.bed"
  echo "$(wc -l < "$PK") peaks"
  rm -f "$XL_C" "$SM"
done

OUT_PK="${PREFIX}_${SUFFIX}_Peaks.bed"
OUT_SM="${PREFIX}_${SUFFIX}_Summits.bed"
sort -k1,1 -k2,2n "$WORK/all_peaks.bed" > "$OUT_PK"
sort -k1,1 -k2,2n "$WORK/all_summits.bed" > "$OUT_SM"

echo
printf "  %10d  %s\n" "$(wc -l < "$OUT_PK")" "$OUT_PK"
printf "  %10d  %s\n" "$(wc -l < "$OUT_SM")" "$OUT_SM"
