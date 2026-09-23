#!/usr/bin/env bash
# GC content for every site in the feature table: the harbouring intron, and
# the 49 nt window around the site that SpliceAI_Inference.py re-scores.
#
# Usage:
#   bash scripts/extract_gc_content.sh <genome.fa> <intron.bed> <site.bed> > results/gc_by_decoy.tsv
# Output (stdout, header):
#   decoyID  intron_gc  site_gc_49nt
# GC fractions are strand-independent, so no -s is needed.

set -euo pipefail

FASTA="${1:?Usage: extract_gc_content.sh <genome.fa> <intron.bed> <site.bed>}"
INTRON_BED="${2:?Usage: extract_gc_content.sh <genome.fa> <intron.bed> <site.bed>}"
SITE_BED="${3:?Usage: extract_gc_content.sh <genome.fa> <intron.bed> <site.bed>}"

for f in "$FASTA" "$FASTA.fai" "$INTRON_BED" "$SITE_BED"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
# Chromosome sizes from the FASTA index, so contig names match the FASTA.
cut -f1,2 "$FASTA.fai" > "$TMP/genome.sizes"

# bedtools nuc appends pct_at then pct_gc after the BED columns, so for BED6
# pct_gc is column 8. Column 4 is decoyID.
bedtools nuc -fi "$FASTA" -bed "$INTRON_BED" | awk 'NR > 1 {print $4 "\t" $8}' | sort -k1,1 > "$TMP/intron.tsv"
bedtools slop -i "$SITE_BED" -g "$TMP/genome.sizes" -b 24 \
  | bedtools nuc -fi "$FASTA" -bed - | awk 'NR > 1 {print $4 "\t" $8}' | sort -k1,1 > "$TMP/site.tsv"

N_IN=$(wc -l < "$INTRON_BED"); N_INTRON=$(wc -l < "$TMP/intron.tsv"); N_SITE=$(wc -l < "$TMP/site.tsv")
if [ "$N_INTRON" -ne "$N_IN" ] || [ "$N_SITE" -ne "$N_IN" ]; then
  echo "ERROR: $N_IN sites in, but $N_INTRON intron and $N_SITE site GC rows out (contig mismatch?)" >&2
  exit 1
fi

printf "decoyID\tintron_gc\tsite_gc_49nt\n"
join -t $'\t' "$TMP/intron.tsv" "$TMP/site.tsv"
