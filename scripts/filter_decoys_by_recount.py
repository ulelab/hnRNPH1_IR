#!/usr/bin/env python3
"""
Filter the predicted decoy splice-site set down to loci with NO support from
Recount/Snaptron annotated junctions.

This is the final filtering step of the decoy pipeline, applied after
`decoy_exon_overlaps.Rmd` produces `Decoys_proteincoding_splicescores.bed`
(15,172 loci). It is the step missing from `figures/decoy_splice_site_flowchart.pdf`.

Rationale
---------
A predicted cryptic splice site that already appears as an *annotated* junction
in Recount (TCGA, SRA, GTEx) is a known, used splice site - not a decoy. Removing
those leaves a set of predicted sites with no evidence of being functional splice
junctions in any of the three databases.

Inputs
------
This script consumes the TSVs written by `query_junctions.py`, which does the
actual SQLite join against the Recount `junctions.sqlite` databases. Those TSVs
already have `WHERE annotated = 1` applied, so every row in them is an annotated
junction match. This script does not re-run that query - see the README for how
to regenerate the TSVs from scratch.

Matching
--------
A locus is keyed as `<gene>_<bed_start>`, matching the `locus_id` column that
`query_junctions.py` emits (after the README's `sed` step). A locus is dropped if
it appears in ANY of the supplied junction TSVs.

Usage
-----
  python3 scripts/filter_decoys_by_recount.py \
      --loci-bed data/Decoys/Decoys_proteincoding_splicescores.bed \
      --junctions results/tcgajunctions.tsv results/srajunctions.tsv results/gtexjunctions.tsv \
      --out-bed data/Decoys/Decoys_proteincoding_recount_filtered.bed \
      --out-removed results/decoys_removed_by_recount.tsv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def load_bed(path):
    """Read the 7-column decoy BED. Returns (rows, key_to_row)."""
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                sys.exit(
                    f"ERROR: {path}:{lineno} has {len(parts)} columns, expected >= 7"
                )
            rows.append(parts)

    keyed = {}
    collisions = defaultdict(int)
    for parts in rows:
        key = f"{parts[3]}_{parts[1]}"  # gene_start, matches query_junctions locus_id
        if key in keyed:
            collisions[key] += 1
        keyed[key] = parts
    if collisions:
        # Would make the filter ambiguous - refuse rather than silently drop loci.
        sys.exit(
            f"ERROR: {len(collisions)} duplicate gene_start keys in {path}; "
            "cannot key the filter unambiguously."
        )
    return rows, keyed


def load_junction_hits(paths, keyed, strict=True):
    """Return {locus_id: set(source labels)} for every annotated junction match."""
    hits = defaultdict(set)
    for path in paths:
        label = os.path.basename(path).replace("junctions.tsv", "").replace(
            ".tsv", ""
        ).upper() or os.path.basename(path)
        n_rows = 0
        unknown = 0
        mismatched = 0
        not_annotated = 0
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                n_rows += 1
                locus_id = row.get("locus_id")
                if locus_id is None:
                    sys.exit(f"ERROR: {path} has no locus_id column")

                # Every row should already be annotated=1 (query_junctions.py
                # filters on it). Verify rather than assume.
                if row.get("annotated") not in (None, "", "1"):
                    not_annotated += 1
                    continue

                bed_row = keyed.get(locus_id)
                if bed_row is None:
                    unknown += 1
                    continue
                if (
                    bed_row[0] != row.get("chrom")
                    or bed_row[5] != row.get("locus_strand")
                    or bed_row[1] != row.get("locus_start")
                ):
                    mismatched += 1
                    continue
                hits[locus_id].add(label)

        print(
            f"  {label:6} rows={n_rows:6}  loci={len(set(r for r in hits if label in hits[r])):6}"
            f"  unknown={unknown}  coord-mismatch={mismatched}  not-annotated={not_annotated}"
        )
        if strict and (unknown or mismatched):
            sys.exit(
                f"ERROR: {path} contains {unknown} loci absent from the BED and "
                f"{mismatched} with mismatched coordinates. The junction TSVs were "
                "likely generated from a different BED - regenerate them before "
                "filtering, or pass --no-strict to ignore."
            )
    return hits


def main():
    ap = argparse.ArgumentParser(
        description="Filter decoy loci to those with no annotated Recount junction support."
    )
    ap.add_argument("--loci-bed", required=True, help="7-column decoy BED")
    ap.add_argument(
        "--junctions",
        required=True,
        nargs="+",
        help="query_junctions.py output TSVs (TCGA/SRA/GTEx)",
    )
    ap.add_argument("--out-bed", required=True, help="Filtered BED (unique decoys)")
    ap.add_argument("--out-removed", help="Optional TSV of removed loci and their sources")
    ap.add_argument(
        "--no-strict",
        dest="strict",
        action="store_false",
        help="Do not abort when junction TSVs disagree with the BED",
    )
    args = ap.parse_args()

    print(f"Loading loci from: {args.loci_bed}")
    rows, keyed = load_bed(args.loci_bed)
    print(f"  {len(rows)} loci, {len(keyed)} unique gene_start keys")

    print("Loading Recount junction matches (annotated = 1):")
    hits = load_junction_hits(args.junctions, keyed, strict=args.strict)
    print(f"  {len(hits)} loci match an annotated junction in at least one database")

    kept = [p for p in rows if f"{p[3]}_{p[1]}" not in hits]
    removed = [p for p in rows if f"{p[3]}_{p[1]}" in hits]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_bed)), exist_ok=True)
    with open(args.out_bed, "w") as fh:
        for p in kept:
            fh.write("\t".join(p) + "\n")

    if args.out_removed:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_removed)), exist_ok=True)
        with open(args.out_removed, "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(
                ["locus_id", "chrom", "start", "end", "gene", "splice_score",
                 "strand", "flag", "recount_sources"]
            )
            for p in removed:
                key = f"{p[3]}_{p[1]}"
                w.writerow(
                    [key, p[0], p[1], p[2], p[3], p[4], p[5], p[6],
                     ",".join(sorted(hits[key]))]
                )

    print()
    print(f"  input loci   : {len(rows)}")
    print(f"  removed      : {len(removed)}  (annotated junction in Recount)")
    print(f"  KEPT (unique): {len(kept)}")
    print(f"  written to   : {args.out_bed}")
    if args.out_removed:
        print(f"  removed list : {args.out_removed}")


if __name__ == "__main__":
    main()
