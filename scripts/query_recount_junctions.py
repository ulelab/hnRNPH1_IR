#!/usr/bin/env python3
"""Match final decoy / supported-cryptic sites to Recount (Snaptron) junctions.

Loads one or more 10-column site BEDs (decoys_final.bed, cryptics_supported_final.bed)
into an in-memory table attached to a read-only Snaptron junctions.sqlite, and reports
every junction whose intron boundary lies within --flank-bp of a site, annotated or not.

Two outputs:
  <out-prefix>_junctions.tsv   one row per (site, junction) pair
  <out-prefix>_summary.tsv     one row per site: junction support at the 5'SS boundary

Standard library only (Python >= 3.6), so it runs on a cluster login node as is.

Usage:
  python3 scripts/query_recount_junctions.py --db junctions.sqlite --compilation tcgav2 \\
      --bed decoy=results/decoys_final.bed --bed cryptic_supported=results/cryptics_supported_final.bed \\
      --out-prefix results/recount_tcgav2
"""
import argparse
import csv
import sqlite3
from collections import defaultdict
from pathlib import Path


def connect(db_path):
    uri = "file:{}?mode=ro".format(Path(db_path).expanduser().resolve())
    conn = sqlite3.connect(uri, uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def load_sites(conn, beds):
    conn.execute("ATTACH DATABASE ':memory:' AS locusdb")
    conn.execute(
        """
        CREATE TABLE locusdb.site (
            site_id TEXT NOT NULL,
            site_class TEXT NOT NULL,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,   -- 1-based site coordinate (BED end)
            strand TEXT NOT NULL,
            gene TEXT,
            spliceai REAL
        )
        """
    )
    conn.execute("CREATE INDEX locusdb.idx_site ON site (chrom, strand, pos)")
    rows = []
    for label, path in beds:
        with open(path) as fh:
            for parts in csv.reader(fh, delimiter="\t"):
                if len(parts) < 6:
                    continue
                start, end = int(parts[1]), int(parts[2])
                rows.append(("{}_{}".format(parts[3], start), label, parts[0], end, parts[5], parts[3], float(parts[4])))
    conn.executemany("INSERT INTO locusdb.site VALUES (?, ?, ?, ?, ?, ?, ?)", rows)
    conn.commit()
    return len(rows)


# Columns the original Recount query relied on; they exist in every Snaptron junctions.sqlite used here.
REQUIRED = ["snaptron_id", "chrom", "start", "end", "strand", "annotated", "left_annotated", "right_annotated",
            "samples_count", "coverage_sum", "coverage_avg", "coverage_median", "source_dataset_id"]
# Included only when the build has them.
OPTIONAL = ["left_motif", "right_motif", "donor", "acceptor"]


def intron_columns(conn):
    cols = [r[1] for r in conn.execute("PRAGMA table_info(intron)")]
    missing = [c for c in REQUIRED if c not in cols]
    if missing:
        raise SystemExit("intron table lacks required columns {}; has {}".format(missing, cols))
    idx = [(r[1], [c[2] for c in conn.execute("PRAGMA index_info('{}')".format(r[1]))])
           for r in conn.execute("PRAGMA index_list(intron)")]
    print("intron columns: {}".format(", ".join(cols)))
    print("intron indexes: {}".format(idx if idx else "none"))
    return [c for c in OPTIONAL if c in cols]


def build_query(optional):
    extra = "".join(", i.{}".format(c) for c in optional)
    # CROSS JOIN fixes the loop order: one sequential pass over the (very large) intron table, probing the
    # small site table through its (chrom, strand, pos) index for either boundary. The ranges are written on
    # s.pos so that index is usable whether or not the intron table has indexes of its own.
    return """
    SELECT
        s.site_id, s.site_class, s.gene, s.chrom, s.pos AS site_pos, s.strand, s.spliceai,
        i.snaptron_id, i.start AS junction_start, i.end AS junction_end, (i.end - i.start + 1) AS junction_length,
        i.annotated, i.left_annotated, i.right_annotated{extra},
        i.samples_count, i.coverage_sum, i.coverage_avg, i.coverage_median, i.source_dataset_id,
        (i.start - s.pos) AS start_offset,
        (i.end   - s.pos) AS end_offset
    FROM intron i CROSS JOIN locusdb.site s
    WHERE s.chrom = i.chrom AND s.strand = i.strand
      AND (s.pos BETWEEN i.start - :flank AND i.start + :flank
        OR s.pos BETWEEN i.end   - :flank AND i.end   + :flank)
    ORDER BY s.site_id, i.samples_count DESC
    """.format(extra=extra)


def junction_type(row, flank):
    # The site is a putative 5' splice site: on + it is the junction start, on - the junction end.
    start_match = abs(row["start_offset"]) <= flank
    end_match = abs(row["end_offset"]) <= flank
    if start_match and end_match:
        return "both"
    if start_match:
        return "5ss" if row["strand"] == "+" else "3ss"
    if end_match:
        return "3ss" if row["strand"] == "+" else "5ss"
    return "none"


def is_unannotated(row):
    # Snaptron: annotated = 1 when the exact junction is in any bundled annotation. A junction
    # with annotated = 0 can still have one annotated end (left_/right_annotated non-empty).
    return int(row["annotated"] or 0) == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="Snaptron junctions.sqlite (read-only)")
    ap.add_argument("--compilation", required=True, help="label written to the output, e.g. tcgav2")
    ap.add_argument("--bed", action="append", required=True, metavar="LABEL=PATH",
                    help="site BED with its class label; repeatable")
    ap.add_argument("--flank-bp", type=int, default=5, help="boundary tolerance in bp [5]")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    beds = []
    for item in a.bed:
        label, _, path = item.partition("=")
        if not path:
            ap.error("--bed expects LABEL=PATH, got {!r}".format(item))
        beds.append((label, path))

    conn = connect(a.db)
    try:
        optional = intron_columns(conn)
        n_sites = load_sites(conn, beds)
        q = build_query(optional)
        for r in conn.execute("EXPLAIN QUERY PLAN " + q, {"flank": a.flank_bp}):
            print("  plan:", r[-1])
        rows = conn.execute(q, {"flank": a.flank_bp}).fetchall()
    finally:
        conn.close()

    junc_path = "{}_junctions.tsv".format(a.out_prefix)
    cols = list(rows[0].keys()) if rows else []
    with open(junc_path, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["compilation"] + cols + ["junction_type"])
        for r in rows:
            w.writerow([a.compilation] + [r[c] for c in cols] + [junction_type(r, a.flank_bp)])

    # Per-site summary over junctions that use the site as a 5' splice site.
    per_site = defaultdict(list)
    for r in rows:
        if junction_type(r, a.flank_bp) in ("5ss", "both"):
            per_site[r["site_id"]].append(r)

    summ_path = "{}_summary.tsv".format(a.out_prefix)
    conn = connect(a.db)  # site list again, so unmatched sites are written too
    try:
        load_sites(conn, beds)
        sites = conn.execute("SELECT site_id, site_class, gene, chrom, pos, strand, spliceai FROM locusdb.site").fetchall()
    finally:
        conn.close()

    n_matched = n_unann = 0
    with open(summ_path, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["compilation", "site_id", "site_class", "gene", "chrom", "site_pos", "strand", "spliceai",
                    "n_junctions_5ss", "n_annotated_5ss", "n_unannotated_5ss",
                    "max_samples_5ss", "max_coverage_sum_5ss",
                    "max_samples_unannotated_5ss", "max_coverage_sum_unannotated_5ss",
                    "best_junction_offset", "best_junction_id"])
        for s in sites:
            js = per_site.get(s["site_id"], [])
            un = [j for j in js if is_unannotated(j)]
            best = max(js, key=lambda j: (j["samples_count"] or 0), default=None)
            if js:
                n_matched += 1
            if un and not [j for j in js if not is_unannotated(j)]:
                n_unann += 1
            w.writerow([
                a.compilation, s["site_id"], s["site_class"], s["gene"], s["chrom"], s["pos"], s["strand"], s["spliceai"],
                len(js), len(js) - len(un), len(un),
                max((j["samples_count"] or 0) for j in js) if js else 0,
                max((j["coverage_sum"] or 0) for j in js) if js else 0,
                max((j["samples_count"] or 0) for j in un) if un else 0,
                max((j["coverage_sum"] or 0) for j in un) if un else 0,
                (best["start_offset"] if s["strand"] == "+" else best["end_offset"]) if best else "",
                best["snaptron_id"] if best else "",
            ])

    print("compilation {}: {} sites, {} (site, junction) rows within +/-{} bp".format(a.compilation, n_sites, len(rows), a.flank_bp))
    print("  sites with a junction using them as 5'SS: {}  (of which unannotated-only: {})".format(n_matched, n_unann))
    print("  wrote {} and {}".format(junc_path, summ_path))


if __name__ == "__main__":
    main()
