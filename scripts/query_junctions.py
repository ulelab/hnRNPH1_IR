#!/usr/bin/env python3
"""
Python 3.6-compatible utility to query Recount junction introns near RBP binding inference loci
"""

import argparse
import csv
import sqlite3
from pathlib import Path


def connect_recount(db_path):
    db_uri = "file:{}?mode=ro".format(str(Path(db_path).expanduser().resolve()))
    conn = sqlite3.connect(db_uri, uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def init_locusdb(conn):
    # Keep recount DB read-only and stage inference loci in attached in-memory DB.
    conn.execute("ATTACH DATABASE ':memory:' AS locusdb")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS locusdb.locus_coordinates (
            locus_id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            gene TEXT,
            splice_score REAL,
            strand TEXT NOT NULL,
            flag INTEGER
        )
        """
    )
    conn.execute(
        """
        CREATE INDEX IF NOT EXISTS locusdb.idx_locus_chr_strand_start_end
        ON locus_coordinates (chrom, strand, start, end)
        """
    )


def load_loci(conn, bed_path):
    conn.execute("DELETE FROM locusdb.locus_coordinates")
    with open(bed_path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = []
        for parts in reader:
            if not parts or len(parts) < 7:
                continue
            rows.append(
                (
                    parts[0],          # chrom
                    int(parts[1]),     # start
                    int(parts[2]),     # end
                    parts[3],          # gene
                    float(parts[4]),   # splice_score
                    parts[5],          # strand
                    int(parts[6]),     # flag
                )
            )
    conn.executemany(
        """
        INSERT INTO locusdb.locus_coordinates (
            chrom, start, end, gene, splice_score, strand, flag
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )
    conn.commit()
    return len(rows)


def query_locus_matches(conn, flank_bp, limit):
    sql = """
    SELECT
        l.locus_id,
        l.gene AS gene,
        l.chrom,
        l.start AS locus_start,
        l.end AS locus_end,
        l.strand AS locus_strand,
        l.splice_score,
        i.snaptron_id,
        i.start AS intron_start,
        i.end AS intron_end,
        i.samples_count,
        i.coverage_sum,
        i.coverage_avg,
        i.coverage_median,
        i.left_annotated,
        i.right_annotated,
        i.source_dataset_id,
        i.annotated
    FROM locusdb.locus_coordinates l
    JOIN intron i
      ON i.chrom = l.chrom
     AND i.strand = l.strand
     AND (
            ABS(i.start - l.start) <= ?
         OR ABS(i.end   - l.end)   <= ?
     )
    WHERE i.annotated = 1
    ORDER BY l.locus_id, i.coverage_sum DESC
    """
    params = [flank_bp, flank_bp]
    if limit is not None:
        sql += "\nLIMIT ?"
        params.append(limit)
    return conn.execute(sql, tuple(params)).fetchall()


def print_rows(rows, max_rows):
    shown = rows[:max_rows]
    print("Showing {} rows".format(len(shown)))
    for idx, row in enumerate(shown, 1):
        print("{:>4}: {}".format(idx, dict(row)))


def write_tsv(rows, out_tsv):
    if not rows:
        with open(out_tsv, "w") as fh:
            fh.write("")
        return
    fieldnames = rows[0].keys()
    with open(out_tsv, "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames)
        for row in rows:
            writer.writerow([row[k] for k in fieldnames])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Query Recount intron junctions near RBP binding inference loci coordinates."
    )
    parser.add_argument(
        "--db",
        required=True,
        help="Path to recount junctions.sqlite",
    )
    parser.add_argument(
        "--loci-bed",
        dest="loci_bed",
        default="Decoys_proteincoding_splicescores.bed",
        help="Path to inference loci BED file (7-column).",
    )
    parser.add_argument(
        "--inference-bed",
        dest="loci_bed",
        help="Alias for --loci-bed.",
    )
    parser.add_argument(
        "--schema-sql",
        default="locus_recount_schema.sql",
        help="Unused: kept for backward compatibility.",
    )
    parser.add_argument(
        "--flank-bp",
        type=int,
        default=5,
        help="Distance tolerance in bp when matching locus to intron boundaries.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional max rows to return from joined query. Default: no limit.",
    )
    parser.add_argument(
        "--print-rows",
        type=int,
        default=20,
        help="How many rows to print to stdout.",
    )
    parser.add_argument(
        "--out-tsv",
        default=None,
        help="Optional output TSV path for all joined rows.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    conn = connect_recount(args.db)
    try:
        init_locusdb(conn)
        n_loci = load_loci(conn, args.loci_bed)
        print("Loaded {} loci".format(n_loci))

        rows = query_locus_matches(conn, args.flank_bp, args.limit)
        print("Joined rows: {}".format(len(rows)))
        print_rows(rows, args.print_rows)
        if args.out_tsv:
            write_tsv(rows, args.out_tsv)
            print("Wrote TSV: {}".format(args.out_tsv))
    finally:
        conn.close()


if __name__ == "__main__":
    main()
