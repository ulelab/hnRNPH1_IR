-- Schema for joining loci coordinates against Recount introns.
-- Use this as an attached in-memory schema from Python.

CREATE TABLE IF NOT EXISTS locus_coordinates (
    locus_id INTEGER PRIMARY KEY AUTOINCREMENT,
    chrom TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    gene TEXT,
    splice_score REAL,
    strand TEXT NOT NULL,
    flag INTEGER
);

CREATE INDEX IF NOT EXISTS idx_locus_chr_strand_start_end
    ON locus_coordinates (chrom, strand, start, end);
