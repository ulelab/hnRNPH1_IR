#!/usr/bin/env python3
"""
Score ONE locus across all species in a MAF alignment with SpliceAI.

Output: one row per species with max1, pos1, max2, pos2 and d, plus the same
values packed as [max1][pos1][max2][pos2][d].

Usage:
  python3 scripts/maf_spliceai_single_locus.py \
      --maf reference/chr5.maf.gz \
      --region chr5:179620038-179620941 --strand - \
      --out results/hnrnph1_179620582_maf_spliceai.tsv
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile

import numpy as np

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--maf", required=True, help="MAF for the chromosome (.maf or .maf.gz)")
    p.add_argument("--region", default="chr5:179620038-179620941",
                   help="chrom:start-end, 1-based inclusive [HNRNPH1 decoy intron+49nt]")
    p.add_argument("--strand", default="-", choices=["+", "-"],
                   help="gene strand; '-' reverse-complements via msa_view -V [-]")
    p.add_argument("--out", required=True)
    p.add_argument("--msa-view", default="msa_view")
    p.add_argument("--mafs-in-region", default="mafsInRegion",
                   help="kent mafsInRegion binary; used to TRIM blocks to the "
                        "region exactly as the original pipeline did")
    p.add_argument("--cache-raw", default=None,
                   help="path to cache the extracted (untrimmed) region MAF; reused "
                        "on re-run so the multi-GB stream happens once")
    p.add_argument("--no-trim", action="store_true",
                   help="use the built-in streaming extractor instead (keeps whole "
                        "overlapping blocks, so positions shift - diagnostics only)")
    p.add_argument("--ref-species", default="hg38",
                   help="assembly whose coordinates --region refers to [hg38]")
    p.add_argument("--decoy-offset", type=int, default=358,
                   help="expected position of the decoy in the extracted window, "
                        "used only for a sanity line in the log [358]")
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--debug", action="store_true")
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def extract_region(maf_path, chrom, start, end, out_maf, ref="hg38", progress=True):
    """Copy MAF blocks whose REFERENCE row overlaps the region.
    """
    kept = 0
    seen = 0
    prefix = ref + "."
    with opener(maf_path) as fin, open(out_maf, "w") as fout:
        fout.write("##maf version=1\n")
        block, ref_ok = [], False
        for line in fin:
            if line.startswith("a"):
                if ref_ok and block:
                    fout.writelines(block)
                    fout.write("\n")   # blank line terminates a MAF block
                    kept += 1
                block, ref_ok = [line], False
                seen += 1
                if progress and seen % 200000 == 0:
                    print(f"    ...scanned {seen} blocks, kept {kept}", file=sys.stderr, flush=True)
            elif line.strip() == "":
                continue
            else:
                if block:
                    block.append(line)
                if line.startswith("s") and not ref_ok:
                    f = line.split()
                    if len(f) >= 6 and f[1].startswith(prefix):
                        c = f[1][len(prefix):]
                        bstart, blen = int(f[2]), int(f[3])
                        if c == chrom and not (bstart + blen <= start or bstart >= end):
                            ref_ok = True
        if ref_ok and block:
            fout.writelines(block)
            fout.write("\n")
            kept += 1
    return kept


def maf_to_fasta(maf_file, msa_view_bin, strand, out_fa, debug=False):
    cmd = [msa_view_bin, "--in-format", "MAF", "--out-format", "FASTA"]
    if strand == "-":
        cmd.append("-V")
    cmd.append(maf_file)
    with open(out_fa, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if r.returncode != 0:
        sys.exit(f"ERROR: msa_view failed: {r.stderr.decode()[:400]}")


def read_fasta(path):
    seqs, hdr, cur = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    seqs[hdr] = "".join(cur).upper()
                h = line[1:].strip()
                hdr = h.split(".")[0] if "." in h else h.split()[0]
                cur = []
            else:
                cur.append(line)
    if hdr is not None:
        seqs[hdr] = "".join(cur).upper()
    return seqs


def find_top2_scores_length_adjusted(donor_probs, hg38_length, current_seq_length):
    """Verbatim port of the v2.5 function (100 nt length-adjusted max1 window)."""
    if len(donor_probs) == 0:
        return (-1.0, -1, -1.0, -1, -1.0)
    donor_probs = np.array(donor_probs)
    d = abs(current_seq_length - hg38_length) / hg38_length if hg38_length > 0 else 0.0

    if current_seq_length <= hg38_length:
        search_end = min(100, len(donor_probs))
    else:
        search_end = min(int(100 * (current_seq_length / hg38_length)), len(donor_probs))
    if search_end == 0:
        return (-1.0, -1, -1.0, -1, d)

    region = donor_probs[:search_end]
    max1 = float(np.max(region))
    pos1 = int(np.argmax(region))

    s2 = pos1 + 20
    if s2 >= len(donor_probs):
        return (max1, pos1, -1.0, -1, d)
    tail = donor_probs[s2:]
    max2 = float(np.max(tail))
    pos2 = int(np.argmax(tail)) + s2
    return (max1, pos1, max2, pos2, d)


def main():
    a = parse_args()
    chrom, rng = a.region.split(":")
    rstart, rend = (int(x.replace(",", "")) for x in rng.split("-"))
    start0 = rstart - 1  # MAF is 0-based

    with tempfile.TemporaryDirectory() as tmp:
        sub_maf = os.path.join(tmp, "region.maf")
        fa = os.path.join(tmp, "region.fa")
        print(f"extracting {chrom}:{rstart}-{rend} from {a.maf} ...", file=sys.stderr)
        raw_maf = a.cache_raw or os.path.join(tmp, "raw.maf")
        if a.cache_raw and os.path.exists(a.cache_raw) and os.path.getsize(a.cache_raw) > 0:
            print(f"  reusing cached extraction: {a.cache_raw}", file=sys.stderr)
        else:
            n = extract_region(a.maf, chrom, start0, rend, raw_maf, ref=a.ref_species)
            print(f"  MAF blocks overlapping region: {n}", file=sys.stderr)
            if n == 0:
                sys.exit("ERROR: no MAF blocks overlapped the region.")

        if a.no_trim:
            sub_maf = raw_maf
        else:
            # Trim to the exact window with kent mafsInRegion
            bed = os.path.join(tmp, "region.bed")
            with open(bed, "w") as fh:
                fh.write(f"{chrom}\t{start0}\t{rend}\n")
            r = subprocess.run([a.mafs_in_region, bed, sub_maf, raw_maf],
                               stderr=subprocess.PIPE)
            if r.returncode != 0 or not os.path.exists(sub_maf):
                sys.exit(f"ERROR: mafsInRegion failed: {r.stderr.decode()[:300]}\n"
                         f"  (pass --no-trim to skip, but positions will shift)")
        maf_to_fasta(sub_maf, a.msa_view, a.strand, fa, a.debug)
        seqs = read_fasta(fa)
    print(f"  species in alignment: {len(seqs)}", file=sys.stderr)
    if "hg38" not in seqs:
        sys.exit("ERROR: hg38 not present in the alignment; cannot set the reference length.")

    from keras.models import load_model
    from pkg_resources import resource_filename
    from spliceai.utils import one_hot_encode

    context = 10000
    models = [load_model(resource_filename("spliceai", f"models/spliceai{i}.h5"), compile=False)
              for i in range(1, 6)]
    pad = "N" * (context // 2)

    hg38_len = len(seqs["hg38"].replace("-", ""))
    print(f"  hg38 ungapped length: {hg38_len}", file=sys.stderr)

    names, cleaned = [], []
    for sp, s in seqs.items():
        c = s.replace("-", "").replace(".", "").replace("*", "").upper()
        names.append(sp)
        cleaned.append(c)

    rows = {}
    order = sorted(range(len(cleaned)), key=lambda i: len(cleaned[i]))
    i = 0
    while i < len(order):
        grp = [order[i]]
        L = len(cleaned[order[i]])
        j = i + 1
        while j < len(order) and len(cleaned[order[j]]) == L and len(grp) < a.batch_size:
            grp.append(order[j]); j += 1
        i = j
        valid = [k for k in grp if L > 0 and cleaned[k].replace("N", "")]
        for k in set(grp) - set(valid):
            rows[names[k]] = (-1.0, -1, -1.0, -1, -1.0)
        if not valid:
            continue
        X = np.stack([one_hot_encode(pad + cleaned[k] + pad) for k in valid]).astype(np.float32)
        Y = np.mean([m(X, training=False).numpy() for m in models], axis=0)
        for bi, k in enumerate(valid):
            y = Y[bi]
            seq_len = len(cleaned[k])
            donor = y[:seq_len, 2] if seq_len < context else y[context // 2: context // 2 + seq_len, 2]
            rows[names[k]] = find_top2_scores_length_adjusted(donor, hg38_len, seq_len)
        print(f"  scored {len(rows)}/{len(names)}", file=sys.stderr, flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("species\tmax1\tpos1\tmax2\tpos2\td\tpacked\n")
        for sp in names:
            m1, p1, m2, p2, d = rows[sp]
            fh.write(f"{sp}\t{m1:.6f}\t{p1}\t{m2:.6f}\t{p2}\t{d:.6f}\t"
                     f"[{m1:.6f}][{p1}][{m2:.6f}][{p2}][{d:.6f}]\n")
    print(f"wrote {a.out}  ({len(names)} species)", file=sys.stderr)

    h = rows.get("hg38")
    if h:
        print(f"  hg38 check -> max1={h[0]:.6f} pos1={h[1]} max2={h[2]:.6f} pos2={h[3]}", file=sys.stderr)
        print(f"  (decoy expected at pos2={a.decoy_offset})", file=sys.stderr)


if __name__ == "__main__":
    main()
