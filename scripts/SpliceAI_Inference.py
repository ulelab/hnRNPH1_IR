import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
from keras.models import load_model  # type: ignore
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode  # type: ignore


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bed", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--batch-size", type=int, default=32,
                   help="sequences scored per model call [32]")
    p.add_argument("--resume", action="store_true",
                   help="reuse scores already written to <out>.scores")
    return p.parse_args()


def load_spliceai_models():
    paths = (f"models/spliceai{x}.h5" for x in range(1, 6))
    return [load_model(resource_filename("spliceai", p), compile=False) for p in paths]


def read_fasta_sequences(fasta_path):
    seqs, current, in_seq = [], [], False
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_seq:
                    seqs.append("".join(current).upper())
                current, in_seq = [], True
                continue
            current.append(line)
        if current:
            seqs.append("".join(current).upper())
        elif in_seq:
            seqs.append("")
    return seqs


def clean(seq):
    return seq.replace("-", "").replace(" ", "").replace("*", "").upper()


def score_all(seqs, models, context=10000, batch_size=32, ckpt=None, done=None):
    """Batched equivalent of the original per-sequence scoring.

    Identical maths: same padding, same 5-model mean, same donor-channel max.
    The only change is that sequences are stacked and passed to each model once
    per batch instead of once per sequence, which removes ~100k individual
    predict() dispatches (and the tf.function retracing they trigger).
    """
    n = len(seqs)
    scores = [None] * n
    done = done or {}

    pending = []
    for i, s in enumerate(seqs):
        if i in done:
            scores[i] = done[i]
            continue
        c = clean(s)
        if len(c) == 0 or len(c.replace("N", "")) == 0:
            scores[i] = -1.0
            continue
        pending.append((i, c))

    # group by length so every batch has a uniform shape (no retracing)
    groups = defaultdict(list)
    for i, c in pending:
        groups[len(c)].append((i, c))

    total = len(pending)
    processed = 0
    pad = "N" * (context // 2)
    fh = open(ckpt, "a", buffering=1) if ckpt else None
    try:
        for L in sorted(groups):
            items = groups[L]
            for b in range(0, len(items), batch_size):
                chunk = items[b:b + batch_size]
                # one_hot_encode returns int64; cast to float32 up front. TF casts
                # internally anyway, and int64 doubles the memory traffic per batch.
                X = np.stack([one_hot_encode(pad + c + pad) for _, c in chunk]).astype(np.float32)
                y = np.mean([m(X, training=False).numpy() for m in models], axis=0)
                for k, (i, c) in enumerate(chunk):
                    yy = y[k]
                    if len(c) >= context:
                        s0 = context // 2
                        region = yy[s0:s0 + len(c), :]
                    else:
                        region = yy[: len(c), :]
                    sc = float(np.max(region[:, 2]))
                    scores[i] = sc
                    if fh:
                        fh.write(f"{i}\t{sc:.6f}\n")
                processed += len(chunk)
                pct = 100.0 * processed / total if total else 100.0
                print(f"  scored {processed}/{total} ({pct:.1f}%)", file=sys.stderr, flush=True)
    finally:
        if fh:
            fh.close()

    for i in range(n):
        if scores[i] is None:
            scores[i] = -1.0
    return scores


def main():
    args = parse_args()
    ckpt = args.out + ".scores"

    done = {}
    if args.resume and os.path.exists(ckpt):
        with open(ckpt) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    done[int(parts[0])] = float(parts[1])
        print(f"resuming: {len(done)} scores already computed", file=sys.stderr)
    elif os.path.exists(ckpt):
        os.remove(ckpt)

    models = load_spliceai_models()

    with tempfile.TemporaryDirectory() as tmpdir:
        slopped_bed = os.path.join(tmpdir, "slop24.bed")
        seq_fa = os.path.join(tmpdir, "slop24.fa")
        subprocess.run(["bedtools", "slop", "-i", args.bed, "-g", args.genome, "-b", "24"],
                       stdout=open(slopped_bed, "w"), check=True)
        subprocess.run(["bedtools", "getfasta", "-fi", args.fasta, "-bed", slopped_bed,
                        "-s", "-name"], stdout=open(seq_fa, "w"), check=True)
        seqs = read_fasta_sequences(seq_fa)
        print(f"sequences: {len(seqs)}  batch size: {args.batch_size}", file=sys.stderr)
        scores = score_all(seqs, models, batch_size=args.batch_size, ckpt=ckpt, done=done)

    with open(args.bed) as fin, open(args.out, "w") as fout:
        for line, score in zip(fin, scores):
            cols = line.rstrip("\n").rstrip("\r").split("\t")
            while len(cols) < 5:
                cols.append(".")
            cols[4] = f"{score:.6f}"
            fout.write("\t".join(cols) + "\n")
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
