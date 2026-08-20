#!/usr/bin/env python3
"""
Download per-sample genome crosslink BEDs (.genome.xl.bed) from Flow for the
samples listed in a selection TSV.

Picks the `<sample>_R1.genome.xl.bed` variant. That matters: PRPF8 samples carry
two crosslink files from two different pipelines -

  <sample>.genome.xl.bed      "hanalysis clipseq"  (NFCORE_CLIPSEQ)
  <sample>_R1.genome.xl.bed   "CLIP-Seq"           (CLIPSEQ:CALC_GENOME_CROSSLINKS)

- while the SNRPB samples only have the CLIP-Seq `_R1` output. Taking `_R1`
everywhere keeps the PRPF8 and SmB categories on the same pipeline, so a
difference between them is biology and not a processing artefact.

Auth: Flow JWT from ~/.config/flow/api-token.

Usage:
  python3 scripts/flow_fetch_crosslinks.py \
      --selection data/flow_prpf8_snrpb_samples_selected.tsv \
      --outdir data/CLIP
"""
import argparse, csv, json, os, sys, time, urllib.parse, urllib.request

BASE = "https://app.flow.bio/api"
TOKEN = open(os.path.expanduser("~/.config/flow/api-token")).read().strip()
HDRS = {"Authorization": f"Bearer {TOKEN}"}


def get(path, **params):
    url = f"{BASE}/{path}" + ("?" + urllib.parse.urlencode(params) if params else "")
    req = urllib.request.Request(url, headers=HDRS)
    for a in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except Exception:
            if a == 3:
                raise
            time.sleep(2 * (a + 1))


def paged(path, key, **params):
    page, seen, total = 1, 0, None
    while True:
        d = get(path, count=100, page=page, **params)
        items = d.get(key, [])
        if total is None:
            total = d.get("count")
        if not items:
            break
        yield from items
        seen += len(items)
        if total is not None and seen >= total:
            break
        page += 1


def download(data_id, filename, dest, expect):
    """Stream a data file to dest; skip if already present at the right size."""
    if os.path.exists(dest) and os.path.getsize(dest) == expect:
        print(f"    already present ({expect:,} bytes) - skipping")
        return True
    url = f"{BASE}/downloads/{data_id}/{urllib.parse.quote(filename)}"
    req = urllib.request.Request(url, headers=HDRS)
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
            got = 0
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                got += len(chunk)
        if expect and got != expect:
            print(f"    SIZE MISMATCH: got {got:,}, expected {expect:,}", file=sys.stderr)
            os.rename(tmp, dest + ".bad")
            return False
        os.rename(tmp, dest)
        print(f"    downloaded {got:,} bytes")
        return True
    except Exception as e:
        print(f"    FAILED: {e}", file=sys.stderr)
        if os.path.exists(tmp):
            os.remove(tmp)
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--suffix", default="_R1.genome.xl.bed")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.selection), delimiter="\t"))
    print(f"{len(rows)} samples in {a.selection}\n")
    ok = fail = 0
    for i, r in enumerate(rows, 1):
        cohort = "PRPF8" if r["cohort"] == "PRPF8" else "SmB"
        d = os.path.join(a.outdir, cohort)
        os.makedirs(d, exist_ok=True)
        print(f"[{i}/{len(rows)}] {r['name']}  ({cohort})")
        files = [f for f in paged(f"samples/{r['sample_id']}/data", "data")
                 if f["filename"].endswith(a.suffix)]
        if not files:
            print(f"    NO FILE ending {a.suffix}", file=sys.stderr)
            fail += 1
            continue
        if len(files) > 1:
            print(f"    {len(files)} candidates, taking largest", file=sys.stderr)
            files.sort(key=lambda f: -f["size"])
        f = files[0]
        dest = os.path.join(d, f["filename"])
        if download(f["id"], f["filename"], dest, f["size"]):
            ok += 1
        else:
            fail += 1

    print(f"\ndone: {ok} ok, {fail} failed")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
