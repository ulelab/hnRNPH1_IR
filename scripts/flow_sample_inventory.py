#!/usr/bin/env python3
"""
Inventory Flow CLIP samples for a given purification target.

Used to scope the re-analysis: which PRPF8 and SNRPB (SmB/B') samples exist on
Flow, and how deep each one is. Writes a TSV of assay, species, cell type,
condition, read counts and crosslink (deduplicated read) counts.

Auth: reads a Flow JWT from ~/.config/flow/api-token (see the Flow API notes in
the lab memory). Token must be current - mint a new one via flowbio if it 401s.

  python3 scripts/flow_sample_inventory.py

Counts come from pipeline logs, not from the .xl.bed files themselves:
  total_reads      STAR Log.final.out  "Number of input reads"
  uniquely_mapped  STAR Log.final.out  "Uniquely mapped reads number"
  dedup_reads      *.genomeUnique.dedup_UMICollapse.log  "Number of reads after deduplicating"
  unique_positions  ... "Number of unique alignment positions"
The source filename for each is recorded in the reads_src / dedup_src columns.
"""
import csv, re, sys
import json, os, time, urllib.parse, urllib.request

import json, os, time, urllib.parse, urllib.request

BASE = "https://app.flow.bio/api"
TOKEN = open(os.path.expanduser("~/.config/flow/api-token")).read().strip()


def get(path, **params):
    url = f"{BASE}/{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Authorization": f"Bearer {TOKEN}"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read())
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))


def paged(path, key, **params):
    """Yield every item across all pages (API caps count at 100)."""
    page, seen, total = 1, 0, None
    while True:
        d = get(path, count=100, page=page, **params)
        items = d.get(key, [])
        if total is None:
            total = d.get("count")
        if not items:
            break
        for it in items:
            yield it
        seen += len(items)
        if total is not None and seen >= total:
            break
        page += 1
        if page > 200:
            break


STAR_RE = re.compile(r"Log\.final\.out$")
DEDUP_RE = re.compile(r"dedup_UMICollapse\.log$")
# Prefer the unique-genome dedup log; these are the two spellings seen.
DEDUP_PREF = ("genomeunique", "unique_genome")


def num(text, label):
    m = re.search(re.escape(label) + r"\s*\|?\s*\t*\s*([\d.]+)", text)
    return int(float(m.group(1))) if m else None


def contents(data_id):
    c = get(f"data/{data_id}/contents")
    if isinstance(c, dict):
        return c.get("contents") or c.get("content") or ""
    return c or ""


def condition_of(name, comments, target):
    """Best-effort condition, from explicit name tokens."""
    n = name.lower()
    bits = []
    if "siprpf8" in n:
        bits.append("siPRPF8 knockdown")
    if "sictrl" in n or "sicontrol" in n or "siluc" in n:
        bits.append("si-control")
    if "stringentlysis" in n:
        bits.append("stringent lysis")
    if "flag" in n:
        bits.append("FLAG-tagged")
    if target == "SMInput":
        bits.append("size-matched input")
    if re.search(r"\bwt\b", n):
        bits.append("WT")
    return "; ".join(bits)


def collect(sample_stub, cohort):
    sid = sample_stub["id"]
    d = get(f"samples/{sid}")
    m = d.get("metadata", {})
    g = lambda k: (m.get(k) or {}).get("value", "") or ""

    row = {
        "cohort": cohort,
        "sample_id": sid,
        "name": d.get("name", ""),
        "purification_target": g("purification_target"),
        "assay": g("experimental_method") or (
            d.get("sample_type") if isinstance(d.get("sample_type"), str) else ""),
        "species": (d.get("organism") or {}).get("name", ""),
        "cell_type": g("source"),
        "condition": condition_of(d.get("name", ""), g("comments"), g("purification_target")),
        "organisation": g("organisation"),
        "pi": g("pi"),
        "geo": g("geo"),
        "project": (d.get("project") or {}).get("name", ""),
        "paired": sample_stub.get("is_paired", ""),
        "total_reads": "", "uniquely_mapped": "", "dedup_reads": "",
        "unique_positions": "", "reads_src": "", "dedup_src": "",
        "n_data_files": 0, "processed": "no",
    }

    try:
        files = list(paged(f"samples/{sid}/data", "data"))
    except Exception as e:
        row["processed"] = f"error: {e}"
        return row
    row["n_data_files"] = len(files)
    if not files:
        return row
    row["processed"] = "yes"

    star = [f for f in files if STAR_RE.search(f["filename"])]
    if star:
        star.sort(key=lambda f: f["filename"])
        f = star[0]
        txt = contents(f["id"])
        row["total_reads"] = num(txt, "Number of input reads") or ""
        row["uniquely_mapped"] = num(txt, "Uniquely mapped reads number") or ""
        row["reads_src"] = f["filename"]

    dd = [f for f in files if DEDUP_RE.search(f["filename"])]
    if dd:
        pref = [f for f in dd if any(p in f["filename"].lower() for p in DEDUP_PREF)]
        f = (pref or dd)[0]
        txt = contents(f["id"])
        row["dedup_reads"] = num(txt, "Number of reads after deduplicating") or ""
        row["unique_positions"] = num(txt, "Number of unique alignment positions") or ""
        row["dedup_src"] = f["filename"]
    return row


def main():
    cohorts = [
        ("PRPF8", dict(purification_target="PRPF8")),
        ("SNRPB (SmB/B')", dict(purification_target="SNRPB")),
    ]
    stubs = {}
    for label, params in cohorts:
        for s in paged("samples/search", "samples", **params):
            stubs.setdefault(s["id"], (s, label))
        print(f"{label}: cumulative {len(stubs)}", file=sys.stderr)

    # SMInput controls explicitly paired to PRPF8 eCLIP (not PRPF8 pulldowns)
    for s in paged("samples/search", "samples", name="PRPF8"):
        if s["id"] not in stubs and s["name"].startswith("PRPF8_Control"):
            stubs[s["id"]] = (s, "SMInput control (PRPF8 eCLIP)")

    print(f"total samples to profile: {len(stubs)}", file=sys.stderr)
    rows = []
    for i, (sid, (stub, label)) in enumerate(sorted(stubs.items()), 1):
        try:
            rows.append(collect(stub, label))
        except Exception as e:
            print(f"  FAILED {stub['name']}: {e}", file=sys.stderr)
            rows.append({"cohort": label, "sample_id": sid,
                         "name": stub.get("name", ""), "processed": f"error: {e}"})
        if i % 5 == 0:
            print(f"  {i}/{len(stubs)}", file=sys.stderr)

    cols = ["cohort", "name", "sample_id", "purification_target", "assay", "species",
            "cell_type", "condition", "total_reads", "uniquely_mapped", "dedup_reads",
            "unique_positions", "paired", "processed", "n_data_files", "organisation",
            "pi", "geo", "project", "reads_src", "dedup_src"]
    out = "flow_prpf8_snrpb_samples.tsv"
    with open(out, "w", newline="") as fh:
        # lineterminator="\n": csv.writer defaults to CRLF, which rides along on the
        # last column and breaks downstream numeric parsing.
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r.get("cohort", ""), r.get("name", ""))):
            w.writerow(r)
    print(f"wrote {out} ({len(rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
