#!/usr/bin/env python3
"""Sashimi-style junction plot and decoy donor usage for the HNRNPH1 intron-4 decoy.

Junction counts come from the Snaptron REST API (Recount3 compilations srav3h, gtexv2, tcgav2),
cached under results/snaptron/. The gene model is the HNRNPH1 MANE Select transcript from
reference/gencode.v49.annotation.gtf.gz.

Outputs:
  figures/hnrnph1_intron4_sashimi.pdf/.png      arc plot, srav3h read counts
  figures/hnrnph1_decoy_usage.pdf/.png          donor usage per compilation and per GTEx tissue: the intron-4
                                                decoy (top) and the annotated alternative 5'SS HsaALTD0003092-2 (bottom)
  results/hnrnph1_donor_usage_gtex_tissue.tsv   per-tissue table behind the second figure

Usage (repository root, rbpnet-env):  python scripts/plot_hnrnph1_junctions.py
"""
import gzip
import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Rectangle

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "results" / "snaptron"
FIG = ROOT / "figures"
CHROM, STRAND = "chr5", "-"
WIN = (179616800, 179621600)
DECOY = 179620582            # last intronic base 3' of the decoy donor = junction end on the minus strand
CANON_DONOR = 179620891      # intron 4 canonical donor (junction end)
CANON_ACCEPTOR = 179619408   # intron 4 canonical acceptor (junction start)
# Annotated alternative 5'SS (Vast-DB HsaALTD0003092-2, donor chr5:179,623,595) and the event's reference
# donor (HsaALTD0003092-1, chr5:179,623,636), as junction ends on the minus strand (donor base - 1).
ALT_WIN = (179622900, 179624200)
ALT_DONOR = 179623594
ALT_REF = 179623635
TEAL = "#1B7837"
COMPS = ["srav3h", "gtexv2", "tcgav2"]
MIN_READS_DECOY = 100        # decoy-donor junctions drawn at >= this many reads
MIN_READS_UNANNOTATED = 20000  # other unannotated junctions drawn only when this strong
RED, GREY, DARK = "#B2182B", "#8C8C8C", "#2D3142"
COLS = ["type", "snaptron_id", "chrom", "start", "end", "length", "strand", "annotated", "left_motif", "right_motif",
        "left_annotated", "right_annotated", "samples", "samples_count", "coverage_sum", "coverage_avg",
        "coverage_median", "source_dataset_id"]


def fetch(url, path):
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=600) as r:
            path.write_bytes(r.read())
    return path


def junctions(comp, win=WIN, tag="hnrnph1"):
    url = f"https://snaptron.cs.jhu.edu/{comp}/snaptron?regions={CHROM}:{win[0]}-{win[1]}&rfilter=strand:{STRAND}"
    df = pd.read_csv(fetch(url, CACHE / f"{tag}_{comp}.tsv"), sep="\t", header=0, names=COLS)
    return df[(df.end >= win[0]) & (df.end <= win[1])].copy()


def donor_reads(df, ends):
    return {e: int(df.loc[df.end == e, "coverage_sum"].sum()) for e in ends}


def tissue_usage(g, meta, donor, ref):
    # Per GTEx tissue: donor reads / (donor + reference donor reads), summed over the tissue's samples.
    ps = per_sample(g[g.end.isin([donor, ref])]).merge(g[["snaptron_id", "end"]], on="snaptron_id")
    ps["kind"] = np.where(ps.end == donor, "site", "reference")
    per = ps.pivot_table(index="rail_id", columns="kind", values="reads", aggfunc="sum", fill_value=0).reset_index()
    for k in ("site", "reference"):
        if k not in per:
            per[k] = 0
    per = per.merge(meta, on="rail_id", how="left")
    tis = per.groupby("SMTS").agg(samples=("rail_id", "size"), samples_with_site=("site", lambda x: int((x > 0).sum())),
                                  site_reads=("site", "sum"), reference_reads=("reference", "sum")).reset_index()
    tis["pct"] = 100 * tis.site_reads / (tis.site_reads + tis.reference_reads)
    return tis


def per_sample(df):
    rows = []
    for jid, s in zip(df.snaptron_id, df.samples):
        for tok in str(s).strip(",").split(","):
            if tok:
                sid, n = tok.split(":")
                rows.append((jid, int(sid), int(n)))
    return pd.DataFrame(rows, columns=["snaptron_id", "rail_id", "reads"])


def exons():
    keep = []
    with gzip.open(ROOT / "reference" / "gencode.v49.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if 'gene_name "HNRNPH1"' in line and "\texon\t" in line and "MANE_Select" in line:
                f = line.split("\t")
                num = int(f[8].split('exon_number ')[1].split(";")[0])
                keep.append((int(f[3]), int(f[4]), num))
    return pd.DataFrame(keep, columns=["start", "end", "exon"])


def usage(df):
    decoy = df[df.end == DECOY].coverage_sum.sum()
    canon = df[df.end == CANON_DONOR].coverage_sum.sum()
    return decoy, canon


def sashimi(df, ex):
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.set_xlim(WIN[1], WIN[0])        # minus strand: transcript runs left to right
    ax.set_ylim(-1.15, 1.25)
    for r in ex.itertuples():
        if r.end >= WIN[0] and r.start <= WIN[1]:
            ax.add_patch(Rectangle((r.start, -0.06), r.end - r.start, 0.12, color=DARK, zorder=3))
            ax.text((r.start + r.end) / 2, 0.1, f"E{r.exon}", ha="center", va="bottom", fontsize=8, color=DARK)
    ax.plot(WIN, [0, 0], color=DARK, lw=1, zorder=2)
    is_d = df.end == DECOY
    show = df[(is_d & (df.coverage_sum >= MIN_READS_DECOY)) | (~is_d & (df.annotated == 1) & (df.coverage_sum >= MIN_READS_DECOY))
              | (~is_d & (df.annotated == 0) & (df.coverage_sum >= MIN_READS_UNANNOTATED))].sort_values("coverage_sum")
    top = np.log10(show.coverage_sum.max())
    for r in show.itertuples():
        is_decoy = r.end == DECOY
        h = 0.25 + 0.7 * np.log10(r.coverage_sum) / top
        sign = -1 if is_decoy else 1
        col = RED if is_decoy else (DARK if r.annotated == 1 else GREY)
        lw = 0.6 + 2.6 * np.log10(r.coverage_sum) / top
        ax.add_patch(Arc(((r.start + r.end) / 2, 0), r.end - r.start, 2 * h, theta1=0 if sign > 0 else 180,
                         theta2=180 if sign > 0 else 360, color=col, lw=lw, alpha=0.9, zorder=1))
        if (is_decoy and r.coverage_sum >= 400) or (r.annotated == 1 and r.coverage_sum >= 1e6):
            ax.text((r.start + r.end) / 2, sign * (h + 0.04), f"{r.coverage_sum:,}", ha="center",
                    va="bottom" if sign > 0 else "top", fontsize=8, color=col,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
    ax.axvline(DECOY, color=RED, lw=1, ls="--", zorder=0)
    ax.text(DECOY - 40, 1.18, "decoy donor\nchr5:179,620,583", color=RED, ha="left", va="top", fontsize=8)
    ax.axvline(CANON_DONOR, color=DARK, lw=0.8, ls=":", zorder=0)
    ax.text(CANON_DONOR + 40, 1.18, "intron 4\ndonor", color=DARK, ha="right", va="top", fontsize=8)
    ax.set_yticks([])
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.set_xlabel(f"{CHROM} (minus strand, transcript 5'→3' left to right)")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title("HNRNPH1 junctions across all SRA samples (Recount3 srav3h), read counts on arcs\n"
                 f"Above: annotated junctions (black), unannotated >= {MIN_READS_UNANNOTATED:,} reads (grey). "
                 f"Below: junctions using the decoy as donor, >= {MIN_READS_DECOY} reads (red)", fontsize=10)
    return fig


def main():
    ex = exons()
    data = {c: junctions(c) for c in COMPS}
    fig = sashimi(data["srav3h"], ex)
    fig.savefig(FIG / "hnrnph1_intron4_sashimi.pdf", bbox_inches="tight")
    fig.savefig(FIG / "hnrnph1_intron4_sashimi.png", dpi=170, bbox_inches="tight")
    print("decoy-donor junctions (srav3h):")
    print(data["srav3h"][data["srav3h"].end == DECOY][["start", "end", "length", "annotated", "samples_count", "coverage_sum"]]
          .sort_values("coverage_sum", ascending=False).to_string(index=False))

    alt = {c: junctions(c, ALT_WIN, "hnrnph1_alt5ss") for c in COMPS}
    meta = pd.read_csv(fetch("https://snaptron.cs.jhu.edu/data/gtexv2/samples.tsv", CACHE / "gtexv2_samples.tsv"),
                       sep="\t", usecols=lambda c: c in {"rail_id", "SMTS", "SMTSD"}, low_memory=False)
    sites = [
        ("decoy", "Intron-4 decoy donor (chr5:179,620,583) vs canonical intron-4 donor", data, DECOY, CANON_DONOR, RED),
        ("alt5ss", "Annotated alternative 5'SS HsaALTD0003092-2 (chr5:179,623,595) vs the event's reference donor", alt, ALT_DONOR, ALT_REF, TEAL),
    ]
    comp_rows, tis_all = [], []
    for key, _, dat, donor, ref, _ in sites:
        for c in COMPS:
            r = donor_reads(dat[c], [donor, ref])
            comp_rows.append({"site": key, "compilation": c, "site_reads": r[donor], "reference_reads": r[ref],
                              "pct": 100 * r[donor] / (r[donor] + r[ref])})
        t = tissue_usage(dat["gtexv2"], meta, donor, ref); t.insert(0, "site", key); tis_all.append(t)
    comp = pd.DataFrame(comp_rows)
    tis = pd.concat(tis_all)
    tis = tis[tis.groupby("SMTS").samples.transform("min") >= 20]
    order = tis[tis.site == "decoy"].sort_values("pct", ascending=False).SMTS.tolist()
    tis.to_csv(ROOT / "results" / "hnrnph1_donor_usage_gtex_tissue.tsv", sep="\t", index=False)
    print(comp.to_string(index=False))

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), gridspec_kw={"width_ratios": [1, 3.2]})
    for row, (key, title, _, _, _, col) in enumerate(sites):
        c = comp[comp.site == key].reset_index(drop=True)
        ax = axes[row, 0]
        ax.bar(c.compilation, c.pct, color=col, alpha=0.85)
        for i, r in c.iterrows():
            ax.text(i, r.pct, f"{r.pct:.3g}%\n{r.site_reads:,} /\n{r.reference_reads:,}", ha="center", va="bottom", fontsize=7)
        ax.set_ylim(0, c.pct.max() * 1.35)
        ax.set_ylabel("% of donor events at this site")
        ax.set_title("By compilation (site reads / reference reads)", fontsize=9)
        t = tis[tis.site == key].set_index("SMTS").loc[order].reset_index()
        ax = axes[row, 1]
        ax.bar(range(len(t)), t.pct, color=col, alpha=0.85)
        for i, r in t.iterrows():
            ax.text(i, r.pct, f"{r.samples_with_site}/{r.samples}", ha="center", va="bottom", fontsize=6, rotation=90)
        ax.set_ylim(0, t.pct.max() * 1.18)
        ax.set_xticks(range(len(t)), t.SMTS if row == 1 else [""] * len(t), rotation=60, ha="right", fontsize=8)
        ax.set_title(f"{title}\nGTEx tissues, ordered by decoy use (label: samples with a junction at the site / samples)", fontsize=9)
    fig.suptitle("HNRNPH1 donor usage: the intron-4 decoy vs an annotated alternative 5' splice site (Recount3)", y=1.0)
    fig.tight_layout()
    fig.savefig(FIG / "hnrnph1_decoy_usage.pdf", bbox_inches="tight")
    fig.savefig(FIG / "hnrnph1_decoy_usage.png", dpi=170, bbox_inches="tight")
    print(tis.pivot(index="SMTS", columns="site", values="pct").loc[order].round(4).head(8).to_string())
    x = tis.pivot(index="SMTS", columns="site", values="pct")
    from scipy.stats import spearmanr
    rho, pv = spearmanr(x.decoy, x.alt5ss)
    print(f"Spearman rho across tissues, decoy vs alt 5'SS use: {rho:.2f} (p = {pv:.2g})")


if __name__ == "__main__":
    main()
