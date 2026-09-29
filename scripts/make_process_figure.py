"""Process figure for dataset generation: decoys, supported cryptics, and the two
feature tables. Counts are those of the current run (see README, Workflow).

Usage: python scripts/make_process_figure.py   -> figures/dataset_generation_process.{pdf,png}
"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures" / "dataset_generation_process"

RED, BLUE, PURPLE, GREY, DARK = "#B2182B", "#2166AC", "#762A83", "#8C8C8C", "#333333"
FILL_RED, FILL_BLUE, FILL_PURPLE, FILL_GREY = "#F6DCDD", "#DCE8F4", "#E6DAEB", "#F0F0F0"

fig, ax = plt.subplots(figsize=(11, 15.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 148); ax.axis("off")
ax.invert_yaxis()

def box(x, y, w, h, title, body, tag=None, edge=DARK, fill="white", title_color=None, fs=8.6, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                linewidth=lw, edgecolor=edge, facecolor=fill, zorder=2))
    if tag:
        ax.text(x + 1.6, y - 0.9, tag, fontsize=7.2, color=edge, fontweight="bold", ha="left", va="bottom",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=edge, linewidth=0.8), zorder=3)
    ax.text(x + w / 2, y + 2.2, title, fontsize=fs + 1.2, fontweight="bold", ha="center", va="top",
            color=title_color or edge, zorder=3)
    ax.text(x + w / 2, y + 5.6, body, fontsize=fs, ha="center", va="top", color=DARK, zorder=3, linespacing=1.35)

def arrow(x1, y1, x2, y2, color=DARK, label=None, lx=1.2, ly=0, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=13,
                                 linewidth=1.3, color=color, linestyle=ls, zorder=1,
                                 shrinkA=0, shrinkB=0))
    if label:
        ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, fontsize=7.2, color=color, ha="left", va="center", zorder=3,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none"))

# ---- trunk -----------------------------------------------------------------
box(20, 3, 60, 11, "SpliceAI whole-intron inference",
    "Donor score at every intronic nucleotide of all Vast-DB HsaIN introns (400 nt exonic flanks)\n"
    "filtered to score >= 0.05\nSplice_All.filtered.05min.bed   794,192 sites", tag="SpliceAI")
arrow(50, 14, 50, 18.5)

box(20, 18.5, 60, 10, "Protein-coding filter",
    "decoy_exon_overlaps.Rmd Part 1: keep sites overlapping a protein_coding gene feature, same strand\n"
    "(gencode.v49.annotation.gtf)\nspliceai_05min_proteincoding.bed   550,282 sites", tag="GenomicRanges")
arrow(50, 28.5, 50, 33)

box(20, 33, 60, 11, "Canonical splice-site split",
    "split_canonical_exonic.sh: intersect with Wide_canonical_splice_sites.bed\n"
    "(192,965 windows, +/-50 nt around each canonical splice site)\n"
    "deep_intronic_splice_sites.bed   367,089 sites", tag="bedtools intersect -v -s")
# exonic side branch
box(84, 35, 15, 7.5, "Exonic borders", "exonic_splice_sites.bed\n183,193 sites", edge=GREY, fill=FILL_GREY, fs=7.4)
arrow(80, 38.75, 84, 38.75, color=GREY, style="-|>")
arrow(50, 44, 50, 48.5)

box(20, 48.5, 60, 10, "Exon-overlap removal",
    "decoy_exon_overlaps.Rmd Part 2: remove sites overlapping any GENCODE v49 exon, either strand\n"
    "(including retained_intron transcripts)\ndeep_intronic_noexon_splice_sites.bed   289,640 sites", tag="GenomicRanges")
arrow(50, 58.5, 50, 63)

box(20, 63, 60, 12, "CLIP support at the site",
    "intersect_spliceai_support.sh -w 0\n"
    "RBPnet PRPF8 predictions  and/or  PRPF8 eCLIP peaks  and/or  SmB iCLIP peaks\n"
    "(Clippy peaks on merged crosslink tracks from Flow)", tag="bedtools window")

# ---- two branches ----------------------------------------------------------
arrow(35, 75, 22, 82, color=RED, label="signal at site", lx=-14, ly=-1.5)
arrow(65, 75, 78, 82, color=BLUE, label="no signal at site", lx=1.5, ly=-1.5)

box(2, 82, 40, 9.5, "Decoys", "decoys.bed   19,942 sites\nRBPnet 9,526 / PRPF8 6,805 / SmB 5,668", edge=RED, fill=FILL_RED)
box(58, 82, 40, 9.5, "Cryptic sites", "cryptic_sites.bed   269,698 sites", edge=BLUE, fill=FILL_BLUE)

arrow(78, 91.5, 78, 95.5, color=BLUE)
box(58, 95.5, 40, 15, "Canonical 5'SS support",
    "compile_decoy_intron_data.Rmd Parts 1-2\n"
    "site -> Vast-DB intron -> canonical donor window\n(Canonical_splice_sites.bed, +/-5 nt)\n"
    "intersect: RBPnet and/or PRPF8 and/or SmB\ncryptics_supported.bed   191,904 sites",
    tag="GenomicRanges + bedtools", edge=BLUE, fill=FILL_BLUE)

# decoy branch: long arrow down to rescore
arrow(22, 91.5, 22, 115, color=RED)
arrow(78, 110.5, 78, 115, color=BLUE)

# ---- rescore + filters (both branches) --------------------------------------
box(2, 115, 96, 12.5, "Local re-score and final filters",
    "SpliceAI_Inference.py: re-score each site over a 49 nt window (SLURM slices for the cryptics)\n"
    "keep local donor score >= 0.1  ->  remove annotated alternative 5' splice sites (Vast-DB HsaALTD donors)\n"
    "decoys_final.bed   6,506            cryptics_supported_final.bed   17,721\n"
    "MaxEntScan subsets: site MaxEnt >= 8 / 9 / 10   ->   5,420 / 4,027 / 2,450   |   11,954 / 7,585 / 3,698",
    tag="SpliceAI + MaxEntScan", edge=PURPLE, fill=FILL_PURPLE)

# ---- feature tables --------------------------------------------------------
arrow(30, 127.5, 30, 132, color=PURPLE)
arrow(70, 127.5, 70, 132, color=PURPLE)

box(2, 132, 46, 14, "Site feature table",
    "compile_decoy_intron_data.Rmd Parts 3-4\n"
    "decoy_intron_features_final.tsv   24,218 sites, site_class = decoy / cryptic_supported\n"
    "distance to canonical 5'SS, MaxEntScan (site, canonical), phastCons 100/470,\n"
    "intron + 49 nt GC, tissues with PSI > 10 (145 tissues, 26 brain), mESC IRFinder",
    edge=PURPLE, fill="white", fs=7.8)
box(52, 132, 46, 14, "Intron feature table",
    "intron_summary.Rmd\n"
    "intron_summary_table.tsv   192,965 Vast-DB introns\n"
    "class: none 172,373 / decoy 5,430 / cryptic 14,377 / both 785; site IDs and counts per intron,\n"
    "length, tissues with PSI > 10, intron and flanking-exon GC, GC ratio",
    edge=PURPLE, fill="white", fs=7.8)
arrow(48, 139, 52, 139, color=PURPLE, style="<|-|>")

fig.savefig(f"{OUT}.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}.png", dpi=170, bbox_inches="tight")
print("wrote", OUT.with_suffix(".pdf"), "and .png")
