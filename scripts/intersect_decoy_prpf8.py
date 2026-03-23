import argparse
import os
import subprocess
from pathlib import Path

from pybedtools import BedTool


# Example:
# python3 intersect_decoy_prpf8.py --crosslinks data/PRPF8 --decoys data/Decoys/DecoySpliceSites_proteincoding.bed --min-sum 5 --output PRPF8_decoys.filtered.bed
DEFAULT_GENOME = "/home/mikej10/advbfx/reference/genomes/Gencode49/genome.sizes"


def parse_args():
    parser = argparse.ArgumentParser(description="Sum PRPF8 crosslink signal around decoys")
    parser.add_argument("--crosslinks", required=True, help="Path to crosslinks BED directory")
    parser.add_argument("--decoys", required=True, help="Path to decoy BED")
    parser.add_argument(
        "--window",
        type=int,
        default=50,
        help="Half-window size in bp (total added width is +/- this value)",
    )
    parser.add_argument("--min-sum", type=float, default=5, help="Minimum summed signal threshold")
    parser.add_argument("--output", required=True, help="Output BED path")
    parser.add_argument(
        "--genome",
        default=DEFAULT_GENOME,
        help="Genome sizes file used by slop (default: local Gencode49 genome.sizes)",
    )
    return parser.parse_args()


def merge_crosslinks(crosslinks_dir: str) -> BedTool:
    """
    Merge all *.genome.xl.bed files into one BED6 with summed score per locus+strand.
    """
    crosslinks_path = Path(crosslinks_dir)
    xl_files = sorted(crosslinks_path.glob("*genome.xl.bed"))
    if not xl_files:
        raise FileNotFoundError(f"No *.genome.xl.bed files found in {crosslinks_dir}")

    merged_raw = crosslinks_path / "PRPF8_merged_raw.bed"
    merged_fixed = crosslinks_path / "PRPF8_merged.bed"

    # groupby output is: chr start end strand sum; convert to BED6: chr start end . sum strand
    # and normalize chromosome names to chr-prefixed style for decoy compatibility.
    merge_cmd = (
        "cat *.genome.xl.bed | "
        "sort -k1,1 -k2,2n -k6,6 | "
        "bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum "
        f"> \"{merged_raw.name}\""
    )
    fix_cmd = (
        "awk 'BEGIN{OFS=\"\\t\"} "
        "{chrom=$1; if (chrom !~ /^chr/) chrom=\"chr\" chrom; print chrom, $2, $3, \".\", $5, $4}' "
        f"\"{merged_raw.name}\" > \"{merged_fixed.name}\""
    )
    subprocess.run(merge_cmd, cwd=str(crosslinks_path), shell=True, check=True)
    subprocess.run(fix_cmd, cwd=str(crosslinks_path), shell=True, check=True)
    return BedTool(str(merged_fixed))


def decoy_window(decoys: BedTool, window: int, genome: str) -> BedTool:
    return decoys.slop(b=window, g=genome)


def sum_crosslinks_by_window(crosslinks: BedTool, windows: BedTool) -> BedTool:
    # Map summed crosslink score (column 5 in crosslinks BED6) to each decoy window.
    return windows.map(b=crosslinks, c=5, o="sum", null="0")


def filter_by_min_sum(mapped_windows: BedTool, min_sum: float, output_path: str) -> int:
    kept = 0
    with open(output_path, "w", encoding="utf-8") as out:
        for interval in mapped_windows:
            total_signal = float(interval.fields[-1])
            if total_signal >= min_sum:
                out.write(str(interval))
                kept += 1
    return kept


def main():
    args = parse_args()
    if not os.path.exists(args.genome):
        raise FileNotFoundError(f"Genome sizes file not found: {args.genome}")

    decoys = BedTool(args.decoys)
    merged_crosslinks = merge_crosslinks(args.crosslinks)
    windows = decoy_window(decoys, args.window, args.genome)
    mapped = sum_crosslinks_by_window(merged_crosslinks, windows)
    kept = filter_by_min_sum(mapped, args.min_sum, args.output)
    print(f"Wrote {kept} decoy windows with summed signal >= {args.min_sum} to {args.output}")


if __name__ == "__main__":
    main()