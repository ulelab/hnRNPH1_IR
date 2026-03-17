import argparse
import pybedtools
from pybedtools import BedTool

#usage: python3 intersect_decoy_prpf8.py --crosslink /path/to/crosslink.bed --decoy decoy.bed --min-sum 25 --output PRPF8decoys.bed

def parse_args():
    parser = argparse.ArgumentParser(description="Sum RBP signal around decoys")
    parser.add_argument("--crosslinks", required=True, help="Path to crosslinks BED")
    parser.add_argument("--decoys", required=True, help="Path to decoy BED")
    parser.add_argument("--window", type=int, default=50,  help="Half-window size (total window is 2x this value)")
    parser.add_argument("--min-sum", type=float, default=5, help="Minimum signal threshold")
    parser.add_argument("--output", required=True, help="Output BED path")
    parser.add_argument("--signal-col", default=5, help="1-based index for the crosslink column signal")

    return parser.parse_args()

genome = "/scratch/prj/ppn_rnp_networks/shared/references/genomes/homo_sapiens/GRCh38-EnsembleRelease109/Homo_sapiens.GRCh38.109.gtf"

def decoy_window(decoys: BedTool, window, genome): 
    return decoys.slop(b=window, genome=genome) 


def intersection(crosslinks, windows):
    Bedtool.intersect()
    crosslinks_with_windows = crosslinks.intersect(windows, u=True)
    crosslinks_with_windows.head()

def main():
    args = parse_args()
    windows = decoy_window(decoys, args.window, genome)
    PRPF8_decoys = intersection(crosslinks, windows)
    print("{PRPF8_decoys}") 