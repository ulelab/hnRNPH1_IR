#!/usr/bin/env python3
"""BEDs of every Vast-DB alternative 5' splice site donor (HsaALTD events), for query_recount_junctions.py.

Donors are read from FullCO, the same source as the alternative-5'SS filter in compile_decoy_intron_data.Rmd,
so the BEDs hold exactly the donors that filter removes sites at. FullCO lists an event's donors in option
order (+ strand: chr:acceptor-D1+D2..,next; - strand: chr:D1+D2..-acceptor,next); the first listed donor is
the event's reference option. Some events list more donors than they have PSI rows, so the rows alone would
miss donors. Donor coordinates are the last exonic base, the convention of the site BEDs (BED end = donor).

Outputs (10-column, the layout of decoys_final.bed, so the query and downstream code read them unchanged):
  results/altd_donors_reference.bed     first-listed donor of each event
  results/altd_donors_alternative.bed   every other listed donor
Column 4 is GENE|EVENT-k/n, so the site ID written by the query is GENE|EVENT-k/n_start.
A donor that is the reference of any event is written only to the reference BED.

Usage (repository root):  python scripts/make_altd_donor_beds.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

psi = pd.read_csv(ROOT / "reference" / "PSI_TABLE-hg38.tab.gz", sep="\t", usecols=["GENE", "EVENT", "LENGTH", "FullCO"])
a = psi[psi.EVENT.str.startswith("HsaALTD")].copy()
a["event"] = a.EVENT.str.extract(r"^(HsaALTD\d+)-\d+/\d+$")[0]

rows = []
for (event, fullco), d in a.groupby(["event", "FullCO"], sort=False):
    chrom, rest = fullco.split(":", 1)
    left, right = rest.split(",", 1)[0].split("-", 1)
    minus = "+" in left
    donors = [int(x) for x in (left if minus else right).split("+")]
    gene = d.GENE.iloc[0]
    for k, pos in enumerate(donors, start=1):
        rows.append((chrom, pos, "-" if minus else "+", f"{gene}|{event}-{k}/{len(donors)}", k == 1, event, k))
don = pd.DataFrame(rows, columns=["chr", "donor", "strand", "name", "is_ref", "event", "k"])

# Check the first-listed = reference convention against rows that carry LENGTH.
a["k"] = a.EVENT.str.extract(r"-(\d+)/\d+$")[0].astype(int)
agree = ((a.k == 1) == (a.LENGTH == 0)).mean()
print(f"rows where (option 1) == (LENGTH 0): {100 * agree:.2f}%")

pos_key = ["chr", "donor", "strand"]
ref_pos = don[don.is_ref].drop_duplicates(pos_key)
alt_pos = don[~don.is_ref].drop_duplicates(pos_key).merge(ref_pos[pos_key], how="left", indicator=True)
n_shared = int((alt_pos["_merge"] == "both").sum())
alt_pos = alt_pos[alt_pos["_merge"] == "left_only"].drop(columns="_merge")

layout = ["chr", "start", "donor", "name", "score", "strand", "rbpnet", "prpf8", "smb", "n_categories"]
for df, path in [(ref_pos, "altd_donors_reference.bed"), (alt_pos, "altd_donors_alternative.bed")]:
    (df.assign(start=df.donor - 1, score=0, rbpnet=0, prpf8=0, smb=0, n_categories=0)
       .sort_values(["chr", "start"])[layout]
       .to_csv(ROOT / "results" / path, sep="\t", header=False, index=False, lineterminator="\n"))

print(f"HsaALTD rows: {len(a):,}  events: {a.event.nunique():,}  genes: {a.GENE.nunique():,}")
print(f"donors listed in FullCO: {len(don):,}  unique positions: {don[pos_key].drop_duplicates().shape[0]:,}")
print(f"reference donors:   {len(ref_pos):,} -> results/altd_donors_reference.bed")
print(f"alternative donors: {len(alt_pos):,} -> results/altd_donors_alternative.bed "
      f"({n_shared:,} that are also another event's reference kept in the reference BED only)")
