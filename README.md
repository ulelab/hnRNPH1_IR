## hnRNPH1 Intron 4 Decoy Analysis

#### This repository contains data, scripts, and outputs for the hnRNPH1 intron 4 retention decoy analysis workflow.

## Decoy Dataset Generation

#### The database of predicted decoy loci was generated as a BED file using the bedtools intersect workflow described in `decoy_splice_site_flowchart.pdf`

#### Briefly, all human intron coordinates were collected from Vast-DB `PSI-TABLE-hg38.tab.gz` as any EVENT with an ID beginning with 'HsaIN'. These coordinates were fed in a strandwise fashion to bedtools getfasta, then SpliceAI to predict splice donor scores for every intronic nucleotide, using 400 nt exonic flanks for internal normalization by the canonical 5' splice site. Scores below 0.01 were filtered out, then all coordinates were intersected against PRPF8-binding datasets of RBPnet predictions generated with the same scaled workflow of all introns, and eCLIP and iCLIP signal from PRPF8 experimental data. Predictions within 100nt of a canonical splice site were removed from the dataset. Finally this dataset `DecoySpliceSites_f50.bed` was filtered to remove exonic loci and keep only 'protein-coding' introns by overlapping with `Gencode.v49.annotation.gtf` using GenomicRanges in the `decoy_exon_overlaps.Rmd` script to output `DecoySpliceSites_proteincoding.bed`.

### bedtools intersect commands

#### `$ bedtools intersect -a Splice_All.filtered.01min.bed -b rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8_iCLIP_HepG2_xlinks.bed PRPF8_eCLIP_HepG2_rollmean10_minHeightAdjust1. 0_minPromAdjust1.0_minGeneCount5_Peaks.bed -c -s | awk '$NF > 0' > SupportedSpliceSites_f50.bed`

#### `$ bedtools intersect -a SupportedSpliceSites_f50.bed -b Wide_Canonical_splice_sites.bed -v -s > DecoySpliceSites_f50.bed`

## Calculate peak SpliceAI score around decoy site

#### The final `DecoySpliceSites_proteincoding.bed` dataset was used as input for `SpliceAI_Inference.py` which slops the coordinates 24 nt wider on either side of the proposed decoy site, then uses strand-aware extraction of the fasta sequence with bedtools getfasta. Nucleotide sequences of the 49 nt region surrounding the decoy site is passed to SpliceAI. Max splice donor score is recorded for the locus and an updated BED file is output. 

### `SpliceAI_Inference.py` args

#### `$ python3 ../scripts/SpliceAI_Inference.py --bed Decoys/DecoySpliceSites_proteincoding.bed --fasta ../../../reference/genomes/Gencode49/GRCh38.primary_assembly.genome.fa --genome ../../../reference/genomes/Gencode49/genome.sizes --out Decoys_proteincoding_splicescores.bed`

## Recount Validation

#### Recount validation is performed with `scripts/query_junctions.py` using inference loci from `Decoys_proteincoding_splicescores.bed` (7-column BED). The script loads loci into an attached in-memory SQLite schema `sql/locus_recount_schema.sql` and joins them to the Recount intron table from `junctions.sqlite` (run separately for TCGA, SRA, and GTEx junction databases) to generate:
#### - `results/tcgajunctions.tsv`
#### - `results/srajunctions.tsv`
#### - `results/gtexjunctions.tsv`

#### - Defines the temporary input table structure used for overlap queries (`chrom`, `start`, `end`, `gene`, `splice_score`, `strand`, `flag`).
#### - Enables consistent coordinate loading before intersecting against Recount intron junction records.

### Script: `scripts/query_junctions.py`
#### - Connects to Recount SQLite in read-only mode, loads inference loci into an attached in-memory database, matches loci to introns on `chrom` + `strand`, with boundary proximity (`flank_bp`) on intron start/end and filters to reference-supported records (`annotated = 1`). Finally writes full joined outputs as TSV (`--out-tsv`) and optional terminal previews.

### Added/returned columns in output TSVs include:
#### - Locus fields: `locus_id`, `gene`, `chrom`, `locus_start`, `locus_end`, `locus_strand`, `splice_score`
#### - Junction identity/coordinates: `snaptron_id`, `intron_start`, `intron_end`
#### - Junction support metrics: `samples_count`, `coverage_sum`, `coverage_avg`, `coverage_median`
#### - Annotation/source metadata: `left_annotated`, `right_annotated`, `source_dataset_id`, `annotated`

#### Unique IDs added back to *junctions.tsv files with `sed -Ei 's/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/\2_\4\t\2\t\3\t\4/' tcgajunctions.tsv`

### Defined 3' or 5' splice sites from Recount by determining if the loci of interest match the junction intron start or end using `annotate_recount.sh`
#### if start_match:
#### strand + -> junction_type = 5ss
#### strand - -> junction_type = 3ss
#### if end_match:
#### strand + -> junction_type = 3ss
#### strand - -> junction_type = 5ss

### Filtering decoys to the Recount-unique set

#### The three `*junctions.tsv` files above are the *matches* - decoy loci that coincide (within `--flank-bp`, default 5) with an **annotated** junction boundary in Recount. A predicted cryptic site that already appears as an annotated junction is a known, used splice site, not a decoy, so these are removed to leave the final unique set.

#### `scripts/filter_decoys_by_recount.py` performs this subtraction. A locus is keyed as `<gene>_<bed_start>` (matching the `locus_id` column) and is dropped if it appears in **any** of the three databases.

#### `$ python3 scripts/filter_decoys_by_recount.py --loci-bed data/Decoys/Decoys_proteincoding_splicescores.bed --junctions results/tcgajunctions.tsv results/srajunctions.tsv results/gtexjunctions.tsv --out-bed data/Decoys/Decoys_proteincoding_recount_filtered.bed --out-removed results/decoys_removed_by_recount.tsv`

#### Result: **15,172 input loci -> 1,259 removed -> 13,913 retained** as `data/Decoys/Decoys_proteincoding_recount_filtered.bed`. The removed loci, with the databases that supported each, are listed in `results/decoys_removed_by_recount.tsv`.

#### The script aborts if any locus in a junction TSV is absent from the BED or has mismatched coordinates, so a stale junction file cannot silently produce a wrong filter (override with `--no-strict`).

#### **Note: this step is not yet drawn in `figures/decoy_splice_site_flowchart.pdf`, which currently terminates at the 15,172-site box.** The flowchart needs a final node: *"Remove loci matching annotated Recount junctions (TCGA/SRA/GTEx, +/-5 bp) -> 13,913 Predicted Cryptic Splice Sites"*.

#### QC of the subtraction:
#### - 1,230 of the 1,259 removed loci (97.7%) were found in all three databases, 23 in two, 6 in one - the overlap is highly consistent, as expected for genuinely annotated junctions.
#### - Removed loci carry a modestly higher SpliceAI score than retained ones (median 0.101 vs 0.072), consistent with them being real splice sites, though the separation is not large.
#### - Removals occur on all 24 chromosomes at broadly similar rates (4-14%), confirming the query results are genome-wide and not truncated.
#### - `HNRNPH1_179620582`, the locus used for the conservation figure, is retained.

## Decoy Feature Table Generation

#### The database of decoys in 'protein-coding' introns is uploaded as `Decoys_proteincoding_recount_filtered.bed` (the Recount-filtered 13,913-locus set; set `splicescores_bed` back to `Decoys_proteincoding_splicescores.bed` in the `paths-v3` chunk to use the unfiltered 15,172) in the `compile_decoy_intron_data.Rmd` to integrate intron retention quantification data from `PSI-TABLE-hg38.tab.gz`. The Rmd file uses Genomic Ranges to integrate intron coordinates, unique identifiers `EVENT` and intron retention PSI values in 145 cell and tissue types with SpliceAI inference. Decoy distance from canonical splice site is calculated with strandwise logic. After overlapping the decoy database with introns, MaxEntScan is used to calculate the strength of decoy predicted splice sites and the canonical 5' splice site for the intron harboring the decoy with the scripts `run_maxentscan_decoy.sh` `run_maxentscan_canonical.sh`.Average phastCons 100-way and 470-way scoring across the intron harboring the decoy is calculated with `extract_phastcons_scores.sh`. Part 2 of the R markdown file reloads the results from MaxEntScan and phastCons and merges into the final feature table.

##### Ten predicted decoys are dropped in feature table generation. Dropped decoyIDs:
##### "TNNI2_1839212" "SLC7A6_68264187" "OAZ1_2270281" "ITPA_3221726" "ARHGAP40_38626901" "DHX35_38962112"        
##### "PLCG1_41162940" "SS18L1_62163387" "RP4-583P15.14_63738551" "BHLHB9_102745917"
##### These loci are largely within protein-coding introns, with the exception of SLC7A6_68264187 and DHX35_38962112 which appear in the 5' UTR/intergenic space. 6 out of the 10 missing loci are in chromosome 20.

## Figure 1 R Markdown (`scripts/hnRNPH1_figure1.rmd`)

#### This R Markdown document generates Figure 1: a comparative multiple-sequence alignment view of the hnRNPH1 intron 4 decoy region (`chr5:179620560-179620600`, hg38). It starts from an extracted UCSC multiz MAF alignment (converted to FASTA), plots the raw alignment, then creates a manuscript-ready alignment by renaming taxa, converting DNA bases from `T` to `U`, removing selected outlier species, and dropping columns that are gaps/missing across all taxa.

The final plot highlights the proposed decoy site (positions 27-33 in the processed alignment; labeled as genomic interval `179,620,576-179,620,582`) and is intended for direct use in manuscript figure generation.

### Inputs used by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.fa`
- `data/tree_to_clade_mapping.tsv`

### Output written by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.processed.fa`


## Cross-species SpliceAI scan across whole-intron MAF alignments

#### This is the workflow behind `results/all_spliceai_scores_v3.1.tsv`, the table used for the HNRNPH1 positional-shift figure. It asks, for every decoy locus and every species in the UCSC 100-way alignment, where the strongest SpliceAI donor signal sits inside the intron relative to the canonical 5' splice site.

#### **Full provenance, including which script version produced the table and the evidence for it, is in [`PROVENANCE.md`](PROVENANCE.md). Read that before citing this pipeline in a manuscript.**

### Pipeline

#### 1. `scripts/split_chromosomes.sh` splits the decoy BED into ~50-locus parts, producing the `chrparts.txt` manifest (109 parts).
#### 2. `scripts/run_mafsInRegion_array.sbatch` (SLURM array 1-109) runs kent `mafsInRegion` to cut a per-locus MAF out of each chromosome's multiz alignment. `scripts/run_ctrl_mafsInRegion_array.sbatch` does the same for the matched control set.
#### 3. `scripts/process_maf_spliceai_v2.5.py`, driven by `scripts/process_maf_spliceai_array.sbatch`, converts each MAF to FASTA with `msa_view` (`-V` to reverse-complement on the minus strand), then runs the 5-model SpliceAI ensemble over each species' ungapped sequence. For every species it records the canonical donor peak (`max1`/`pos1`) and the strongest downstream cryptic peak (`max2`/`pos2`), searching from `pos1 + 20` onward.
#### 4. `scripts/combine_spliceai_tsv.py` concatenates the 109 per-part TSVs into the single wide table.

### Output format

#### One row per decoy locus, one column per species (100 species from `data/tree_to_clade_mapping.tsv`). Each cell is `[max1][pos1][max2][pos2][d]`:

#### - `max1` / `pos1` - peak SpliceAI donor score in the canonical 5' splice-site search window, and its position in nt from the intron start
#### - `max2` / `pos2` - peak donor score downstream of `pos1 + 20` (the proposed cryptic regulatory element), and its position
#### - `d` - normalized intron-length divergence from hg38, `d = abs(n - v) / v`, where `v` is the hg38 ungapped intron length and `n` is this species'. This is a length ratio, **not** a phylogenetic distance; it is what widens the `max1` search window for species with longer introns.
#### - `-1` in any field means the score could not be computed: species absent from the alignment block, empty or invalid MAF, or an all-gap sequence

### The v3 re-scoring pass

#### `scripts/process_maf_spliceai_v3.py` is a **later, separate** pass that takes an existing v3.1 TSV via `--v31_tsv`, re-extracts a local 50-nt window around each recorded `pos1`/`pos2`, and rewrites only `max1`/`max2` while copying `pos1`, `pos2` and `d` through unchanged. It did not produce `all_spliceai_scores_v3.1.tsv` - it consumes it. Its driver `scripts/process_maf_spliceai_array_v3.sbatch` still carries an unfilled `V31_TSV` placeholder and no output of this pass exists in the repository. It is kept for the record; see `PROVENANCE.md`.

### Single-locus streaming

#### `scripts/stream_maf_region.py` and `scripts/test_stream_maf_hnrnph1.sh` pull the alignment for one region without going through the full array, which is how the HNRNPH1 intron 4 locus was inspected interactively.

## Positional-shift figure (`scripts/hnrnph1decoyposition.Rmd`)

#### Reads `results/all_spliceai_scores_v3.1.tsv` and `data/tree_to_clade_mapping.tsv`, parses the bracketed cells, filters to `decoyID == "HNRNPH1_179620582"`, and plots `pos2` against `max2` for each species.

#### Species are collapsed into four plotted groups: **placental mammals** (Afrotheria, Euarchontoglires, Laurasiatheria, Primate, plus Armadillo, which the mapping file files under the catch-all `Mammal` clade despite being a xenarthran), **other mammals** (opossum, Tasmanian devil, wallaby, platypus), **Aves**, and **Fish**. Sarcopterygii is dropped from this figure. Nine species are labelled: cat, pig, human, rat, cow, gorilla, chimp, baboon, mouse.

#### Group is encoded redundantly by both fill colour and point shape. The palette was checked for colour-vision deficiency across all pairs (OKLab dE under Machado 2009 severity-1.0 simulation); the mammal red/orange pair and the red/green pair sit in warning bands, so the shape scale is load-bearing and must not be removed.


### Recount filtering applied upstream of the figure

#### The `recount-filter` chunk restricts the score table to the Recount-filtered decoy set before parsing. Of the 4,886 decoyIDs in `all_spliceai_scores_v3.1.tsv`, 271 match an annotated Recount junction and are removed.

#### A further 309 decoyIDs are **not present in `Decoys_proteincoding_splicescores.bed` at all** - the MSA/SpliceAI scan was run on `DecoysInConservedGenes50thpercentile.bed`, a different decoy list, so these were never queried against Recount and cannot be judged by this filter. The `DROP_UNASSESSED` flag at the top of the chunk decides their fate:

#### - `TRUE` (default) - treat the filtered BED strictly as a whitelist, leaving **4,306** decoyIDs
#### - `FALSE` - remove only the 271 Recount-flagged decoys and keep the 309 unassessed ones, leaving **4,615**

#### Both counts are printed when the chunk runs, so the 309 are never dropped silently. The chunk aborts if `HNRNPH1_179620582` does not survive, since the figure would otherwise be empty.

### Outputs

- `data/Decoys/hnrnph1_spliceaispecies.tsv`
- `figures/hnrnph1_max2_pos2_scatterplot.pdf`
- `figures/hnrnph1_max2_pos2_scatterplot.png`
