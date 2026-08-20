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

#### **Applies to the v1 dataset only.** The v2 workflow below does not currently include a Recount filter; the sections here are retained because they document how `Decoys_proteincoding_recount_filtered.bed` was produced and remain re-runnable. Note also that `query_junctions.py` filters on `annotated = 1`, so it only removes sites matching *annotated* junctions - unannotated junctions with real read support are not considered.

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

## v2 dataset: merged CLIP tracks, Clippy peak calling, and the cross-CLIP intersect

#### This supersedes the original single-track workflow above. The `-b` side of the intersect is now three peak sets rather than a mix of peaks and raw crosslinks, and the CLIP data is merged across samples pulled from Flow.

### Merged CLIP crosslink tracks

#### Samples were selected from Flow with `scripts/flow_sample_inventory.py` (writes `data/flow_prpf8_snrpb_samples.tsv`; the retained subset is `data/flow_prpf8_snrpb_samples_selected.tsv`, exclusions and their reasons in `..._excluded.tsv`). Crosslinks were downloaded with `scripts/flow_fetch_crosslinks.py` and merged per target with `scripts/merge_crosslinks.sh`.

#### **Per-sample provenance for the merged tracks is in `results/supplementary_merged_clip_inputs.tsv`** (12 rows: Flow sample name and ID, purification target, assay, cell type, condition, source filename, read and crosslink counts, GEO accession). This is the supplementary table for the paper.

| Track | Samples | Assay | Cell lines | Crosslink reads | Unique positions |
|---|---|---|---|---|---|
| `PRPF8_merged.xl.bed` | 4 | eCLIP (ENCODE) | HepG2, K562 | 23,626,946 | 16,052,984 |
| `SmB_merged.xl.bed` | 8 | iCLIP (mild lysis) | HEK293, HepG2, K562 | 16,450,562 | 12,855,363 |

#### Selection criteria: siRNA-treated samples, cell-cycle-phase samples, size-matched inputs, and non-mild-lysis SmB samples were all excluded. SmB is filed on Flow under the gene symbol **SNRPB** - `purification_target=SmB` returns nothing.

#### **Contig naming:** Flow `.genome.xl.bed` files use Ensembl contigs (`1`, `MT`); every other track here is UCSC (`chr1`, `chrM`). `merge_crosslinks.sh` converts them and sums cDNA counts at shared positions. Intersecting the two conventions returns zero overlaps *with exit code 0*, so this is silent if missed.

### Clippy peak calling

#### Peaks are called with the pinned biocontainer `quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0` via `scripts/run_clippy_chunked.sh`, one chromosome at a time (a genome-wide invocation is OOM-killed at ~5.6 GB RSS). Clippy calls peaks per gene and no gene spans two chromosomes, so chunking is equivalent to a whole-genome run, not an approximation - verified by the chr5 subset of a genome-wide file matching a standalone chr5 run exactly.

#### Parameters were tuned separately per target against the crosslink bigWigs in Clippy's interactive mode, because the two tracks differ in assay and depth:

| Track | Parameters | Peaks | Mean width | Median width |
|---|---|---|---|---|
| PRPF8 | `-n 80 -w 0.5 -x 3.0 -mx 3.0 -mg 5` | 227,651 | 88.2 nt | 82 nt |
| SmB | `-n 40 -w 0.5 -x 5.0 -mx 8.0 -mg 5` | 152,983 | 44.8 nt | 42 nt |

#### `-n` rolling-mean window, `-w` width, `-x` min prominence adjust, `-mx` min height adjust, `-mg` min gene counts. Both thresholds are multiples of each gene's mean smoothed coverage, so they normalise to local coverage rather than applying an absolute floor.

#### Two behaviours worth knowing when re-tuning:
#### - **`-mx` has no effect unless it exceeds `-x`.** Prominence can never exceed a peak's absolute height, so any peak passing the prominence test automatically passes an equal-or-lower height test.
#### - **`--width` is not encoded in Clippy's output filename** (only rollmean, minHeightAdjust, minPromAdjust, minGeneCount). Runs differing only in width overwrite each other unless the width is put in the output prefix, which is why the peak files carry `w0.5`.

#### PRPF8 peaks are ~2x the width of SmB peaks, a direct consequence of the 80 vs 40 nt window. Raw support counts are therefore not directly comparable between the two categories - PRPF8 presents roughly twice the genomic target per peak.

#### Crosslink bigWigs for IGV are generated by `scripts/xl_to_bedgraph.sh` (stranded, minus-strand values negative).

### Cross-CLIP support intersect

#### `scripts/intersect_spliceai_support.sh` retains a SpliceAI inference if it falls within `-w` nt of a feature in **any** supplied category, on the same strand, and reports per-category counts alongside the union.

#### `$ scripts/intersect_spliceai_support.sh -a data/Decoys/Splice_All.filtered.05min.bed -o results/supported_splice_sites -w 0 RBPNET=data/Decoys/rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8=data/CLIP/PRPF8_clippy_w0.5_rollmean80_minHeightAdjust3.0_minPromAdjust3.0_minGeneCount5_Peaks.bed SmB=data/CLIP/SmB_clippy_n40_w0.5_rollmean40_minHeightAdjust8.0_minPromAdjust5.0_minGeneCount5_Peaks.bed`

#### The SpliceAI inference threshold is now **0.05** (not 0.01): `data/Decoys/Splice_All.filtered.05min.bed`, 794,192 sites, derived from `Splice_All.filtered.01min.bed` with `awk -F'\t' '$5>=0.05'`. A 0.1 threshold was tested and rejected - the HNRNPH1 target scores 0.073 in the raw whole-intron inference and is lost above 0.08.

#### **A zero-nt window is used**, because the peak widths above (88 and 45 nt mean) already supply the positional tolerance that a `-w` window previously provided.

| Category | Sites supported |
|---|---|
| RBPnet | 250,963 |
| PRPF8 | 136,732 |
| SmB | 28,621 |
| **Union** | **302,025** of 794,192 (38.0%) |
| Supported by >= 2 categories | 103,889 |

### Removing canonical splice sites

#### The script writes `results/supported_splice_sites_w0.bed`; that file is kept as **`results/supported_splice_sites.bed`** (the output of the first intersect) and canonical sites are then removed from it:

#### `$ bedtools intersect -a results/supported_splice_sites.bed -b data/Decoys/Wide_Canonical_splice_sites.bed -v -s > results/noncanonical_predicted_splice_sites.bed`

#### **302,025 -> 36,240** sites (265,785 removed as canonical). Of the survivors, 349 are supported by all three categories, 3,233 by exactly two, and 32,658 by one. `HNRNPH1_179620582` is retained with support from all three.

#### Output `results/noncanonical_predicted_splice_sites.bed` is 10 columns: `chr, start, end, gene, spliceai_score, strand, rbpnet, prpf8, smb, n_categories`. The three count columns and `n_categories` are carried through every downstream step, so a tier can be selected at any point with e.g. `awk -F'\t' '$10==3'`.

### Downstream: exon removal, re-scoring, local-score filter

#### 1. **`scripts/decoy_exon_overlaps.Rmd`** - input `results/noncanonical_predicted_splice_sites.bed`. Removes exon overlaps and keeps only protein-coding introns using `reference/gencode.v49.annotation.gtf.gz` (GenomicRanges). Output `data/Decoys/noncanonical_proteincoding.bed`.

#### 2. **`scripts/SpliceAI_Inference.py`** - re-scores each surviving site. Slops 24 nt either side, extracts the strand-aware sequence with `bedtools getfasta`, and runs SpliceAI over the 49 nt window to obtain a **local** donor score. Only column 5 is rewritten; all other columns pass through, so the support counts survive.

#### `$ python3 scripts/SpliceAI_Inference.py --bed data/Decoys/noncanonical_proteincoding.bed --fasta ../../../reference/genomes/Gencode49/GRCh38.primary_assembly.genome.fa --genome ../../../reference/genomes/Gencode49/genome.sizes --out data/Decoys/noncanonical_proteincoding_splicescores.bed`

#### 3. **`scripts/compile_decoy_intron_data.Rmd`** - loads the re-scored BED and applies a **local SpliceAI score filter of >= 0.1** (`MIN_LOCAL_SPLICEAI` in the `paths-v3` chunk) before the intron overlap. This threshold is applied to the local 49 nt score and is **not** comparable to the 0.05 threshold used upstream on the whole-intron inference - the two are different measurements.

## Decoy Feature Table Generation

#### The database of predicted splice sites in 'protein-coding' introns is uploaded as `data/Decoys/noncanonical_proteincoding_splicescores.bed` - the re-scored, local-SpliceAI-filtered set produced by the v2 workflow above - in the `compile_decoy_intron_data.Rmd` to integrate intron retention quantification data from `PSI-TABLE-hg38.tab.gz`. The Rmd file uses Genomic Ranges to integrate intron coordinates, unique identifiers `EVENT` and intron retention PSI values in 145 cell and tissue types with SpliceAI inference. Decoy distance from canonical splice site is calculated with strandwise logic. After overlapping the decoy database with introns, MaxEntScan is used to calculate the strength of decoy predicted splice sites and the canonical 5' splice site for the intron harboring the decoy with the scripts `run_maxentscan_decoy.sh` `run_maxentscan_canonical.sh`.Average phastCons 100-way and 470-way scoring across the intron harboring the decoy is calculated with `extract_phastcons_scores.sh`. Part 2 of the R markdown file reloads the results from MaxEntScan and phastCons and merges into the final feature table.

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
