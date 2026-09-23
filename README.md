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

## Running Clippy without a container

#### `env/clippy_environment.yml` builds a conda environment that runs Clippy 1.5.0 natively. Every version in it is pinned deliberately - Clippy 1.5.0 predates NumPy 2, pandas 2, plotly 6 and modern bedtools, and an unpinned solve installs all four.

#### `$ conda env create -f env/clippy_environment.yml`
#### `$ conda activate clippy15`
#### `$ pip install -e /path/to/ulelab/clippy --no-deps`

#### `--no-deps` is required; without it pip re-resolves the dependencies and undoes the pinning.

#### The four failure modes the pins prevent, all of which were hit while building this analysis:

| Package | Unpinned result | Symptom |
|---|---|---|
| numpy | >= 2.x | numexpr/pandas ABI error on `import clip` |
| pandas | >= 2.x | peak-calling internals change behaviour |
| **bedtools** | > 2.26 | `bedtools merge ... -c 11,6 ... only has fields 1 - 0` - the broad-peak merge gets an empty file. **Persists even with the Python stack pinned; this was the hardest one to find.** |
| plotly | >= 6.x | interactive mode draws the height threshold as `y = position` instead of a horizontal line |

#### The upstream `environment.yml` shipped with Clippy pins only `bedtools`, `dash`, `dash-bootstrap-components` and `werkzeug`; it leaves numpy, pandas, scipy and pybedtools unbounded and omits plotly entirely, so it no longer produces a working environment.

#### **Interactive mode under WSL2:** Clippy hardcodes `host="127.0.0.1"` in `clip/interaction.py`, which binds only inside the WSL VM and is unreachable from a Windows browser. Changing that call to read host/port from the environment (defaulting to `0.0.0.0`) makes `clippy ... -int` reachable at `http://localhost:8050`. Running with `debug=True` also starts Flask's reloader as a second process, which re-binds the port and produces `OSError: [Errno 98] Address already in use`.

## Line endings

#### `.gitattributes` forces LF on all text files. Knitting the Rmds from **RStudio on Windows against the WSL filesystem** writes CRLF, because Windows R uses `\r\n` for both `write.table` and `data.table::fwrite`. This does not raise an error - the trailing `\r` attaches to the **last column**, so `n_categories` becomes `"1\r"` and every numeric test on it silently fails. `scripts/decoy_exon_overlaps.Rmd` also sets `eol = "\n"` explicitly on its BED output.

## v2 dataset: merged CLIP tracks, Clippy peak calling, and the cross-CLIP intersect

#### This supersedes the original single-track workflow above. The `-b` side of the intersect is now three peak sets rather than a mix of peaks and raw crosslinks, and the CLIP data is merged across samples pulled from Flow.

### Merged CLIP crosslink tracks

#### Samples were selected from Flow, and their crosslink files were downloaded and merged per target.

#### **Per-sample provenance for the merged tracks is in `results/supplementary_merged_clip_inputs.tsv`** (12 rows: Flow sample name and ID, purification target, assay, cell type, condition, source filename, read and crosslink counts, GEO accession). This is the supplementary table for the paper.

| Track | Samples | Assay | Cell lines | Crosslink reads | Unique positions |
|---|---|---|---|---|---|
| `PRPF8_merged.xl.bed` | 4 | eCLIP (ENCODE) | HepG2, K562 | 23,626,946 | 16,052,984 |
| `SmB_merged.xl.bed` | 8 | iCLIP (mild lysis) | HEK293, HepG2, K562 | 16,450,562 | 12,855,363 |

#### Selection criteria: siRNA-treated samples, cell-cycle-phase samples, size-matched inputs, and non-mild-lysis SmB samples were all excluded. SmB is filed on Flow under the gene symbol **SNRPB** - `purification_target=SmB` returns nothing.

#### **Contig naming:** Flow `.genome.xl.bed` files use Ensembl contigs (`1`, `MT`); every other track here is UCSC (`chr1`, `chrM`). They were converted to UCSC names before merging, and cDNA counts at shared positions were summed. Intersecting the two conventions returns zero overlaps *with exit code 0*, so this is silent if missed.

### Clippy peak calling

#### Peaks are called with the pinned biocontainer `quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0`, one chromosome at a time (a genome-wide invocation is OOM-killed at ~5.6 GB RSS). Clippy calls peaks per gene and no gene spans two chromosomes, so chunking is equivalent to a whole-genome run, not an approximation - verified by the chr5 subset of a genome-wide file matching a standalone chr5 run exactly.

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

### Ordering (v3): protein-coding base set first

#### The filtering order was changed so that the protein-coding set is the base dataset, and CLIP support is the LAST split rather than the first. The previous order intersected against CLIP before anything else, which meant every downstream file was conditioned on CLIP support and there was no way to ask which SpliceAI sites in protein-coding introns have *no* support. Exon removal now runs only on the deep-intronic branch, after the canonical split, so the canonical and alternative 5' splice sites still reach the exonic branch.

```
data/Decoys/Splice_All.filtered.05min.bed              794,192  BED6
  |  [1]  scripts/decoy_exon_overlaps.Rmd Part 1   protein_coding keep (GRanges, stranded)
  v
data/Decoys/spliceai_05min_proteincoding.bed                    BED6
  |  [2]  scripts/split_canonical_exonic.sh   vs Wide_canonical_splice_sites.bed
  |-- -v -s --> results/deep_intronic_splice_sites.bed          BED6
  \-- -c -s --> results/proteincoding_splice_sites.bed          BED6 + canon_count
                   \- awk '$NF>0' -> results/exonic_splice_sites.bed
  |  [2b] scripts/decoy_exon_overlaps.Rmd Part 2   exon-overlap removal (unstranded)
  v
data/Decoys/deep_intronic_noexon_splice_sites.bed               BED6
  |  [3]  scripts/intersect_spliceai_support.sh -w 0   RBPnet / PRPF8 / SmB at the site
  |-- hits > 0  --> results/decoys.bed          FINAL decoy set   BED6 + 3 counts + hits
  \-- hits == 0 --> results/cryptic_sites.bed   no support at the site
                       |  [4] scripts/compile_decoy_intron_data.Rmd Parts 1-2
                       |      cryptic -> Vast-DB intron -> canonical 5'SS window
                       |      -> intersect with RBPnet / PRPF8 / SmB peaks
                       v
                     results/cryptics_supported.bed   canonical 5'SS supported by any track
  |  [5]  scripts/SpliceAI_Inference.py   local re-score of decoys.bed + cryptics_supported.bed
  v
  [6]  scripts/compile_decoy_intron_data.Rmd Parts 3-4   local SpliceAI >= 0.1 -> remove HsaALTD donors
       -> results/decoys_final.bed, results/cryptics_supported_final.bed -> feature table
```

### Step 1 - protein-coding base set

#### Part 1 of `scripts/decoy_exon_overlaps.Rmd` keeps SpliceAI inferences (>= 0.05) that overlap a feature of a `gene_type == "protein_coding"` gene, on the same strand, using `reference/gencode.v49.annotation.gtf.gz`. Exonic sites are **retained** here - the exonic/intronic split is step 2.

#### The protein-coding keep is **row-level** (`queryHits`). The previous version selected on gene-symbol membership, which kept a site whenever any other site sharing its gene symbol overlapped; the document prints both counts so the difference is visible.

### Step 2 - split into deep intronic and exonic borders

#### `$ bash scripts/split_canonical_exonic.sh`

#### `Wide_canonical_splice_sites.bed` is BED6 with 192,965 intervals, each exactly 100 bp (+/-50 nt around a canonical splice site). "Exonic" therefore means *within 50 nt of an annotated splice site*, not *anywhere in an exon*.

| Output | How | Contents |
|---|---|---|
| `results/deep_intronic_splice_sites.bed` | `intersect -v -s` | no canonical-window overlap |
| `results/proteincoding_splice_sites.bed` | `intersect -c -s` | full set + canonical-overlap count |
| `results/exonic_splice_sites.bed` | `awk '$NF>0'` | canonical + alternative 5'SS |

#### The script asserts `intronic + exonic == input` and fails loudly otherwise, so a truncated dataset cannot flow downstream unnoticed.

### Step 2b - exon-overlap removal on the deep-intronic branch

#### Knit `scripts/decoy_exon_overlaps.Rmd` again after step 2. Part 2 reloads `results/deep_intronic_splice_sites.bed`, removes every site that overlaps any GENCODE v49 exon on either strand, and writes `data/Decoys/deep_intronic_noexon_splice_sites.bed`. The canonical split only removes sites within 50 nt of a canonical splice site. This step also removes sites inside alternative or internal exons, and inside exons of overlapping genes. The exonic branch is not touched.

#### The exon set includes `retained_intron` transcripts, whose exons span introns, so a site inside an annotated retained intron is removed. The Rmd prints how many sites are removed only because of those transcripts.

#### On the current data, 77,449 of the 367,089 deep-intronic sites are removed, 17,432 of them only because of `retained_intron` transcripts. That leaves **289,640** in `deep_intronic_noexon_splice_sites.bed`. The HNRNPH1 intron-4 decoy is not inside any exon, and the Rmd stops if it is ever removed.

### Step 3 - CLIP support split

#### `$ bash scripts/intersect_spliceai_support.sh -a data/Decoys/deep_intronic_noexon_splice_sites.bed -o results/clip_support -w 0 RBPNET=data/Decoys/rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8=data/CLIP/PRPF8_clippy_w0.5_rollmean80_minHeightAdjust3.0_minPromAdjust3.0_minGeneCount5_Peaks.bed SmB=data/CLIP/SmB_clippy_n40_w0.5_rollmean40_minHeightAdjust8.0_minPromAdjust5.0_minGeneCount5_Peaks.bed`
#### `$ mv results/clip_support_w0.bed results/decoys.bed`
#### `$ mv results/clip_support_w0_nohits.bed results/cryptic_sites.bed`

#### One pass writes **both** branches:
#### - `results/decoys.bed` - supported by at least one category at the site. This is the **final decoy set**.
#### - `results/cryptic_sites.bed` - supported by none. Filtered further in step 4.

#### On the current data, the 289,640 sites split into **19,942 decoys** and **269,698 cryptic sites**. Among the decoys, RBPnet supports 9,526, PRPF8 6,805 and SmB 5,668; 1,865 are supported by two or more.

#### Both carry the same layout (`chr, start, end, gene, spliceai_score, strand, rbpnet, prpf8, smb, hits`), so a tier is selectable at any point with e.g. `awk -F'\t' '$10==3'`. The script asserts `supported + unsupported == input`.

#### **A zero-nt window is used**, because the peak widths (88 nt mean for PRPF8, 45 nt for SmB) already supply the positional tolerance a `-w` window previously provided.

#### The SpliceAI threshold is **0.05**: `data/Decoys/Splice_All.filtered.05min.bed`, 794,192 sites, derived from `Splice_All.filtered.01min.bed` with `awk -F'\t' '$5>=0.05'`. A 0.1 threshold was tested and rejected - the HNRNPH1 target scores 0.073 in the raw whole-intron inference and is lost above 0.08.

### Step 4 - supported cryptics: canonical 5'SS supported by RBPnet, PRPF8 or SmB

#### A cryptic site has no support at the site itself. It is kept only if the canonical 5' splice site of its intron is supported by RBPnet, PRPF8 or SmB - the same three tracks as step 3. Part 1 of `scripts/compile_decoy_intron_data.Rmd` overlaps each cryptic site with the Vast-DB introns in `PSI_TABLE-hg38.tab.gz`, takes each intron's donor (`+`: intron start, `-`: intron end), and matches it to a `Canonical_splice_sites.bed` window. It writes the matched windows to `results/cryptic_canonical_sites.bed`, with the cryptic ID (`GENE_start`, e.g. `HNRNPH1_179620582`) in column 4. After the intersect below, Part 2 keeps the cryptic sites that have at least one supported window, and writes `results/cryptics_supported.bed` in the same 10-column layout as `decoys.bed`.

#### `$ bash scripts/intersect_spliceai_support.sh -a results/cryptic_canonical_sites.bed -o results/canonical_support -w 0 RBPNET=data/Decoys/rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8=data/CLIP/PRPF8_clippy_w0.5_rollmean80_minHeightAdjust3.0_minPromAdjust3.0_minGeneCount5_Peaks.bed SmB=data/CLIP/SmB_clippy_n40_w0.5_rollmean40_minHeightAdjust8.0_minPromAdjust5.0_minGeneCount5_Peaks.bed`
#### `$ mv results/canonical_support_w0.bed results/canonical_supported_sites.bed`

#### On the current data, 269,698 cryptic sites → 267,465 inside a Vast-DB intron → 267,309 with a canonical donor for that intron (276,032 site × window rows) → **191,904 supported cryptics** (72%). By track at the donor: RBPnet 135,874, PRPF8 119,032, SmB 25,196; 66,651 are supported by RBPnet alone.

#### RBPnet hits 134,621 of the 192,965 canonical windows (70%), and PRPF8 or SmB 78,626 (41%), so with all three tracks most cryptics with a canonical donor pass. The selectivity of this step comes mainly from requiring a Vast-DB intron with an annotated canonical donor.

#### **`Canonical_splice_sites.bed` is built by `scripts/CreateCanonicalSpliceSiteBed.sh`** from the central ±5 nt of each `Wide_canonical_splice_sites.bed` window. An earlier build read an `introns.bed` that still carried the 400 nt SpliceAI exonic flanks. That shifted every window 400 nt into the flanking exon (+ strand -400, - strand +400), and PRPF8/SmB support at canonical donors came out at 6% instead of 41%. The Rmd stops if no window covers the HNRNPH1 intron-4 donor (chr5:179,620,891).

### Steps 5-6 - re-score and local-score filter

#### `scripts/SpliceAI_Inference.py` re-scores both sets: it slops 24 nt either side, extracts the strand-aware sequence with `bedtools getfasta`, and runs SpliceAI over the 49 nt window to obtain a **local** donor score. Only column 5 is rewritten, so the support counts survive.

#### `$ conda activate spliceai-env`
#### `$ python3 scripts/SpliceAI_Inference.py --bed results/decoys.bed --fasta <GRCh38.primary_assembly.genome.fa> --genome <genome.sizes> --batch-size 64 --out data/Decoys/decoys_splicescores.bed`
#### `$ python3 scripts/SpliceAI_Inference.py --bed results/cryptics_supported.bed --fasta <GRCh38.primary_assembly.genome.fa> --genome <genome.sizes> --batch-size 64 --out data/Decoys/cryptics_supported_splicescores.bed`

#### It batches sequences and checkpoints to `<out>.scores`; pass `--resume` to continue an interrupted run. Roughly 0.65 s/site on 12 CPU cores.

#### Part 3 of `scripts/compile_decoy_intron_data.Rmd` loads both re-scored files as one table with a `site_class` column (`decoy` / `cryptic_supported`) and applies the **local SpliceAI >= 0.1** filter (`MIN_LOCAL_SPLICEAI`) before the intron overlap. This threshold applies to the local 49 nt score and is **not** comparable to the 0.05 used upstream on the whole-intron inference - they are different measurements.

#### Part 3 then removes **annotated alternative 5' splice sites**: any site whose position (BED start + 1) and strand match a donor of a Vast-DB `HsaALTD` (Alt5) event in `PSI_TABLE-hg38.tab.gz`. Every donor is parsed from the event's `FullCO` (133,342 unique donors). The match is exact: on `exonic_splice_sites.bed`, 45,088 sites match at offset 0, against at most 730 at ±1–2 nt. Applied before the local-score filter, it would remove 1,543 of 19,942 decoys and 3,989 of 125,253 supported cryptics. The final sets are written to `results/decoys_final.bed` and `results/cryptics_supported_final.bed`.

#### Parts 3-4 compute every feature for both classes in one pass: the phastCons, MaxEntScan and GC (`scripts/extract_gc_content.sh`: intron GC and GC of the 49 nt re-score window) commands run once on BEDs that hold decoys and supported cryptics together, and Part 4 merges the results back by ID. The final feature table `results/decoy_intron_features_final.tsv` therefore has identical columns for both classes, with `site_class` telling them apart, so they can be compared as two groups on tissue counts with PSI > 10, intron length, GC, local SpliceAI, MaxEnt and phastCons. `results/decoy_vs_cryptic_class_summary.tsv` gives n, median and IQR of each feature per class, and `figures/decoy_vs_cryptic_class_comparison.pdf` shows the distributions side by side. The existing decoy-only figures are unchanged; two figures plot the supported cryptics' phastCons 100-way and 470-way against the number of tissues with PSI > 10.

### Known gap: the alternative 5' splice site

#### The exonic branch is intended to hold both the canonical and the alternative 5' splice sites, but the HNRNPH1 alternative 5'SS at **chr5:179,623,595** (VastDB `HsaALTD0003092-2`) is currently captured by neither branch: it is absent from `Splice_All.filtered.05min.bed` (nearest site >= 0.05 is 298 nt away, at 179,623,297) and no `Wide_canonical_splice_sites.bed` window covers it. Its SpliceAI score is presumably below 0.05; this cannot be confirmed from the repo because `Splice_All.filtered.01min.bed` is not checked in. Resolving this needs either a lower inference threshold or an annotation source that includes annotated alternative donors.

## Decoy Feature Table Generation

#### The re-scored decoys (`data/Decoys/decoys_splicescores.bed`) and supported cryptics (`data/Decoys/cryptics_supported_splicescores.bed`) are loaded as one table, with a `site_class` column, in Part 3 of `compile_decoy_intron_data.Rmd` to integrate intron retention quantification data from `PSI-TABLE-hg38.tab.gz`. The Rmd file uses Genomic Ranges to integrate intron coordinates, unique identifiers `EVENT` and intron retention PSI values in 145 cell and tissue types with SpliceAI inference. Decoy distance from canonical splice site is calculated with strandwise logic. After overlapping the decoy database with introns, MaxEntScan is used to calculate the strength of decoy predicted splice sites and the canonical 5' splice site for the intron harboring the decoy with the scripts `run_maxentscan_decoy.sh` `run_maxentscan_canonical.sh`.Average phastCons 100-way and 470-way scoring across the intron harboring the decoy is calculated with `extract_phastcons_scores.sh`. Part 4 of the R markdown file reloads the results from MaxEntScan and phastCons and merges into the final feature table.

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

## Decoys vs supported cryptics: exploration notebook

#### `notebooks/decoy_vs_cryptic_features.ipynb` (kernel `rbpnet-env`; needs pandas, statsmodels, seaborn) loads `results/decoy_intron_features_final.tsv` and compares the two classes with proportions and effect sizes rather than raw counts, since there are 2.7x more cryptics than decoys. HNRNPH1 is marked on every plot. Figures go to `figures/notebook/`.

#### - **A.** Retention propensity: proportion of sites with PSI > 10 in at least k tissues against k, with a logistic fit on log(1+k) per class.
#### - **B.** Tissue count against local SpliceAI, MaxEnt and intron GC, with negative-binomial regression per class (the count is over-dispersed: variance ~25x the mean).
#### - **C.** GC content (intron, 49 nt site window) against tissue count as 2D histograms normalised to percent of each class, plus the per-bin difference.
#### - **D.** Statistics. ANOVA is avoided: with two groups it is a t-test, the features are skewed or bounded, and at n = 6.5k vs 17.7k every difference is "significant". Instead: Mann-Whitney with Cliff's delta and BH correction per feature (`figures/notebook/D_effect_sizes.tsv`), a logistic regression of class on standardised features (with VIFs), and a negative-binomial regression of tissue count on class adjusted for the features that differ between classes.
#### - **E.** ECDFs, site competitiveness (MaxEnt site minus canonical) and intron position against retention, decoys by support tier, and Spearman correlation heatmaps per class.

#### Regenerate with `jupyter nbconvert --to notebook --execute --inplace notebooks/decoy_vs_cryptic_features.ipynb` from the repo root.

## HNRNPH1 cross-species SpliceAI scores (`scripts/maf_spliceai_single_locus.py`)

#### Scores the HNRNPH1 intron-4 decoy locus in every species of the UCSC 100-way multiz alignment. The script trims `reference/chr5.maf.gz` to the region with kent `mafsInRegion`, converts the alignment to FASTA with PHAST `msa_view` (`-V` reverse-complements, because HNRNPH1 is on the minus strand), and runs SpliceAI on each species' ungapped sequence. It needs `spliceai-env`, plus `mafsInRegion` and `msa_view` on `PATH` (or pass `--mafs-in-region` / `--msa-view`).

#### `$ python3 scripts/maf_spliceai_single_locus.py --maf reference/chr5.maf.gz --out results/hnrnph1_179620582_maf_spliceai.tsv`

#### The default region is `chr5:179620038-179620941` (1-based, inclusive). The output has one row per species (100 rows). `max1`/`pos1` is the canonical 5' splice-site peak and its position in the window. `max2`/`pos2` is the strongest donor peak further into the intron; the decoy is at `pos2 = 358` in hg38. `d` is the intron-length divergence from hg38, `abs(n - v) / v`. `scripts/hnrnph1decoyposition.Rmd` plots `max2` against `pos2` from this file.
