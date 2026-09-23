## hnRNPH1 Intron 4 Decoy Analysis

#### This repository contains data, scripts, and outputs for the hnRNPH1 intron 4 retention decoy analysis workflow.

## Decoy Dataset Generation

#### The database of predicted decoy loci was generated as a BED file using the workflow described in `figures/decoy_splice_site_flowchart.pdf`.

#### Briefly, all human intron coordinates were collected from Vast-DB `PSI_TABLE-hg38.tab.gz` as any EVENT with an ID beginning with `HsaIN`. These coordinates were fed in a strandwise fashion to bedtools getfasta, then SpliceAI to predict splice donor scores for every intronic nucleotide, using 400 nt exonic flanks for internal normalization by the canonical 5' splice site. Scores below 0.05 were filtered out (`data/Decoys/Splice_All.filtered.05min.bed`, 794,192 sites). Sites were restricted to protein-coding genes, sites within 50 nt of a canonical splice site were separated as an exonic set, and sites overlapping any annotated exon were removed, using `gencode.v49.annotation.gtf.gz` and GenomicRanges. The remaining deep-intronic sites were intersected against PRPF8 RBPnet predictions and PRPF8 and SmB CLIP peaks: sites with signal from any track form the decoy set (`decoys.bed`); sites with no signal whose intron's canonical 5' splice site carries signal from any track form the supported cryptic set (`cryptics_supported.bed`). Both sets were re-scored with SpliceAI over a local 49 nt window, filtered to a local donor score of at least 0.1, and cleared of annotated alternative 5' splice sites.

### Merged CLIP crosslink tracks

#### Samples were selected from Flow, and their crosslink files were downloaded and merged per target. Per-sample provenance for the merged tracks is in `results/supplementary_merged_clip_inputs.tsv` (Flow sample name and ID, purification target, assay, cell type, condition, source filename, read and crosslink counts, GEO accession).

| Track                 | Samples | Assay              | Cell lines          | Crosslink reads | Unique positions |
| --------------------- | ------- | ------------------ | ------------------- | --------------- | ---------------- |
| `PRPF8_merged.xl.bed` | 4       | eCLIP (ENCODE)     | HepG2, K562         | 23,626,946      | 16,052,984       |
| `SmB_merged.xl.bed`   | 8       | iCLIP (mild lysis) | HEK293, HepG2, K562 | 16,450,562      | 12,855,363       |

### Clippy peak calling

#### Peaks were called with Clippy 1.5.0 (`quay.io/biocontainers/clippy:1.5.0--pyhdfd78af_0`). Parameters were tuned separately per target against the crosslink bigWigs in Clippy's interactive mode, because the two tracks differ in assay and depth. `-n` rolling-mean window, `-w` width, `-x` minimum prominence adjust, `-mx` minimum height adjust, `-mg` minimum gene counts.

| Track | Parameters                          | Peaks   | Mean width | Median width |
| ----- | ----------------------------------- | ------- | ---------- | ------------ |
| PRPF8 | `-n 80 -w 0.5 -x 3.0 -mx 3.0 -mg 5` | 227,651 | 88.2 nt    | 82 nt        |
| SmB   | `-n 40 -w 0.5 -x 5.0 -mx 8.0 -mg 5` | 152,983 | 44.8 nt    | 42 nt        |

### Workflow

```
data/Decoys/Splice_All.filtered.05min.bed              794,192  BED6
  |  [1]  scripts/decoy_exon_overlaps.Rmd Part 1   protein_coding keep (GRanges, stranded)
  v
data/Decoys/spliceai_05min_proteincoding.bed           550,282  BED6
  |  [2]  scripts/split_canonical_exonic.sh   vs Wide_canonical_splice_sites.bed
  |-- -v -s --> results/deep_intronic_splice_sites.bed  367,089  BED6
  \-- -c -s --> results/proteincoding_splice_sites.bed           BED6 + canon_count
                   \- awk '$NF>0' -> results/exonic_splice_sites.bed  183,193
  |  [2b] scripts/decoy_exon_overlaps.Rmd Part 2   exon-overlap removal (unstranded)
  v
data/Decoys/deep_intronic_noexon_splice_sites.bed      289,640  BED6
  |  [3]  scripts/intersect_spliceai_support.sh -w 0   RBPnet / PRPF8 / SmB at the site
  |-- hits > 0  --> results/decoys.bed           19,942   BED6 + 3 counts + hits
  \-- hits == 0 --> results/cryptic_sites.bed   269,698   no signal at the site
                       |  [4] scripts/compile_decoy_intron_data.Rmd Parts 1-2
                       |      cryptic -> Vast-DB intron -> canonical 5'SS window
                       |      -> intersect with RBPnet / PRPF8 / SmB peaks
                       v
                     results/cryptics_supported.bed   191,904   canonical 5'SS with signal
  |  [5]  scripts/SpliceAI_Inference.py   local re-score of decoys.bed and cryptics_supported.bed
  v
  [6]  scripts/compile_decoy_intron_data.Rmd Parts 3-4   local SpliceAI >= 0.1 -> remove HsaALTD donors
       -> results/decoys_final.bed (6,506), results/cryptics_supported_final.bed (17,721) -> feature table
```

### Step 1 - protein-coding base set

#### Part 1 of `scripts/decoy_exon_overlaps.Rmd` keeps SpliceAI inferences (>= 0.05) that overlap a feature of a `gene_type == "protein_coding"` gene, on the same strand, using `reference/gencode.v49.annotation.gtf.gz`. Exonic sites are retained at this step; the exonic/intronic split is step 2.

### Step 2 - split into deep intronic and exonic borders

#### `$ bash scripts/split_canonical_exonic.sh`

#### `Wide_canonical_splice_sites.bed` is BED6 with 192,965 intervals, each 100 bp (+/-50 nt around a canonical splice site). "Exonic" therefore means within 50 nt of an annotated splice site.

| Output                                   | How               | Contents                           |
| ---------------------------------------- | ----------------- | ---------------------------------- |
| `results/deep_intronic_splice_sites.bed` | `intersect -v -s` | no canonical-window overlap        |
| `results/proteincoding_splice_sites.bed` | `intersect -c -s` | full set + canonical-overlap count |
| `results/exonic_splice_sites.bed`        | `awk '$NF>0'`     | canonical + alternative 5'SS       |

#### The script asserts `intronic + exonic == input`.

### Step 2b - exon-overlap removal on the deep-intronic branch

#### Part 2 of `scripts/decoy_exon_overlaps.Rmd` reloads `results/deep_intronic_splice_sites.bed`, removes every site that overlaps any GENCODE v49 exon on either strand, and writes `data/Decoys/deep_intronic_noexon_splice_sites.bed`. This removes sites inside alternative or internal exons and inside exons of overlapping genes. The exon set includes `retained_intron` transcripts, so a site inside an annotated retained intron is removed. 77,449 of 367,089 deep-intronic sites are removed, 17,432 of them only because of `retained_intron` transcripts, leaving 289,640.

### Step 3 - CLIP support split

#### `$ bash scripts/intersect_spliceai_support.sh -a data/Decoys/deep_intronic_noexon_splice_sites.bed -o results/clip_support -w 0 RBPNET=data/Decoys/rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8=data/CLIP/PRPF8_clippy_w0.5_rollmean80_minHeightAdjust3.0_minPromAdjust3.0_minGeneCount5_Peaks.bed SmB=data/CLIP/SmB_clippy_n40_w0.5_rollmean40_minHeightAdjust8.0_minPromAdjust5.0_minGeneCount5_Peaks.bed`

#### `$ mv results/clip_support_w0.bed results/decoys.bed`

#### `$ mv results/clip_support_w0_nohits.bed results/cryptic_sites.bed`

#### One pass writes both branches. `results/decoys.bed` holds sites supported by at least one track at the site (19,942: RBPnet 9,526, PRPF8 6,805, SmB 5,668; 1,865 supported by two or more). `results/cryptic_sites.bed` holds sites supported by none (269,698). Both carry the layout `chr, start, end, gene, spliceai_score, strand, rbpnet, prpf8, smb, hits`, so a support tier is selectable with e.g. `awk -F'\t' '$10==3'`. The script asserts `supported + unsupported == input`.

#### A zero-nt window is used because the peak widths (88 nt mean for PRPF8, 45 nt for SmB) supply the positional tolerance.

### Step 4 - supported cryptics: canonical 5'SS supported by RBPnet, PRPF8 or SmB

#### A cryptic site has no signal at the site itself and is kept only if the canonical 5' splice site of its intron is supported by RBPnet, PRPF8 or SmB. Part 1 of `scripts/compile_decoy_intron_data.Rmd` overlaps each cryptic site with the Vast-DB introns in `PSI_TABLE-hg38.tab.gz`, takes each intron's donor (`+`: intron start, `-`: intron end), and matches it to a `Canonical_splice_sites.bed` window. It writes the matched windows to `results/cryptic_canonical_sites.bed`, with the cryptic ID (`GENE_start`, e.g. `HNRNPH1_179620582`) in column 4. After the intersect below, Part 2 keeps the cryptic sites that have at least one supported window and writes `results/cryptics_supported.bed` in the same 10-column layout as `decoys.bed`.

#### `$ bash scripts/intersect_spliceai_support.sh -a results/cryptic_canonical_sites.bed -o results/canonical_support -w 0 RBPNET=data/Decoys/rbpnet_clippy_f50_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed PRPF8=data/CLIP/PRPF8_clippy_w0.5_rollmean80_minHeightAdjust3.0_minPromAdjust3.0_minGeneCount5_Peaks.bed SmB=data/CLIP/SmB_clippy_n40_w0.5_rollmean40_minHeightAdjust8.0_minPromAdjust5.0_minGeneCount5_Peaks.bed`

#### `$ mv results/canonical_support_w0.bed results/canonical_supported_sites.bed`

#### 269,698 cryptic sites → 267,465 inside a Vast-DB intron → 267,309 with a canonical donor for that intron → 191,904 supported cryptics. By track at the donor: RBPnet 135,874, PRPF8 119,032, SmB 25,196.

#### `Canonical_splice_sites.bed` is built by `scripts/CreateCanonicalSpliceSiteBed.sh` from the central +/-5 nt of each `Wide_canonical_splice_sites.bed` window.

### Step 5 - local SpliceAI re-score

#### `scripts/SpliceAI_Inference.py` re-scores both sets: it slops 24 nt either side, extracts the strand-aware sequence with `bedtools getfasta`, and runs SpliceAI over the 49 nt window to obtain a local donor score. Only column 5 is rewritten, so the support counts are preserved. The environment is `env/spliceai_environment.yml`.

#### `$ conda activate spliceai-env`

#### `$ python3 scripts/SpliceAI_Inference.py --bed results/decoys.bed --fasta <GRCh38.primary_assembly.genome.fa> --genome <genome.sizes> --batch-size 64 --out data/Decoys/decoys_splicescores.bed`

#### `$ python3 scripts/SpliceAI_Inference.py --bed results/cryptics_supported.bed --fasta <GRCh38.primary_assembly.genome.fa> --genome <genome.sizes> --batch-size 64 --out data/Decoys/cryptics_supported_splicescores.bed`

#### The script checkpoints to `<out>.scores`; `--resume` continues an interrupted run. `scripts/slurm_rescore_cryptics.sh` runs the same script on SLURM as parallel slices of the input and merges the output in input order.

### Step 6 - local-score filter and alternative 5' splice-site removal

#### Part 3 of `scripts/compile_decoy_intron_data.Rmd` loads both re-scored files as one table with a `site_class` column (`decoy` / `cryptic_supported`) and applies the local SpliceAI >= 0.1 filter (`MIN_LOCAL_SPLICEAI`): 19,942 → 7,021 decoys and 191,904 → 18,808 supported cryptics. This threshold applies to the local 49 nt score and is not comparable to the 0.05 applied to the whole-intron inference.

#### Annotated alternative 5' splice sites are then removed: any site whose position (BED start + 1) and strand match a donor of a Vast-DB `HsaALTD` (Alt5) event in `PSI_TABLE-hg38.tab.gz`, with every donor parsed from the event's `FullCO` field (133,342 unique donors). This removes 515 decoys and 1,087 supported cryptics. The final sets are `results/decoys_final.bed` (6,506) and `results/cryptics_supported_final.bed` (17,721).

#### The HNRNPH1 alternative 5' splice site at chr5:179,623,595 (Vast-DB `HsaALTD0003092-2`) is not present in `Splice_All.filtered.05min.bed` and is therefore absent from both the exonic and intronic branches.

## Decoy Feature Table Generation

#### Part 3 of `scripts/compile_decoy_intron_data.Rmd` overlaps the final decoy and supported cryptic sites with Vast-DB intron coordinates using GenomicRanges, integrating the unique identifier `EVENT` and intron retention PSI values in 145 cell and tissue types. The shortest overlapping intron is retained per site. Distance from the canonical 5' splice site is calculated with strandwise logic. Part 3 writes `results/decoy_intron_overlap_step1.tsv` and the BED files for the following steps, which run once on both classes together:

#### `$ bash scripts/extract_phastcons_scores.sh reference/hg38.phastCons100way.bw results/intron_segments_for_phastcons.bed > results/phastcons100_by_decoy.tsv`

#### `$ bash scripts/extract_phastcons_scores.sh reference/hg38.phastCons470way.bw results/intron_segments_for_phastcons.bed > results/phastcons470_by_decoy.tsv`

#### `$ bash scripts/run_maxentscan_decoys.sh results/decoy_coords_for_maxent.bed results/maxent_decoy.bed <GRCh38.primary_assembly.genome.fa>`

#### `$ bash scripts/run_maxentscan_canonical.sh results/canonical_coords_for_maxent.bed results/maxent_canonical.bed <GRCh38.primary_assembly.genome.fa>`

#### `$ bash scripts/extract_gc_content.sh <GRCh38.primary_assembly.genome.fa> results/intron_segments_for_phastcons.bed results/decoy_coords_for_maxent.bed > results/gc_by_decoy.tsv`

#### MaxEntScan scores the strength of each site and of the canonical 5' splice site of its intron. phastCons 100-way and 470-way scores are averaged across the intron harboring each site. GC content is calculated for the intron and for the 49 nt window around the site. Part 4 reloads the step-1 table, merges these results by site ID, counts the tissues with PSI >= 10 (all tissues, and Brain tissues from `data/vastdb_Sample_Groups.csv`), adds mESC IRFinder retention for introns lifted to mm10, and writes the final feature table `results/decoy_intron_features_final.tsv`. Both classes have identical columns, with `site_class` separating them. `results/decoy_vs_cryptic_class_summary.tsv` gives n, median and interquartile range of each feature per class.

##### Nine decoys have no overlapping Vast-DB intron and are dropped at the overlap step:

##### "FCGR2A_161510691" "PRR36_7873710" "ZNF44_12276147" "RP13-152O15.5_64057209" "PBRM1_52679911" "CYP3A5_99665358" "PRAG1_8386528" "VAV2_133780190" "GCNA_71597874"

##### All sites in these genes are written to `results/dropped_overlap_genes_all_splicescores_decoys.bed` for inspection.

## Figure 1 R Markdown (`scripts/hnRNPH1_figure1.rmd`)

#### This R Markdown document generates Figure 1: a comparative multiple-sequence alignment view of the hnRNPH1 intron 4 decoy region (`chr5:179620560-179620600`, hg38). It starts from an extracted UCSC multiz MAF alignment (converted to FASTA), plots the raw alignment, then creates a manuscript-ready alignment by renaming taxa, converting DNA bases from `T` to `U`, removing selected outlier species, and dropping columns that are gaps/missing across all taxa.

The final plot highlights the proposed decoy site (positions 27-33 in the processed alignment; labeled as genomic interval `179,620,576-179,620,582`) and is intended for direct use in manuscript figure generation.

### Inputs used by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.fa`
- `data/tree_to_clade_mapping.tsv`

### Output written by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.processed.fa`

## HNRNPH1 cross-species SpliceAI scores (`scripts/maf_spliceai_single_locus.py`)

#### Scores the HNRNPH1 intron-4 decoy locus in every species of the UCSC 100-way multiz alignment. The script trims `reference/chr5.maf.gz` to the region with kent `mafsInRegion`, converts the alignment to FASTA with PHAST `msa_view` (`-V` reverse-complements, because HNRNPH1 is on the minus strand), and runs SpliceAI on each species' ungapped sequence. It requires `spliceai-env`, with `mafsInRegion` and `msa_view` on `PATH` (or passed with `--mafs-in-region` / `--msa-view`).

#### `$ python3 scripts/maf_spliceai_single_locus.py --maf reference/chr5.maf.gz --out results/hnrnph1_179620582_maf_spliceai.tsv`

#### The default region is `chr5:179620038-179620941` (1-based, inclusive). The output has one row per species (100 rows). `max1`/`pos1` is the canonical 5' splice-site peak and its position in the window. `max2`/`pos2` is the strongest donor peak further into the intron; the decoy is at `pos2 = 358` in hg38. `d` is the intron-length divergence from hg38, `abs(n - v) / v`. `scripts/hnrnph1decoyposition.Rmd` plots `max2` against `pos2` from this file.
