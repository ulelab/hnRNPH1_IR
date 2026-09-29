## hnRNPH1 Intron 4 Decoy Analysis

#### This repository contains data, scripts, and outputs for the hnRNPH1 intron 4 retention decoy analysis workflow.

## Decoy Dataset Generation

#### The database of predicted decoy loci was generated as a BED file using the workflow described in `figures/decoy_splice_site_flowchart.pdf`.

#### Briefly, all human intron coordinates were collected from Vast-DB `PSI_TABLE-hg38.tab.gz` as any EVENT with an ID beginning with `HsaIN`. These coordinates were fed in a strandwise fashion to bedtools getfasta, then SpliceAI to predict splice donor scores for every intronic nucleotide, using 400 nt exonic flanks for internal normalization by the canonical 5' splice site. Scores below 0.05 were filtered out (`data/Decoys/Splice_All.filtered.05min.bed`, 794,192 sites). Sites were restricted to protein-coding genes, sites within 50 nt of a canonical splice site were separated as an exonic set, and sites overlapping any annotated exon were removed, using `gencode.v49.annotation.gtf.gz` and GenomicRanges. The remaining deep-intronic sites were intersected against PRPF8 RBPnet predictions and PRPF8 and SmB CLIP peaks: sites with signal from any track form the decoy set (`decoys.bed`); sites with no signal whose intron's canonical 5' splice site carries signal from any track form the supported cryptic set (`cryptics_supported.bed`). Both sets were re-scored with SpliceAI over a local 49 nt window, filtered to a local donor score of at least 0.1, cleared of annotated alternative 5' splice sites, and cleared of sites within 100 nt upstream of a canonical 3' splice site.

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
                     results/cryptics_supported.bed   191,902   canonical 5'SS with signal
  |  [5]  scripts/SpliceAI_Inference.py   local re-score of decoys.bed and cryptics_supported.bed
  v
  [6]  scripts/compile_decoy_intron_data.Rmd Parts 3-4   local SpliceAI >= 0.1 -> remove HsaALTD donors
       -> remove sites <= 100 nt upstream of a 3'SS
       -> results/decoys_final.bed (5,689), results/cryptics_supported_final.bed (16,295) -> feature table
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

#### 269,698 cryptic sites → 267,465 inside a Vast-DB intron → 267,309 with a canonical donor for that intron → 191,902 supported cryptics. By track at the donor: RBPnet 135,872, PRPF8 119,032, SmB 25,196. Sites are matched only to introns on their own strand (strand taken from `FullCO`).

#### `Canonical_splice_sites.bed` is built by `scripts/CreateCanonicalSpliceSiteBed.sh` from the central +/-5 nt of each `Wide_canonical_splice_sites.bed` window.

### Step 5 - local SpliceAI re-score

#### `scripts/SpliceAI_Inference.py` re-scores both sets: it slops 24 nt either side, extracts the strand-aware sequence with `bedtools getfasta`, and runs SpliceAI over the 49 nt window to obtain a local donor score. Only column 5 is rewritten, so the support counts are preserved. The environment is `env/spliceai_environment.yml`.

#### `$ conda activate spliceai-env`

#### `$ python3 scripts/SpliceAI_Inference.py --bed results/decoys.bed --fasta <GRCh38.primary_assembly.genome.fa> --genome <genome.sizes> --batch-size 64 --out data/Decoys/decoys_splicescores.bed`

#### `$ python3 scripts/SpliceAI_Inference.py --bed results/cryptics_supported.bed --fasta <GRCh38.primary_assembly.genome.fa> --genome <genome.sizes> --batch-size 64 --out data/Decoys/cryptics_supported_splicescores.bed`

#### The script checkpoints to `<out>.scores`; `--resume` continues an interrupted run. `scripts/slurm_rescore_cryptics.sh` runs the same script on SLURM as parallel slices of the input and merges the output in input order.

### Step 6 - local-score filter and alternative 5' splice-site removal

#### Part 3 of `scripts/compile_decoy_intron_data.Rmd` loads both re-scored files as one table with a `site_class` column (`decoy` / `cryptic_supported`) and applies the local SpliceAI >= 0.1 filter (`MIN_LOCAL_SPLICEAI`): 19,942 → 7,021 decoys and 191,902 → 18,808 supported cryptics. This threshold applies to the local 49 nt score and is not comparable to the 0.05 applied to the whole-intron inference.

#### Annotated alternative 5' splice sites are then removed: any site whose position (BED start + 1) and strand match a donor of a Vast-DB `HsaALTD` (Alt5) event in `PSI_TABLE-hg38.tab.gz`, with every donor parsed from the event's `FullCO` field (133,342 unique donors). This removes 515 decoys and 1,087 supported cryptics.

#### Sites 0–100 nt upstream of a canonical 3' splice site of any same-strand Vast-DB intron are removed next (`NEAR_3SS_NT`). This stretch holds the branch point and polypyrimidine tract, where SmB and PRPF8 crosslink as part of the spliceosome: within 31–90 nt of the 3' splice site, 66.5% of decoys were SmB-supported against 26.4% beyond 200 nt. It also holds `AG|GT` sites at the last intronic base, where the downstream exon begins with a donor-like `GT`. This removes 817 decoys and 1,426 supported cryptics. The final sets are `results/decoys_final.bed` (5,689) and `results/cryptics_supported_final.bed` (16,295).

#### The HNRNPH1 alternative 5' splice site at chr5:179,623,595 (Vast-DB `HsaALTD0003092-2`) is not present in `Splice_All.filtered.05min.bed` and is therefore absent from both the exonic and intronic branches.

## Decoy Feature Table Generation

#### Part 3 of `scripts/compile_decoy_intron_data.Rmd` overlaps the final decoy and supported cryptic sites with Vast-DB intron coordinates using GenomicRanges, integrating the unique identifier `EVENT` and intron retention PSI values in 145 cell and tissue types. The overlap is stranded, with intron strand taken from `FullCO`, and the shortest overlapping intron is retained per site. Distance from the canonical 5' splice site is calculated with strandwise logic. Part 3 writes `results/decoy_intron_overlap_step1.tsv` and the BED files for the following steps, which run once on both classes together:

#### `$ bash scripts/extract_phastcons_scores.sh reference/hg38.phastCons100way.bw results/intron_segments_for_phastcons.bed > results/phastcons100_by_decoy.tsv`

#### `$ bash scripts/extract_phastcons_scores.sh reference/hg38.phastCons470way.bw results/intron_segments_for_phastcons.bed > results/phastcons470_by_decoy.tsv`

#### `$ bash scripts/run_maxentscan_decoys.sh results/decoy_coords_for_maxent.bed results/maxent_decoy.bed <GRCh38.primary_assembly.genome.fa>`

#### `$ bash scripts/run_maxentscan_canonical.sh results/canonical_coords_for_maxent.bed results/maxent_canonical.bed <GRCh38.primary_assembly.genome.fa>`

#### `$ bash scripts/extract_gc_content.sh <GRCh38.primary_assembly.genome.fa> results/intron_segments_for_phastcons.bed results/decoy_coords_for_maxent.bed > results/gc_by_decoy.tsv`

#### MaxEntScan scores the strength of each site and of the canonical 5' splice site of its intron. phastCons 100-way and 470-way scores are averaged across the intron harboring each site. GC content is calculated for the intron and for the 49 nt window around the site. Part 4 reloads the step-1 table, merges these results by site ID, counts the tissues with PSI >= 10 (all tissues, and Brain tissues from `data/vastdb_Sample_Groups.csv`), adds mESC IRFinder retention for introns lifted to mm10, and writes the final feature table `results/decoy_intron_features_final.tsv`. Both classes have identical columns, with `site_class` separating them. `results/decoy_vs_cryptic_class_summary.tsv` gives n, median and interquartile range of each feature per class.

#### Part 4 also redraws the class comparison with both classes filtered to a site MaxEntScan score of at least 8, 9 or 10 (`figures/decoy_vs_cryptic_class_comparison_maxent{8,9,10}.pdf`), with the counts in `results/maxent_threshold_counts.tsv`.

##### Ten decoys have no overlapping same-strand Vast-DB intron and are dropped at the overlap step. Nine lie 164-339 nt outside their scored `HsaIN` intron, within the 400 nt inference flank but beyond a flanking exon shorter than 400 nt, and therefore in an adjacent intron that has no Vast-DB intron retention event. The tenth, `MSTO1_155610087`, lies only inside an intron of the antisense gene `RP11-29H23.4`:

##### "FCGR2A_161510691" "PRR36_7873710" "ZNF44_12276147" "RP13-152O15.5_64057209" "PBRM1_52679911" "CYP3A5_99665358" "PRAG1_8386528" "VAV2_133780190" "GCNA_71597874" "MSTO1_155610087"

##### All sites in these genes are written to `results/dropped_overlap_genes_all_splicescores_decoys.bed` for inspection.

## Recount Junction Support

#### `scripts/query_recount_junctions.py` matches each final site to junctions in the hg38 Snaptron/Recount compilations `tcgav2`, `srav3h` and `gtexv2` (`junctions.sqlite` from `https://snaptron.cs.jhu.edu/data/<compilation>/`). Sites are loaded into an in-memory table attached to the read-only database and joined to the `intron` table on chromosome and strand, with either junction boundary within 5 bp of the site. All matching junctions are returned, annotated or not; a junction is typed `5ss` when its boundary at the site is the donor for the site's strand (junction start on `+`, junction end on `-`). `scripts/slurm_query_recount.sh` runs the three compilations as a SLURM array and downloads the databases when absent.

#### `$ bash scripts/slurm_query_recount.sh`

#### Outputs per compilation: `results/recount_<compilation>_junctions.tsv` (one row per site-junction pair with `snaptron_id`, boundary offsets, `annotated`, `samples_count`, `coverage_sum`) and `results/recount_<compilation>_summary.tsv` (one row per site: numbers of annotated and unannotated 5'SS junctions and their maximum sample and read support). Part 5 of `scripts/compile_decoy_intron_data.Rmd` reads the junction tables and counts a junction as using a site as its 5' splice site when the junction's donor boundary is the first intronic base after the site (junction start = site + 1 on `+`, junction end = site − 1 on `−`; 250,758 of the 270,878 junction rows with a donor-side boundary within ±5 bp) and it has at least 5 reads summed over the compilation's samples. Each site gets `recount_support`: `annotated` (an annotated junction uses the site as its donor), `unannotated_only`, or `no_junction`, plus the number of compilations supporting it and the sample and read counts of its best unannotated junction.

| Site class | Annotated | Unannotated only | No junction |
| --- | --- | --- | --- |
| Decoys (5,679) | 534 (9.4%) | 4,879 (85.9%) | 266 (4.7%) |
| Supported cryptics (16,295) | 936 (5.7%) | 13,223 (81.1%) | 2,136 (13.1%) |

#### The HNRNPH1 decoy is used as a donor only in unannotated junctions, in all three compilations; its junction to the canonical intron-4 acceptor (1,175 nt) is found in 3,845 SRA samples (5,382 reads).

## Intron Summary Table

#### `scripts/intron_summary.Rmd` builds one row per Vast-DB `HsaIN` intron (192,965) from `PSI_TABLE-hg38.tab.gz` and marks the introns that hold decoys and supported cryptics from the final feature table, giving four intron classes: `none`, `decoy`, `cryptic` and `both`. Site IDs (`GENE_start`) and site counts per intron are carried over. Flanking exon coordinates are parsed from `FullCO`, and `bedtools nuc` gives the GC fraction of the intron and of the two flanking exons together; `gc_ratio` is intron GC divided by flanking-exon GC. Tissue counts with PSI >= 10 are computed for all 145 tissues and for the 26 Brain tissues. The genome FASTA path is the `fasta` parameter.

#### `$ Rscript -e 'rmarkdown::render("scripts/intron_summary.Rmd", params = list(fasta = "<GRCh38.primary_assembly.genome.fa>"))'`

#### Outputs: `results/intron_summary_table.tsv`, `results/intron_class_summary.tsv` (n, median and interquartile range of each feature per class), `figures/intron_class_comparison.pdf` and `figures/intron_class_retention.pdf`. Of the 192,965 introns, 5,430 hold decoys only, 14,377 supported cryptics only and 785 both; 172,373 hold neither.

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
