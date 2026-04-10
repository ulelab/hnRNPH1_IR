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

## Decoy Feature Table Generation

#### The database of decoys in 'protein-coding' introns is uploaded as `Decoys_proteincoding_splicescores.bed` in the `compile_decoy_intron_data.Rmd` to integrate intron retention quantification data from `PSI-TABLE-hg38.tab.gz`. The Rmd file uses Genomic Ranges to integrate intron coordinates, unique identifiers `EVENT` and intron retention PSI values in 145 cell and tissue types with SpliceAI inference. Decoy distance from canonical splice site is calculated with strandwise logic. After overlapping the decoy database with introns, MaxEntScan is used to calculate the strength of decoy predicted splice sites and the canonical 5' splice site for the intron harboring the decoy with the scripts `run_maxentscan_decoy.sh` `run_maxentscan_canonical.sh`.Average phastCons 100-way and 470-way scoring across the intron harboring the decoy is calculated with `extract_phastcons_scores.sh`. Part 2 of the R markdown file reloads the results from MaxEntScan and phastCons and merges into the final feature table.

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

