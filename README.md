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

#### `$ python3 ../scripts/SpliceAI_Inference.py --bed Decoys/DecoySpliceSites_proteincoding.bed --fasta ../../../reference/genomes/Gencode49/GRCh38.primary_assembly.genome.fa --genome ../../../reference/genomes/Gencode49/genome.sizes --out Decoys_proteincoding_slop24.bed`

### Decoy Feature Table Generation

#### The database of decoys in 'protein-coding' introns is uploaded as `Decoys_proteincoding_slop24.bed` in the `compile_decoy_intron_data.Rmd` to integrate intron retention quantification data from `PSI-TABLE-hg38.tab.gz`. The Rmd file uses Genomic Ranges to combine PSI values in 145 cell and tissue types with SpliceAI inferences, MaxEntScan scoring for canonical and decoy 5' splice sites, calculate nucleotide distance of decoys from the canonical splice site, and integrate average phastCons 100-way and 470-way scoring across the intron harboring the decoy. 

## Figure 1 R Markdown (`scripts/hnRNPH1_figure1.rmd`)

#### This R Markdown document generates Figure 1: a comparative multiple-sequence alignment view of the hnRNPH1 intron 4 decoy region (`chr5:179620560-179620600`, hg38). It starts from an extracted UCSC multiz MAF alignment (converted to FASTA), plots the raw alignment, then creates a manuscript-ready alignment by renaming taxa, converting DNA bases from `T` to `U`, removing selected outlier species, and dropping columns that are gaps/missing across all taxa.

The final plot highlights the proposed decoy site (positions 27-33 in the processed alignment; labeled as genomic interval `179,620,576-179,620,582`) and is intended for direct use in manuscript figure generation.

### Inputs used by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.fa`
- `data/tree_to_clade_mapping.tsv`

### Output written by the Figure 1 workflow

- `data/hnRNPH1_intron4decoyMSA.processed.fa`

