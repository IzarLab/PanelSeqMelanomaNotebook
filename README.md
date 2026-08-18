# PanelSeq Melanoma Reproduction Code

This repository contains the code used to generate the staged inputs and analysis notebooks for the melanoma panel-seq figures.

It is meant to stand on its own apart from raw cohort data, external references, and one preformatted clinical metadata file used by the notebooks in `analysis/`.

## Repository layout

`01_fastq_to_bam/`

- BAM-generation workflow and helper wrappers

`02_somatic_vcf/`

- somatic calling
- post-calling filtering
- per-sample OncoKB annotation
- merged OncoKB summary generation

`03_baf_and_coverage/`

- BAF extraction from pair VCFs
- coverage-matrix generation for ClinCNV

`04_clincnv_runs/`

- ClinCNV run scripts
- helper scripts for building the combined panel cohort
- aneuploidy-score generation for Figure 1

`analysis/`

- figure notebooks
- minimal PyClone-VI wrapper for the Figure 4B clone-dominance inputs

## Analysis notebooks

`analysis/figure_1_scripts.ipynb`

- Figure 1C-F

`analysis/figure_2_scripts.ipynb`

- Figure 2 panels from the shared clinical metadata file plus staged VCF, OncoKB, and ClinCNV inputs

`analysis/figure_3_scripts.ipynb`

- paired primary/metastasis plots and ternary summaries

`analysis/figure_4_scripts.ipynb`

- Figure 4B clone-dominance score plot

`analysis/cox_analysis_scripts.ipynb`

- Cox tables and forest plots rebuilt from the shared clinical metadata file plus staged VCF, OncoKB, and ClinCNV inputs

## Notes

- Some filenames retain their original working names. In this repo they should be read simply as the exact commands or helpers used in the analysis chain.
- The notebooks are intended to consume already-staged inputs, not to do raw metadata cleanup.
- Clinical metadata is expected as one already-prepared external file rather than being rebuilt inside the notebooks.
- Some stages depend on external infrastructure being available, especially:
  - AWS/S3 access
  - local references under `/home/ubuntu/reference`
  - Conda environment `mamba-env`
  - `BedCoverage`, `bedtools`, `samtools`, `bwa`, and GATK in the expected locations
- This repo is code-first: the analysis notebooks are thin consumers of staged outputs rather than all-in-one pipeline drivers.
