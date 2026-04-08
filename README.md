# FASTQ To ClinCNV Transparency Bundle

This folder collects the local scripts and notebooks used to go from raw paired-end FASTQs to somatic ClinCNV calls in the WES BedCoverage workflow.

The goal is transparency rather than abstraction: these are copies of the working scripts we used, grouped in pipeline order, plus a stitched notebook that explains how the outputs of one step become the inputs of the next.

## Overall Pipeline

1. `01_fastq_to_bam/`
   GATK-style preprocessing from paired FASTQs to deduplicated BAMs and BAIs.

2. `02_somatic_vcf/`
   Tumor-normal somatic calling with `Mutect2`, followed by `LearnReadOrientationModel`, `CalculateContamination`, and `FilterMutectCalls`.

3. `03_baf_and_coverage/`
   Creation of the cleaned per-sample BAF folder from pair VCFs, generation of the WES off-target BED, and creation of the BedCoverage matrices used by ClinCNV.

4. `04_clincnv_runs/`
   Final ClinCNV run scripts for the WES BedCoverage analysis:
   one with `--bafFolder` and one without.

## What Gets Computed At Each Stage

### Stage 1: FASTQs to BAMs

Starting inputs:

- paired FASTQ files per sample, stored in S3
- GRCh38 reference FASTA, FAI, and dictionary
- capture BED used for the panel/WES setup

Computed outputs:

- one deduplicated BAM per sample
- one BAM index per sample
- per-sample logs and metrics

How the inputs are assembled:

- wrapper scripts read a CSV with sample IDs and S3 locations
- FASTQs are discovered under the sample S3 prefix, or read from an explicit FASTQ list column
- if multiple FASTQ members exist for a mate, they are concatenated into one merged `R1` and one merged `R2`
- a per-sample `config.yaml` is written and passed into the Snakemake workflow

How the BAMs are computed:

- `FastqToSam` creates an unmapped BAM
- `RevertSam` normalizes the unmapped BAM
- `MarkIlluminaAdapters` marks adapters
- `SamToFastq | bwa mem | MergeBamAlignment` performs alignment and merges aligned reads back with the unmapped BAM metadata
- `MarkDuplicates` creates the final deduplicated BAM and BAI

### Stage 2: BAMs to Somatic VCFs

Starting inputs:

- tumor and normal deduplicated BAMs and BAIs
- GRCh38 reference FASTA, FAI, and dictionary
- target BED
- panel of normals VCF and index
- pair list describing tumor-normal matches

Computed outputs:

- unfiltered paired somatic VCFs
- Mutect2 `.stats` files
- F1R2 tarballs for orientation-bias modeling
- filtered VCFs after post-processing
- contamination tables and segmentation outputs

How the inputs are assembled:

- `vcf_gen_v3.sh` reads a `pairs.txt` file containing `tumor,normal`
- for each pair, the script downloads the tumor and normal BAM/BAI files from the BAM S3 prefix
- `filter_vcfs.sh` then scans the VCF output prefix for pair folders and downloads those local results again for post-processing

How the VCFs are computed:

- `Mutect2` is run per tumor-normal pair with the panel BED and a panel of normals
- `LearnReadOrientationModel` is run from the Mutect2 F1R2 tarball
- `GetPileupSummaries` is run for tumor and normal against a common-SNP resource
- `CalculateContamination` estimates contamination and segmentation
- `FilterMutectCalls` uses the orientation model and contamination estimates to emit the filtered somatic VCF

### Stage 3A: Pair VCFs to BAF Folder

Starting inputs:

- one pair VCF per tumor-normal pair
- pair list CSV for all included samples

Computed outputs:

- one TSV per sample in the BAF folder
- a skipped-pairs log
- a README describing the generated folder

How the inputs are assembled:

- `baf_from_pair_vcfs.sh` reads `pairs_df_filtered.csv`
- for each pair it expects a file named `unfiltered.<tumor>__<normal>.vcf.gz`
- by default it reads from `/mnt/myvolume/panel_seq/new_bed_analysis/vcfs`

How the BAF values are computed:

- only biallelic SNVs are retained
- each sample record must have non-missing `GT`, `AD`, and `DP`
- a minimum depth threshold of `DP >= 10` is applied
- BAF is computed as `alt_count / DP`
- output columns are `chr start end chr_pos baf depth`
- if the same sample appears in more than one pair, duplicate positions are collapsed by keeping the row with the highest depth

### Stage 3B: BAMs to BedCoverage Matrices

Starting inputs:

- sample BAMs from the BAM S3 prefix
- target BED
- reference FASTA and FAI
- WES off-target BED definition

Computed outputs:

- on-target coverage matrix for tumors
- on-target coverage matrix for normals
- off-target coverage matrix for tumors
- off-target coverage matrix for normals
- cached per-sample BedCoverage outputs
- logs for matrix construction

How the inputs are assembled:

- `gen_wes_offtarget_bed_100kb.sh` derives the WES off-target BED from the genome minus merged target regions
- the script bins off-target space into 100 kb windows and removes windows smaller than 50 kb
- GC is then computed per off-target window using `bedtools nuc`
- `build_wes_bedcoverage_matrices.sh` downloads each sample BAM from S3 and runs `BedCoverage`

How the coverage matrices are computed:

- on-target bins use `BedCoverage` with `min_mapq 0`
- off-target bins use `BedCoverage` with `min_mapq 10`
- per-sample outputs are cached and later merged into matrix-style `.cov` files
- tumor and normal sample order is reconstructed from header source coverage files

### Stage 4: BedCoverage Matrices to ClinCNV Calls

Starting inputs:

- tumor and normal on-target coverage matrices
- tumor and normal off-target coverage matrices
- on-target BED with GC column
- off-target BED with GC column
- tumor-normal pair CSV
- optionally, the per-sample BAF folder

Computed outputs:

- ClinCNV cohort reanalysis output
- final pair-level somatic ClinCNV calls
- ClinCNV logs and summary tables in the run directory

How the inputs are assembled:

- both run scripts point to the same BedCoverage matrices and BED files
- both use the same pair list and the same ClinCNV R script
- the only branch point is whether `--bafFolder` is supplied

How the final calls are computed:

- ClinCNV is run in WES mode with off-target support
- key parameters are:
  - `--colNum 4`
  - `--lengthS 9`
  - `--filterStep 2`
  - `--scoreS 200`
  - `--reanalyseCohort`

## File-By-File Description

### `01_fastq_to_bam/`

`snakemake_wf.py`

- Main GATK-style preprocessing workflow from paired FASTQs to final deduplicated BAM/BAI outputs.
- Expects per-sample FASTQs named `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` in the sample input folder.
- Reads tool and reference locations from `config.yaml`.
- Produces rule-specific logs and metrics for every sample.

`config.yaml`

- Template configuration for the Snakemake workflow.
- Defines input and output roots, sample IDs, GRCh38 reference files, target BED, padding, and tool names.
- In the wrappers, a per-sample config is generated dynamically from this structure.

`run_all.py`

- General wrapper for samples with one or more FASTQ members per mate under an S3 prefix.
- Reads a CSV with `sample` and `s3_dir`.
- Discovers all FASTQs under the prefix, merges all R1 files together and all R2 files together, writes a per-sample config, runs Snakemake, uploads final BAM/BAI files to S3, uploads non-FASTQ run artifacts, and deletes the local work directory.

`double_csv.py`

- Variant wrapper for samples whose FASTQs may be defined more explicitly in a CSV.
- Can either read an explicit FASTQ list column or discover FASTQs under an S3 prefix.
- Handles cases where several R1 or R2 files should be concatenated before running the Snakemake workflow.

`double_fastq.py`

- More restrictive wrapper for cases where the CSV explicitly provides exactly one R1 and one R2 URI per sample.
- Resolves FASTQ paths, downloads the two mates, writes a per-sample config, runs Snakemake, and preserves logs in a dedicated local logs directory.
- Useful when the FASTQ layout is irregular and direct discovery under a prefix is not reliable.

### `02_somatic_vcf/`

`vcf_gen_v3.sh`

- Paired `Mutect2` caller.
- Reads `pairs.txt` to define tumor-normal matches.
- Downloads the deduplicated BAMs and BAIs for each pair from the BAM S3 prefix.
- Runs `Mutect2` with the target BED, panel of normals, and GRCh38 reference.
- Emits unfiltered VCFs, `.stats`, F1R2 tarballs, and logs, then uploads them back to the VCF S3 prefix.

`filter_vcfs.sh`

- Post-processing script for the paired Mutect2 outputs.
- Scans the VCF S3 prefix for pair folders, downloads the unfiltered VCFs and F1R2 tarballs, fetches the corresponding BAMs, and runs the downstream GATK filtering steps.
- Produces artifact priors, pileup summaries, contamination tables, contamination segments, and filtered VCFs.

### `03_baf_and_coverage/`

`baf_from_pair_vcfs.sh`

- Converts pair VCFs into a ClinCNV-compatible BAF folder.
- Reads each pair VCF, extracts per-sample SNV allele fractions from `AD` and `DP`, and writes one TSV per sample.
- Also writes a skipped-pairs report and a README describing the generated folder.

`gen_wes_offtarget_bed_100kb.sh`

- Builds the WES off-target BED used by ClinCNV.
- Starts from the genome, subtracts the merged target intervals, windows the remaining space into 100 kb bins, removes bins smaller than 50 kb, and computes GC content for each bin.
- Produces the GC-annotated off-target BED consumed by the BedCoverage and ClinCNV steps.

`build_wes_bedcoverage_matrices.sh`

- Generates WES coverage matrices from BAMs using `BedCoverage`.
- Creates sorted versions of the target and off-target BEDs, downloads each sample BAM, computes on-target and off-target depth profiles, caches the per-sample outputs, and merges them into final tumor and normal matrix files.
- Uses MAPQ 0 on-target and MAPQ 10 off-target.

`bedcoverage_wes_full_transparency.ipynb`

- Existing explanatory notebook focused on the WES BedCoverage and ClinCNV portion of the workflow.
- Documents how the cleaned BAF folder, off-target BED, coverage matrices, and final ClinCNV runs were created.

### `04_clincnv_runs/`

`run_clincnv_baf_bedcoverage_wes.sh`

- Final ClinCNV command for the WES BedCoverage run with BAF support.
- Uses the on-target and off-target BedCoverage matrices together with the cleaned BAF folder.
- Writes outputs to `/mnt/myvolume/panel_seq/new_bed_analysis/new_clincnv_runs/full_wes_s200_l9_f2_baf_bedcoverage_wes`.

`run_clincnv_no_baf_bedcoverage_wes.sh`

- Final ClinCNV command for the matching WES BedCoverage run without BAF support.
- Uses the same matrices, BED files, pair list, and core ClinCNV parameters, but omits `--bafFolder`.
- Writes outputs to `/mnt/myvolume/panel_seq/new_bed_analysis/new_clincnv_runs/full_wes_s200_l9_f2_no_baf_bedcoverage_wes`.

## Main Bundle Files

`MANIFEST.tsv`

- Maps each bundled file to its original source path and gives a one-line purpose.

`fastq_to_clincnv_end_to_end_transparency.ipynb`

- Narrative notebook that ties the copied scripts together and spells out the full flow from FASTQs to ClinCNV.

## Important Split

The two final ClinCNV runs are identical except for BAF usage:

- `04_clincnv_runs/run_clincnv_baf_bedcoverage_wes.sh`
  Includes `--bafFolder /mnt/myvolume/panel_seq/new_bed_analysis/baf_from_pair_vcfs_clean`
- `04_clincnv_runs/run_clincnv_no_baf_bedcoverage_wes.sh`
  Omits `--bafFolder`

## Included Caveats

- These scripts point at the real local paths, references, S3 locations, and working directories that were used here.
- Some stages depend on external infrastructure being available, especially:
  - AWS/S3 access
  - local references under `/home/ubuntu/reference`
  - Conda environment `mamba-env`
  - `BedCoverage`, `bedtools`, `samtools`, `bwa`, and GATK in the expected locations
- The notebook is intentionally explanatory and reproducibility-oriented; it is not meant to rerun the whole pipeline automatically inside a single notebook session.
