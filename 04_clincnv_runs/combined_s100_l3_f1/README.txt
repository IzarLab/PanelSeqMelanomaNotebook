combined_s100_l3_f1

This folder contains the helper scripts used to assemble the combined ClinCNV cohort used by the figure notebooks:
/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/combined_bed_split_s100_l3_f1

Workflow:
1. build_pair_subsets.py splits pairs_df_filtered.csv into panel-specific pair CSVs using MEL_Annotation.xlsx.
2. build_subset_coverages.py writes panel-specific coverage matrices under coverage_subsets/.
3. launch_split_runs.sh starts the ssSC_v4LI and ssSC_v5 ClinCNV runs.
4. build_combined_run.py copies the pair-level outputs back into one combined somatic tree and injects provenance headers.
5. build_aneuploidy_scores.py converts the combined pair-level ClinCNV CNA tables into autosomal arm-level aneuploidy inputs and score summaries used by the Figure 1 metric plots.

Expected counts:
- ssSC_v4LI: 232 pairs
- ssSC_v5: 108 pairs
- combined output: 340 pairs
