# Decomposition of somatic mutation profiles into mutational signatures

## Project Overview
## Data Description 
## Analysis Workflow
## Results Summary 

## Repository Structure
```text
Internship_project_Mutational_signatures_analysis/
├── Internship.pdf
├── InternshipFinal.qmd
├── README.md
├── Zhuravleva et al. (2025).pdf
├── booksIntern.bib
└── Code_Outputs/
    ├── maf_from_excel_sheet4.tsv
    ├── mut_matrix_96.csv
    ├── mut_matrix_96_normalized_hg19.csv
    ├── Clustering_robustness/
    │   ├── robustness_ARI_heatmap.pdf
    │   └── robustness_summary.csv
    ├── Clusters&Mutational_profile_heatmap&KM_curve/
    │   ├── cluster_assignments.csv
    │   ├── kaplan_meier_by_cluster.pdf
    │   └── mutational_profiles_heatmap.pdf
    ├── Dominant_signature_per_cluster/
    │   ├── cluster_majority_dominant_counts.csv
    │   ├── cluster_mean_exposures.csv
    │   ├── heatmap_cluster_mean_exposures.pdf
    │   ├── samples_with_dominant_signature.csv
    │   └── stacked_by_sample_sorted_within_cluster.pdf
    ├── Final_k3_signatures/
    │   ├── cosmic_similarity_wide.csv
    │   ├── exposures.csv
    │   ├── signature_exposures_cosmic.pdf
    │   ├── signature_profiles.pdf
    │   ├── signature_profiles_96bins.csv
    │   ├── top3_cosine.csv
    │   └── top3_peak_sensitive.csv
    ├── Mutational_Landscape_plots/
    │   ├── Figure_3B_analog.pdf
    │   ├── Maf_summary.pdf
    │   ├── QC_Rainfall_plots.pdf
    │   ├── QC_TiTv_plot.pdf
    │   ├── SNP_counts_histogram.pdf
    │   ├── SNP_per_sample_all.pdf
    │   └── oncoplot_top25.pdf
    └── Signature_count_selection/
        ├── k_sweep_NMF_estimateRank_overview.pdf
        ├── k_sweep_mean_best_cosines.csv
        ├── k_sweep_mean_best_cosines.pdf
        ├── k_2/
        │   ├── cosmic_similarity_wide.csv
        │   ├── exposures.csv
        │   ├── signature_profiles_96bins.csv
        │   ├── top3_cosine.csv
        │   └── top3_peak_sensitive.csv
        ├── k_3/ (analogous to k_2 folder contents)
        ├── k_4/ (analogous to k_2 folder contents)
        ├── k_5/ (analogous to k_2 folder contents)
        ├── k_6/ (analogous to k_2 folder contents)
        └── k_7/ (analogous to k_2 folder contents)
```
## Target Audience

**References:** Full citations and discussion of methods are provided in the report and accompanying bibliography.
