# Internship project: Decomposition of somatic mutation profiles into mutational signatures

**Author:** Anastasiia Deviataieva (Student Number: 24100519)

This repository contains the Quarto document `InternshipFinal.qmd` which reads a MAF-like table from Excel, performs mutation quality control, extracts mutational signatures (SBS-96), evaluates the number of signatures, clusters samples by mutational profiles, and generates clinical and survival visualizations. All artifacts are written under a single output directory.

## Reproduce the analysis

1. **Install prerequisites**
   - [R](https://cran.r-project.org/) ≥ 4.2
   - [Quarto](https://quarto.org/) ≥ 1.4
   - Required R packages (installed automatically in the document if missing):  
     `AnnotationDbi, BSgenome.Hsapiens.UCSC.hg19, ComplexHeatmap, IRanges, NMF, RColorBrewer, SomaticSignatures, VariantAnnotation, circlize, clue, cluster, dendextend, dplyr, dynamicTreeCut, ggplot2, ggplotify, maftools, org.Hs.eg.db, pheatmap, readxl, reshape2, scales, sigminer, survival, survminer, tidyr`

2. **Configure paths (optional)**
   - Edit `InternshipFinal.qmd` near the top:
     - `out_dir <- "~/Desktop/Internship_outputs_2"` — directory where all results will be written.
     - `excel_path <- "/path/to/your.xlsx"` and `excel_sheet <- 4` — input data location.

3. **Render**
   ```bash
   quarto render InternshipFinal.qmd --to pdf
   ```

All output files will be created under **`~/Desktop/Internship_outputs_2`** (subfolders are created as needed).

## Key outputs

> The table below lists every file the script writes and what it represents, in plain English.

| File | What it means |
|---|---|
| `maf_from_excel_sheet4.tsv` | MAF-like TSV exported from Excel sheet 4 after cell-type normalization. |
| `Maf_summary.pdf` | Cohort-level MAF summary plots (variant classifications, sample/mutation counts). |
| `oncoplot_top25.pdf` | Oncoplot of the 25 most frequently mutated genes across samples. |
| `QC_TiTv_plot.pdf` | Transition/Transversion (Ti/Tv) quality-control plot per sample. |
| `QC_Rainfall_plots.pdf` | Rainfall plots for up to three most mutated samples showing inter-mutation distances (kataegis detection). |
| `SNP_per_sample_all.pdf` | Bar chart of total SNP counts per sample (mutational load by sample). |
| `SNP_counts_histogram.pdf` | Histogram of SNP count distribution across the cohort. |
| `mut_matrix_96.csv` | Raw SBS-96 mutation count matrix (rows: samples, columns: 96 trinucleotide contexts). |
| `mut_matrix_96_normalized_hg19.csv` | Opportunity-normalized SBS-96 matrix (per-sample frequencies; hg19 context). |
| `k_{2..7}/signature_profiles_96bins.csv` | Extracted mutational signature profiles for each k (96-context spectra). |
| `k_{2..7}/cosmic_similarity_wide.csv` | Similarity of extracted signatures to COSMIC references (cosine & peak-sensitive). |
| `k_{2..7}/top3_cosine.csv` | Top-3 COSMIC matches by cosine similarity for each extracted signature. |
| `k_{2..7}/top3_peak_sensitive.csv` | Top-3 COSMIC matches by peak-sensitive cosine for each extracted signature. |
| `k_{2..7}/exposures.csv` | Signature exposures per sample for each k (mutation counts per signature). |
| `k_sweep_NMF_estimateRank_overview.pdf` | NMF diagnostic multi-panel plot across k=2..7 (cophenetic, RSS, dispersion, etc.). |
| `k_sweep_mean_best_cosines.csv` | Table of the best COSMIC match similarity per k (cosine & peak-sensitive). |
| `k_sweep_mean_best_cosines.pdf` | Line plot of mean best COSMIC similarity vs k (two metrics). |
| `finalK/k_3/signature_profiles.pdf` | Plots of the final signatures’ 96-context spectra for k=3. |
| `finalK/k_3/signature_exposures_cosmic.pdf` | Plots of sample exposures to the final signatures (k=3). |
| `mutational_profiles_heatmap.pdf` | Heatmap of per-sample normalized SBS-96 profiles with hierarchical clustering. |
| `cluster_assignments.csv` | Mapping from sample to cluster label (from mutational profile clustering). |
| `robustness_summary.csv` | Clustering robustness summary over distance/linkage choices (one row per setting). |
| `robustness_ARI_heatmap.pdf` | Adjusted Rand Index heatmap comparing clusterings across methods. |
| `samples_with_dominant_signature.csv` | Per-sample dominant signature and its fraction, along with cluster assignment. |
| `cluster_mean_exposures.csv` | Mean signature fractions per cluster (rows: clusters; columns: signatures). |
| `cluster_majority_dominant_counts.csv` | Counts of samples per cluster where each signature is dominant (majority vote). |
| `heatmap_cluster_mean_exposures.pdf` | Heatmap of cluster-level mean signature fractions (no clustering). |
| `stacked_by_sample_sorted_within_cluster.pdf` | Stacked bar chart of per-sample signature fractions ordered within clusters. |
| `kaplan_meier_by_cluster.pdf` | Kaplan–Meier survival curves stratified by cluster with statistics. |
| `Figure_3B_analog.pdf` | Composite clinicogenomic heatmap (clinical annotations + z-scored numeric features). |


## Directory layout (after render)

```
~/Desktop/Internship_outputs_2/
├── maf_from_excel_sheet4.tsv
├── Maf_summary.pdf
├── oncoplot_top25.pdf
├── QC_TiTv_plot.pdf
├── QC_Rainfall_plots.pdf
├── SNP_per_sample_all.pdf
├── SNP_counts_histogram.pdf
├── mut_matrix_96.csv
├── mut_matrix_96_normalized_hg19.csv
├── k_2/ … k_7/
│   ├── signature_profiles_96bins.csv
│   ├── cosmic_similarity_wide.csv
│   ├── top3_cosine.csv
│   ├── top3_peak_sensitive.csv
│   └── exposures.csv
├── finalK/
│   └── k_3/
│       ├── signature_profiles.pdf
│       └── signature_exposures_cosmic.pdf
├── mutational_profiles_heatmap.pdf
├── cluster_assignments.csv
├── robustness_summary.csv
├── robustness_ARI_heatmap.pdf
├── samples_with_dominant_signature.csv
├── cluster_mean_exposures.csv
├── cluster_majority_dominant_counts.csv
├── heatmap_cluster_mean_exposures.pdf
├── stacked_by_sample_sorted_within_cluster.pdf
├── kaplan_meier_by_cluster.pdf
└── Figure_3B_analog.pdf
```

## Notes

- The analysis uses **hg19** contexts (`BSgenome.Hsapiens.UCSC.hg19`) for SBS-96 normalization.
- Random seeds are set (`set.seed(42)`) to improve reproducibility, but NMF has inherent stochasticity.
- The document prints informative `INFO` messages during filtering and QC; these are not saved as files.

---

*Generated automatically from the code in `InternshipFinal.qmd`.*

