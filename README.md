# Decomposition of somatic mutation profiles into mutational signatures

## Project Overview
Ampullary carcinoma (AMPAC) is a rare malignancy arising at the ampulla of Vater—the junction of the biliary, pancreatic, and duodenal epithelia. While surgical resection can be curative, a substantial fraction of patients present with unresectable disease, and reported five-year survival remains low (≈30%). Traditional pathology-based classification into intestinal vs pancreatobiliary subtypes has limited prognostic utility and does not reliably predict therapeutic response. At the molecular level, recurrent driver alterations appear across both histologic categories, and many tumors do not fall into a clear molecular class—motivating the need for a more robust taxonomy that reflects underlying mutational processes and can better support precision treatment strategies.

A principled route to refined stratification is mutational signature analysis, which models somatic mutation patterns as the cumulative footprint of endogenous and exogenous mutagenic processes. In resent study, [Zhuravleva et al. (2025)](https://gut.bmj.com/content/74/5/804) decomposed whole-exome sequencing data from 170 AMPAC tumors using non-negative matrix factorization (NMF) and then clustered tumors by their signature exposures. They reported three patient groups: an MMR-deficient/CpG-deamination cluster (SBS1 and SBS6), a cluster linked to transcription-coupled nucleotide excision repair (SBS40/SBS5), and a polymerase-η–driven hypermutation cluster (SBS9/SBS5). These subtypes were associated with Wnt alterations, DNA methylation patterns, and immune features, and the polymerase-η cluster showed improved survival—supporting the idea that signature-based classification could help identify patients most likely to benefit from Wnt-targeted approaches or immunotherapy.

The project sought to reproduce the mutational signature-based sub-typing of AMPAC from  [Zhuravleva et al. (2025)](https://gut.bmj.com/content/74/5/804) using publicly available data, and explore whether similar biological subgroups can be identified from scratch.

## Data Description 
**Source:** The primary data for analysis were obtained from the electronic supplementary material to the article by [Zhuravleva et al. (2025)](https://gut.bmj.com/content/74/5/804). The table **gutjnl-2024-333368supp011.xlsx** (journal online supplement) was used. The workbook was imported as provided; no manual edits were introduced. During import, Excel auto-formatting artefacts were programmatically reversed. 

**Data used from the workbook:**
- **Somatic SNVs (MAF-like table):** used to build the cleaned mutation table exported as `Code_Outputs/maf_from_excel_sheet4.tsv`. *(excel sheet 4)*
- **Clinical annotations (per sample):** used for the clinicogenomic heatmap and interpretation, including `Subtype_Morphology`, `gender`, `msi_mss`, `Group_Wnt`. *(excel sheet 3)*
- **Survival variables:** `SURVIVAL_TIME_60months` and `patient.vital_status_60months` were used for Kaplan–Meier analysis by cluster. *(excel sheet 3)*

## Analysis Workflow
All analyses were run in R using `InternshipFinal.qmd`.
1. **Data loading and recovery (Excel → MAF):** Imported the SNV table from the supplementary Excel file and performed value recovery to undo Excel auto-conversions. The raw table was heavily distorted by Excel’s auto-formatting (e.g. strings like 3/8 or MARCH1 were misinterpreted as dates). I therefore loaded cells with raw types and systematically recovered original values. Gene symbols were restored using their Entrez IDs (via `org.Hs.eg.db`). Numeric coordinates (`Start_Position`, `End_Position`) that Excel had turned into dates were converted back by subtracting Excel’s base date offset. Columns like `Exon_Number`, `cDNA_position`, etc., were similarly parsed to correct partial “MM/YYYY” or “M/D” date formats. The recovered table was converted into a MAF object and exported as `maf_from_excel_sheet4.tsv`.
2. **Variant/sample filtering:** Restricted the dataset to SNVs (`Variant_Type == "SNP"`), retained only high-confidence calls (`FILTER == "PASS"`), and excluded samples with <15 SNVs to improve robustness of downstream NMF-based signature extraction.
3. **Quality control (QC) plots:** 
    - Mutational burden: Plot overall mutation burden and variant distribution per sample. This includes a `maftools` summary PDF (`Maf_summary.pdf`)
   and a bar chart of SNP counts (`SNP_per_sample_all.pdf`, `SNP_counts_histogram.pdf`).
    - Transition/transversion ratio: Compute Ti/Tv per sample and plot as `QC_TiTv_plot.pdf`.
    - Rainfall plots: For the top 2–3 samples by mutation count, generate rainfall plots to check for hypermutation foci (`QC_Rainfall_plots.pdf`).
 
4. **SBS-96 matrix construction (hg19/GRCh37):** Built 96-channel trinucleotide spectra using two complementary routes:
    - `sigminer::sig_tally()` to obtain an SBS-96 matrix (saved as `mut_matrix_96.csv`) directly from the cleaned MAF; (*raw counts, used for COSMIC matching*)
    - `SomaticSignatures` (`VRanges` → `mutationContext()` → `motifMatrix()`) to build an alternative SBS-96 representation (saved as `mut_matrix_96_normalized_hg19.csv`), allowing alternative opportunity-based normalization. (*used for clustering*)

      There are two normalization schemes: 
        1. Per-sample frequency (convert counts to fractions so each sample sums to 1), which preserves spectrum shape. *(better with COSMIC)*
        2. Opportunity normalization (divide counts by trinucleotide background frequency then renormalize).  *(closest to the [Zhuravleva et al.'s (2025)](https://gut.bmj.com/content/74/5/804) method)*.
      
       Ultimately I used normalization (a) for clustering because it retains compatibility with COSMIC scales. Context channels were aligned to the standard SBS-96 naming/order (e.g., CA A.A to A[C>A]A) for compatibility across tools. 
  
> [!TIP]
> Use normalization (a), because it preserves the relative shape of the 96-channel spectrum on the same scale in which COSMIC signatures are specified (distributions with sum = 1). With normalization (b) we
> additionally divide by contexts by the frequencies of 3-mers from the reference set; this changes the ratios between channels and especially
> relabels rare contexts. Also, due to the use of “genomic” norms, the profile can additionally be shifted because of differences in the composition
> of triplets, and the cosine similarity with COSMIC decreases. Moreover, estimates via `kmerFrequency` can introduce variability, and dividing
> by very rare triplets amplifies noise with a small number of mutations. As a result, the extracted signatures turn out to be farther from the COSMIC standards. If the priority is stable matching with COSMIC, use normalization (a).
           
5. **Signature number selection (k-sweep):** Perform non-negative matrix factorization (NMF) over k = 2:7 signatures. We first estimate optimal rank via `NMF::nmfEstimateRank`(examined cophenetic correlation (stability of clustering of solutions), residual errors (RSS/Residuals), proportion of explained variability, cluster silhouette, etc.), saving `k_sweep_NMF_estimateRank_overview.pdf`. For each k, fit NMF (500 restarts), compute COSMIC similarity, and save outputs in subfolders `k_2/, k_3/, …, k_7/` under `Code_Outputs/Signature_extraction/`. Each `k_*` folder contains:
    - `signature_profiles_96bins.csv` (96×k signature profiles),
    - `cosmic_similarity_wide.csv` (wide table of cosine and peak-sensitive similarities to all COSMIC signatures),
    - `top3_cosine.csv` and `top3_peak_sensitive.csv` (top-3 COSMIC matches by each metric),
    - `exposures.csv` (sample exposures to each signature).
      
    I also compiled the mean top-COSMIC similarity vs. k into `k_sweep_mean_best_cosines.csv` and plotted it as `k_sweep_mean_best_cosines.pdf`.

6. **Final signature extraction (k=3):** Based on the sweep results and rank estimation, fixed k = 3 and ran NMF with 1000 random restarts for stability. Results are in `Code_Outputs/Signature_extraction/Final_k3_signatures/k_3/`, including the same CSVs as above. In addition, I saved visualizations of the final signatures: `signature_profiles.pdf` (SBS plots) and `signature_exposures_cosmic.pdf` (stacked exposure barplots).
   
7. **COSMIC matching:** Compared the extracted signatures to COSMIC v3.3 SBS references using two metrics: standard cosine similarity (all 96 channels) and a peak-sensitive cosine that focuses on dominant mutation peaks (ignoring minor noise). Top matches were recorded.

8. **Hierarchical clustering:** Clustered the sample exposures (z-score normalized by signature) using complete linkage on 1−Pearson distance. The resulting heatmap of mutational profiles (`mutational_profiles_heatmap.pdf`) shows clusters annotated by color. Cluster assignments (sample → cluster) were saved as `cluster_assignments.csv`.
 
9. **Cluster interpretation:** Analyzed signature exposures by cluster. Generated:
    - `cluster_mean_exposures.csv`: the mean exposure fraction of each signature in each cluster.
    - `cluster_majority_dominant_counts.csv`: for each cluster, the most frequent “dominant” signature (the signature with highest fraction in most samples).
    - `samples_with_dominant_signature.csv`: for each sample, its dominant signature and fraction.

    Also plotted a heatmap of cluster-mean exposures (`heatmap_cluster_mean_exposures.pdf`).

10. **Cluster stability:** To assess robustness, applied the auto-clustering pipeline across multiple distance/linkage combinations. Summarized the chosen number of clusters and silhouette scores per metric in `robustness_summary.csv` and visualized agreement between solutions via an ARI (Adjusted Rand Index) heatmap (`robustness_ARI_heatmap.pdf`).
    
11. **Kaplan–Meier survival:** Merged cluster labels with patient survival data and plotted Kaplan–Meier curves (log-rank test) for overall survival by cluster and saved as `kaplan_meier_by_cluster.pdf`. Pairwise log-rank tests were also computed.

12. **Clinicogenomic heatmap:** Generated a composite heatmap showing sample annotations (morphology, gender, MSI status, etc.) and numeric tracks (mutations, copy-number events) alongside the clustered mutation profiles. (saved as `Figure_3B_analog.pdf`).

Each step’s code output (plots, tables, heatmaps) is saved under `Code_Outputs/`.

## Results Summary 
This project reproduced a mutational signature–based classification of AMPAC. Identified three stable signatures. Sig1 matched COSMIC SBS1/SBS6/SBS15, consistent with spontaneous deamination of methylated cytosines and defective mismatch repair. Sig2 corresponded primarily to SBS4 and SBS29 (tobacco smoking and aflatoxin exposure) and included SBS95, a recently described signature of tobacco carcinogens. Sig3 matched SBS3, a hallmark of homologous-recombination deficiency, and was reminiscent of the polymerase-η–associated SBS9 in the original study. Comparing signature exposures across samples revealed three clusters with distinct aetiological themes:

- **Cluster 3 (Sig1‐dominant)** – characterised by high levels of SBS1/SBS6/SBS15, this cluster reflects a hypermutational, mismatch‐repair–deficient phenotype. The enrichment of SBS6 and SBS15 indicates microsatellite instability, and the accumulation of unrepaired point mutations and indels is analogous to the MSI‐high subtype described in the [Zhuravleva et al.(2025)](https://gut.bmj.com/content/74/5/804) study. In other cancers, tumours with MMR deficiency respond well to PD‐1/PD‐L1 blockade; whether this applies to AMPAC remains to be tested.
- **Cluster 2 (Sig2‐dominant)** –enriched for SBS4, SBS29 and SBS95, this cluster likely represents tumours exposed to exogenous carcinogens, most plausibly tobacco smoke. SBS4 arises from bulky DNA adducts produced by benzoapyrene in tobacco, while SBS29 and SBS95 appear to capture variations of this smoking signature. The correspondence between this cluster and the TC‐NER‐associated C2/C3 clusters in the original article suggests that environmental mutagens are important drivers of AMPAC in a subset of patients.
- **Cluster 1 (Sig3‐dominant)** – dominated by SBS3, this cluster points to homologous‐recombination deficiency (HRD), a defect in DNA double‐strand‐break repair often caused by BRCA1/2 mutations. HRD can lead to an ultra‐mutated phenotype and recruitment of translesion polymerases such as polymerase η, linking this cluster to the SBS9‐driven C3 group reported by [Zhuravleva et al.(2025)](https://gut.bmj.com/content/74/5/804). Tumours with HRD signatures may be sensitive to platinum drugs and PARP inhibitors; thus, identifying this subgroup has therapeutic implications.

While our signatures broadly matched Zhuravleva et al. (2025), cluster correspondence shifted (Sig3-dominant ≈ their C3; Sig1-dominant ≈ their C1).  The presence of SBS95 in the smoking-related cluster suggests either previously unrecognised variation within tobacco-associated processes or technical artefacts. Clustering showed low silhouette and poor reproducibility across distance metrics, consistent with a continuum rather than discrete groups. Together with the modest cohort (n=71), this likely contributed to the lack of survival differences (log-rank p≈0.71). Despite these limitations, the study underscores the potential value of mutational signatures in understanding AMPAC biology. Showed that distinct mutational processes – MMR deficiency, HRD and exposure to tobacco carcinogens – can shape the AMPAC mutational landscape. Future work should integrate larger cohorts, whole-genome sequencing and additional omics layers (transcriptomics, methylation, immune profiling) to refine classification and validate therapeutic associations.

## Software and Requirements
**Runtime**
- **R:** v4.5.1
- **Execution format:** The full pipeline is implemented in Quarto (`InternshipFinal.qmd`) and can be run in RStudio (Render/Knitr) or via the Quarto CLI.

**Dependencies:** All required packages are installed and loaded at the beginning of `InternshipFinal.qmd` (BLOCK 0) using a helper function (`install_if_missing()`) that installs missing packages from `CRAN` or `Bioconductor` with a fixed CRAN mirror.

**Core packages by role:**
- Mutational signatures / SBS-96 construction: `sigminer`, `SomaticSignatures`
- MAF processing and mutation QC (Ti/Tv, oncoplots, rainfall): `maftools`
- Reference genome & genomic infrastructure (hg19/GRCh37): `BSgenome.Hsapiens.UCSC.hg19`, `VariantAnnotation`, `IRanges`
- Gene annotation / ID mapping: `org.Hs.eg.db`, `AnnotationDbi`
- NMF and clustering: `NMF`, `dynamicTreeCut`, `cluster`, `clue`, `dendextend`
- Heatmaps and complex annotations: `ComplexHeatmap`, `pheatmap`, `circlize`, `RColorBrewer`, `ggplotify`
- Data wrangling and plotting: `dplyr`, `tidyr`, `reshape2`, `ggplot2`, `scales`
- Survival analysis: `survival`, `survminer`
- Input I/O: `readxl`

**Input requirements:** The file `gutjnl-2024-333368supp011.xlsx` must be located in the **same directory** as `InternshipFinal.qmd` (project root), as the script reads it by relative path.

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
## References
A complete description of the methodology and all cited sources are documented in `Internship.pdf`. The corresponding BibTeX bibliography used in the report is provided in `booksIntern.bib`.

