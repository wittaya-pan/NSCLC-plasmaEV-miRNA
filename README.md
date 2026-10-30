# Code for: Plasma Extracellular Vesicle-Derived MicroRNAs Associated with Overall Survival in Advanced Non-Targeted NSCLC: A Pilot Study

Manuscript CTS-2026-0503, *Clinical and Translational Science*.
Panvongsa W, Oranratnachai S, Trachu N, Khiewngam K, Monnamo N, Reungwetwattana T, Chairoungdua A.

This repository contains the code for every analysis in the manuscript and its Supporting Information,
organized by manuscript section. **No patient data are included**; the plasma EV miRNA data and de-identified clinical characteristics are in
GEO under accession [GSE352464](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352464); see `DATA_README.md`
for the inputs.

## Structure

| Folder | Analyses | Manuscript items |
|---|---|---|
| `1_discovery/` | Analysis metadata; differential expression, univariate Cox, random survival forest, LASSO, final 3-miRNA model, risk score | Figure 3, Figure 4A, 4C, 4D; Tables S1–S4 |
| `2_internal_validation/` | Bootstrap optimism correction and selection frequencies; final model check; time-dependent ROC | Figure 4B; Tables S5–S6 |
| `3_clinical_factors/` | Univariate, bivariate and pre-specified multivariable Cox models; power; forest plot | Figure 4E; Table S8 |
| `4_sensitivity/` | 12- and 18-month landmark analyses; stage IV-only and stage-stratified analyses | Table S7; Figure S1 |
| `5_TCGA_exploratory/` | TCGA survival, tumor vs normal, immune deconvolution, single-miRNA expression, smoking | Figure 5A–D; Figures S2, S3, S5, S6; Table S9 |
| `6_pathway/` | Over-representation analysis of validated target genes | Figure 5E–F; Figure S4; Tables S10–S13 |

Scripts are numbered in the order they should be run (e.g. `1.1`, `1.2`, ...).

## How to run

Software: R 4.4.0 (Windows); Python 3 with `python-docx` for `4.4` (table formatting only). Package versions are listed below.

**R scripts (`.R`, `.py`)** are run from the repository root. They read inputs from `data/` and write to `results/`:

```
Rscript 2_internal_validation/2.1_bootstrap_optimism_selection_frequency.R 1000 rounded   # about 25 min
Rscript 4_sensitivity/4.1_landmark_analysis.R 12
Rscript 4_sensitivity/4.1_landmark_analysis.R 18
```

**R Markdown files (`.Rmd`)** read and write files in their working directory. Knit each one with the working
directory set to a folder that contains its inputs:

| Script | Inputs in the working directory | Main outputs used later |
|---|---|---|
| `1.2_discovery_pipeline.Rmd` | `Count_Matrix_removed_unknown.csv`, `Metadata_Limma_dead_time_all.csv` | `Final_Cox_betas_for_risk_score_borda_zscore.csv` (used by 5.4) |
| `5.1_TCGA_download_normal_tumor_miRNA.R` | none (downloads from GDC) | `TCGA_*_miRNA_NormalTumor_SE.rds` |
| `5.2_TCGA_all_cohort.Rmd` | none (downloads from GDC) | `TCGA_*_miRNA_df.rds`, `TCGA_*_clinical.rds`, `TCGA_*_risk_scores_3mir.csv` |
| `5.3_TCGA_non_targeted_cohort.Rmd` | outputs of 5.2; cBioPortal `luad_tcga_pan_can_atlas_2018/` and `lusc_tcga_pan_can_atlas_2018/` (`data_sv.txt`) | `*_nontargeted_cases.csv` |
| `5.4_TCGA_tumor_vs_normal_score.Rmd` | outputs of 1.2, 5.1 and 5.3 | |
| `5.5_TCGA_CIBERSORT.Rmd` | outputs of 5.2 (downloads mRNA-seq from GDC; LM22 from the `IOBR` package) | |
| `5.6_TCGA_single_miRNA_tumor_normal.Rmd`, `5.7_TCGA_smoking.Rmd` | outputs of 5.2 | |
| `6.1_target_gene_pathway_ORA.Rmd` | MSigDB files in `MSigDB/` (see below); targets are queried with `multiMiR` | |

`5.1` and the `5_TCGA` and `6_pathway` scripts can share one working directory, as in the original analysis.

Notes on running:

- The random survival forest in `1.2`, `2.1` and `4.1` runs on a single core (`options(rf.cores = 1)`); with
  multiple cores its results change slightly from run to run. With this setting `1.2` reproduces Table S3
  exactly. The published bootstrap (`2.1`) was run with multiple cores (940 valid resamples, mean optimism
  0.117, optimism-corrected C-index 0.639, selection frequencies 23.7%, 9.6% and 4.6% for miR-377-3p,
  miR-509-3p and miR-378h). A single-core run of the deposited script gives 943 valid resamples, mean optimism
  0.117, optimism-corrected C-index 0.639 and selection frequencies of 23.8%, 9.5% and 4.2%, with identical
  conclusions.
- TCGA data are downloaded from GDC on the first run of `5.1`, `5.2`, `5.3` and `5.5` and cached as `.rds`
  files; later runs reuse the cache. GDC data releases are updated over time, so sample numbers may differ
  slightly from the manuscript if the data are downloaded at a later date. `5.5` and `6.1` also need internet
  access (LM22 signature; multiMiR database).
- `5.3` adds two log2 expression columns to the `TCGA_*_risk_scores_3mir.csv` files written by `5.2` if they
  are missing. In the non-targeted subsets (Figure S3), the score uses miR-509-3p and miR-377-3p only, with
  z-scores recomputed within each subset.
- The TCGA scripts were written before the analysis was framed as exploratory; "validation" in their comments,
  function names and output file names refers to applying the plasma model to TCGA data, not to external
  validation. Several scripts also contain additional exploratory analyses that are not reported in the
  manuscript (for example, treatment response in TCGA); the reported items are listed in the table below.

## Where each result comes from

| Manuscript item | Script |
|---|---|
| Table 1 (patient characteristics) | descriptive statistics of the metadata from `1.1_build_metadata.R` |
| Figure 1 (study design), Figure 2 (EV characterization) | schematic and laboratory data, no code |
| Figure 3A (volcano), Figure 3B and Table S2 (univariate Cox), Figure 3C–E and Table S3 (random survival forest) | `1.2_discovery_pipeline.Rmd` |
| Figure 4A and Table S4 (final multivariable Cox), Figure 4C (risk score distribution), Figure 4D (Kaplan-Meier) | `1.2_discovery_pipeline.Rmd`; coefficients checked by `2.2_final_model_check.R` |
| Table S1 (29 differentially expressed miRNAs) | `1.3_TableS1_differential_expression.R` |
| Table S5 (optimism-corrected C-index), Table S6 (bootstrap selection frequencies) | `2.1_bootstrap_optimism_selection_frequency.R` |
| Figure 4B (time-dependent ROC) | `2.3_Figure4B_timeROC.R` |
| Table S8 (clinical Cox models) | `3.1_clinical_cox_models.R` |
| Minimum detectable hazard ratios | `3.2_power_clinical_factors.R` |
| Figure 4E (forest plot) | `3.3_Figure4E_forest.R` (after 3.1) |
| Table S7, Figure S1 (sensitivity analyses) | `4.1` (run with 12 and 18), `4.2`, then `4.3`; `4.4` formats the Word table and figure |
| Figures S2–S3 (TCGA Kaplan-Meier) | `5.2_TCGA_all_cohort.Rmd`, `5.3_TCGA_non_targeted_cohort.Rmd` |
| Figure 5A–B (risk score, tumor vs normal) | `5.4_TCGA_tumor_vs_normal_score.Rmd` |
| Figure 5C–D, Table S9 (CIBERSORT) | `5.5_TCGA_CIBERSORT.Rmd` |
| Figure S5 (tissue expression of the three miRNAs) | `5.6_TCGA_single_miRNA_tumor_normal.Rmd` |
| Figure S6 (expression and smoking) | `5.7_TCGA_smoking.Rmd` |
| Figure 5E–F, Figure S4, Tables S10–S13 (target-gene over-representation) | `6.1_target_gene_pathway_ORA.Rmd` |

## Note on the bootstrap C-index

In `2.1`, every model is scored on log2 expression standardized with parameters from the batch-corrected
training data, but without applying ComBat to the scored data (ComBat here protects survival status, so
applying it would use outcome information). The apparent C-index in this frame is 0.756; the final
coefficients scored on ComBat-corrected expression give 0.806. All other reported performance metrics use
the batch-corrected score.

## Relation to the scripts that were run

The scripts are copies of the versions used for the manuscript. Only file names, input/output paths and the
labels listed below were changed; no computation was altered.

| Script here | Original file name |
|---|---|
| `1.1_build_metadata.R` | `12_build_rounded_metadata.R` |
| `1.2_discovery_pipeline.Rmd` | `new_34_z-score.Rmd` |
| `1.3_TableS1_differential_expression.R` | `09_A6_regenerate_TableS1.R` |
| `2.1_bootstrap_optimism_selection_frequency.R` | `02_A1_nested_bootstrap.R` |
| `2.2_final_model_check.R` | `05_locked_3miRNA_model.R` |
| `2.3_Figure4B_timeROC.R` | `15_Figure4B_timeROC.R` |
| `3.1_clinical_cox_models.R` | `07_A2_A3_adjusted_clinical_cox.R` |
| `3.2_power_clinical_factors.R` | `13_power_clinical_factors.R` |
| `3.3_Figure4E_forest.R` | `14_Figure4E_forest.R` |
| `4.1_landmark_analysis.R` | `16_R1_landmark12_DE.R` |
| `4.2_stage_analysis.R` | `17_R1_stageIV_sensitivity.R` |
| `4.3_TableS7_FigureS1.R` | `18_R1_sensitivity_table_forest.R` |
| `4.4_format_TableS7_FigureS1.py` | `make_table_R1_sensitivity.py` |
| `5.1_TCGA_download_normal_tumor_miRNA.R` | download section of `new34_TCGA.Rmd` (copied unchanged) |
| `5.2_TCGA_all_cohort.Rmd` | `Script_TCGA_all_cohort.txt` |
| `5.3_TCGA_non_targeted_cohort.Rmd` | `Script_TCGA_non_targeted__cohort.txt` |
| `5.4_TCGA_tumor_vs_normal_score.Rmd` | `norm_cancer_TCGA_zscore.Rmd` |
| `5.5_TCGA_CIBERSORT.Rmd` | `Script_TCGA_CIBERSORT.txt` |
| `5.6_TCGA_single_miRNA_tumor_normal.Rmd` | `single_miRNA_normal_tumor.Rmd` |
| `5.7_TCGA_smoking.Rmd` | `Smoking_status.Rmd` |
| `6.1_target_gene_pathway_ORA.Rmd` | `Script_pathway_of_target_gene.txt` |

Other changes:

1. `.R` scripts: input and output folders renamed to `data/` and `results/`; the default metadata file in
   `2.1` and the metadata path in `1.3` point to the final whole-month file used for the manuscript.
2. `6.1`: an output file name and plot titles read "miR-337-3p", a typographical error for miR-377-3p; target
   genes were always queried for `hsa-miR-377-3p`. The MSigDB file paths point to `MSigDB/`.
3. Scripts stored as `.txt` are R Markdown and are provided as `.Rmd`.
4. Changes needed for the files to run from start to finish in a fresh R session (R 4.4, current package
   versions); none changes a computation, and all reported outputs were re-generated and checked against the
   original results:
   - `1.2`: removed three diagnostic lines that used objects before they were created and a final chunk that
     converted the original file to a script; set the random survival forest to a single core (also in `2.1`
     and `4.1`); added a missing
     `y` aesthetic in an unreported clinical forest plot (required by current `ggplot2`).
   - `5.2`, `5.3`, `5.5`: downloads and the log2-column step run automatically when cached files are absent
     (previously manual switches); `5.5` declares `base` set functions in `conflicted` and loads LM22 in a way
     that works with current `IOBR`.
   - `5.7`: corrected the pack-years column name in a diagnostic chunk and a `print()` call.
   - `1.1`: a comment that referred to an individual patient was reworded.

Superseded script versions are not included.

## Gene set files

`6.1` reads MSigDB v2026.1 human gene set files (GMT, gene symbols), e.g. `h.all.v2026.1.Hs.symbols.gmt`, from
a folder `MSigDB/` in its working directory. Download them from https://www.gsea-msigdb.org/gsea/msigdb
(free registration).

## Licence

MIT (see `LICENSE`).

## R package versions (as installed when the analyses were run or re-run)

| Package | Version |
|---|---|
| BiocManager | 1.30.27 |
| DOSE | 4.0.1 |
| IOBR | 2.2.0 |
| RColorBrewer | 1.1.3 |
| ReactomePA | 1.50.0 |
| STRINGdb | 2.18.0 |
| SummarizedExperiment | 1.36.0 |
| TCGAbiolinks | 2.34.1 |
| UpSetR | 1.4.0 |
| clusterProfiler | 4.14.6 |
| conflicted | 1.2.0 |
| dplyr | 1.2.1 |
| dunn.test | 1.3.7 |
| enrichplot | 1.26.6 |
| forestplot | 3.2.0 |
| ggplot2 | 4.0.3 |
| ggpubr | 0.6.3 |
| ggraph | 2.2.2 |
| ggrepel | 0.9.8 |
| ggtext | 0.1.2 |
| glmnet | 4.1.10 |
| igraph | 2.2.3 |
| knitr | 1.51 |
| limma | 3.62.2 |
| maftools | 2.22.0 |
| multiMiR | 1.28.0 |
| org.Hs.eg.db | 3.20.0 |
| pROC | 1.19.0.1 |
| patchwork | 1.3.2 |
| pheatmap | 1.0.13 |
| plotROC | 2.3.3 |
| randomForestSRC | 3.6.2 |
| readxl | 1.4.5 |
| rms | 8.1.1 |
| scales | 1.4.0 |
| survival | 3.8.6 |
| survminer | 0.5.2 |
| sva | 3.54.0 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| tidyverse | 2.0.0 |
| timeROC | 0.4.1 |
| msigdbr | 26.1.0 |
