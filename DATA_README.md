# Data needed to run this code

No data files are included in this repository.

| Input | Place it at | Source |
|---|---|---|
| NanoString nCounter miRNA counts, nSolver-normalized (798 miRNAs x 34 patients) | `data/Count_Matrix_removed_unknown.csv` (R scripts) and `Count_Matrix_removed_unknown.csv` in the working directory of `1.2` | GEO **GSE352464** (series matrix; raw RCC files are also deposited, one per sample) |
| De-identified clinical metadata | `data/Metadata_FINAL_rounded_34pt.csv` (R scripts) and `Metadata_Limma_dead_time_all.csv` in the working directory of `1.2` | GEO GSE352464 sample characteristics; see below |
| TCGA-LUAD and TCGA-LUSC miRNA-seq, RNA-seq, mutation and clinical data | downloaded automatically by `5.1`, `5.2`, `5.3` and `5.5` on first run | public, via `TCGAbiolinks` (GDC) |
| TCGA PanCancer Atlas structural variant files | `luad_tcga_pan_can_atlas_2018/data_sv.txt` and `lusc_tcga_pan_can_atlas_2018/data_sv.txt` in the TCGA working directory | public, cBioPortal (https://www.cbioportal.org/datasets) |
| CIBERSORT LM22 signature | downloaded by `5.5` from the IOBR GitHub repository on first run | IOBR |
| Validated miRNA targets | queried by `6.1` | `multiMiR` (online database) |
| Gene set collections | `MSigDB/*.gmt` in the working directory of `6.1` | MSigDB v2026.1, https://www.gsea-msigdb.org |

## File formats

**Count matrix** (`Count_Matrix_removed_unknown.csv`): comma-separated, **one row per patient**, first column
`SampleID` (pt01, pt02, ...), then one column per miRNA named as in the deposited file (e.g. `hsa-let-7a-5p`).
The GEO series matrix (`GSE352464_series_matrix.txt.gz`, or `GEOquery::getGEO("GSE352464")`) has miRNAs in rows
(`ID_REF`) and one column per sample (GSM accession; sample title `NSCLC_EV_pt01` etc.). Transpose it and use the
sample codes (characteristic `patient code`: pt01, pt02, ...) as `SampleID`.

**Metadata** (one row per patient; `Metadata_Limma_dead_time_all.csv` for `1.2` needs only the columns from
`SampleID` to `Age`, the other file needs all of them). Build it from the GEO sample characteristics
(`pData()` of the GEOquery object):

| Column | Values | GEO characteristic |
|---|---|---|
| `SampleID` | pt01, pt02, ... | `patient code` |
| `dead.status` | 1 = deceased, 0 = alive at last follow-up | `vital status` (deceased / alive) |
| `Time.Dead` | overall survival in whole months | `overall survival (months)` |
| `Batch` | 1, 2, 3 (nCounter cartridge) | `batch` (nCounter cartridge 1 to 3) |
| `sex` | `Male`, `Female` | `Sex` (male / female) |
| `Age` | `<60`, `>=60` | `age group` (<60 years / >=60 years) |
| `smoking` | `never smoker`, `former smoker`, `current smoker` | `smoking status` |
| `histo_detail` | e.g. `Adenocarcinoma`, `Squamous cell carcinoma` | `histology` (lower case) |
| `stage` | 3 or 4 | first part of `stage (AJCC 8th edition)` (III or IV) |
| `ajcc8_stage` | e.g. `IIIB`, `IVA` | `stage (AJCC 8th edition)` |
| `ECOG_gr` | `0-1`, `2`, `unknown` | `ECOG performance status` (0, 1, 2, not available) |
| `PDL1_gr` | `<1`, `1-49`, `≥50`, `unknown` | `PD-L1 tumor proportion score` (<1%, 1-49%, >=50%, not available) |
| `firstline_category` | `platinum-based doublet chemotherapy`, `immunotherapy`, `single-agent chemotherapy`, `chemotherapy followed by immunotherapy` (used by `1.2`) | `first-line treatment category` |

`1.2` also uses `T`, `N`, `M` and `firstline_detail` (regimen) in descriptive clinical summaries that are not
reported in the manuscript; these are not shared, and those summaries cannot be reproduced from the public data.
The treatment grouping in `1.2` (Chemo / Immuno / Chemo+Immuno) uses `firstline_category` only.
In the metadata file used by `1.2`, the age group of one patient was coded `>=60` instead of `<60`; this affects
only the unreported clinical summaries in `1.2`. The corrected coding (`1.1`, GEO) is used for all
reported clinical analyses (`3.1`, `3.2`).

**Note on clinical metadata.** `1.1_build_metadata.R` builds the analysis metadata from the hospital clinical
record, which contains dates and is not shared (ethics approval COA.MURA2021/956). The GEO sample
characteristics provide the variables above for each patient. Exact ages, TNM categories and treatment regimens
are not shared.
