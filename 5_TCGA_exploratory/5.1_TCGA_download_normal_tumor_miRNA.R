# Downloads TCGA-LUAD and TCGA-LUSC miRNA-seq (primary tumor + solid tissue normal) from GDC and saves
# TCGA_LUAD_miRNA_NormalTumor_SE.rds / TCGA_LUSC_miRNA_NormalTumor_SE.rds, used by 5.4.
# Extracted unchanged from the original working script (new34_TCGA.Rmd).

library(TCGAbiolinks)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)


# ==============================================================================
# PART 1 — Download Normal + Tumor miRNA-seq (run once, then load from RDS)
# ==============================================================================

download_mirna <- function(project) {
  query <- GDCquery(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling",
    sample.type   = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(query, method = "api", files.per.chunk = 20)
  GDCprepare(query)
}

luad_se <- download_mirna("TCGA-LUAD")
lusc_se <- download_mirna("TCGA-LUSC")

saveRDS(luad_se, "TCGA_LUAD_miRNA_NormalTumor_SE.rds")
saveRDS(lusc_se, "TCGA_LUSC_miRNA_NormalTumor_SE.rds")
