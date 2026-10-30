# =====================================================================
#  A6 — Regenerate Table S1
#
#  PROBLEM: the published Table S1 lists 31 miRNAs, but it was exported
#  with the P < 0.05 filter ONLY -- the |logFC| > 0.5 criterion stated in
#  Methods was never applied. Two rows fail it:
#      hsa-miR-30e-3p   logFC = -0.4632
#      hsa-miR-513b-5p  logFC = -0.3915
#  (hsa-miR-514a-3p at -0.5047 passes; it only looks borderline at 2 dp.)
#
#  Correct count: 29 miRNAs = 22 up + 7 down.
#  The manuscript states 30 (abstract, section 3.4) and 31 with "23 up and
#  8 down" (section 3.3). All three figures are wrong. (Reviewer 1, Minor 1.)
#
#  NOTE: differential expression contrasts Dead vs Alive using dead.status
#  only -- survival TIME is not used -- so Table S1 is identical under the
#  floored and exact time codings. The A4 correction does not affect it.
#
#  Downstream impact: none. Neither failing miRNA reached univariate Cox
#  significance, so the 15 -> RSF -> LASSO chain is unchanged.
# =====================================================================

suppressPackageStartupMessages({library(limma); library(sva); library(dplyr)})

OUT <- "results/A6_DE"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

counts <- read.csv("data/Count_Matrix_removed_unknown.csv", row.names = 1)
meta   <- read.csv("data/Metadata_FINAL_rounded_34pt.csv",
                   stringsAsFactors = FALSE)
expr   <- t(as.matrix(counts)); colnames(expr) <- trimws(colnames(expr))
meta   <- meta[match(colnames(expr), trimws(meta$SampleID)), ]

status <- as.integer(meta$dead.status)
sf <- factor(status, levels = c(0, 1), labels = c("Alive", "Dead"))
lg <- log2(expr + 1)
lg <- lg[(apply(lg[, status == 1], 1, median) > 1) |
         (apply(lg[, status == 0], 1, median) > 1), ]
cb <- ComBat(lg, batch = factor(meta$Batch), mod = model.matrix(~ sf),
             par.prior = TRUE, prior.plots = FALSE)
des <- model.matrix(~ 0 + sf); colnames(des) <- levels(sf)
fit <- eBayes(contrasts.fit(lmFit(cb, des), makeContrasts(Dead - Alive, levels = des)),
              trend = TRUE, robust = TRUE)
de  <- topTable(fit, adjust = "fdr", number = Inf)
de$miRNA <- gsub("\\.", "-", rownames(de))

pub <- de[de$P.Value < 0.05, ]                                   # what was exported: 31
cor <- de[de$P.Value < 0.05 & abs(de$logFC) > 0.5, ]             # Methods rule:     29

cat("=================== A6: Table S1 ===================\n")
cat(sprintf("P < 0.05 only                    : %d  <- what the published Table S1 lists\n", nrow(pub)))
cat(sprintf("P < 0.05 AND |logFC| > 0.5       : %d  <- the stated Methods criterion\n", nrow(cor)))
cat(sprintf("  up-regulated   (logFC > 0)     : %d\n", sum(cor$logFC > 0)))
cat(sprintf("  down-regulated (logFC < 0)     : %d\n", sum(cor$logFC < 0)))

cat("\n--- rows removed by applying the stated filter ---\n")
drop <- pub[!(rownames(pub) %in% rownames(cor)), c("miRNA", "logFC", "P.Value")]
print(data.frame(miRNA = drop$miRNA, logFC = round(drop$logFC, 4),
                 P = signif(drop$P.Value, 3)), row.names = FALSE)

cat("\n--- manuscript statements vs truth ---\n")
cat(sprintf("  Abstract      : 'Thirty'            -> %d\n", nrow(cor)))
cat(sprintf("  Section 3.3   : '31 ... 23 up, 8 down' -> %d ... %d up, %d down\n",
            nrow(cor), sum(cor$logFC > 0), sum(cor$logFC < 0)))
cat(sprintf("  Section 3.4   : 'the 30'            -> %d\n", nrow(cor)))

cat("\n--- downstream check ---\n")
sig15 <- c("hsa-miR-377-3p","hsa-miR-509-3p","hsa-miR-1268b","hsa-miR-1276",
           "hsa-miR-548ak","hsa-miR-30b-5p","hsa-miR-378h","hsa-miR-575",
           "hsa-miR-541-3p","hsa-miR-585-3p","hsa-miR-18b-5p","hsa-miR-508-5p",
           "hsa-miR-514a-3p","hsa-miR-620","hsa-miR-487b-5p")
cat(sprintf("  all 15 univariate-significant miRNAs still present: %s\n",
            all(sig15 %in% cor$miRNA)))
cat(sprintf("  either removed miRNA among the 15: %s\n",
            any(drop$miRNA %in% sig15)))
cat("  -> selection chain (29 -> 15 -> RSF -> LASSO) unchanged\n")

out <- cor %>%
  arrange(P.Value) %>%
  transmute(miRNA,
            logFC   = round(logFC, 2),
            P.Value = round(P.Value, 3),
            adj.P.Val = signif(adj.P.Val, 3),
            Direction = ifelse(logFC > 0, "Up", "Down"))
write.csv(out, file.path(OUT, "Table_S1_corrected_29miRNA.csv"), row.names = FALSE)
write.csv(de[, c("miRNA","logFC","AveExpr","t","P.Value","adj.P.Val","B")],
          file.path(OUT, "DE_full_results.csv"), row.names = FALSE)
cat(sprintf("\nSaved Table_S1_corrected_29miRNA.csv (%d rows)\n", nrow(out)))
print(head(out, 8), row.names = FALSE)
