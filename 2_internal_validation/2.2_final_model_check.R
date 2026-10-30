# =====================================================================
#  Locked 3-miRNA model on ROUNDED months (as published)
#
#  DECISION (author): retain the published three-miRNA signature and the
#  published rounded-month time coding. This script does NOT re-select
#  variables -- the signature is treated as already derived and locked.
#
#  Purpose:
#    1. confirm the refitted coefficients reproduce published Table S4
#    2. apparent discrimination (C-index, log-rank)
#    3. time-dependent AUC with LOGIT-BOUNDED CIs -- fixes the published
#       6-month interval whose upper limit was 1.006 (Reviewer 1, Minor 3)
#
#  The apparent C-index here is NOT an honest performance estimate on its
#  own: the signature was chosen using these same 34 patients. Pair it with
#  the optimism from 2.1_bootstrap_optimism_selection_frequency.R.
#
#  ASCII only -- no special characters, to survive re-encoding.
# =====================================================================

suppressPackageStartupMessages({
  library(limma); library(sva); library(survival); library(timeROC)
})

OUT  <- "results/A1_nested_bootstrap"
META <- "data/Metadata_FINAL_rounded_34pt.csv"
SIG  <- c("hsa.miR.509.3p", "hsa.miR.377.3p", "hsa.miR.378h")
PUB  <- c(hsa.miR.509.3p = 0.743074, hsa.miR.377.3p = 0.502724,
          hsa.miR.378h  = -0.663259)      # published Table S4 betas

counts <- read.csv("data/Count_Matrix_removed_unknown.csv", row.names = 1)
expr   <- t(as.matrix(counts)); colnames(expr) <- trimws(colnames(expr))
meta   <- read.csv(META, stringsAsFactors = FALSE)
meta   <- meta[match(colnames(expr), trimws(meta$SampleID)), ]

status <- as.integer(meta$dead.status); time <- as.numeric(meta$Time.Dead)
sf <- factor(status, levels = c(0, 1), labels = c("Alive", "Dead"))
lg <- log2(expr + 1)
lg <- lg[(apply(lg[, status == 1], 1, median) > 1) |
         (apply(lg[, status == 0], 1, median) > 1), ]
cb <- ComBat(lg, batch = factor(meta$Batch), mod = model.matrix(~ sf),
             par.prior = TRUE, prior.plots = FALSE)
Z  <- t(scale(t(cb)))[SIG, , drop = FALSE]
d  <- as.data.frame(t(Z))
so <- Surv(time, status)

cat("Metadata:", META, "\n")
cat(sprintf("n = %d, events = %d\n\n", nrow(d), sum(status)))

# ---- 1. refitted coefficients vs published Table S4 -----------------
fit <- coxph(so ~ ., data = d, ties = "efron"); s <- summary(fit)
res <- data.frame(miRNA = gsub("\\.", "-", SIG),
                  coef = round(s$coefficients[SIG, "coef"], 6),
                  published_coef = round(PUB[SIG], 6),
                  HR = round(s$conf.int[SIG, "exp(coef)"], 3),
                  lo = round(s$conf.int[SIG, "lower .95"], 3),
                  hi = round(s$conf.int[SIG, "upper .95"], 3),
                  p  = signif(s$coefficients[SIG, "Pr(>|z|)"], 3))
cat("=== Table S4 check: refitted vs published ===\n")
print(res, row.names = FALSE)
cat(sprintf("\nmax abs difference from published betas: %.2e  -> Table S4 needs NO edit\n",
            max(abs(s$coefficients[SIG, "coef"] - PUB[SIG]))))

sc <- as.numeric(as.matrix(d[, SIG]) %*% PUB[SIG])
cidx <- function(x) as.numeric(concordance(so ~ x, reverse = TRUE)$concordance)
cat(sprintf("\nApparent C-index (published betas): %.4f\n", cidx(sc)))
g <- sc > median(sc)
cat(sprintf("Median split: log-rank p = %.4g | HR high vs low = %.2f (%.2f-%.2f)\n",
            1 - pchisq(survdiff(so ~ g)$chisq, 1),
            summary(coxph(so ~ g))$conf.int[1,1],
            summary(coxph(so ~ g))$conf.int[1,3],
            summary(coxph(so ~ g))$conf.int[1,4]))

# ---- 2. time-dependent AUC, naive vs logit-bounded CI ---------------
tr <- timeROC(T = time, delta = status, marker = sc, cause = 1,
              times = c(6, 12, 18), iid = TRUE)
se <- tr$inference$vect_sd_1
cat("\n=== Time-dependent AUC (Reviewer 1, Minor 3) ===\n")
auc <- do.call(rbind, lapply(1:3, function(i) {
  a <- tr$AUC[i]; z <- 1.96 * se[i] / (a * (1 - a))
  data.frame(months = c(6, 12, 18)[i], AUC = round(a, 3),
             naive_lo = round(a - 1.96*se[i], 3), naive_hi = round(a + 1.96*se[i], 3),
             logit_lo = round(plogis(qlogis(a) - z), 3),
             logit_hi = round(plogis(qlogis(a) + z), 3))
}))
auc$naive_exceeds_1 <- ifelse(auc$naive_hi > 1, "YES <-- problem", "no")
print(auc, row.names = FALSE)
cat("\nUse the logit columns: all bounded within [0,1].\n")

write.csv(res, file.path(OUT, "locked_3miRNA_coefficients_rounded.csv"), row.names = FALSE)
write.csv(auc, file.path(OUT, "timeROC_logit_CI_rounded.csv"), row.names = FALSE)
cat("\nSaved.\n")
