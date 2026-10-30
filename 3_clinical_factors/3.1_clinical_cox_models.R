# =====================================================================
#  A2 + A3 -- Clinical Cox, done properly        [ROUNDED-MONTH VERSION]
#
#  A2. THE PROBLEM. Table S6 / Figure 4E report the risk score as an
#      "independent prognostic factor, HR = 2.72". That HR is
#      2.71828182887954 = e, because the score IS the linear predictor of a
#      Cox model fitted on these same 34 patients: regressing it on itself
#      returns beta = 1.000000000 by construction. And per section 3.6 no
#      clinical covariate cleared the P<0.10 screen, so the "multivariable"
#      model contained the score ALONE -- nothing was adjusted for, despite
#      the abstract claiming independence from five clinical factors.
#      (Reviewer 1, Major comment 3.)
#
#      FIX: enter the score and clinical covariates TOGETHER. With 18 events
#      a 6-covariate model is overfitted (EPV = 3), so the primary analysis
#      is a set of bivariate models (score + one covariate) plus one
#      pre-specified multivariable model.
#
#  A3. ECOG PS and PD-L1 TPS are in Table 1 but were absent from every
#      analysis file. Added here. (Reviewer 1, Minor comment 2.)
#
#  Survival time: ROUNDED MONTHS, exactly as published (author decision).
#  Model: LOCKED to the published 3-miRNA signature with published betas.
#  NOTE: ASCII only -- no special characters, to survive re-encoding.
# =====================================================================

suppressPackageStartupMessages({library(limma); library(sva); library(survival)})

OUT  <- "results/A2_A3_clinical"
META <- "data/Metadata_FINAL_rounded_34pt.csv"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SIG <- c("hsa.miR.509.3p", "hsa.miR.377.3p", "hsa.miR.378h")
PUB <- c(hsa.miR.509.3p = 0.743074, hsa.miR.377.3p = 0.502724,
         hsa.miR.378h  = -0.663259)

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

d <- meta
d$score  <- as.numeric(t(Z) %*% PUB[SIG])
d$scoreZ <- as.numeric(scale(d$score))
so <- Surv(time, status)

d$Age_bin   <- factor(d$Age, levels = c("<60", ">=60"))
d$Sex       <- factor(d$sex, levels = c("Female", "Male"))
d$Smoking   <- factor(ifelse(d$smoking == "never smoker", "Never", "Ever"),
                      levels = c("Never", "Ever"))
d$Stage     <- factor(ifelse(d$stage == 4, "IV", "III"), levels = c("III", "IV"))
d$Histology <- factor(ifelse(d$histo_detail == "Adenocarcinoma", "Adeno", "Other"),
                      levels = c("Other", "Adeno"))
d$ECOG      <- factor(ifelse(d$ECOG_gr == "unknown", NA, d$ECOG_gr),
                      levels = c("0-1", "2"))
# PD-L1: dichotomised at the clinically used 50% cut-off.
# (The 3-level version gives two contrasts -- TPS1-49 vs <1 HR 2.60 p=0.256 and
#  TPS50+ vs <1 HR 4.28 p=0.067, overall LRT p=0.12 -- which cannot be shown as
#  a single forest-plot row. Dichotomising keeps one interpretable estimate.)
pl <- trimws(as.character(d$PDL1_gr))
d$PDL1 <- factor(ifelse(pl == "unknown", NA,
                        ifelse(pl %in% c("<1", "1-49"), "TPS<50", "TPS>=50")),
                 levels = c("TPS<50", "TPS>=50"))

cat("Metadata:", META, "\n")
cat("n =", nrow(d), " events =", sum(status),
    " | ECOG missing:", sum(is.na(d$ECOG)), " PD-L1 missing:", sum(is.na(d$PDL1)), "\n")

# ---------------------------------------------------------------- 1
cat("\n===== 1. The tautology, demonstrated =====\n")
f_alone <- coxph(so ~ score, data = d, ties = "efron")
b <- coef(f_alone)[["score"]]
cat(sprintf("Score alone : beta = %.12f\n", b))
cat(sprintf("              HR   = %.14f\n", exp(b)))
cat(sprintf("exp(1) = e         %.14f\n", exp(1)))
cat(sprintf("difference from e: %.3e\n", abs(exp(b) - exp(1))))
cat("Published Table S6 reported HR = 2.71828182887954, i.e. exactly e.\n")
cat("beta = 1 by construction: the score is this model's own linear\n")
cat("predictor, so regressing it on itself must return 1. No information.\n")

f_sd <- coxph(so ~ scoreZ, data = d, ties = "efron")
cat(sprintf("\nPer 1 SD (SD = %.3f): HR = %.3f (%.2f-%.2f), p = %.3g\n",
            sd(d$score), summary(f_sd)$conf.int[1,1], summary(f_sd)$conf.int[1,3],
            summary(f_sd)$conf.int[1,4], summary(f_sd)$coefficients[1,5]))
cat("  (Table S5 labelled 2.808 'per SD'; that is per 1 unit, not per SD.)\n")

# ---------------------------------------------------------------- 2
covs <- c(Age_bin = "Age (>=60 vs <60)", Sex = "Sex (Male vs Female)",
          Smoking = "Smoking (Ever vs Never)", Stage = "Stage (IV vs III)",
          Histology = "Histology (Adeno vs Other)", ECOG = "ECOG PS (2 vs 0-1)",
          PDL1 = "PD-L1 TPS (>=50 vs <50)")
cat("\n===== 2. Univariate Cox, clinical factors (incl. ECOG, PD-L1) =====\n")
uni <- do.call(rbind, lapply(names(covs), function(v) {
  dd <- d[!is.na(d[[v]]), ]
  s <- summary(coxph(Surv(Time.Dead, dead.status) ~ dd[[v]], data = dd, ties = "efron"))
  data.frame(variable = covs[[v]], n = nrow(dd), events = sum(dd$dead.status),
             HR = round(s$conf.int[1,1], 3), lo = round(s$conf.int[1,3], 3),
             hi = round(s$conf.int[1,4], 3), p = signif(s$coefficients[1,5], 3))
}))
print(uni, row.names = FALSE)

# ---------------------------------------------------------------- 3
cat("\n===== 3. PRIMARY: bivariate models (score + one covariate) =====\n")
biv <- do.call(rbind, lapply(names(covs), function(v) {
  dd <- d[!is.na(d[[v]]), ]
  f <- try(coxph(Surv(Time.Dead, dead.status) ~ scoreZ + dd[[v]], data = dd,
                 ties = "efron"), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  s <- summary(f)
  data.frame(adjusted_for = covs[[v]], n = nrow(dd),
             score_HR = round(s$conf.int["scoreZ",1], 3),
             lo = round(s$conf.int["scoreZ",3], 3),
             hi = round(s$conf.int["scoreZ",4], 3),
             p  = signif(s$coefficients["scoreZ",5], 3))
}))
print(biv, row.names = FALSE)

# ---------------------------------------------------------------- 4
cat("\n===== 4. Pre-specified multivariable model =====\n")
cat("score + ECOG PS + Stage\n\n")
dm <- d[!is.na(d$ECOG), ]
sm <- summary(coxph(Surv(Time.Dead, dead.status) ~ scoreZ + ECOG + Stage,
                    data = dm, ties = "efron"))
mv <- data.frame(term = rownames(sm$conf.int),
                 HR = round(sm$conf.int[,1], 3), lo = round(sm$conf.int[,3], 3),
                 hi = round(sm$conf.int[,4], 3), p = signif(sm$coefficients[,5], 3))
print(mv, row.names = FALSE)
cat(sprintf("\nn = %d, events = %d (EPV = %.1f)\n",
            nrow(dm), sum(dm$dead.status), sum(dm$dead.status)/3))

# ---------------------------------------------------------------- 5
cat("\n===== 5. Signature vs clinical variables (same patients) =====\n")
C <- function(f) as.numeric(summary(f)$concordance["C"])
base <- c("Age_bin","Sex","Smoking","Stage","Histology")
comp <- do.call(rbind, lapply(list(list(base, "All 34 (5 clinical, as published)"),
                                   list(c(base,"ECOG"), "n=31 (+ ECOG PS)")),
  function(z) {
    cl <- z[[1]]; dd <- d[complete.cases(d[, cl]), ]
    dd$scoreZ <- as.numeric(scale(dd$score))
    fc <- coxph(as.formula(paste("Surv(Time.Dead,dead.status)~", paste(cl, collapse="+"))),
                data = dd, ties = "efron")
    fs <- coxph(Surv(Time.Dead, dead.status) ~ scoreZ, data = dd, ties = "efron")
    fb <- coxph(as.formula(paste("Surv(Time.Dead,dead.status)~scoreZ+", paste(cl, collapse="+"))),
                data = dd, ties = "efron")
    lr <- anova(fc, fb)
    data.frame(cohort = z[[2]], n = nrow(dd), events = sum(dd$dead.status),
               C_clinical = round(C(fc),4), C_signature = round(C(fs),4),
               C_combined = round(C(fb),4), LRT_p = signif(lr$`Pr(>|Chi|)`[2], 3))
  }))
print(comp, row.names = FALSE)

write.csv(uni,  file.path(OUT, "univariate_clinical.csv"), row.names = FALSE)
write.csv(biv,  file.path(OUT, "bivariate_score_adjusted.csv"), row.names = FALSE)
write.csv(mv,   file.path(OUT, "multivariable_prespecified.csv"), row.names = FALSE)
write.csv(comp, file.path(OUT, "model_comparison_cindex.csv"), row.names = FALSE)
cat("\nSaved to", OUT, "\n")

