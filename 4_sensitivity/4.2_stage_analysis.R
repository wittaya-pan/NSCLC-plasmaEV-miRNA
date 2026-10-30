# =====================================================================
#  17 - Reviewer 1, point 2: Stage IV-only sensitivity analysis
#
#  The cohort pools Stage III (n = 4) and Stage IV (n = 30). Reviewer 1 asks
#  whether model performance is consistent when Stage III is removed, or
#  when stage is accounted for.
#
#  The 3-miRNA model is applied LOCKED, exactly as in 2.2_final_model_check.R:
#  ComBat and z-scoring on the full 34-patient cohort, published Table S4
#  betas, no re-selection and no re-fitting of the score. Performance is then
#  evaluated in the Stage IV subset, and stage is handled in the full cohort
#  by stratification.
#
#  All C-indices here are APPARENT (the signature was derived in these
#  patients); interpret alongside the optimism from 2.1_bootstrap_optimism_selection_frequency.R.
#
#  Run from project root:  Rscript 4_sensitivity/4.2_stage_analysis.R
#  ASCII only.
# =====================================================================
suppressPackageStartupMessages({
  library(limma); library(sva); library(survival); library(timeROC)
})
OUT  <- "results/R1_stageIV"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
META <- "data/Metadata_FINAL_rounded_34pt.csv"
SIG  <- c("hsa.miR.509.3p", "hsa.miR.377.3p", "hsa.miR.378h")
PUB  <- c(hsa.miR.509.3p = 0.743074, hsa.miR.377.3p = 0.502724,
          hsa.miR.378h  = -0.663259)

# ---- locked score, built exactly as script 05 -------------------------
counts <- read.csv("data/Count_Matrix_removed_unknown.csv", row.names = 1)
expr   <- t(as.matrix(counts)); colnames(expr) <- trimws(colnames(expr))
meta   <- read.csv(META, stringsAsFactors = FALSE)
meta   <- meta[match(colnames(expr), trimws(meta$SampleID)), ]
status <- as.integer(meta$dead.status); time <- as.numeric(meta$Time.Dead)
stage  <- as.integer(meta$stage)
sf <- factor(status, levels = c(0, 1), labels = c("Alive", "Dead"))
lg <- log2(expr + 1)
lg <- lg[(apply(lg[, status == 1], 1, median) > 1) |
         (apply(lg[, status == 0], 1, median) > 1), ]
cb <- ComBat(lg, batch = factor(meta$Batch), mod = model.matrix(~ sf),
             par.prior = TRUE, prior.plots = FALSE)
Z  <- t(scale(t(cb)))[SIG, , drop = FALSE]
d  <- as.data.frame(t(Z))
sc     <- as.numeric(as.matrix(d[, SIG]) %*% PUB[SIG])
sc_sd  <- (sc - mean(sc)) / sd(sc)            # per-SD unit fixed on the full cohort
cut    <- median(sc)                           # locked full-cohort median cut-off
high   <- sc > cut

cidx <- function(t, s, x) as.numeric(concordance(Surv(t, s) ~ x, reverse = TRUE)$concordance)
hr   <- function(fit, term) { s <- summary(fit)
  c(HR = s$conf.int[term, 1], lo = s$conf.int[term, 3], hi = s$conf.int[term, 4],
    p = s$coefficients[term, "Pr(>|z|)"]) }
lr_p <- function(t, s, g) 1 - pchisq(survdiff(Surv(t, s) ~ g)$chisq, length(unique(g)) - 1)
auc_tab <- function(t, s, x, times) {
  tr <- timeROC(T = t, delta = s, marker = x, cause = 1, times = times, iid = TRUE)
  se <- tr$inference$vect_sd_1
  do.call(rbind, lapply(seq_along(times), function(i) {
    a <- tr$AUC[i]; z <- 1.96 * se[i] / (a * (1 - a))
    data.frame(months = times[i], AUC = a,
               logit_lo = plogis(qlogis(a) - z), logit_hi = plogis(qlogis(a) + z))
  }))
}

summarise <- function(idx, label) {
  t <- time[idx]; s <- status[idx]; x <- sc[idx]; xs <- sc_sd[idx]; g <- high[idx]
  f_sd <- coxph(Surv(t, s) ~ xs, ties = "efron")
  f_g  <- coxph(Surv(t, s) ~ g,  ties = "efron")
  h1 <- hr(f_sd, "xs"); h2 <- hr(f_g, "gTRUE")
  data.frame(cohort = label, n = sum(idx), events = sum(s),
             n_high = sum(g), n_low = sum(!g),
             C_index_apparent = cidx(t, s, x),
             HR_perSD = h1["HR"], HR_perSD_lo = h1["lo"], HR_perSD_hi = h1["hi"], HR_perSD_p = h1["p"],
             HR_high_vs_low = h2["HR"], HR_hl_lo = h2["lo"], HR_hl_hi = h2["hi"],
             logrank_p = lr_p(t, s, g), row.names = NULL)
}

all_idx <- rep(TRUE, length(sc)); iv_idx <- stage == 4
perf <- rbind(summarise(all_idx, "All patients (published)"),
              summarise(iv_idx,  "Stage IV only"))

# ---- stage handled in the full cohort -----------------------------------
f_strat <- coxph(Surv(time, status) ~ sc_sd + strata(stage), ties = "efron")
f_adj   <- coxph(Surv(time, status) ~ sc_sd + factor(stage), ties = "efron")
strat <- rbind(
  data.frame(model = "Stratified by stage (III / IV)", t(hr(f_strat, "sc_sd")), row.names = NULL),
  data.frame(model = "Adjusted for stage (covariate)", t(hr(f_adj, "sc_sd")), row.names = NULL))

# ---- time-dependent AUC --------------------------------------------------
auc <- rbind(cbind(cohort = "All patients", auc_tab(time, status, sc, c(6, 12, 18))),
             cbind(cohort = "Stage IV only", auc_tab(time[iv_idx], status[iv_idx], sc[iv_idx], c(6, 12, 18))))

# ---- per-miRNA and refitted 3-miRNA Cox within Stage IV ------------------
per_mir <- do.call(rbind, lapply(SIG, function(m) {
  f <- coxph(Surv(time[iv_idx], status[iv_idx]) ~ d[iv_idx, m], ties = "efron")
  s <- summary(f)
  data.frame(miRNA = gsub("\\.", "-", m), HR_perSD = s$conf.int[1, 1],
             lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5])
}))
refit <- summary(coxph(Surv(time[iv_idx], status[iv_idx]) ~ ., data = d[iv_idx, SIG], ties = "efron"))
refit_tab <- data.frame(miRNA = gsub("\\.", "-", SIG),
                        published_beta = PUB[SIG],
                        stageIV_beta = refit$coefficients[SIG, "coef"],
                        stageIV_p = refit$coefficients[SIG, "Pr(>|z|)"], row.names = NULL)

# ---- Stage III patients, descriptive only --------------------------------
st3 <- data.frame(SampleID = meta$SampleID[stage == 3], Time.Dead = time[stage == 3],
                  dead.status = status[stage == 3], risk_group = ifelse(high[stage == 3], "High", "Low"),
                  score = round(sc[stage == 3], 3))

# ---- save BEFORE printing ------------------------------------------------
write.csv(perf,      file.path(OUT, "stageIV_performance.csv"), row.names = FALSE)
write.csv(strat,     file.path(OUT, "stage_stratified_cox.csv"), row.names = FALSE)
write.csv(auc,       file.path(OUT, "stageIV_timeROC_logitCI.csv"), row.names = FALSE)
write.csv(per_mir,   file.path(OUT, "stageIV_per_miRNA_cox.csv"), row.names = FALSE)
write.csv(refit_tab, file.path(OUT, "stageIV_refitted_betas.csv"), row.names = FALSE)
write.csv(st3,       file.path(OUT, "stageIII_patients.csv"), row.names = FALSE)

options(width = 140)
cat("=== Performance of the locked score ===\n");          print(perf, digits = 3)
cat("\n=== Stage accounted for in the full cohort ===\n");  print(strat, digits = 3)
cat("\n=== Time-dependent AUC (logit-bounded CI) ===\n");   print(auc, digits = 3, row.names = FALSE)
cat("\n=== Per-miRNA Cox, Stage IV only ===\n");            print(per_mir, digits = 3, row.names = FALSE)
cat("\n=== Refitted 3-miRNA Cox, Stage IV only ===\n");     print(refit_tab, digits = 3, row.names = FALSE)
cat("\n=== Stage III patients ===\n");                      print(st3, row.names = FALSE)
