# =====================================================================
#  16 - Reviewer 1, point 1: corrected 12-month landmark analysis
#
#  Published DE grouped patients by status at LAST FOLLOW-UP, so survivors
#  followed 4-5 months were pooled with survivors followed 23-28 months.
#
#  Here the DE grouping uses 12-month landmark status:
#     died at <= 12 months                    -> 1
#     alive (or died later) at >= 12 months    -> 0
#     alive with follow-up < 12 months         -> EXCLUDED (status unknown)
#  The time-to-event steps (univariate Cox, RSF, LASSO, final Cox) keep the
#  full OS time and event status, which handle censoring correctly.
#
#  The old Other_version/Fix_time_to12month run is NOT reused: it counted
#  unknown-status patients as survivors.
#
#  Step 0 reproduces the published model with the original grouping to
#  verify the engine before the landmark run.
#  Run from project root:  Rscript 4_sensitivity/4.1_landmark_analysis.R [months]
#  (default 12; e.g. 18 for the 18-month landmark)
# =====================================================================
suppressPackageStartupMessages({
  library(limma); library(sva); library(survival)
  library(randomForestSRC); library(glmnet)
})
options(rf.cores = 1, mc.cores = 1)  # single core: random survival forest results are reproducible
args <- commandArgs(trailingOnly = TRUE)
L    <- if (length(args) >= 1) as.numeric(args[1]) else 12
TAG  <- paste0("landmark", L, "_")
DATA <- "data"
META <- "data/Metadata_FINAL_rounded_34pt.csv"
OUT  <- paste0("results/R1_landmark", L)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PUB  <- c("hsa.miR.509.3p", "hsa.miR.377.3p", "hsa.miR.378h")

counts_in <- read.csv(file.path(DATA, "Count_Matrix_removed_unknown.csv"), row.names = 1)
meta_in   <- read.csv(META, stringsAsFactors = FALSE)
expr_all  <- t(as.matrix(counts_in)); colnames(expr_all) <- trimws(colnames(expr_all))
meta_in$SampleID <- trimws(meta_in$SampleID)
common   <- intersect(colnames(expr_all), meta_in$SampleID)
expr_all <- expr_all[, common, drop = FALSE]
meta_all <- meta_in[match(common, meta_in$SampleID), ]
stopifnot(ncol(expr_all) == 34)

# ---- engine: grp drives filter/ComBat/DE; (time, status) drive survival steps
run <- function(expr, grp, time, status, batch, seed = 42) {
  lg <- log2(expr + 1)
  md <- apply(lg[, grp == 1, drop = FALSE], 1, median)
  ma <- apply(lg[, grp == 0, drop = FALSE], 1, median)
  lg <- lg[(md > 1) | (ma > 1), , drop = FALSE]
  sf <- factor(grp, levels = c(0, 1), labels = c("Alive", "Dead"))
  cb <- ComBat(dat = lg, batch = factor(batch), mod = model.matrix(~ sf),
               par.prior = TRUE, prior.plots = FALSE)
  des <- model.matrix(~ 0 + sf); colnames(des) <- levels(sf)
  fit <- lmFit(cb, des)
  fit <- contrasts.fit(fit, makeContrasts(Dead - Alive, levels = des))
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  de  <- topTable(fit, adjust = "fdr", number = Inf)
  cand <- rownames(de)[de$P.Value < 0.05 & abs(de$logFC) > 0.5]

  zm <- rowMeans(cb[cand, , drop = FALSE]); zs <- apply(cb[cand, , drop = FALSE], 1, sd)
  Z  <- t((cb[cand, , drop = FALSE] - zm) / zs)
  so <- Surv(time, status)
  pv <- vapply(seq_along(cand), function(j)
          summary(coxph(so ~ Z[, j], ties = "efron"))$coefficients[1, 5], numeric(1))
  ok  <- !is.na(pv) & pv < 0.05
  sig <- cand[ok][order(pv[ok])]
  res <- list(de = de, cb = cb, cand = cand, sig = sig, sel = character(0),
              top = character(0), betas = NULL, fin = NULL)
  if (length(sig) < 2) return(res)

  rf_df  <- data.frame(Time.Dead = time, dead.status = status, Z[, sig, drop = FALSE],
                       check.names = TRUE)
  rf_ord <- rf_df[order(rf_df$dead.status, rf_df$Time.Dead), , drop = FALSE]
  vsig   <- setdiff(colnames(rf_df), c("Time.Dead", "dead.status"))
  set.seed(seed)
  rf <- rfsrc(Surv(Time.Dead, dead.status) ~ ., data = rf_ord, ntree = 1000,
              nodesize = 10, mtry = ceiling(sqrt(length(vsig))),
              importance = TRUE, na.action = "na.impute")
  res$leaf_mean <- mean(rf$leaf.count)
  if (res$leaf_mean <= 1) { res$rsf_note <- "RSF degenerate: no tree made any split"; return(res) }
  vi <- rf$importance[vsig]
  set.seed(seed)
  mdep  <- max.subtree(rf, conservative = FALSE)$order[, 1][vsig]
  borda <- (rank(-vi, ties.method = "min") + rank(mdep, ties.method = "min")) / 2
  top   <- vsig[vi > 0 & borda <= 5]
  res$top <- top
  res$rsf <- data.frame(var = vsig, VIMP = vi, min_depth = mdep, avg_rank = borda)
  if (length(top) < 1) return(res)
  sel <- top
  if (length(top) >= 2) {
    set.seed(seed)
    cv <- cv.glmnet(x = as.matrix(rf_df[, top, drop = FALSE]), y = so,
                    family = "cox", alpha = 1, nfolds = 10, maxit = 10000)
    cf <- as.matrix(coef(cv, s = "lambda.1se")); s1 <- rownames(cf)[cf[, 1] != 0]
    if (length(s1) == 0) {
      cf <- as.matrix(coef(cv, s = "lambda.min")); s1 <- rownames(cf)[cf[, 1] != 0]
    }
    if (length(s1) > 0) sel <- s1
  }
  fin <- coxph(so ~ ., data = rf_df[, sel, drop = FALSE], ties = "efron")
  res$sel <- sel; res$betas <- coef(fin); res$fin <- summary(fin)$coefficients
  res
}

# ---- step 0: reproduce published model (original grouping, n = 34)
o <- run(expr_all, meta_all$dead.status, meta_all$Time.Dead, meta_all$dead.status,
         meta_all$Batch)

# ---- landmark coding
tt <- meta_all$Time.Dead; dd <- meta_all$dead.status
lm <- ifelse(dd == 1 & tt <= L, 1L, ifelse(tt >= L, 0L, NA_integer_))
grp_tab <- data.frame(SampleID = meta_all$SampleID, Time.Dead = tt, dead.status = dd,
                      stage = meta_all$stage, landmark = lm)
write.csv(grp_tab, file.path(OUT, paste0(TAG, "group_assignment.csv")), row.names = FALSE)
keep <- !is.na(lm)
m2   <- meta_all[keep, ]; e2 <- expr_all[, keep, drop = FALSE]
lmk  <- run(e2, lm[keep], m2$Time.Dead, m2$dead.status, m2$Batch)

# ---- per-miRNA comparison for the 3 published candidates
uni <- function(cb, ids, time, status) t(sapply(ids, function(id) {
  if (!id %in% rownames(cb)) return(c(HR = NA, lo = NA, hi = NA, p = NA))
  z <- as.numeric(scale(cb[id, ]))
  s <- summary(coxph(Surv(time, status) ~ z, ties = "efron"))
  c(HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4],
    p = s$coefficients[1, 5])
}))
pick <- function(de, ids) t(sapply(ids, function(id)
  if (id %in% rownames(de)) unlist(de[id, c("logFC", "P.Value")])
  else c(logFC = NA, P.Value = NA)))
po <- pick(o$de, PUB); pl <- pick(lmk$de, PUB)
cmp <- data.frame(miRNA = gsub("\\.", "-", PUB),
                  pub_logFC = po[, 1], pub_P = po[, 2],
                  lm_logFC = pl[, 1], lm_P = pl[, 2],
                  lm_DEcandidate = PUB %in% lmk$cand,
                  lm_univarCox_sig = PUB %in% lmk$sig,
                  lm_RSF_top = make.names(PUB) %in% lmk$top,
                  lm_final = make.names(PUB) %in% lmk$sel)
u <- uni(lmk$cb, PUB, m2$Time.Dead, m2$dead.status)
cmp$lm_coxHR_perSD <- u[, "HR"]; cmp$lm_cox_lo <- u[, "lo"]
cmp$lm_cox_hi      <- u[, "hi"]; cmp$lm_cox_p  <- u[, "p"]

# ---- save BEFORE printing
write.csv(lmk$de, file.path(OUT, paste0(TAG, "DE_full.csv")))
write.csv(cmp, file.path(OUT, paste0(TAG, "published_miRNA_comparison.csv")), row.names = FALSE)
if (!is.null(lmk$fin)) write.csv(lmk$fin, file.path(OUT, paste0(TAG, "final_cox.csv")))
if (!is.null(lmk[["rsf"]])) write.csv(lmk[["rsf"]][order(lmk[["rsf"]]$avg_rank), ], file.path(OUT, paste0(TAG, "RSF_ranking.csv")), row.names = FALSE)
summ <- c(
  sprintf("REPRODUCTION (original grouping, n=34): DE cand=%s, univar sig=%s, RSF top=%s, final=%s",
          length(o$cand), length(o$sig), paste(o$top, collapse = ","),
          paste(o$sel, collapse = ",")),
  sprintf("LANDMARK %s months: died<=L=%s, alive>=L=%s, excluded=%s (n analysed=%s)",
          L, sum(lm == 1, na.rm = TRUE), sum(lm == 0, na.rm = TRUE), sum(is.na(lm)), sum(keep)),
  sprintf("Excluded: %s", paste(meta_all$SampleID[is.na(lm)], collapse = ", ")),
  sprintf("Events within analysed n (full OS): %s", sum(m2$dead.status)),
  sprintf("LANDMARK DE cand=%s (overlap with published %s: %s)", length(lmk$cand),
          length(o$cand), length(intersect(lmk$cand, o$cand))),
  sprintf("LANDMARK univar Cox sig=%s: %s", length(lmk$sig), paste(lmk$sig, collapse = ", ")),
  sprintf("RSF mean leaves per tree: reproduction=%s, landmark=%s", round(o$leaf_mean, 2), round(lmk$leaf_mean, 2)),
  sprintf("LANDMARK RSF note: %s", ifelse(is.null(lmk$rsf_note), "-", lmk$rsf_note)),
  sprintf("LANDMARK RSF top: %s", paste(lmk$top, collapse = ", ")),
  sprintf("LANDMARK final model: %s", paste(lmk$sel, collapse = ", ")))
writeLines(summ, file.path(OUT, paste0(TAG, "summary.txt")))
cat(summ, sep = "\n"); cat("\n"); print(cmp, digits = 3)
if (!is.null(lmk$fin)) { cat("\nLandmark final Cox:\n"); print(round(lmk$fin, 4)) }
