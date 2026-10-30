# =====================================================================
#  A1 — Nested bootstrap optimism correction
#
#  WHY: the original Script_boost_strap.txt resampled patients but passed
#  the FIXED final betas into every replicate, so the only thing varying
#  was the z-scoring. It never re-ran DE / univariate Cox / RSF / LASSO,
#  which is where nearly all the optimism lives. Reported optimism was
#  0.0037 -- effectively zero, and not a real estimate.
#
#  THIS SCRIPT re-executes the ENTIRE selection pipeline inside every
#  bootstrap replicate:
#     log2 -> expression filter -> ComBat -> limma DE -> candidate filter
#     -> z-score -> univariate Cox -> RSF (VIMP + minimal depth, Borda)
#     -> LASSO Cox (lambda.1se) -> unpenalised multivariate Cox
#
#  Harrell optimism:  optimism = mean(C_boot - C_orig)
#                     C_corrected = C_apparent - optimism
#
#  Usage:  Rscript 2_internal_validation/2.1_bootstrap_optimism_selection_frequency.R [B] [tag]
# =====================================================================

suppressPackageStartupMessages({
  library(limma); library(sva); library(survival)
  library(randomForestSRC); library(glmnet); library(dplyr)
})
options(rf.cores = 1, mc.cores = 1)  # single core: random survival forest results are reproducible

args <- commandArgs(trailingOnly = TRUE)
B    <- if (length(args) >= 1) as.integer(args[1]) else 1000
TAG  <- if (length(args) >= 2) args[2] else "full"
META <- if (length(args) >= 3) args[3] else
          "data/Metadata_FINAL_rounded_34pt.csv"

DATA <- "data"
OUT  <- "results/A1_nested_bootstrap"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- inputs
# NOTE: read exactly as the published script does -- no check.names = FALSE,
# so "hsa-miR-377-3p" becomes "hsa.miR.377.3p", matching their saved files.
counts_in <- read.csv(file.path(DATA, "Count_Matrix_removed_unknown.csv"),
                      row.names = 1)
meta_in   <- read.csv(META, stringsAsFactors = FALSE)
if (!"Time.Dead" %in% names(meta_in)) stop("metadata lacks Time.Dead")
cat("Metadata:", META, "\n")

expr_all <- t(as.matrix(counts_in))          # miRNAs x samples
colnames(expr_all) <- trimws(colnames(expr_all))
meta_in$SampleID   <- trimws(meta_in$SampleID)

common   <- intersect(colnames(expr_all), meta_in$SampleID)
expr_all <- expr_all[, common, drop = FALSE]
meta_all <- meta_in[match(common, meta_in$SampleID), ]

stopifnot(nrow(meta_all) == 34, ncol(expr_all) == 34)
cat(sprintf("Loaded: %d miRNAs x %d samples | events = %d\n",
            nrow(expr_all), ncol(expr_all), sum(meta_all$dead.status)))

# ------------------------------------------------------------- utilities
# Direction-fixed C-index: higher risk score = shorter survival.
# (Deliberately NOT re-fitting a coxph to get C -- refitting lets the
#  coefficient flip sign on a bad model and silently reports C > 0.5,
#  which would understate optimism.)
cindex <- function(time, status, score) {
  if (length(unique(score)) < 2) return(NA_real_)
  cc <- tryCatch(concordance(Surv(time, status) ~ score, reverse = TRUE),
                 error = function(e) NULL)
  if (is.null(cc)) NA_real_ else as.numeric(cc$concordance)
}

# ---------------------------------------------------------------- engine
# Runs the complete selection pipeline on one dataset.
# Returns betas, the z-scoring parameters, and the selected miRNAs -- or
# NULL with a reason if any stage collapses on this resample.
run_pipeline <- function(expr, meta, seed) {

  fail <- function(w) list(ok = FALSE, why = w)

  status <- as.integer(meta$dead.status)
  time   <- as.numeric(meta$Time.Dead)
  batch  <- factor(meta$Batch)

  if (length(unique(status)) < 2)      return(fail("one_status_only"))
  if (sum(status) < 3)                 return(fail("too_few_events"))
  if (any(table(batch) < 2))           return(fail("degenerate_batch"))
  if (any(table(batch, status) == 0))  return(fail("batch_status_confounded"))

  # duplicated sample names break ComBat/limma column handling
  colnames(expr) <- paste0("s", seq_len(ncol(expr)))

  # 1-2. log2 + expression filter (median log2 > 1 in dead OR alive)
  lg <- log2(expr + 1)
  md <- apply(lg[, status == 1, drop = FALSE], 1, median)
  ma <- apply(lg[, status == 0, drop = FALSE], 1, median)
  lg <- lg[(md > 1) | (ma > 1), , drop = FALSE]
  if (nrow(lg) < 10) return(fail("too_few_mirnas"))

  # 3. ComBat, protecting survival status
  sf  <- factor(status, levels = c(0, 1), labels = c("Alive", "Dead"))
  cb  <- tryCatch(
    ComBat(dat = lg, batch = batch, mod = model.matrix(~ sf),
           par.prior = TRUE, prior.plots = FALSE),
    error = function(e) NULL)
  if (is.null(cb)) return(fail("combat_failed"))

  # 4. limma DE, Dead vs Alive
  des <- model.matrix(~ 0 + sf); colnames(des) <- levels(sf)
  de <- tryCatch({
    fit <- lmFit(cb, des)
    fit <- contrasts.fit(fit, makeContrasts(Dead - Alive, levels = des))
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    topTable(fit, adjust = "fdr", number = Inf)
  }, error = function(e) NULL)
  if (is.null(de)) return(fail("limma_failed"))

  # 5. candidates
  cand <- rownames(de)[de$P.Value < 0.05 & abs(de$logFC) > 0.5]
  if (length(cand) < 2) return(fail("no_de_candidates"))

  # 6. z-score across samples (keep params to apply to the original data)
  zm <- rowMeans(cb[cand, , drop = FALSE])
  zs <- apply(cb[cand, , drop = FALSE], 1, sd)
  keep <- zs > 0; cand <- cand[keep]; zm <- zm[keep]; zs <- zs[keep]
  if (length(cand) < 2) return(fail("no_variance"))
  Z <- t((cb[cand, , drop = FALSE] - zm) / zs)     # samples x miRNAs

  # 7. univariate Cox screen
  so <- Surv(time, status)
  pv <- vapply(seq_along(cand), function(j) {
    f <- tryCatch(coxph(so ~ Z[, j], ties = "efron"), error = function(e) NULL)
    if (is.null(f)) NA_real_ else summary(f)$coefficients[1, "Pr(>|z|)"]
  }, numeric(1))
  # Order by univariate p-value, matching the column order of the published
  # Merged_sigCox_zscore_with_metadata.csv. This matters: rfsrc() samples
  # mtry variables per split, so column order changes the forest for a
  # given seed. Verified identical to their ordering on the full data.
  ok  <- !is.na(pv) & pv < 0.05
  sig <- cand[ok][order(pv[ok])]
  if (length(sig) < 2) return(fail("no_univariate_sig"))

  # 8-9. RSF -> VIMP + minimal depth -> Borda score < 6, VIMP > 0
  rf_df <- data.frame(Time.Dead = time, dead.status = status,
                      Z[, sig, drop = FALSE], check.names = TRUE)
  # Row order matters to rfsrc for a given seed. The published RSF input
  # (Merged_sigCox_zscore_with_metadata.csv) is ordered by status then
  # survival time; reproduce that convention so the apparent model matches
  # the published signature. Verified: identical data in SampleID order
  # instead selects miR-509-3p in only 3% of seeds vs 100% in this order.
  # NOTE: applies to the RSF call ONLY -- `so`, LASSO and the final Cox all
  # use rf_df in its original row order, or rows would desync from `so`.
  rf_ord <- rf_df[order(rf_df$dead.status, rf_df$Time.Dead), , drop = FALSE]
  vsig  <- setdiff(colnames(rf_df), c("Time.Dead", "dead.status"))
  set.seed(seed)
  rf <- tryCatch(
    rfsrc(Surv(Time.Dead, dead.status) ~ ., data = rf_ord, ntree = 1000,
          nodesize = 10, mtry = ceiling(sqrt(length(vsig))),
          importance = TRUE, na.action = "na.impute"),
    error = function(e) NULL)
  if (is.null(rf)) return(fail("rsf_failed"))

  vi <- rf$importance[vsig]
  set.seed(seed)          # published script re-seeds before max.subtree
  ms <- tryCatch(max.subtree(rf, conservative = FALSE), error = function(e) NULL)
  if (is.null(ms)) return(fail("maxsubtree_failed"))
  mdep <- ms$order[, 1][vsig]

  # Selection rule (author-confirmed): average of the VIMP rank and the
  # minimal-depth rank <= 5, with positive VIMP.
  #   Note: the original script coded this as `Borda_Score < 6`. On the real
  #   data the two are identical -- scores jump 4 -> 6, nothing lies between --
  #   so the published four candidates are unaffected. They differ only when
  #   the average is a half-integer (e.g. 5.5 passes `< 6` but fails `<= 5`),
  #   which does occur across bootstrap replicates. `<= 5` is used here.
  borda <- (rank(-vi, ties.method = "min") + rank(mdep, ties.method = "min")) / 2
  top   <- vsig[vi > 0 & borda <= 5]
  if (length(top) < 1) return(fail("no_borda_selected"))

  # 10. LASSO Cox at lambda.1se (needs >= 2 columns; else carry forward)
  if (length(top) >= 2) {
    set.seed(seed)
    cv <- tryCatch(cv.glmnet(x = as.matrix(rf_df[, top, drop = FALSE]),
                             y = so, family = "cox", alpha = 1,
                             nfolds = 10, maxit = 10000),
                   error = function(e) NULL)
    if (!is.null(cv)) {
      cf  <- as.matrix(coef(cv, s = "lambda.1se"))
      sel <- rownames(cf)[cf[, 1] != 0]
      if (length(sel) == 0) {
        cf  <- as.matrix(coef(cv, s = "lambda.min"))
        sel <- rownames(cf)[cf[, 1] != 0]
      }
      if (length(sel) > 0) top <- sel
    }
  }
  if (length(top) < 1) return(fail("lasso_empty"))

  # 11. unpenalised multivariate Cox -> final betas
  fin <- tryCatch(coxph(so ~ ., data = rf_df[, top, drop = FALSE], ties = "efron"),
                  error = function(e) NULL)
  if (is.null(fin) || any(is.na(coef(fin)))) return(fail("final_cox_failed"))

  # map the sanitised RSF names back to the original miRNA ids
  orig <- setNames(sig, make.names(sig))
  list(ok = TRUE, betas = coef(fin), vars = top,
       mirnas = unname(orig[top]),
       zmean = zm, zsd = zs, ncand = length(cand), nsig = length(sig))
}

# score original data using a model fitted on a bootstrap replicate
score_with <- function(fit, expr, meta) {
  lg <- log2(expr + 1)
  mi <- fit$mirnas
  if (!all(mi %in% rownames(lg))) return(NULL)
  z <- t((lg[mi, , drop = FALSE] - fit$zmean[mi]) / fit$zsd[mi])
  as.numeric(z %*% fit$betas)
}

# ------------------------------------------------------------- apparent
cat("\n--- Apparent model (full data) ---\n")
app <- run_pipeline(expr_all, meta_all, seed = 42)
if (!isTRUE(app$ok)) stop("Apparent model failed: ", app$why)
cat("Selected:", paste(gsub("\\.", "-", app$mirnas), collapse = ", "), "\n")
print(round(app$betas, 4))
app_score <- score_with(app, expr_all, meta_all)
C_app <- cindex(meta_all$Time.Dead, meta_all$dead.status, app_score)
cat(sprintf("Apparent C-index = %.4f   (candidates %d -> univariate-sig %d)\n",
            C_app, app$ncand, app$nsig))

# ------------------------------------------------------------ bootstrap
if (B == 0) { cat("\nB = 0: apparent model only, stopping here.\n"); quit(save = "no") }

cat(sprintf("\n--- Nested bootstrap: B = %d ---\n", B))
n   <- ncol(expr_all)
rec <- data.frame(b = integer(), C_boot = numeric(), C_orig = numeric(),
                  optimism = numeric(), k = integer(), stringsAsFactors = FALSE)
picks   <- list()
reasons <- character()

t0 <- Sys.time()
set.seed(20260812)
seeds <- sample.int(1e6, B)

for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  fit <- run_pipeline(expr_all[, idx, drop = FALSE], meta_all[idx, ], seeds[b])
  if (!isTRUE(fit$ok)) { reasons <- c(reasons, fit$why); next }

  s_boot <- score_with(fit, expr_all[, idx, drop = FALSE], meta_all[idx, ])
  s_orig <- score_with(fit, expr_all, meta_all)
  if (is.null(s_boot) || is.null(s_orig)) { reasons <- c(reasons, "score_failed"); next }

  cb_ <- cindex(meta_all$Time.Dead[idx], meta_all$dead.status[idx], s_boot)
  co_ <- cindex(meta_all$Time.Dead,      meta_all$dead.status,      s_orig)
  if (is.na(cb_) || is.na(co_)) { reasons <- c(reasons, "cindex_na"); next }

  rec <- rbind(rec, data.frame(b = b, C_boot = cb_, C_orig = co_,
                               optimism = cb_ - co_, k = length(fit$mirnas)))
  picks[[length(picks) + 1]] <- fit$mirnas

  if (b %% 50 == 0) cat(sprintf("  %d/%d  valid=%d  elapsed=%.1f min\n",
        b, B, nrow(rec), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

# --------------------------------------------------------------- results
opt   <- mean(rec$optimism)
C_cor <- C_app - opt
ci    <- quantile(C_app - rec$optimism, c(0.025, 0.975), na.rm = TRUE)

# ---- SAVE FIRST -----------------------------------------------------
# Everything below this point is formatting. Persist the results before
# printing so a formatting error can never discard a 23-minute run again.
freq <- sort(table(unlist(picks)), decreasing = TRUE) / length(picks) * 100
ft   <- data.frame(miRNA = gsub("\\.", "-", names(freq)),
                   pct   = round(as.numeric(freq), 1))

saveRDS(list(rec = rec, picks = picks, C_app = C_app, opt = opt,
             C_cor = C_cor, ci = ci, apparent = app, reasons = reasons),
        file.path(OUT, sprintf("bootstrap_raw_%s.rds", TAG)))
write.csv(rec, file.path(OUT, sprintf("bootstrap_iterations_%s.csv", TAG)), row.names = FALSE)
write.csv(ft,  file.path(OUT, sprintf("selection_frequency_%s.csv", TAG)), row.names = FALSE)
write.csv(data.frame(
  Metric = c("n", "events", "B_requested", "B_valid", "Apparent C-index",
             "Mean optimism", "Optimism-corrected C-index",
             "95% CI lower", "95% CI upper", "Runtime (min)",
             "Apparent model"),
  Value  = c(n, sum(meta_all$dead.status), B, nrow(rec), round(C_app, 4),
             round(opt, 4), round(C_cor, 4), round(ci[1], 4), round(ci[2], 4),
             round(el, 1), paste(gsub("\\.", "-", app$mirnas), collapse = " + "))),
  file.path(OUT, sprintf("bootstrap_summary_%s.csv", TAG)), row.names = FALSE)
cat("\n[saved to", OUT, "]\n")

cat(sprintf("\n================ A1 RESULTS (%s) ================\n", TAG))
cat(sprintf("Runtime                     : %.1f min\n", el))
cat(sprintf("Valid replicates            : %d / %d (%.1f%% failed)\n",
            nrow(rec), B, 100 * (B - nrow(rec)) / B))
if (length(reasons)) { cat("Failure reasons:\n"); print(table(reasons)) }
cat(sprintf("\nApparent C-index            : %.4f\n", C_app))
cat(sprintf("Mean optimism               : %.4f\n", opt))
cat(sprintf("OPTIMISM-CORRECTED C-index  : %.4f  (95%% CI %.4f - %.4f)\n",
            C_cor, ci[1], ci[2]))
cat("\nFor reference, the published values were 0.8031 / 0.0037 / 0.7994\n")
cat(sprintf("Model size across replicates: median %.1f miRNAs (range %d-%d)\n",
            median(rec$k), min(rec$k), max(rec$k)))

cat(sprintf("\n--- Selection frequency across %d valid replicates ---\n", length(picks)))
print(head(ft, 20), row.names = FALSE)
cat("\nThe three published miRNAs:\n")
print(ft[ft$miRNA %in% c("hsa-miR-509-3p", "hsa-miR-377-3p", "hsa-miR-378h"), ],
      row.names = FALSE)

cat("\nDone. Results in", OUT, "\n")
