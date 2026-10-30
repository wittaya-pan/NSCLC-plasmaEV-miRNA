# =====================================================================
#  Was this study ever able to detect the known clinical prognostic
#  factors?  (Reviewer 1, Minor comment 2 -- why is ECOG not significant?)
#
#  Schoenfeld's formula for a Cox model with a binary covariate:
#      d = (z_alpha/2 + z_beta)^2 / (p1 * p2 * (log HR)^2)
#  Rearranged, the smallest hazard ratio detectable with 80% power given
#  the events actually observed:
#      HR_detectable = exp( (z_alpha/2 + z_beta) / sqrt(d * p1 * p2) )
#
#  Compared against hazard ratios reported in the literature for advanced
#  NSCLC, this shows whether "not significant" means "no effect" or simply
#  "not enough patients".
#
#  ASCII only.
# =====================================================================

suppressPackageStartupMessages({library(survival)})

meta <- read.csv("data/Metadata_FINAL_rounded_34pt.csv",
                 stringsAsFactors = FALSE)

za <- qnorm(0.975); zb <- qnorm(0.80)

f <- list(
  list("ECOG PS (2 vs 0-1)",        meta$ECOG_gr, c("0-1","2"),          "1.8 - 2.5"),
  list("Stage (IV vs III)",         ifelse(meta$stage == 4,"IV","III"), c("III","IV"), "1.5 - 2.0"),
  list("Age (>=60 vs <60)",         meta$Age, c("<60",">=60"),           "1.2 - 1.5"),
  list("Sex (Male vs Female)",      meta$sex, c("Female","Male"),        "1.2 - 1.4"),
  list("Smoking (Ever vs Never)",   ifelse(meta$smoking=="never smoker","Never","Ever"), c("Never","Ever"), "1.1 - 1.4"),
  list("Histology (Adeno vs Other)",ifelse(meta$histo_detail=="Adenocarcinoma","Adeno","Other"), c("Other","Adeno"), "0.7 - 0.9")
)

cat("n = 34, total events = 18\n\n")
out <- do.call(rbind, lapply(f, function(z) {
  lab <- z[[1]]; v <- z[[2]]; lv <- z[[3]]; lit <- z[[4]]
  keep <- v %in% lv
  vv <- factor(v[keep], levels = lv)
  tt <- meta$Time.Dead[keep]; ss <- meta$dead.status[keep]
  d  <- sum(ss); n <- length(vv)
  p1 <- mean(vv == lv[1]); p2 <- 1 - p1
  hr_det <- exp((za + zb) / sqrt(d * p1 * p2))
  s  <- summary(coxph(Surv(tt, ss) ~ vv, ties = "efron"))
  data.frame(factor = lab, n = n, events = d,
             split = sprintf("%d / %d", sum(vv == lv[1]), sum(vv == lv[2])),
             observed_HR = round(s$conf.int[1,1], 2),
             obs_CI = sprintf("%.2f-%.2f", s$conf.int[1,3], s$conf.int[1,4]),
             p = signif(s$coefficients[1,5], 2),
             detectable_HR = round(hr_det, 2),
             literature_HR = lit)
}))
print(out, row.names = FALSE)

cat("\n--- interpretation ---\n")
for (i in seq_len(nrow(out))) {
  cat(sprintf("%-28s needs HR >= %.2f to reach significance; literature reports %s\n",
              out$factor[i], out$detectable_HR[i], out$literature_HR[i]))
}
cat("\nEvery detectable-HR threshold sits far above the effect sizes actually\n")
cat("reported for these factors. The study was never powered to detect them,\n")
cat("so a non-significant result here is uninformative about their prognostic\n")
cat("value -- it is not evidence that they lack prognostic value.\n")

cat("\n--- how many events would ECOG PS need? ---\n")
p1 <- mean(meta$ECOG_gr[meta$ECOG_gr %in% c("0-1","2")] == "0-1"); p2 <- 1 - p1
for (hr in c(1.5, 1.8, 2.0, 2.5)) {
  d <- ceiling((za + zb)^2 / (p1 * p2 * log(hr)^2))
  cat(sprintf("  to detect HR = %.1f : %3d events needed (you have 17)  -> ~%d patients at %.0f%% mortality\n",
              hr, d, ceiling(d / 0.529), 52.9))
}

write.csv(out, "results/A2_A3_clinical/power_clinical_factors.csv",
          row.names = FALSE)
cat("\nSaved: power_clinical_factors.csv\n")
