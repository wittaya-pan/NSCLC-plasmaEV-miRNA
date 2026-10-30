# =====================================================================
#  18 - Reviewer 1: combined sensitivity table and forest plot
#
#  Collects the published analysis, the 12- and 18-month landmark analyses
#  (script 16) and the Stage IV-only analysis (script 17) into:
#    Panel A - the three candidate miRNAs: DE log2FC / P and univariate Cox
#              HR per SD (95% CI) / P under each grouping
#    Panel B - the 3-miRNA risk score (fixed coefficients): apparent C-index,
#              HR per SD, high vs low risk, and the stage-stratified HR
#  plus a forest plot of the per-miRNA hazard ratios.
#
#  Requires outputs of scripts 16 (run with 12 and 18) and 17.
#  The "Published" per-miRNA Cox row is recomputed here and checked against
#  published Table S2 before anything is written.
#
#  Run from project root:  Rscript 4_sensitivity/4.3_TableS7_FigureS1.R
#  ASCII only.
# =====================================================================
suppressPackageStartupMessages({
  library(limma); library(sva); library(survival); library(ggplot2)
})
RES  <- "results"
OUT  <- file.path(RES, "R1_sensitivity_combined")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
META <- "data/Metadata_FINAL_rounded_34pt.csv"
SIG  <- c("hsa.miR.509.3p", "hsa.miR.377.3p", "hsa.miR.378h")
LAB  <- c(hsa.miR.509.3p = "miR-509-3p", hsa.miR.377.3p = "miR-377-3p",
          hsa.miR.378h = "miR-378h")
# published Table S2 (univariate Cox, HR per SD) - used as a check only
S2   <- c(hsa.miR.509.3p = 2.109, hsa.miR.377.3p = 2.005, hsa.miR.378h = 0.510)

# ---- published per-miRNA Cox, full cohort (as Table S2) -------------------
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
pub_cox <- do.call(rbind, lapply(SIG, function(m) {
  s <- summary(coxph(Surv(time, status) ~ Z[m, ], ties = "efron"))
  data.frame(id = m, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3],
             hi = s$conf.int[1, 4], p = s$coefficients[1, 5])
}))
chk <- abs(round(pub_cox$HR, 3) - S2[pub_cox$id])
if (any(chk > 0.001)) stop("Published per-miRNA Cox does not reproduce Table S2: ",
                           paste(round(pub_cox$HR, 3), collapse = ", "))

# ---- landmark results (script 16) ------------------------------------------
lmk <- function(L) {
  d <- read.csv(file.path(RES, sprintf("R1_landmark%s/landmark%s_published_miRNA_comparison.csv", L, L)))
  g <- read.csv(file.path(RES, sprintf("R1_landmark%s/landmark%s_group_assignment.csv", L, L)))
  d$id <- gsub("-", ".", d$miRNA)
  list(d = d, n = sum(!is.na(g$landmark)), died = sum(g$landmark == 1, na.rm = TRUE),
       alive = sum(g$landmark == 0, na.rm = TRUE),
       events = sum(g$dead.status[!is.na(g$landmark)]))
}
l12 <- lmk(12); l18 <- lmk(18)

# ---- Stage IV results (script 17) ------------------------------------------
s4m  <- read.csv(file.path(RES, "R1_stageIV/stageIV_per_miRNA_cox.csv"))
s4m$id <- gsub("-", ".", sub("^hsa-", "hsa.", s4m$miRNA)); s4m$id <- gsub("-", ".", s4m$id)
perf <- read.csv(file.path(RES, "R1_stageIV/stageIV_performance.csv"))
strt <- read.csv(file.path(RES, "R1_stageIV/stage_stratified_cox.csv"))

# ---- Panel A: long table ----------------------------------------------------
rowA <- function(id, analysis, ngrp, nev, lfc, dp, hr, lo, hi, p)
  data.frame(miRNA = LAB[[id]], analysis = analysis, n_groups = ngrp, n_events = nev,
             DE_log2FC = lfc, DE_P = dp, HR_perSD = hr, HR_lo = lo, HR_hi = hi, Cox_P = p)
A <- do.call(rbind, lapply(SIG, function(id) {
  pc <- pub_cox[pub_cox$id == id, ]
  a  <- l12$d[l12$d$id == id, ]; b <- l18$d[l18$d$id == id, ]; s <- s4m[s4m$id == id, ]
  rbind(
    rowA(id, "Primary: status at last follow-up", "34 (18 / 16)", sum(status),
         a$pub_logFC, a$pub_P, pc$HR, pc$lo, pc$hi, pc$p),
    rowA(id, "12-month landmark", sprintf("%d (%d / %d)", l12$n, l12$died, l12$alive), l12$events,
         a$lm_logFC, a$lm_P, a$lm_coxHR_perSD, a$lm_cox_lo, a$lm_cox_hi, a$lm_cox_p),
    rowA(id, "18-month landmark", sprintf("%d (%d / %d)", l18$n, l18$died, l18$alive), l18$events,
         b$lm_logFC, b$lm_P, b$lm_coxHR_perSD, b$lm_cox_lo, b$lm_cox_hi, b$lm_cox_p),
    rowA(id, "Stage IV only", "30", perf$events[perf$cohort == "Stage IV only"],
         NA, NA, s$HR_perSD, s$lo, s$hi, s$p))
}))
rownames(A) <- NULL

# ---- Panel B: locked risk score ---------------------------------------------
pr <- function(cohort, label) { r <- perf[perf$cohort == cohort, ]
  data.frame(analysis = label, n = r$n, events = r$events, C_index = r$C_index_apparent,
             HR_perSD = r$HR_perSD, HR_lo = r$HR_perSD_lo, HR_hi = r$HR_perSD_hi, HR_P = r$HR_perSD_p,
             HR_high_low = r$HR_high_vs_low, HL_lo = r$HR_hl_lo, HL_hi = r$HR_hl_hi,
             logrank_P = r$logrank_p) }
sr <- function(model, label) { r <- strt[strt$model == model, ]
  data.frame(analysis = label, n = 34, events = sum(status), C_index = NA,
             HR_perSD = r$HR, HR_lo = r$lo, HR_hi = r$hi, HR_P = r$p,
             HR_high_low = NA, HL_lo = NA, HL_hi = NA, logrank_P = NA) }
B <- rbind(pr("All patients (published)", "All patients"),
           pr("Stage IV only", "Stage IV only"),
           sr("Stratified by stage (III / IV)", "All patients, stratified by stage"))

# ---- save tables BEFORE plotting --------------------------------------------
write.csv(A, file.path(OUT, "R1_sensitivity_panelA_miRNAs.csv"), row.names = FALSE)
write.csv(B, file.path(OUT, "R1_sensitivity_panelB_riskscore.csv"), row.names = FALSE)

# ---- forest plot ----------------------------------------------------------------
ord <- c("Primary: status at last follow-up", "12-month landmark",
         "18-month landmark", "Stage IV only")
fp <- A
fp$analysis <- factor(fp$analysis, levels = rev(ord))
fp$miRNA    <- factor(fp$miRNA, levels = c("miR-377-3p", "miR-378h", "miR-509-3p"))
cols <- c("Primary: status at last follow-up" = "#000000",   # Okabe-Ito, colour-blind safe
          "12-month landmark" = "#0072B2", "18-month landmark" = "#56B4E9",
          "Stage IV only" = "#E69F00")
fp$lab <- sprintf("%.2f (%.2f-%.2f)", fp$HR_perSD, fp$HR_lo, fp$HR_hi)
g <- ggplot(fp, aes(x = HR_perSD, y = analysis, colour = analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = HR_lo, xmax = HR_hi), orientation = "y", width = 0.25,
                linewidth = 0.6) +
  geom_point(size = 2.4, shape = 15) +
  geom_text(aes(x = 5.6, label = lab), hjust = 0, size = 2.8, colour = "black") +
  geom_text(data = data.frame(miRNA = factor("miR-377-3p", levels = levels(fp$miRNA)),
                              x = 5.6, y = 4.85, label = "HR per SD (95% CI)"),
            aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0,
            size = 2.8, fontface = "bold") +
  facet_grid(miRNA ~ ., switch = "y") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
  coord_cartesian(xlim = c(0.2, 4.6), clip = "off") +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "Hazard ratio per SD of expression (log scale)", y = NULL) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
        panel.spacing = unit(4, "mm"), axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        plot.margin = margin(9, 32, 4, 4, "mm"))
ggsave(file.path(OUT, "FigureS_R1_forest.pdf"), g, width = 178, height = 110, units = "mm",
       device = cairo_pdf)
ggsave(file.path(OUT, "FigureS_R1_forest.png"), g, width = 178, height = 110, units = "mm",
       dpi = 600)

options(width = 160)
cat("Published per-miRNA Cox reproduces Table S2 (max |dHR| = ",
    format(max(chk), digits = 2), ")\n\n", sep = "")
print(A, digits = 3, row.names = FALSE); cat("\n"); print(B, digits = 3, row.names = FALSE)
cat("\nSaved to", OUT, "\n")
