# =====================================================================
#  Figure 4B -- time-dependent ROC, REGENERATED with logit-bounded CIs
#
#  PROBLEM (Reviewer 1, Minor comment 3): the published 6-month interval
#  was 0.892 (0.778-1.006). An AUC cannot exceed 1. The cause is at
#  1.2_discovery_pipeline.Rmd lines 2019-2022, which build a Wald interval on the
#  raw AUC scale:
#      AUC +/- 1.96 * SE
#  Nothing constrains that to [0,1].
#
#  FIX: compute the interval on the logit scale and back-transform.
#  SE is carried across by the delta method: se_logit = se / (AUC*(1-AUC)).
#  AUC point estimates are unchanged; only the intervals move.
#
#  Style matched to the published figure: colours #B22222 / #E88B00 /
#  #1874CD, legend inside at (0.62, 0.18), aspect.ratio = 1, theme_bw.
#  Survival time: ROUNDED MONTHS, as published.
# =====================================================================

suppressPackageStartupMessages({
  library(limma); library(sva); library(survival); library(timeROC); library(ggplot2); library(grid)
})

OUT <- "results/Updated_figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SIG <- c("hsa.miR.509.3p", "hsa.miR.377.3p", "hsa.miR.378h")
PUB <- c(hsa.miR.509.3p = 0.743074, hsa.miR.377.3p = 0.502724,
         hsa.miR.378h  = -0.663259)

counts <- read.csv("data/Count_Matrix_removed_unknown.csv", row.names = 1)
expr   <- t(as.matrix(counts)); colnames(expr) <- trimws(colnames(expr))
meta   <- read.csv("data/Metadata_FINAL_rounded_34pt.csv",
                   stringsAsFactors = FALSE)
meta   <- meta[match(colnames(expr), trimws(meta$SampleID)), ]

status <- as.integer(meta$dead.status); time <- as.numeric(meta$Time.Dead)
sf <- factor(status, levels = c(0, 1), labels = c("Alive", "Dead"))
lg <- log2(expr + 1)
lg <- lg[(apply(lg[, status == 1], 1, median) > 1) |
         (apply(lg[, status == 0], 1, median) > 1), ]
cb <- ComBat(lg, batch = factor(meta$Batch), mod = model.matrix(~ sf),
             par.prior = TRUE, prior.plots = FALSE)
score <- as.numeric(t(t(scale(t(cb)))[SIG, , drop = FALSE]) %*% PUB[SIG])

troc <- timeROC(T = time, delta = status, marker = score,
                cause = 1, times = c(6, 12, 18), iid = TRUE)

# ---- the fix -------------------------------------------------------
auc_raw <- as.numeric(troc$AUC)
auc_se  <- as.numeric(troc$inference$vect_sd_1)
z_logit <- 1.96 * auc_se / (auc_raw * (1 - auc_raw))

auc_df <- data.frame(
  Timepoint = c(6, 12, 18),
  AUC       = auc_raw,
  AUC_lower = plogis(qlogis(auc_raw) - z_logit),
  AUC_upper = plogis(qlogis(auc_raw) + z_logit),
  Wald_lower = auc_raw - 1.96 * auc_se,
  Wald_upper = auc_raw + 1.96 * auc_se
)
cat("--- AUC: published Wald vs corrected logit ---\n")
print(round(auc_df, 4))
cat(sprintf("\nWald upper bounds exceeding 1: %d\n", sum(auc_df$Wald_upper > 1)))

# ---- plot, matching the published style -----------------------------
roc_df <- data.frame(
  FPR  = c(troc$FP[, 1], troc$FP[, 2], troc$FP[, 3]),
  TPR  = c(troc$TP[, 1], troc$TP[, 2], troc$TP[, 3]),
  Time = factor(rep(paste0(c(6, 12, 18), " months"), each = nrow(troc$FP)),
                levels = paste0(c(6, 12, 18), " months"))
)

auc_labels <- paste0(c(6, 12, 18), " months  AUC = ", sprintf("%.3f", auc_df$AUC),
                     " (", sprintf("%.3f", auc_df$AUC_lower), "–",
                     sprintf("%.3f", auc_df$AUC_upper), ")")
names(auc_labels) <- paste0(c(6, 12, 18), " months")

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Time)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("6 months"  = "#B22222",
                                "12 months" = "#E88B00",
                                "18 months" = "#1874CD"),
                     labels = auc_labels, name = NULL) +
  scale_x_continuous(name = "1-Specificity", limits = c(0, 1), expand = c(0.01, 0)) +
  scale_y_continuous(name = "Sensitivity",   limits = c(0, 1), expand = c(0.01, 0)) +
  theme_bw(base_size = 11) +
  theme(legend.position   = c(0.60, 0.16),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text       = element_text(size = 12),
        legend.key.height = unit(1.1, "lines"),
        panel.grid.minor  = element_blank(),
        aspect.ratio      = 1)

ggsave(file.path(OUT, "Figure4B_timeROC_logitCI.pdf"), p_roc, width = 6.5, height = 6.5)
ggsave(file.path(OUT, "Figure4B_timeROC_logitCI.png"), p_roc, width = 6.5, height = 6.5, dpi = 600)
ggsave(file.path(OUT, "Figure4B_timeROC_logitCI.svg"), p_roc, width = 6.5, height = 6.5)

# version without the AUC legend, for hand-labelling
p_nolab <- p_roc + scale_color_manual(
  values = c("6 months" = "#B22222", "12 months" = "#E88B00", "18 months" = "#1874CD"),
  name = NULL)
ggsave(file.path(OUT, "Figure4B_timeROC_nolabel.pdf"), p_nolab, width = 6.5, height = 6.5)
ggsave(file.path(OUT, "Figure4B_timeROC_nolabel.png"), p_nolab, width = 6.5, height = 6.5, dpi = 600)

write.csv(auc_df, file.path(OUT, "Figure4B_AUC_values.csv"), row.names = FALSE)
cat("\nSaved to", OUT, "\n")


