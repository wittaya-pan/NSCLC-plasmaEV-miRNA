# =====================================================================
#  Figure 4E -- univariate Cox forest plot, REGENERATED
#
#  Changes vs the published figure:
#    1. ECOG PS added        (Reviewer 1, Minor comment 2)
#    2. PD-L1 TPS added      (same situation as ECOG; avoids "why one and
#                             not the other?")
#    3. Risk score row is now per 1 SD (HR 3.58). The published row showed
#       HR 2.72 = e, an artefact of regressing the score on itself.
#
#  NO HR/CI TEXT ANNOTATIONS -- author will add these in Photoshop.
#  Style matched to the original: colours #B22222 / #1874CD, diamond for
#  the risk score, circle for clinical variables, log10 x-axis, theme_bw.
#
#  Outputs PDF (vector, best for Photoshop) + PNG 600 dpi + SVG.
# =====================================================================

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

OUT <- "results/Updated_figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

uni <- read.csv("results/A2_A3_clinical/univariate_clinical.csv",
                stringsAsFactors = FALSE)

# risk score, per 1 SD -- from 07_A2_A3 section 1
score_row <- data.frame(variable = "3-miRNA risk score (per 1 SD)",
                        n = 34, events = 18, HR = 3.583, lo = 1.990,
                        hi = 6.455, p = 2.16e-05)
fd <- rbind(score_row, uni)

fd <- fd %>%
  mutate(Type = ifelse(grepl("risk score", variable),
                       "3-miRNA Risk Score", "Clinical Variable"),
         Direction = ifelse(HR > 1, "Risk (HR > 1)", "Protective (HR < 1)"),
         label = factor(variable, levels = rev(variable)))

color_vals <- c("Risk (HR > 1)" = "#B22222", "Protective (HR < 1)" = "#1874CD")
ci_cap <- quantile(fd$hi, 0.95, na.rm = TRUE) * 1.5

p <- ggplot(fd, aes(x = HR, y = label, color = Direction)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = pmax(lo, 0.01), xmax = pmin(hi, ci_cap)),
                 height = 0.3, linewidth = 0.7) +
  geom_point(aes(shape = Type), size = 3.5) +
  scale_color_manual(values = color_vals, name = NULL) +
  scale_shape_manual(values = c("3-miRNA Risk Score" = 18,
                                "Clinical Variable"  = 16), name = NULL) +
  scale_x_continuous(name = "Hazard Ratio (95% CI) - log scale",
                     trans = "log10",
                     limits = c(0.1, ci_cap * 3), expand = c(0, 0)) +
  labs(y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        axis.text.y      = element_text(size = 9))

h <- max(4, 1.8 + nrow(fd) * 0.4)
ggsave(file.path(OUT, "Figure4E_forest_nolabel.pdf"), p, width = 8, height = h)
ggsave(file.path(OUT, "Figure4E_forest_nolabel.png"), p, width = 8, height = h, dpi = 600)
ggsave(file.path(OUT, "Figure4E_forest_nolabel.svg"), p, width = 8, height = h)

# plot-only version: no axis text, no legend -- pure geometry for overlaying
p_bare <- p + theme(axis.text = element_blank(), axis.title = element_blank(),
                    legend.position = "none", axis.ticks = element_blank())
ggsave(file.path(OUT, "Figure4E_forest_bare.pdf"), p_bare, width = 8, height = h)
ggsave(file.path(OUT, "Figure4E_forest_bare.png"), p_bare, width = 8, height = h, dpi = 600)

# ---- LABELLED version, matching the original figure's annotation style ----
# original: geom_text(aes(x = pmin(hi, ci_cap) + 0.05, label = annot),
#                     hjust = 0, size = 2.6, color = "grey20")
fmt_p <- function(x) ifelse(x < 0.001, "p<0.001", sprintf("p=%.3f", x))
fd$annot <- paste0("HR=", sprintf("%.2f", fd$HR),
                   " (", sprintf("%.2f", fd$lo), "–", sprintf("%.2f", fd$hi), ")",
                   "  ", fmt_p(fd$p),
                   ifelse(fd$n != 34, paste0("  (n=", fd$n, ")"), ""))

# pass data = fd explicitly: `p` was built before the annot column existed
p_lab <- p +
  geom_text(data = fd, aes(x = pmin(hi, ci_cap) + 0.05, label = annot),
            hjust = 0, size = 2.6, color = "grey20", show.legend = FALSE)

ggsave(file.path(OUT, "Figure4E_forest_labelled.pdf"), p_lab, width = 10, height = h)
ggsave(file.path(OUT, "Figure4E_forest_labelled.png"), p_lab, width = 10, height = h, dpi = 600)
ggsave(file.path(OUT, "Figure4E_forest_labelled.svg"), p_lab, width = 10, height = h)

# the values, in plot order, for hand-labelling
lab <- fd %>%
  transmute(row_top_to_bottom = row_number(), variable, n, events,
            HR = round(HR, 2), CI = sprintf("%.2f-%.2f", lo, hi),
            p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
write.csv(lab, file.path(OUT, "Figure4E_values_for_labelling.csv"), row.names = FALSE)
print(lab, row.names = FALSE)
cat("\nSaved to", OUT, "\n")
