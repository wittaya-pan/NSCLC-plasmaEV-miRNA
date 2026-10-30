# =====================================================================
#  Build the FINAL analysis metadata — AUTHOR DECISION: rounded months
#
#  Time.Dead is kept EXACTLY as published (whole months). This keeps the
#  whole pipeline internally consistent and reproducible from Methods, and
#  leaves Tables S2, S3 and S4 at their published values -- no edits.
#
#  Only the genuine errors are corrected:
#     1. one patient aged just under 60 was coded ">=60"  ->  "<60"
#     2. ECOG PS and PD-L1 TPS added (present in Table 1, absent from every
#        analysis file, which is why Reviewer 1 asked why ECOG was omitted)
#
#  Exact times are carried alongside as reference columns only
#  (Time.Dead.exact, Time.Dead.days). NOTHING downstream uses them.
# =====================================================================

suppressPackageStartupMessages({library(readxl); library(dplyr); library(survival)})

OUT <- "data"

old <- read.csv("data/Metadata_OLD_floored.csv", stringsAsFactors = FALSE)
xl  <- as.data.frame(read_excel("data/Metadata34.xlsx", col_types = "text"))
norm <- function(x) sub("^pt0*", "pt", trimws(tolower(x)))
old$.k <- norm(old$SampleID); xl$.k <- norm(xl$SampleID)
stopifnot(setequal(old$.k, xl$.k), nrow(old) == 34)

as_date <- function(v) {
  v <- trimws(as.character(v)); o <- as.Date(rep(NA_real_, length(v)), origin = "1970-01-01")
  n <- suppressWarnings(as.numeric(v)); s <- !is.na(n) & n > 20000 & n < 60000
  o[s] <- as.Date(n[s], origin = "1899-12-30")
  i <- is.na(o) & grepl("^\\d{4}-", v); o[i] <- as.Date(substr(v[i], 1, 10)); o
}
st <- setNames(old$dead.status, old$.k)
xl$start_d <- as_date(xl$start1)
xl$end_d   <- as.Date(ifelse(st[xl$.k] == 1, as_date(xl$deaddate), as_date(xl$lastFUdate)),
                      origin = "1970-01-01")
stopifnot(!any(is.na(xl$start_d)), !any(is.na(xl$end_d)))

meta <- old %>%
  select(-Age) %>%
  left_join(xl %>% transmute(.k,
              Time.Dead.days  = as.numeric(end_d - start_d),
              Time.Dead.exact = as.numeric(end_d - start_d) / 30.4375,
              Age.years = as.numeric(age),
              ECOG      = trimws(as.character(ECOG)),
              PDL1_TPS  = suppressWarnings(as.numeric(`PDL1_22c3 TPS`)),
              PDL1_gr   = trimws(as.character(PDL1_gr))), by = ".k") %>%
  mutate(Age     = ifelse(Age.years >= 60, ">=60", "<60"),        # age group recomputed from exact age
         ECOG_gr = case_when(ECOG %in% c("0","1") ~ "0-1",
                             ECOG == "2"          ~ "2", TRUE ~ "unknown"),
         PDL1_gr = ifelse(PDL1_gr == "0", "<1", PDL1_gr)) %>%
  select(SampleID, sex, smoking, histo_detail, stage, ajcc8_stage, T, N, M,
         firstline_detail, dead.status, Time.Dead, Batch, Age, Age.years,
         ECOG, ECOG_gr, PDL1_TPS, PDL1_gr, Time.Dead.exact, Time.Dead.days)

# Time.Dead must be byte-identical to the published file
stopifnot(identical(meta$Time.Dead[match(old$SampleID, meta$SampleID)], old$Time.Dead))
cat("Time.Dead identical to published file: TRUE\n")

write.csv(meta, file.path(OUT, "Metadata_FINAL_rounded_34pt.csv"), row.names = FALSE)

cat(sprintf("\nn = %d, events = %d\n", nrow(meta), sum(meta$dead.status)))
f  <- survfit(Surv(Time.Dead, dead.status) ~ 1, data = meta, conf.type = "log-log")
tb <- summary(f)$table
cat(sprintf("Median OS        = %.1f months (95%% CI %.1f - %s)   [MS says 14.3]\n",
            tb["median"], tb["0.95LCL"],
            ifelse(is.na(tb["0.95UCL"]), "NR", sprintf("%.1f", tb["0.95UCL"]))))
rev <- survfit(Surv(Time.Dead, 1 - dead.status) ~ 1, data = meta)
cat(sprintf("Median follow-up = %.1f months (reverse KM)\n", summary(rev)$table["median"]))
cat(sprintf("Follow-up range  = %g - %g months\n", min(meta$Time.Dead), max(meta$Time.Dead)))
cat("\nAge group:\n"); print(table(meta$Age))
cat("\nECOG:\n");      print(table(meta$ECOG_gr))
cat("\nPD-L1:\n");     print(table(meta$PDL1_gr))
cat("\nSaved: Metadata_FINAL_rounded_34pt.csv\n")
