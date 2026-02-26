# =============================================================================
# R/00_make_k_priors.R
# Purpose:
#   Build parasite-specific priors for the aggregation parameter k.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

dir.create("data/parameters", recursive = TRUE, showWarnings = FALSE)

# Helpers
is_bad_num <- function(x) is.na(x) | !is.finite(x)

# Committed parasite label set (MUST match downstream exactly)
expected <- c("Schistosomiasis", "Ascariasis", "Trichuriasis", "Hookworm disease")

# =============================================================================
# 1) SUPPLIED PRIORS (FILL THESE)
# =============================================================================
# You MUST replace NA with literature-derived values.
# Guidance:
#   - k0 must be > 0
#   - sigma_k must be > 0
#   - source_note must cite the paper/table/estimate you used
k_priors <- tribble(
  ~parasite,            ~k0,      ~sigma_k,      ~source_note,
  "Schistosomiasis",     0.2385,   0.9109642,    "Chan et al. 1995 Table 1 constant k=0.2385; Kura et al. 2019 states k=0.04 (low) and k=0.24 (moderate/high)",
  "Ascariasis",          0.59,     0.2520503,    "Truscott et al. 2015 Table 1 k=0.90; Anderson et al. 2013 Table 1 Ascaris k values/ranges incl 0.36–0.54, 0.80, 0.81, 0.59; Martin et al. 1983 k=0.44",
  "Trichuriasis",        0.335,    0.5681889,    "Truscott et al. 2015 Table 1 k=0.38; Anderson et al. 2013 Table 1 Trichuris k range 0.11–0.65; Bundy et al. 1985 k=0.29",
  "Hookworm disease",    0.345,    0.2907727,    "Truscott et al. 2015 Table 1 k=0.35; Anderson et al. 2013 Table 1 Necator k=0.34 and range 0.33–0.61; Lwambo et al. 1992 k=0.34"
) %>%
  mutate(
    parasite = as.character(parasite),
    k0 = as.numeric(k0),
    sigma_k = as.numeric(sigma_k),
    source_note = as.character(source_note)
  )

# =============================================================================
# 2) QC: label contract + 1 row per parasite
# =============================================================================
bad_labels <- setdiff(unique(k_priors$parasite), expected)
missing_labels <- setdiff(expected, unique(k_priors$parasite))
if (length(bad_labels) > 0 || length(missing_labels) > 0) {
  stop(
    paste0(
      "Parasite labels must match downstream exactly.\n",
      "Unexpected: ", paste(bad_labels, collapse = ", "), "\n",
      "Missing: ", paste(missing_labels, collapse = ", ")
    ),
    call. = FALSE
  )
}

dup_any <- k_priors %>% count(parasite, name = "n") %>% filter(n != 1)
if (nrow(dup_any) > 0) {
  write_csv(dup_any, "data/parameters/qc_k_priors_dup_or_missing.csv")
  stop(
    "Committed k prior table must contain exactly 1 row per parasite. ",
    "QC written: data/parameters/qc_k_priors_dup_or_missing.csv",
    call. = FALSE
  )
}

# =============================================================================
# 3) QC: values must be present and valid
# =============================================================================
qc <- k_priors %>%
  mutate(
    flag_k0_missing = is_bad_num(k0),
    flag_sigma_missing = is_bad_num(sigma_k),
    flag_k0_nonpositive = !flag_k0_missing & k0 <= 0,
    flag_sigma_nonpositive = !flag_sigma_missing & sigma_k <= 0,
    flag_source_missing = is.na(source_note) | trimws(source_note) == ""
  )

write_csv(qc, "data/parameters/qc_k_priors_validation.csv")

if (any(qc$flag_k0_missing | qc$flag_sigma_missing |
        qc$flag_k0_nonpositive | qc$flag_sigma_nonpositive |
        qc$flag_source_missing, na.rm = TRUE)) {
  stop(
    paste0(
      "k prior QC failed: you must fill k0 (>0), sigma_k (>0), and source_note (non-empty).\n",
      "See: data/parameters/qc_k_priors_validation.csv"
    ),
    call. = FALSE
  )
}

# =============================================================================
# 4) Write outputs
# =============================================================================
write_csv(k_priors, "data/parameters/k_priors_default.csv")
write_csv(k_priors %>% select(parasite, k0, sigma_k),
          "data/parameters/k_priors_default_minimal.csv")

message("Wrote: data/parameters/k_priors_default.csv")
message("Wrote: data/parameters/k_priors_default_minimal.csv")
message("DONE.")