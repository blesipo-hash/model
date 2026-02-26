# Purpose:
# Build Beta(alpha,beta) priors for drug efficacy expressed as ERR (egg reduction rate),
# using mean + 95% CI from literature.


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
})

dir.create("data/parameters", recursive = TRUE, showWarnings = FALSE)

# This script intentionally has no WHO raw-file parsing.
# WHO PC inputs are produced in 01b_build_who_pc_inputs.R and consumed downstream.

clip01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

fit_beta_from_mean_interval <- function(m, lo, hi, p_lo = 0.025, p_hi = 0.975, eps = 1e-6) {
  m <- clip01(as.numeric(m), eps)
  lo <- clip01(as.numeric(min(lo, hi)), eps)
  hi <- clip01(as.numeric(max(lo, hi)), eps)
  
  if (!(lo < m && m < hi)) {
    stop(
      sprintf("ERR interval must straddle mean: lo=%.6f mean=%.6f hi=%.6f", lo, m, hi),
      call. = FALSE
    )
  }
  
  obj <- function(log_kappa) {
    kappa <- exp(log_kappa)
    a <- m * kappa
    b <- (1 - m) * kappa
    qlo <- qbeta(p_lo, a, b)
    qhi <- qbeta(p_hi, a, b)
    (qlo - lo)^2 + (qhi - hi)^2
  }
  
  opt <- optim(
    par = log(20),
    fn = obj,
    method = "Brent",
    lower = log(1e-2),
    upper = log(1e6)
  )
  
  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning("optim() did not converge for Beta fit; results may be unreliable.", call. = FALSE)
  }
  
  kappa <- exp(opt$par)
  alpha <- max(m * kappa, 1e-6)
  beta <- max((1 - m) * kappa, 1e-6)
  
  list(
    alpha = alpha,
    beta = beta,
    kappa = kappa,
    qlo_hat = qbeta(p_lo, alpha, beta),
    qhi_hat = qbeta(p_hi, alpha, beta),
    opt_value = opt$value,
    opt_convergence = opt$convergence
  )
}

# Committed regimen set: exactly one row per parasite label used downstream.
drug_err <- tribble(
  ~parasite,            ~drug,          ~dose,      ~err_mean, ~err_min, ~err_max, ~source_note,
  "Schistosomiasis",    "Praziquantel", "40 mg/kg",  0.95,      0.84,     0.97,     "Fukushige 2021 PLoS NTD",
  "Ascariasis",         "Albendazole",  "400 mg",    0.985,     0.949,    1.000,    "Moser et al. 2017 BMJ, Table 1 (ERR, 95% CI)",
  "Trichuriasis",       "Albendazole",  "400 mg",    0.499,     0.390,    0.606,    "Moser et al. 2017 BMJ, Table 1 (ERR, 95% CI)",
  "Hookworm disease",   "Albendazole",  "400 mg",    0.896,     0.819,    0.973,    "Moser et al. 2017 BMJ, Table 1 (ERR, 95% CI)"
) %>%
  mutate(
    err_mean_raw = as.numeric(err_mean),
    err_min_raw = as.numeric(err_min),
    err_max_raw = as.numeric(err_max),
    err_mean = clip01(err_mean_raw),
    err_min = clip01(err_min_raw),
    err_max = clip01(err_max_raw)
  )

expected <- c("Schistosomiasis", "Ascariasis", "Trichuriasis", "Hookworm disease")
bad_labels <- setdiff(unique(drug_err$parasite), expected)
missing_labels <- setdiff(expected, unique(drug_err$parasite))
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

dup_any <- drug_err %>% count(parasite) %>% filter(n != 1)
if (nrow(dup_any) > 0) {
  write_csv(dup_any, "data/parameters/qc_err_committed_dup_or_missing.csv")
  stop(
    "Committed regimen table must contain exactly 1 row per parasite. QC written to qc_err_committed_dup_or_missing.csv",
    call. = FALSE
  )
}

input_qc <- drug_err %>%
  transmute(
    parasite, drug, dose,
    err_mean_raw, err_min_raw, err_max_raw,
    err_mean, err_lo = pmin(err_min, err_max), err_hi = pmax(err_min, err_max),
    flag_mean_not_in_ci = !(err_lo < err_mean & err_mean < err_hi),
    flag_raw_outside_0_1 = err_mean_raw < 0 | err_mean_raw > 1 |
      err_min_raw < 0 | err_min_raw > 1 |
      err_max_raw < 0 | err_max_raw > 1
  )
write_csv(input_qc, "data/parameters/qc_drug_err_input_validation.csv")

if (any(input_qc$flag_mean_not_in_ci, na.rm = TRUE)) {
  stop("ERR input QC failed: some means are not strictly inside their CI bounds. See qc_drug_err_input_validation.csv", call. = FALSE)
}

write_csv(drug_err, "data/parameters/drug_err_full.csv")

priors <- drug_err %>%
  mutate(
    fit = pmap(list(err_mean, err_min, err_max), ~ fit_beta_from_mean_interval(..1, ..2, ..3)),
    alpha = map_dbl(fit, "alpha"),
    beta = map_dbl(fit, "beta"),
    kappa = map_dbl(fit, "kappa"),
    qlo_hat = map_dbl(fit, "qlo_hat"),
    qhi_hat = map_dbl(fit, "qhi_hat"),
    opt_value = map_dbl(fit, "opt_value"),
    opt_convergence = map_int(fit, "opt_convergence")
  ) %>%
  select(
    parasite, drug, dose,
    err_mean_raw, err_min_raw, err_max_raw,
    err_mean, err_lo = err_min, err_hi = err_max,
    alpha, beta, kappa,
    qlo_hat, qhi_hat,
    opt_value, opt_convergence,
    source_note
  )

qc <- priors %>%
  mutate(
    abs_lo_err = abs(qlo_hat - err_lo),
    abs_hi_err = abs(qhi_hat - err_hi),
    flag_nonfinite_alpha_beta = !is.finite(alpha) | !is.finite(beta),
    flag_alpha_beta_nonpositive = alpha <= 0 | beta <= 0
  )

write_csv(qc, "data/parameters/qc_drug_err_committed_fit.csv")

if (any(qc$flag_nonfinite_alpha_beta | qc$flag_alpha_beta_nonpositive, na.rm = TRUE)) {
  stop("ERR prior fitting produced invalid alpha/beta. See qc_drug_err_committed_fit.csv", call. = FALSE)
}

write_csv(priors, "data/parameters/drug_err_default.csv")
write_csv(priors %>% select(parasite, alpha, beta), "data/parameters/drug_err_default_minimal.csv")