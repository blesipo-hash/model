# =============================================================================
# R/02_prepare_stan_inputs.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tibble)
})

# ---------------------------
# 0) Settings (edit here)
# ---------------------------
SEED <- 1L
set.seed(SEED)

EPS <- 1e-6

PANEL_IN      <- "data/processed/model_input/sac_panel_merged.csv"
DALY_IN       <- "data/processed/gbd/daly_sac_by_iso3_year_parasite.csv"
DRUG_ERR_IN   <- "data/parameters/drug_err_default.csv"
K_PRIORS_IN   <- "data/parameters/k_priors_default_minimal.csv"
ADHERENCE_IN  <- "data/parameters/adherence_scenarios.csv"
HYPERPRIOR_IN <- "data/parameters/hyperprior_scenarios.csv"
ADJ_EDGES     <- "data/input/adjacency_edges.csv"        # optional: iso3_from, iso3_to

# --- Scenario switches (PRIMARY vs SENSITIVITY) ---
PC_MISSING_POLICY_SCENARIO <- "locf"  # "locf" | "zero"
OBS_POST_PC_SCENARIO <- 0L            # 0L | 1L

ADHERENCE_SCENARIO  <- "base"         # base | avg_gap | moderate | stress
HYPERPRIOR_SCENARIO <- "base"         # base | tight | loose

PC_MISSING_POLICY <- PC_MISSING_POLICY_SCENARIO
OBS_POST_PC <- as.integer(OBS_POST_PC_SCENARIO)

if (!PC_MISSING_POLICY %in% c("locf", "zero")) stop("PC_MISSING_POLICY must be 'locf' or 'zero'.", call. = FALSE)
if (!OBS_POST_PC %in% c(0L, 1L)) stop("OBS_POST_PC must be 0L or 1L.", call. = FALSE)

# Rounds logic (truth-first)
ROUNDS_CONSTANT <- 1L
if (!is.integer(ROUNDS_CONSTANT) || ROUNDS_CONSTANT < 1) stop("ROUNDS_CONSTANT must be integer >= 1.", call. = FALSE)

# ERR prior default (only used if DRUG_ERR_IN missing or incomplete—should not happen)
ERR_A_DEFAULT <- 8
ERR_B_DEFAULT <- 2

# Hyperpriors (WILL be overwritten by hyperprior scenario loader below)
BETA0_SD        <- NA_real_
SIGMA_PROC_RATE <- NA_real_
SIGMA_V_RATE    <- NA_real_
SIGMA_U_RATE    <- NA_real_
SIGMA_INIT_RATE <- NA_real_
SIGMA_OBS_RATE  <- NA_real_

# Output paths (scenario-stamped so runs don't overwrite)
tag <- sprintf(
  "adh_%s_hyp_%s_obs%d_pc_%s",
  ADHERENCE_SCENARIO,
  HYPERPRIOR_SCENARIO,
  OBS_POST_PC,
  PC_MISSING_POLICY
)

OUT_PANEL  <- sprintf("data/processed/model_input/sac_panel_model_%s.csv", tag)
OUT_STAN   <- sprintf("data/processed/stan_data_sac_%s.rds", tag)
OUT_LEVEL  <- sprintf("data/processed/levels_sac_%s.rds", tag)
OUT_QC_DIR <- sprintf("data/processed/model_input/qc_02_%s", tag)

# Backward-compatible aliases expected by some downstream scripts/runners
OUT_PANEL_ALIAS <- "data/processed/model_input/sac_panel_model.csv"
OUT_STAN_ALIAS <- "data/processed/stan_data_sac.rds"

dir.create(dirname(OUT_PANEL), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_STAN),  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# 1) Helpers (defensive)
# ---------------------------
clip01 <- function(p, eps = EPS) pmin(pmax(p, eps), 1 - eps)
logit  <- function(p) log(p / (1 - p))

need_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(name, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

assert_no_dups <- function(df, keys, name = "data") {
  dups <- df %>%
    count(across(all_of(keys)), name = "n") %>%
    filter(n > 1)
  if (nrow(dups) > 0) {
    write_csv(dups, file.path(OUT_QC_DIR, paste0("dup_keys_", name, ".csv")))
    stop(name, " has duplicate keys for: ", paste(keys, collapse = ", "),
         ". Wrote duplicates QC to: ", file.path(OUT_QC_DIR, paste0("dup_keys_", name, ".csv")),
         call. = FALSE)
  }
}

log_join_result <- function(before_df, after_df, join_name) {
  n_before <- nrow(before_df)
  n_after <- nrow(after_df)
  dropped <- n_before - n_after
  pct <- if (n_before > 0) round(100 * dropped / n_before, 2) else 0
  message(sprintf("[%s] rows: before=%s after=%s dropped=%s (%.2f%%)",
                  join_name, n_before, n_after, dropped, pct))
}

stop_if_outside_01 <- function(x, nm, allow_na = TRUE) {
  bad <- which(!is.na(x) & (x < 0 | x > 1))
  if (length(bad) > 0) {
    stop(nm, " has values outside [0,1]. First few indices: ",
         paste(head(bad, 10), collapse = ", "),
         call. = FALSE)
  }
  if (!allow_na && any(is.na(x))) stop(nm, " contains NA but allow_na=FALSE.", call. = FALSE)
}

fill_locf_vec <- function(x) {
  for (i in seq_along(x)) if (is.na(x[i]) && i > 1) x[i] <- x[i - 1]
  x
}

beta_from_mean_ci <- function(m, lo, hi, pl = 0.025, ph = 0.975, eps = 1e-6) {
  m  <- pmin(pmax(m,  eps), 1 - eps)
  lo <- pmin(pmax(min(lo, hi), eps), 1 - eps)
  hi <- pmin(pmax(max(lo, hi), eps), 1 - eps)
  
  z <- 1.96
  sd0 <- max((hi - lo) / (2 * z), 1e-3)
  v0 <- sd0^2
  v0 <- min(v0, m * (1 - m) - 1e-8)
  kappa0 <- max(m * (1 - m) / v0 - 1, 1e-2)
  a0 <- max(m * kappa0, 1e-2)
  b0 <- max((1 - m) * kappa0, 1e-2)
  
  obj <- function(par) {
    a <- exp(par[1]); b <- exp(par[2])
    mu <- a / (a + b)
    ql <- qbeta(pl, a, b)
    qh <- qbeta(ph, a, b)
    (mu - m)^2 + 5 * (ql - lo)^2 + 5 * (qh - hi)^2
  }
  
  fit <- tryCatch(
    optim(par = log(c(a0, b0)), fn = obj, method = "L-BFGS-B",
          lower = log(c(1e-6, 1e-6)), upper = log(c(1e6, 1e6))),
    error = function(e) NULL
  )
  
  if (is.null(fit) || !is.finite(fit$value)) return(c(alpha = a0, beta = b0))
  a <- exp(fit$par[1]); b <- exp(fit$par[2])
  c(alpha = a, beta = b)
}

# ---------------------------
# 1a) Write scenario audit (always)
# ---------------------------
write_csv(
  tibble(
    tag = tag,
    adherence_scenario = ADHERENCE_SCENARIO,
    hyperprior_scenario = HYPERPRIOR_SCENARIO,
    obs_post_pc = OBS_POST_PC,
    pc_missing_policy = PC_MISSING_POLICY,
    rounds_constant = ROUNDS_CONSTANT
  ),
  file.path(OUT_QC_DIR, "scenario_audit.csv")
)

# ---------------------------
# 1b) Load adherence scenario
# ---------------------------
if (!file.exists(ADHERENCE_IN)) {
  dir.create(dirname(ADHERENCE_IN), recursive = TRUE, showWarnings = FALSE)
  write_csv(
    tibble(
      scenario = c("base","avg_gap","moderate","stress"),
      adherence = c(1.0, 0.95, 0.85, 0.70),
      source_note = c(
        "Assume programmatic coverage approximates effective treatment on average",
        "Mean coverage–compliance gap sensitivity (see compliance systematic review in Methods)",
        "Moderate adherence shortfall sensitivity",
        "Stress test adherence shortfall sensitivity"
      )
    ),
    ADHERENCE_IN
  )
  stop("Missing ADHERENCE_IN. Template written to: ", ADHERENCE_IN,
       "\nReview/edit it, then rerun.", call. = FALSE)
}

adh_tbl <- read_csv(ADHERENCE_IN, show_col_types = FALSE) %>%
  mutate(
    scenario = as.character(scenario),
    adherence = as.numeric(adherence),
    source_note = as.character(source_note)
  )

need_cols(adh_tbl, c("scenario","adherence","source_note"), basename(ADHERENCE_IN))
assert_no_dups(adh_tbl, c("scenario"), "adherence_scenarios")

if (!ADHERENCE_SCENARIO %in% adh_tbl$scenario) {
  stop(
    "ADHERENCE_SCENARIO not found in ", ADHERENCE_IN,
    ". Got: ", ADHERENCE_SCENARIO,
    ". Allowed: ", paste(sort(unique(adh_tbl$scenario)), collapse = ", "),
    call. = FALSE
  )
}

ADHERENCE <- adh_tbl %>%
  filter(scenario == ADHERENCE_SCENARIO) %>%
  pull(adherence) %>%
  .[[1]]

if (!is.finite(ADHERENCE) || ADHERENCE <= 0 || ADHERENCE > 1) {
  stop("ADHERENCE must be in (0,1]. Got: ", ADHERENCE, call. = FALSE)
}

write_csv(
  adh_tbl %>% filter(scenario == ADHERENCE_SCENARIO),
  file.path(OUT_QC_DIR, "adherence_used.csv")
)

# ---------------------------
# 1c) Load hyperprior scenario
# ---------------------------
if (!file.exists(HYPERPRIOR_IN)) {
  dir.create(dirname(HYPERPRIOR_IN), recursive = TRUE, showWarnings = FALSE)
  write_csv(
    tibble(
      scenario = c("base","tight","loose"),
      beta0_sd = c(2.5, 2.5, 2.5),
      sigma_proc_rate = c(3, 6, 1.5),
      sigma_v_rate    = c(3, 6, 1.5),
      sigma_u_rate    = c(3, 6, 1.5),
      sigma_init_rate = c(2, 4, 1),
      sigma_obs_rate  = c(5, 10, 2.5),
      source_note = c(
        "Weakly-informative regularization; prior-predictive + sensitivity",
        "Tighter shrinkage (rates x2)",
        "Looser shrinkage (rates x0.5)"
      )
    ),
    HYPERPRIOR_IN
  )
  stop("Missing HYPERPRIOR_IN. Template written to: ", HYPERPRIOR_IN,
       "\nReview/edit it, then rerun.", call. = FALSE)
}

hyp_tbl <- read_csv(HYPERPRIOR_IN, show_col_types = FALSE) %>%
  mutate(
    scenario = as.character(scenario),
    beta0_sd = as.numeric(beta0_sd),
    sigma_proc_rate = as.numeric(sigma_proc_rate),
    sigma_v_rate    = as.numeric(sigma_v_rate),
    sigma_u_rate    = as.numeric(sigma_u_rate),
    sigma_init_rate = as.numeric(sigma_init_rate),
    sigma_obs_rate  = as.numeric(sigma_obs_rate),
    source_note = as.character(source_note)
  )

need_cols(
  hyp_tbl,
  c("scenario","beta0_sd","sigma_proc_rate","sigma_v_rate","sigma_u_rate",
    "sigma_init_rate","sigma_obs_rate","source_note"),
  basename(HYPERPRIOR_IN)
)
assert_no_dups(hyp_tbl, c("scenario"), "hyperprior_scenarios")

if (!HYPERPRIOR_SCENARIO %in% hyp_tbl$scenario) {
  stop(
    "HYPERPRIOR_SCENARIO not found in ", HYPERPRIOR_IN,
    ". Got: ", HYPERPRIOR_SCENARIO,
    ". Allowed: ", paste(sort(unique(hyp_tbl$scenario)), collapse = ", "),
    call. = FALSE
  )
}

hyp_row <- hyp_tbl %>% filter(scenario == HYPERPRIOR_SCENARIO)

BETA0_SD        <- hyp_row$beta0_sd[[1]]
SIGMA_PROC_RATE <- hyp_row$sigma_proc_rate[[1]]
SIGMA_V_RATE    <- hyp_row$sigma_v_rate[[1]]
SIGMA_U_RATE    <- hyp_row$sigma_u_rate[[1]]
SIGMA_INIT_RATE <- hyp_row$sigma_init_rate[[1]]
SIGMA_OBS_RATE  <- hyp_row$sigma_obs_rate[[1]]

vals <- c(BETA0_SD, SIGMA_PROC_RATE, SIGMA_V_RATE, SIGMA_U_RATE, SIGMA_INIT_RATE, SIGMA_OBS_RATE)
if (any(!is.finite(vals)) || any(vals <= 0)) stop("Hyperprior values must be finite and > 0.", call. = FALSE)

write_csv(hyp_row, file.path(OUT_QC_DIR, "hyperpriors_used.csv"))

# ---------------------------
# 2) Load panel + validate contract
# ---------------------------
if (!file.exists(PANEL_IN)) stop("Missing input panel. Run script 01 first: ", PANEL_IN, call. = FALSE)

panel <- read_csv(PANEL_IN, show_col_types = FALSE) %>%
  mutate(
    iso3 = toupper(as.character(iso3)),
    year = as.integer(year),
    parasite = as.character(parasite)
  )

need_cols(panel, c("iso3","year","parasite","prev_obs","pop_sac","pc_cov","pc_status"), "sac_panel_merged.csv")
assert_no_dups(panel, c("iso3","year","parasite"), "sac_panel_merged")

if (any(!is.na(panel$prev_obs) & (panel$prev_obs <= 0 | panel$prev_obs >= 1))) {
  bad <- panel %>% filter(!is.na(prev_obs) & (prev_obs <= 0 | prev_obs >= 1)) %>% slice_head(n = 200)
  write_csv(bad, file.path(OUT_QC_DIR, "bad_prev_obs_outside_open_interval.csv"))
  stop("prev_obs must be strictly between 0 and 1 before clamp. QC: bad_prev_obs_outside_open_interval.csv",
       call. = FALSE)
}

stop_if_outside_01(panel$pc_cov, "pc_cov", allow_na = TRUE)
if (any(!is.na(panel$pop_sac) & panel$pop_sac < 0)) stop("pop_sac has negative values.", call. = FALSE)

panel <- panel %>%
  mutate(
    prev_obs = clip01(prev_obs, EPS),
    prev_lo  = if ("prev_lo" %in% names(.)) if_else(!is.na(prev_lo), clip01(prev_lo, EPS), NA_real_) else NA_real_,
    prev_hi  = if ("prev_hi" %in% names(.)) if_else(!is.na(prev_hi), clip01(prev_hi, EPS), NA_real_) else NA_real_,
    y_obs    = logit(prev_obs)
  )

# ---------------------------
# 3) Add DALYs (required)
# ---------------------------
if (!file.exists(DALY_IN)) {
  template <- panel %>%
    distinct(iso3, year, parasite) %>%
    arrange(iso3, year, parasite) %>%
    mutate(daly_gbd = NA_real_)
  dir.create(dirname(DALY_IN), recursive = TRUE, showWarnings = FALSE)
  write_csv(template, DALY_IN)
  stop("DALY file missing. Template written to: ", DALY_IN,
       "\nFill daly_gbd (>=0) and rerun.",
       call. = FALSE)
}

daly <- read_csv(DALY_IN, show_col_types = FALSE) %>%
  mutate(
    iso3 = toupper(as.character(iso3)),
    year = as.integer(year),
    parasite = as.character(parasite),
    daly_gbd = as.numeric(daly_gbd)
  )

need_cols(daly, c("iso3","year","parasite","daly_gbd"), basename(DALY_IN))
assert_no_dups(daly, c("iso3","year","parasite"), basename(DALY_IN))

if (any(is.na(daly$daly_gbd))) {
  write_csv(daly %>% filter(is.na(daly_gbd)),
            file.path(OUT_QC_DIR, "missing_daly_gbd_rows.csv"))
  stop("daly_gbd contains NA. QC: missing_daly_gbd_rows.csv", call. = FALSE)
}
if (any(daly$daly_gbd < 0, na.rm = TRUE)) stop("daly_gbd contains negative values.", call. = FALSE)

panel_before_daly_join <- panel
panel <- panel %>% left_join(daly, by = c("iso3","year","parasite"))
log_join_result(panel_before_daly_join, panel %>% dplyr::filter(!is.na(daly_gbd)), "panel_to_daly_join")

miss_d <- panel %>% filter(is.na(daly_gbd)) %>% distinct(iso3, year, parasite)
if (nrow(miss_d) > 0) {
  write_csv(miss_d, file.path(OUT_QC_DIR, "panel_keys_missing_daly.csv"))
  stop("Some panel keys missing DALYs. QC: panel_keys_missing_daly.csv", call. = FALSE)
}

# ---------------------------
# 4) Observation error SD (logit scale)
# ---------------------------
panel <- panel %>%
  mutate(
    lo2 = pmin(prev_lo, prev_hi, na.rm = TRUE),
    hi2 = pmax(prev_lo, prev_hi, na.rm = TRUE),
    has_ui = !is.na(lo2) & !is.na(hi2) & (hi2 > lo2),
    sigma_ui = if_else(
      has_ui,
      (logit(clip01(hi2, EPS)) - logit(clip01(lo2, EPS))) / 3.92,
      NA_real_
    ),
    sigma_fallback = if ("sigma_logit" %in% names(.)) as.numeric(sigma_logit) else NA_real_,
    sigma_obs_raw  = coalesce(sigma_ui, sigma_fallback),
    has_sigma_obs  = if_else(is.finite(sigma_obs_raw) & sigma_obs_raw > 0, 1L, 0L),
    sigma_obs_data = if_else(has_sigma_obs == 1L, sigma_obs_raw, 1.0)
  ) %>%
  select(-lo2, -hi2, -sigma_obs_raw)

write_csv(
  panel %>%
    summarise(
      n = n(),
      frac_has_ui = mean(has_ui, na.rm = TRUE),
      frac_has_sigma = mean(has_sigma_obs == 1L, na.rm = TRUE),
      sigma_ui_median = median(sigma_ui, na.rm = TRUE),
      sigma_fallback_median = median(sigma_fallback, na.rm = TRUE),
      sigma_obs_median = median(sigma_obs_data, na.rm = TRUE)
    ),
  file.path(OUT_QC_DIR, "sigma_obs_summary.csv")
)

# ---------------------------
# 5) PC coverage fill (panel years only)
# ---------------------------
panel <- panel %>%
  group_by(iso3, parasite) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    pc_cov_filled = case_when(
      pc_status == "no_pc_required" ~ 0,
      TRUE ~ as.numeric(pc_cov)
    ),
    pc_cov_filled = if (PC_MISSING_POLICY == "zero") {
      if_else(is.na(pc_cov_filled), 0, pc_cov_filled)
    } else {
      x <- fill_locf_vec(pc_cov_filled)
      if_else(is.na(x), 0, x)
    },
    pc_cov_filled = pmin(pmax(pc_cov_filled, 0), 1)
  ) %>%
  ungroup()

write_csv(
  panel %>%
    summarise(
      frac_pc_missing_before = mean(is.na(pc_cov), na.rm = TRUE),
      frac_pc_zero_after = mean(pc_cov_filled == 0, na.rm = TRUE),
      frac_pc_one_after  = mean(pc_cov_filled == 1, na.rm = TRUE)
    ),
  file.path(OUT_QC_DIR, "pc_fill_summary.csv")
)

# ---------------------------
# 6) Build full year grid + indices
# ---------------------------
years_in <- sort(unique(panel$year))
if (length(years_in) == 0) stop("No years found in panel.", call. = FALSE)

years_full <- seq(min(years_in), max(years_in))
missing_years <- setdiff(years_full, years_in)
if (length(missing_years) > 0) {
  write_csv(tibble(missing_year = missing_years), file.path(OUT_QC_DIR, "missing_years_in_input_panel.csv"))
  message("NOTE: Panel missing some years; grid includes them. QC: ",
          file.path(OUT_QC_DIR, "missing_years_in_input_panel.csv"))
}

countries <- sort(unique(panel$iso3))
parasites <- sort(unique(panel$parasite))

C <- length(countries)
P <- length(parasites)
T <- length(years_full)

# --- Load parasite-specific k priors aligned to parasite order ---
if (!file.exists(K_PRIORS_IN)) {
  template <- tibble(parasite = parasites, k0 = NA_real_, sigma_k = NA_real_)
  dir.create(dirname(K_PRIORS_IN), recursive = TRUE, showWarnings = FALSE)
  write_csv(template, K_PRIORS_IN)
  stop(
    "K priors file missing. Template written to: ", K_PRIORS_IN,
    "\nFill k0 and sigma_k (>0) for each parasite and rerun.",
    call. = FALSE
  )
}

kpri <- read_csv(K_PRIORS_IN, show_col_types = FALSE) %>%
  mutate(
    parasite = as.character(parasite),
    k0 = as.numeric(k0),
    sigma_k = as.numeric(sigma_k)
  )

need_cols(kpri, c("parasite","k0","sigma_k"), basename(K_PRIORS_IN))

kpri2 <- tibble(parasite = parasites) %>% left_join(kpri, by = "parasite")

if (any(is.na(kpri2$k0)) || any(is.na(kpri2$sigma_k))) {
  write_csv(kpri2, file.path(OUT_QC_DIR, "missing_k_priors_after_join.csv"))
  stop("Some parasites are missing k priors. QC: missing_k_priors_after_join.csv", call. = FALSE)
}
if (any(!is.finite(kpri2$k0) | kpri2$k0 <= 0) ||
    any(!is.finite(kpri2$sigma_k) | kpri2$sigma_k <= 0)) {
  write_csv(kpri2, file.path(OUT_QC_DIR, "bad_k_priors_nonpositive.csv"))
  stop("Invalid k priors (non-finite or <=0). QC: bad_k_priors_nonpositive.csv", call. = FALSE)
}

k0_vec <- kpri2$k0
sigma_k_vec <- kpri2$sigma_k

write_csv(
  tibble(parasite = parasites, k0 = k0_vec, sigma_k = sigma_k_vec),
  file.path(OUT_QC_DIR, "k_priors_used_by_stan.csv")
)

# --- series map ---
series_map <- panel %>%
  distinct(iso3, parasite) %>%
  arrange(iso3, parasite) %>%
  mutate(
    series_id   = row_number(),
    country_id  = match(iso3, countries),
    parasite_id = match(parasite, parasites)
  )

S <- nrow(series_map)

grid <- series_map %>%
  select(iso3, parasite, series_id, country_id, parasite_id) %>%
  crossing(year = years_full) %>%
  left_join(panel, by = c("iso3","parasite","year")) %>%
  mutate(t = match(year, years_full)) %>%
  arrange(series_id, t)

# Apply PC missing policy on FULL grid (including inserted years)
grid <- grid %>%
  group_by(series_id) %>%
  arrange(t, .by_group = TRUE) %>%
  mutate(
    pc_grid = pc_cov_filled,
    pc_grid = if (PC_MISSING_POLICY == "zero") {
      if_else(is.na(pc_grid), 0, pc_grid)
    } else {
      x <- fill_locf_vec(pc_grid)
      if_else(is.na(x), 0, x)
    },
    pc_grid = pmin(pmax(pc_grid, 0), 1)
  ) %>%
  ungroup()

# =============================================================================
# ROUNDS LOGIC (TRUTH-FIRST, CONSISTENT WITH ANNUAL NATIONAL COVERAGE)
# =============================================================================
rounds_col <- NULL
if ("rounds_use" %in% names(grid)) rounds_col <- "rounds_use"
if (is.null(rounds_col) && "rounds" %in% names(grid)) rounds_col <- "rounds"

if (!is.null(rounds_col)) {
  grid <- grid %>%
    group_by(series_id) %>%
    arrange(t, .by_group = TRUE) %>%
    mutate(
      rounds_grid = as.numeric(.data[[rounds_col]]),
      rounds_grid = if (PC_MISSING_POLICY == "zero") {
        if_else(is.na(rounds_grid), 0, rounds_grid)
      } else {
        x <- fill_locf_vec(rounds_grid)
        if_else(is.na(x), 0, x)
      },
      rounds_grid = pmax(rounds_grid, 0),
      rounds_grid = if_else(pc_grid > 0, rounds_grid, 0)
    ) %>%
    ungroup()
  
  if (sum(grid$pc_grid > 0) > 0 && all(grid$rounds_grid[grid$pc_grid > 0] == 0)) {
    grid <- grid %>% mutate(rounds_grid = if_else(pc_grid > 0, as.numeric(ROUNDS_CONSTANT), 0))
  }
} else {
  grid <- grid %>% mutate(rounds_grid = if_else(pc_grid > 0, as.numeric(ROUNDS_CONSTANT), 0))
}

# Observations list (ONLY where prevalence exists)
obs <- grid %>%
  filter(!is.na(prev_obs)) %>%
  transmute(
    s_obs = as.integer(series_id),
    t_obs = as.integer(t),
    y_obs = y_obs,
    sigma_obs_data = sigma_obs_data,
    has_sigma_obs  = as.integer(has_sigma_obs)
  )

if (nrow(obs) == 0) stop("No prevalence observations found (N_obs=0). Stan requires N_obs>=1.", call. = FALSE)

# Matrices S x T for PC and rounds
pc_mat <- grid %>%
  select(series_id, t, pc_grid) %>%
  arrange(series_id, t) %>%
  pivot_wider(names_from = t, values_from = pc_grid) %>%
  arrange(series_id) %>%
  select(-series_id) %>%
  as.matrix()

rounds_mat <- grid %>%
  select(series_id, t, rounds_grid) %>%
  arrange(series_id, t) %>%
  pivot_wider(names_from = t, values_from = rounds_grid, values_fill = 0) %>%
  arrange(series_id) %>%
  select(-series_id) %>%
  as.matrix()

# DALY matrices S x T (Stan requires no NA)
daly_mat <- grid %>%
  mutate(
    has_daly = if_else(!is.na(daly_gbd), 1L, 0L),
    daly0    = if_else(!is.na(daly_gbd), as.numeric(daly_gbd), 0.0)
  ) %>%
  select(series_id, t, daly0) %>%
  pivot_wider(names_from = t, values_from = daly0, values_fill = 0.0) %>%
  arrange(series_id) %>%
  select(-series_id) %>%
  as.matrix()

has_daly_mat <- grid %>%
  mutate(has_daly = if_else(!is.na(daly_gbd), 1L, 0L)) %>%
  select(series_id, t, has_daly) %>%
  pivot_wider(names_from = t, values_from = has_daly, values_fill = 0L) %>%
  arrange(series_id) %>%
  select(-series_id) %>%
  as.matrix()

# Validate matrices
if (any(!is.finite(pc_mat))) stop("pc_mat has non-finite values.", call. = FALSE)
if (any(pc_mat < 0 | pc_mat > 1)) stop("pc_mat outside [0,1].", call. = FALSE)
if (any(!is.finite(rounds_mat)) || any(rounds_mat < 0)) stop("rounds_mat invalid (<0 or non-finite).", call. = FALSE)
if (any(!is.finite(daly_mat)) || any(daly_mat < 0)) stop("daly_mat invalid (<0 or non-finite).", call. = FALSE)
storage.mode(has_daly_mat) <- "integer"

# ---------------------------
# 7) Optional adjacency edges -> ICAR inputs
# ---------------------------
use_car <- 0L
E <- 0L
node1 <- integer(0)
node2 <- integer(0)

if (file.exists(ADJ_EDGES)) {
  edges0 <- read_csv(ADJ_EDGES, show_col_types = FALSE)
  need_cols(edges0, c("iso3_from","iso3_to"), basename(ADJ_EDGES))
  
  edges <- edges0 %>%
    transmute(
      iso3_from = toupper(as.character(iso3_from)),
      iso3_to   = toupper(as.character(iso3_to))
    ) %>%
    filter(!is.na(iso3_from), !is.na(iso3_to), iso3_from != iso3_to) %>%
    distinct() %>%
    mutate(
      i = match(iso3_from, countries),
      j = match(iso3_to, countries)
    )
  
  dropped <- edges %>% filter(is.na(i) | is.na(j))
  if (nrow(dropped) > 0) write_csv(dropped, file.path(OUT_QC_DIR, "adj_edges_dropped_not_in_panel.csv"))
  
  edges <- edges %>% filter(!is.na(i), !is.na(j))
  
  if (nrow(edges) > 0) {
    use_car <- 1L
    
    # ---- FIX: keep UNDIRECTED edges ONCE (standard ICAR scaling) ----
    edges2 <- edges %>%
      transmute(i = pmin(i, j), j = pmax(i, j)) %>%
      distinct()
    
    node1 <- as.integer(edges2$i)
    node2 <- as.integer(edges2$j)
    E <- as.integer(nrow(edges2))
    
    deg <- tibble(node = c(node1, node2)) %>%
      count(node, name = "degree") %>%
      right_join(tibble(node = seq_len(C)), by = "node") %>%
      mutate(degree = replace_na(degree, 0L)) %>%
      mutate(iso3 = countries[node]) %>%
      select(iso3, degree)
    
    write_csv(deg, file.path(OUT_QC_DIR, "adj_degree_by_country.csv"))
  }
}

use_car <- as.integer(use_car)
E <- as.integer(E)

# ---------------------------
# 8) ERR priors (prefer alpha/beta if present)
# ---------------------------
a_ERR <- rep(ERR_A_DEFAULT, P)
b_ERR <- rep(ERR_B_DEFAULT, P)

if (!file.exists(DRUG_ERR_IN)) stop("Missing DRUG_ERR_IN: ", DRUG_ERR_IN, call. = FALSE)

derr <- read_csv(DRUG_ERR_IN, show_col_types = FALSE) %>%
  mutate(parasite = as.character(parasite))

has_ab <- all(c("alpha","beta") %in% names(derr))
has_ci <- all(c("err_mean","err_lo","err_hi") %in% names(derr))
if (!has_ab && !has_ci) {
  stop("DRUG_ERR_IN must contain (alpha,beta) or (err_mean,err_lo,err_hi). File: ", DRUG_ERR_IN, call. = FALSE)
}

err_prior_table <- tibble(parasite = parasites) %>%
  left_join(derr, by = "parasite")

if (has_ab) {
  err_prior_table <- err_prior_table %>%
    mutate(alpha = as.numeric(alpha),
           beta  = as.numeric(beta))
  
  if (any(!is.finite(err_prior_table$alpha) | err_prior_table$alpha <= 0) ||
      any(!is.finite(err_prior_table$beta)  | err_prior_table$beta  <= 0)) {
    write_csv(err_prior_table, file.path(OUT_QC_DIR, "bad_err_alpha_beta.csv"))
    stop("Invalid alpha/beta in DRUG_ERR_IN. QC: bad_err_alpha_beta.csv", call. = FALSE)
  }
  
  a_ERR <- err_prior_table$alpha
  b_ERR <- err_prior_table$beta
} else {
  err_prior_table <- err_prior_table %>%
    mutate(
      err_mean = pmin(pmax(as.numeric(err_mean), EPS), 1 - EPS),
      err_lo   = pmin(pmax(as.numeric(err_lo),   EPS), 1 - EPS),
      err_hi   = pmin(pmax(as.numeric(err_hi),   EPS), 1 - EPS)
    ) %>%
    rowwise() %>%
    mutate(
      pars  = list(beta_from_mean_ci(err_mean, err_lo, err_hi)),
      alpha = pars[["alpha"]],
      beta  = pars[["beta"]]
    ) %>%
    ungroup()
  
  a_ERR <- err_prior_table$alpha
  b_ERR <- err_prior_table$beta
}

write_csv(
  tibble(parasite = parasites, a_ERR = a_ERR, b_ERR = b_ERR),
  file.path(OUT_QC_DIR, "err_priors_used_by_stan.csv")
)

# ---------------------------
# 8b) Force integer storage for Stan int inputs
# ---------------------------
series_map <- series_map %>%
  mutate(
    country_id  = as.integer(country_id),
    parasite_id = as.integer(parasite_id),
    series_id   = as.integer(series_id)
  )

node1 <- as.integer(node1)
node2 <- as.integer(node2)

# ---------------------------
# 9) Save outputs
# ---------------------------
write_csv(grid, OUT_PANEL)
write_csv(grid, OUT_PANEL_ALIAS)
message("Wrote: ", OUT_PANEL)
message("Wrote: ", OUT_PANEL_ALIAS)

stan_data <- list(
  C = as.integer(C), P = as.integer(P), T = as.integer(T), S = as.integer(S),
  
  country_id  = series_map$country_id,
  parasite_id = series_map$parasite_id,
  
  pc = pc_mat,
  rounds = rounds_mat,
  adherence = ADHERENCE,
  
  obs_post_pc = as.integer(OBS_POST_PC),
  
  daly = daly_mat,
  has_daly = has_daly_mat,
  
  N_obs = as.integer(nrow(obs)),
  s_obs = obs$s_obs,
  t_obs = obs$t_obs,
  y_obs = obs$y_obs,
  sigma_obs_data = obs$sigma_obs_data,
  has_sigma_obs  = obs$has_sigma_obs,
  
  use_car = use_car,
  E = E,
  node1 = node1,
  node2 = node2,
  
  eps = EPS,
  
  k0 = k0_vec,
  sigma_k = sigma_k_vec,
  
  a_ERR = a_ERR,
  b_ERR = b_ERR,
  
  beta0_sd = BETA0_SD,
  sigma_proc_rate = SIGMA_PROC_RATE,
  sigma_v_rate = SIGMA_V_RATE,
  sigma_u_rate = SIGMA_U_RATE,
  sigma_init_rate = SIGMA_INIT_RATE,
  sigma_obs_rate = SIGMA_OBS_RATE
)

saveRDS(stan_data, OUT_STAN)
saveRDS(stan_data, OUT_STAN_ALIAS)

saveRDS(
  list(
    years = years_full,
    countries = countries,
    parasites = parasites,
    series_map = series_map,
    tag = tag,
    adherence_scenario = ADHERENCE_SCENARIO,
    hyperprior_scenario = HYPERPRIOR_SCENARIO,
    obs_post_pc = OBS_POST_PC,
    pc_missing_policy = PC_MISSING_POLICY,
    rounds_constant = ROUNDS_CONSTANT
  ),
  OUT_LEVEL
)

writeLines(capture.output(sessionInfo()),
           file.path(OUT_QC_DIR, "sessionInfo_02_prepare_stan_inputs.txt"))

message("Wrote: ", OUT_STAN)
message("Wrote: ", OUT_STAN_ALIAS)
message("Wrote: ", OUT_LEVEL)
message("QC dir: ", OUT_QC_DIR)
message("DONE.")