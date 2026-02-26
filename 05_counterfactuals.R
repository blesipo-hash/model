# =============================================================================
# 05_counterfactuals.R
# Counterfactual simulation (WHO-PC vs PC75 vs No-PC) using posterior draws
# =============================================================================

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ---------------------------
# Config
# ---------------------------
FIT_RDS     <- "data/processed/fit_sac.rds"
LEVELS_RDS  <- "data/processed/levels_sac.rds"
MODEL_PANEL <- "data/processed/model_input/sac_panel_model.csv"

OUT_DIR <- "outputs/counterfactuals"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EPS <- 1e-6
clip01 <- function(p, eps = EPS) pmin(pmax(p, eps), 1 - eps)
need_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(name, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

# How many posterior draws to use for simulation (tradeoff: speed vs precision)
N_DRAWS_SIM <- 400

# Scenarios
SCENARIOS <- c("realized", "scaleup_set_75", "scaleup_floor_75", "scaleup_100", "NO_PC")

# Save a small sample of draws for debugging/plots (set 0 to skip)
SAVE_DRAWS_SAMPLE <- 1
N_SAVE_DRAWS <- 50

set.seed(1)

# ---------------------------
# Load
# ---------------------------
stopifnot(file.exists(FIT_RDS), file.exists(LEVELS_RDS), file.exists(MODEL_PANEL))

fit <- readRDS(FIT_RDS)
lev <- readRDS(LEVELS_RDS)

grid <- readr::read_csv(MODEL_PANEL, show_col_types = FALSE) %>%
  mutate(
    iso3 = toupper(iso3),
    year = as.integer(year),
    parasite = as.character(parasite),
    series_id = as.integer(series_id),
    t = as.integer(t)
  )

need_cols(grid, c("series_id","t","pc_cov_filled","rounds","iso3","year","parasite"), "sac_panel_model")
years <- lev$years
countries <- lev$countries
parasites <- lev$parasites
series_map <- lev$series_map

S <- nrow(series_map)
Tt <- length(years)
P <- length(parasites)

# Build SxT matrices from grid
pc_who <- grid %>%
  select(series_id, t, pc_cov_filled) %>%
  pivot_wider(names_from = t, values_from = pc_cov_filled) %>%
  arrange(series_id) %>%
  select(-series_id) %>%
  as.matrix()

rounds <- grid %>%
  select(series_id, t, rounds) %>%
  pivot_wider(names_from = t, values_from = rounds) %>%
  arrange(series_id) %>%
  select(-series_id) %>%
  as.matrix()

# Scenario pc matrices
pc_list <- list(
  realized = pc_who,
  scaleup_set_75 = matrix(0.75, nrow = S, ncol = Tt),
  scaleup_floor_75 = pmax(pc_who, 0.75),
  scaleup_100 = matrix(1.0, nrow = S, ncol = Tt),
  NO_PC  = matrix(0.0,  nrow = S, ncol = Tt)
)

# Series mappings
country_id  <- series_map$country_id
parasite_id <- series_map$parasite_id

# ---------------------------
# Helper: intensity -> prevalence (vectorized)
# p = 1 - (k/(k+W))^k
# ---------------------------
prev_from_W <- function(W, k) {
  # W and k are numeric vectors of same length
  1 - (k / (k + W))^k
}

# ---------------------------
# Extract posterior draws needed for simulation
# (we don't use latent w from Stan; we simulate forward ourselves)
# Need: beta0[p], v[c,p], u[c,p] or u=0, r[p], ERR[p], k[p], sigma_proc[p], sigma_init, w0_raw[s], eps_raw[s,t]
# ---------------------------

# Identify variables
var_needed <- c(
  "beta0", "r", "ERR", "k", "sigma_proc", "sigma_init",
  "v", "u", "use_car", "adherence"  # adherence is data; may not exist as param
)

# w0_raw and eps_raw are large; pull them selectively
# We'll pull full arrays but subsample draws to keep it manageable
draws_all <- fit$draws()

# subsample draws
nd_total <- posterior::ndraws(draws_all)
keep <- seq_len(min(N_DRAWS_SIM, nd_total))
draws <- posterior::subset_draws(draws_all, draw_indices = keep)

# Pull parameters as matrices/arrays
beta0      <- posterior::as_draws_matrix(draws, variable = "beta0")           # [draw, P]
r_draw     <- posterior::as_draws_matrix(draws, variable = "r")               # [draw, P]
ERR_draw   <- posterior::as_draws_matrix(draws, variable = "ERR")             # [draw, P]
k_draw     <- posterior::as_draws_matrix(draws, variable = "k")               # [draw, P]
sigproc    <- posterior::as_draws_matrix(draws, variable = "sigma_proc")      # [draw, P]
sigma_init <- posterior::as_draws_vector(draws, variable = "sigma_init")      # [draw]

# v[c,p]
v_mat <- posterior::as_draws_matrix(draws, variable = "v")                    # columns like v[1,1]
# u[c,p] (will be 0 in non-spatial mode; still present in generated/transformed)
u_mat <- posterior::as_draws_matrix(draws, variable = "u")

# w0_raw[s]
w0_raw <- posterior::as_draws_matrix(draws, variable = "w0_raw")              # [draw, S]

# eps_raw[s,t] for t=1..T-1
eps_raw <- posterior::as_draws_matrix(draws, variable = "eps_raw")            # [draw, S*(T-1)]

D <- nrow(beta0)

# Helpers to pull v/u into arrays [D, C, P]
parse_cp_array <- function(mat_draws, C, P, prefix) {
  # mat_draws columns are like prefix[i,j]
  out <- array(0, dim = c(nrow(mat_draws), C, P))
  cn <- colnames(mat_draws)
  m <- stringr::str_match(cn, paste0("^", prefix, "\\[(\\d+),(\\d+)\\]$"))
  idx <- which(!is.na(m[,1]))
  for (k0 in idx) {
    i <- as.integer(m[k0,2])
    j <- as.integer(m[k0,3])
    out[, i, j] <- mat_draws[, k0]
  }
  out
}

v_arr <- parse_cp_array(v_mat, C = length(countries), P = P, prefix = "v")
u_arr <- parse_cp_array(u_mat, C = length(countries), P = P, prefix = "u")

# eps_raw into [D, S, T-1]
eps_arr <- array(0, dim = c(D, S, Tt - 1))
cn_eps <- colnames(eps_raw)
m_eps <- stringr::str_match(cn_eps, "^eps_raw\\[(\\d+),(\\d+)\\]$")
idx_eps <- which(!is.na(m_eps[,1]))
for (k0 in idx_eps) {
  s <- as.integer(m_eps[k0,2])
  tt <- as.integer(m_eps[k0,3])
  eps_arr[, s, tt] <- eps_raw[, k0]
}

# ---------------------------
# Simulate forward for each scenario and draw
# ---------------------------

# storage: we’ll summarize on the fly to avoid huge memory,
# but to keep it straightforward, we’ll compute p for all and then summarize.
# Array p: [scenario, draw, series, time]
p_out <- array(NA_real_, dim = c(length(SCENARIOS), D, S, Tt),
               dimnames = list(SCENARIOS, NULL, NULL, NULL))

for (sc_i in seq_along(SCENARIOS)) {
  sc <- SCENARIOS[sc_i]
  pc_mat <- pc_list[[sc]]
  
  for (d in 1:D) {
    # parasite-level params for draw d
    k_p   <- as.numeric(k_draw[d, ])
    r_p   <- as.numeric(r_draw[d, ])
    err_p <- as.numeric(ERR_draw[d, ])
    sp_p  <- as.numeric(sigproc[d, ])
    s_init <- as.numeric(sigma_init[d])
    
    # mu_series[s] = beta0[p] + v[c,p] + u[c,p]
    mu_s <- numeric(S)
    for (s in 1:S) {
      c <- country_id[s]
      p <- parasite_id[s]
      mu_s[s] <- as.numeric(beta0[d, p]) + v_arr[d, c, p] + u_arr[d, c, p]
    }
    
    # latent w[s,t]
    w <- matrix(NA_real_, nrow = S, ncol = Tt)
    w[, 1] <- mu_s + s_init * w0_raw[d, ]
    
    for (t in 1:(Tt - 1)) {
      # apply PC pulse at t
      cov_t <- pmin(pmax(pc_mat[, t], 0), 1)
      mult  <- pmax(1e-9, 1 - err_p[parasite_id] * cov_t)   # elementwise by series' parasite
      w_post <- w[, t] + rounds[, t] * log(mult)
      
      # mean reversion + process noise
      rr <- r_p[parasite_id]
      sp <- sp_p[parasite_id]
      w[, t + 1] <- w_post + rr * (mu_s - w_post) + sp * eps_arr[d, , t]
    }
    
    # convert to prevalence
    W <- exp(w)
    kk <- k_p[parasite_id]
    p_hat <- 1 - (kk / (kk + W))^kk
    p_out[sc_i, d, , ] <- clip01(p_hat, EPS)
  }
  
  message("Sim done: ", sc)
}

# ---------------------------
# Summarize by iso3-year-parasite-scenario
# ---------------------------
series_df <- series_map %>%
  transmute(
    series_id = series_id,
    iso3 = iso3,
    parasite = parasite
  )

# build long summary without exploding too hard:
# summarise per (scenario, series, t) over draws
summ_list <- vector("list", length(SCENARIOS))
names(summ_list) <- SCENARIOS

for (sc_i in seq_along(SCENARIOS)) {
  sc <- SCENARIOS[sc_i]
  # for each series/time, compute q025/q50/q975 across draws
  q025 <- apply(p_out[sc_i, , , , drop = FALSE], c(3,4), quantile, probs = 0.025, na.rm = TRUE)
  q050 <- apply(p_out[sc_i, , , , drop = FALSE], c(3,4), quantile, probs = 0.50,  na.rm = TRUE)
  q975 <- apply(p_out[sc_i, , , , drop = FALSE], c(3,4), quantile, probs = 0.975, na.rm = TRUE)
  
  # q matrices are [S, T]
  df <- as_tibble(expand.grid(
    series_id = seq_len(S),
    t = seq_len(Tt)
  )) %>%
    mutate(
      scenario = sc,
      year = years[t],
      prev_q025 = as.vector(q025),
      prev_q50  = as.vector(q050),
      prev_q975 = as.vector(q975)
    ) %>%
    left_join(series_df, by = "series_id") %>%
    select(iso3, year, parasite, scenario, prev_q50, prev_q025, prev_q975)
  
  summ_list[[sc]] <- df
}

cf_summary <- bind_rows(summ_list) %>%
  arrange(scenario, parasite, iso3, year)

readr::write_csv(cf_summary, file.path(OUT_DIR, "cf_prevalence_summary.csv"))
message("Wrote: ", file.path(OUT_DIR, "cf_prevalence_summary.csv"))

# ---------------------------
# Optional: save small sample of draw-level trajectories for debugging
# ---------------------------
if (SAVE_DRAWS_SAMPLE == 1) {
  keep_d <- seq_len(min(N_SAVE_DRAWS, D))
  sample_long <- as_tibble(expand.grid(
    scenario = SCENARIOS,
    draw = keep_d,
    series_id = seq_len(S),
    t = seq_len(Tt)
  )) %>%
    mutate(
      year = years[t],
      prev = as.vector(p_out[scenario, draw, series_id, t])
    ) %>%
    left_join(series_df, by = "series_id") %>%
    select(scenario, draw, iso3, parasite, year, prev)
  
  readr::write_csv(sample_long, file.path(OUT_DIR, "cf_prevalence_draws_sample.csv"))
  message("Wrote: ", file.path(OUT_DIR, "cf_prevalence_draws_sample.csv"))
}

message("DONE.")