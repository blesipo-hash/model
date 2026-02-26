# =============================================================================
# 04_ppc.R
# Posterior predictive checks for SAC transmission inference model
# Inputs:
#   - data/processed/fit_sac.rds
#   - data/processed/stan_data_sac.rds
#   - data/processed/levels_sac.rds
#   - data/processed/model_input/sac_panel_model.csv
# Outputs:
#   - outputs/ppc/ppc_pred_summary_selected.csv
#   - outputs/ppc/ppc_scatter_obs_vs_pred.png
#   - outputs/ppc/ppc_timeseries_examples.png
#   - outputs/ppc/ppc_residuals_by_bin.png
#   - outputs/ppc/ppc_diagnostics.txt
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

# ---------------------------
# Config
# ---------------------------
FIT_RDS       <- "data/processed/fit_sac.rds"
STAN_DATA_RDS <- "data/processed/stan_data_sac.rds"
LEVELS_RDS    <- "data/processed/levels_sac.rds"
MODEL_PANEL   <- "data/processed/model_input/sac_panel_model.csv"

OUT_DIR <- "outputs/ppc"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# sampling for PPC (keeps memory sane)
SET_SEED <- 1
N_SCATTER <- 2000
TOP_COUNTRIES_PER_PARASITE <- 6

EPS <- 1e-6
clip01 <- function(p, eps = EPS) pmin(pmax(p, eps), 1 - eps)
logit  <- function(p) log(p/(1-p))

# ---------------------------
# Load
# ---------------------------
stopifnot(file.exists(FIT_RDS), file.exists(STAN_DATA_RDS), file.exists(LEVELS_RDS), file.exists(MODEL_PANEL))

fit <- readRDS(FIT_RDS)
stan_data <- readRDS(STAN_DATA_RDS)
lev <- readRDS(LEVELS_RDS)

grid <- readr::read_csv(MODEL_PANEL, show_col_types = FALSE) %>%
  mutate(
    iso3 = toupper(iso3),
    year = as.integer(year),
    parasite = as.character(parasite)
  )

# ---------------------------
# Diagnostics snapshot
# ---------------------------
sum_all <- fit$summary()

n_bad_rhat <- sum(is.finite(sum_all$rhat) & sum_all$rhat > 1.01)
n_bad_ess  <- sum(is.finite(sum_all$ess_bulk) & sum_all$ess_bulk < 200)
n_div <- NA_integer_

# cmdstanr stores sampler diagnostics; simplest: parse from diagnose file if present
diag_path <- "data/processed/fit_sac_cmdstan_diagnose.txt"
if (file.exists(diag_path)) {
  diag_txt <- paste(readLines(diag_path, warn = FALSE), collapse = "\n")
} else {
  diag_txt <- fit$cmdstan_diagnose()
}

writeLines(
  c(
    "PPC diagnostics snapshot",
    paste0("Parameters with R-hat > 1.01: ", n_bad_rhat),
    paste0("Parameters with bulk ESS < 200: ", n_bad_ess),
    "",
    "CmdStan diagnose output:",
    diag_txt
  ),
  file.path(OUT_DIR, "ppc_diagnostics.txt")
)

# ---------------------------
# Map obs index n -> iso3/year/parasite (matches p_hat_obs[n])
# ---------------------------
stopifnot(stan_data$N_obs == length(stan_data$s_obs))

obs_map <- tibble(
  n = seq_len(stan_data$N_obs),
  series_id = stan_data$s_obs,
  t = stan_data$t_obs,
  parasite_id = lev$series_map$parasite_id[series_id]
) %>%
  mutate(
    year = lev$years[t],
    parasite = lev$parasites[parasite_id],
    iso3 = lev$series_map$iso3[series_id]
  ) %>%
  left_join(
    grid %>% select(iso3, parasite, year, prev_obs),
    by = c("iso3","parasite","year")
  ) %>%
  mutate(prev_obs = clip01(prev_obs, EPS))

# ---------------------------
# Choose PPC subsets
# ---------------------------
set.seed(SET_SEED)

# Scatter subset (random obs)
n_scatter <- min(N_SCATTER, nrow(obs_map))
obs_scatter <- obs_map %>% slice_sample(n = n_scatter)

# Time-series subset (top mean prevalence countries per parasite)
top_c <- obs_map %>%
  group_by(parasite, iso3) %>%
  summarise(mean_prev = mean(prev_obs, na.rm = TRUE), .groups = "drop") %>%
  group_by(parasite) %>%
  slice_max(mean_prev, n = TOP_COUNTRIES_PER_PARASITE, with_ties = FALSE) %>%
  ungroup()

obs_ts <- obs_map %>%
  semi_join(top_c, by = c("parasite","iso3")) %>%
  arrange(parasite, iso3, year)

idx <- sort(unique(c(obs_scatter$n, obs_ts$n)))
vars <- paste0("p_hat_obs[", idx, "]")

# ---------------------------
# Extract posterior summaries for selected indices
# ---------------------------
draws_sel <- fit$draws(variables = vars)

summ_sel <- posterior::summarise_draws(
  draws_sel,
  quantiles = c(0.025, 0.5, 0.975)
) %>%
  transmute(
    variable,
    n = as.integer(str_match(variable, "\\[(\\d+)\\]")[,2]),
    p_lo = q2.5,
    p_med = q50,
    p_hi = q97.5
  )

ppc_df <- obs_map %>%
  filter(n %in% idx) %>%
  left_join(summ_sel, by = "n") %>%
  mutate(
    p_lo = clip01(p_lo, EPS),
    p_med = clip01(p_med, EPS),
    p_hi = clip01(p_hi, EPS),
    resid_logit = logit(prev_obs) - logit(p_med)
  )

readr::write_csv(ppc_df, file.path(OUT_DIR, "ppc_pred_summary_selected.csv"))

# ---------------------------
# Plot 1: calibration scatter (obs vs pred median)
# ---------------------------
df_sc <- ppc_df %>% filter(n %in% obs_scatter$n)

p1 <- ggplot(df_sc, aes(x = prev_obs, y = p_med)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "PPC: Observed vs Predicted Prevalence (median)",
    x = "Observed prevalence",
    y = "Predicted prevalence (posterior median)"
  ) +
  theme_bw()

ggsave(filename = file.path(OUT_DIR, "ppc_scatter_obs_vs_pred.png"), plot = p1, width = 7, height = 6, dpi = 300)

# ---------------------------
# Plot 2: time-series overlay for example countries
# ---------------------------
df_ts <- ppc_df %>%
  filter(n %in% obs_ts$n) %>%
  arrange(parasite, iso3, year)

# keep the plot size reasonable: facet by parasite, and show lines per iso3
p2 <- ggplot(df_ts, aes(x = year)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi, group = iso3), alpha = 0.20) +
  geom_line(aes(y = p_med, group = iso3), linewidth = 0.6) +
  geom_point(aes(y = prev_obs, group = iso3), size = 0.8, alpha = 0.7) +
  facet_wrap(~ parasite, scales = "free_y") +
  labs(
    title = "PPC: Example Country Time Series (observed points over posterior median + 95% CrI)",
    y = "Prevalence",
    x = "Year"
  ) +
  theme_bw()

ggsave(filename = file.path(OUT_DIR, "ppc_timeseries_examples.png"), plot = p2, width = 10, height = 6, dpi = 300)

# ---------------------------
# Plot 3: residuals by observed prevalence bin
# ---------------------------
breaks <- c(0, 0.01, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0)
labels <- c("[0,0.01]", "(0.01,0.05]", "(0.05,0.10]", "(0.10,0.20]", "(0.20,0.40]",
            "(0.40,0.60]", "(0.60,0.80]", "(0.80,1]")

df_bins <- df_sc %>%
  mutate(bin = cut(prev_obs, breaks = breaks, include.lowest = TRUE, labels = labels)) %>%
  group_by(parasite, bin) %>%
  summarise(
    n = n(),
    resid_mean = mean(resid_logit, na.rm = TRUE),
    resid_p10  = quantile(resid_logit, 0.10, na.rm = TRUE),
    resid_p90  = quantile(resid_logit, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

p3 <- ggplot(df_bins, aes(x = bin, y = resid_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_errorbar(aes(ymin = resid_p10, ymax = resid_p90), width = 0.2) +
  facet_wrap(~ parasite, scales = "free_x") +
  labs(
    title = "PPC: Residuals on logit scale by observed prevalence bin",
    x = "Observed prevalence bin",
    y = "logit(obs) - logit(pred median)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(OUT_DIR, "ppc_residuals_by_bin.png"), plot = p3, width = 11, height = 6, dpi = 300)

message("Wrote PPC outputs to: ", OUT_DIR)
message("DONE.")