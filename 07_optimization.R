# =============================================================================
# 07_optimization.R
# Discrete constrained optimization of PC policy using CBA outputs
#
# Inputs:
#   - outputs/cba/cba_by_stratum.csv
#   - data/processed/model_input/sac_panel_model.csv
#   - outputs/counterfactuals/cf_prevalence_summary.csv  (for triggers based on prevalence)
#
# Outputs:
#   - outputs/opt/best_policy_by_country.csv
#   - outputs/opt/policy_rankings_sample.csv
#   - outputs/opt/policy_definitions.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# ---------------------------
# Config (EDIT)
# ---------------------------
CBA_IN   <- "outputs/cba/cba_by_stratum.csv"
PANEL_IN <- "data/processed/model_input/sac_panel_model.csv"
CF_IN    <- "outputs/counterfactuals/cf_prevalence_summary.csv"

OUT_DIR <- "outputs/opt"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Choose objective: "DALY" or "NMB"
OBJECTIVE <- "NMB"

# Choose one sensitivity slice to optimize on (you can loop later)
# (These must exist in cba_by_stratum.csv)
ALPHA_USE <- 1.0
D_USE     <- 0.02
V_USE     <- 5000

# Budget constraint (discounted). Set NA to ignore.
# Example: 1e8 means $100M discounted budget over horizon (per country).
BUDGET_TOTAL <- NA_real_

# Policy constraints
COOLDOWN_YEARS <- 3
COVERAGE_CAP   <- 0.90

# Trigger thresholds (based on predicted prevalence under WHO baseline)
THRESH_HIGH <- 0.20
THRESH_LOW  <- 0.10

# Available levers
FREQ_SET <- c(1, 2)          # rounds: annual=1, biannual=2
COV_SET  <- c("WHO", "PC75") # baseline coverage sources we can use without re-simulating
ESCALATE_TO <- "PC75"        # when high prevalence, go to PC75
DEESCALATE_TO <- "WHO"       # when low prevalence, go back to WHO

# ---------------------------
# Helpers
# ---------------------------
discount_factor <- function(year, d, t0 = min(year)) {
  1 / ((1 + d)^(year - t0))
}

log_join_result <- function(before_df, after_df, join_name) {
  n_before <- nrow(before_df)
  n_after <- nrow(after_df)
  dropped <- n_before - n_after
  pct <- if (n_before > 0) round(100 * dropped / n_before, 2) else 0
  message(sprintf("[%s] rows: before=%s after=%s dropped=%s (%.2f%%)", join_name, n_before, n_after, dropped, pct))
}

enforce_cooldown <- function(cov_choice, years, cooldown = 3) {
  # cov_choice: character vector length T with entries like "WHO" / "PC75"
  # enforce cannot change more often than every `cooldown` years
  out <- cov_choice
  last_switch_year <- years[1]
  for (i in 2:length(years)) {
    if (out[i] != out[i-1]) {
      if ((years[i] - last_switch_year) < cooldown) {
        out[i] <- out[i-1] # revert switch
      } else {
        last_switch_year <- years[i]
      }
    }
  }
  out
}

# ---------------------------
# Load
# ---------------------------
stopifnot(file.exists(CBA_IN), file.exists(PANEL_IN), file.exists(CF_IN))

cba <- readr::read_csv(CBA_IN, show_col_types = FALSE)
panel <- readr::read_csv(PANEL_IN, show_col_types = FALSE) %>%
  mutate(iso3 = toupper(iso3), year = as.integer(year), parasite = as.character(parasite))

cf <- readr::read_csv(CF_IN, show_col_types = FALSE) %>%
  mutate(iso3 = toupper(iso3), year = as.integer(year), parasite = as.character(parasite))

# Filter the sensitivity slice we optimize on
cba0 <- cba %>%
  filter(alpha == ALPHA_USE, d == D_USE, V == V_USE) %>%
  select(iso3, year, parasite, scenario, daly_averted, cost_disc, nb_disc)

# Create baseline prevalence signal for trigger rules (use WHO_PC median)
who_prev <- cf %>%
  filter(scenario == "WHO_PC") %>%
  select(iso3, year, parasite, prev_q50)

# ---------------------------
# Define policies
# ---------------------------
# Policy structure: coverage_source_rule + rounds_rule
# coverage_source_rule choices:
#   - "WHO" (use observed WHO coverage)
#   - "PC75" (use 0.75)
#   - "TRIGGER" (switch between WHO and PC75 based on thresholds, with cooldown)
#
# rounds_rule choices:
#   - 1 or 2 fixed (annual/biannual)
#
policies <- tidyr::crossing(
  cov_rule = c("WHO", "PC75", "TRIGGER"),
  rounds_rule = FREQ_SET
) %>%
  mutate(policy_id = sprintf("P%02d", row_number()))

readr::write_csv(policies, file.path(OUT_DIR, "policy_definitions.csv"))

# ---------------------------
# Build country-year-parasite policy paths
# ---------------------------
# Start from panel coverage WHO + implied PC75
cov_base <- panel %>%
  distinct(iso3, year, parasite, pc_cov_filled, rounds) %>%
  mutate(
    cov_WHO  = pmin(pmax(pc_cov_filled, 0), 1),
    cov_PC75 = 0.75
  ) %>%
  select(iso3, year, parasite, cov_WHO, cov_PC75)

years_all <- sort(unique(panel$year))
t0 <- min(years_all)

# For TRIGGER rule: need prevalence series by iso3-parasite-year
trigger_df <- who_prev %>%
  mutate(
    trigger_state = case_when(
      prev_q50 >= THRESH_HIGH ~ "HIGH",
      prev_q50 <= THRESH_LOW  ~ "LOW",
      TRUE ~ "MID"
    )
  )

# Function to create coverage choice path for one (iso3, parasite)
make_cov_choice <- function(df_one) {
  # df_one has year, cov_WHO, cov_PC75, prev_q50, trigger_state
  yrs <- df_one$year
  choice <- rep(DEESCALATE_TO, length(yrs)) # start at WHO by default
  
  # simple hysteresis: escalate on HIGH, de-escalate on LOW, hold otherwise
  for (i in seq_along(yrs)) {
    if (i == 1) {
      choice[i] <- DEESCALATE_TO
    } else {
      choice[i] <- choice[i-1]
      if (df_one$trigger_state[i] == "HIGH") choice[i] <- ESCALATE_TO
      if (df_one$trigger_state[i] == "LOW")  choice[i] <- DEESCALATE_TO
    }
  }
  
  # apply cooldown constraint
  choice <- enforce_cooldown(choice, yrs, cooldown = COOLDOWN_YEARS)
  choice
}

# Build policy-expanded table: iso3-year-parasite-policy -> scenario proxy + rounds
policy_grid <- cov_base %>%
  left_join(trigger_df, by = c("iso3","year","parasite")) %>%
  group_by(iso3, parasite) %>%
  group_modify(~{
    d <- .x %>% arrange(year)
    d$cov_choice_trigger <- make_cov_choice(d)
    d
  }) %>%
  ungroup()

# Expand to policies and map to scenario proxies
policy_expanded <- policy_grid %>%
  tidyr::crossing(policies) %>%
  mutate(
    cov_choice = case_when(
      cov_rule == "WHO" ~ "WHO",
      cov_rule == "PC75" ~ "PC75",
      cov_rule == "TRIGGER" ~ cov_choice_trigger,
      TRUE ~ NA_character_
    ),
    # choose "scenario" to grab DALY_averted from cba0:
    # WHO coverage corresponds to scenario WHO_PC
    # PC75 corresponds to scenario PC75
    scenario_proxy = case_when(
      cov_choice == "WHO"  ~ "WHO_PC",
      cov_choice == "PC75" ~ "PC75",
      TRUE ~ "WHO_PC"
    ),
    coverage = case_when(
      cov_choice == "WHO" ~ cov_WHO,
      cov_choice == "PC75" ~ cov_PC75,
      TRUE ~ cov_WHO
    ),
    coverage = pmin(coverage, COVERAGE_CAP),
    rounds_policy = rounds_rule
  ) %>%
  select(iso3, year, parasite, policy_id, cov_rule, rounds_rule, cov_choice,
         scenario_proxy, coverage, rounds_policy)

# ---------------------------
# Score each policy using CBA slice (DALYs averted + costs/NB)
# ---------------------------
# We take daly_averted from the scenario proxy, but costs need adjustment for rounds changes.
# Approximation:
#   - benefits use scenario proxy (WHO vs PC75) from cf-derived DALYs
#   - costs scale linearly with rounds_policy / (observed rounds in panel)
#
# If you want full fidelity: re-simulate prevalence under changed rounds/cov in Script 5.
base_rounds <- panel %>%
  distinct(iso3, year, parasite, rounds) %>%
  rename(rounds_base = rounds)

policy_before_scoring <- policy_expanded
policy_scored <- policy_expanded %>%
  left_join(base_rounds, by = c("iso3","year","parasite")) %>%
  left_join(
    cba0 %>% select(iso3, year, parasite, scenario, daly_averted, cost_disc, nb_disc),
    by = c("iso3","year","parasite", "scenario_proxy" = "scenario")
  ) %>%
  mutate(
    # scale costs by rounds ratio (discount already applied in cost_disc in CBA)
    rounds_ratio = if_else(!is.na(rounds_base) & rounds_base > 0, rounds_policy / rounds_base, 1),
    cost_disc_adj = cost_disc * rounds_ratio,
    nb_disc_adj   = (nb_disc + cost_disc) - cost_disc_adj,  # keep benefit fixed, replace cost
    # if V=0, NB is just -cost; still fine
    score = case_when(
      OBJECTIVE == "DALY" ~ daly_averted,
      OBJECTIVE == "NMB"  ~ nb_disc_adj,
      TRUE ~ nb_disc_adj
    )
  )

log_join_result(policy_expanded, policy_scored %>% dplyr::filter(!is.na(score)), "policy_to_cba_slice")

# ---------------------------
# Apply budget constraint (optional) and choose best policy per country
# ---------------------------
policy_country <- policy_scored %>%
  group_by(iso3, policy_id, cov_rule, rounds_rule) %>%
  summarise(
    score_total = sum(score, na.rm = TRUE),
    daly_averted_total = sum(daly_averted, na.rm = TRUE),
    cost_disc_total = sum(cost_disc_adj, na.rm = TRUE),
    nb_disc_total = sum(nb_disc_adj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    feasible_budget = if_else(is.na(BUDGET_TOTAL), TRUE, cost_disc_total <= BUDGET_TOTAL)
  )

best_policy <- policy_country %>%
  filter(feasible_budget) %>%
  group_by(iso3) %>%
  slice_max(order_by = score_total, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(score_total))

readr::write_csv(best_policy, file.path(OUT_DIR, "best_policy_by_country.csv"))
message("Wrote: ", file.path(OUT_DIR, "best_policy_by_country.csv"))

# Save a ranked list (sample) for inspection
ranked <- policy_country %>%
  arrange(iso3, desc(score_total)) %>%
  group_by(iso3) %>%
  slice_head(n = 10) %>%
  ungroup()

readr::write_csv(ranked, file.path(OUT_DIR, "policy_rankings_sample.csv"))
message("Wrote: ", file.path(OUT_DIR, "policy_rankings_sample.csv"))

message("DONE.")