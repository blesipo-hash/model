# =============================================================================
# 06b_food_savings_mechanistic.R
# Mechanistic food-savings module for SAC (2010-2023)
# Key design features:
#   - Uses repository-native scenario prevalence (realized/NO_PC/scaleups)
#   - Restricts to WHO-PC-eligible SAC country-years
#   - Isolates MDA effect via NO_PC differencing
#   - Uses ISO3 keys and explicit input validation
#   - Includes placeholders for new external data files
#
# Output:
#   outputs/cba_food/* (or OUT_DIR env var)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(countrycode)
  library(readxl)
})

# ---------------------------
# Inputs from existing repo
# ---------------------------
CF_IN <- "outputs/counterfactuals/cf_prevalence_summary.csv"
ELIGIBLE_CY_IN <- "data/processed/who_pc/country_year_sac_pc_required_2010_2023.csv"
WHO_PC_MODEL_IN <- "data/processed/who_pc/pc_sac_for_model.csv"
PANEL_IN <- "data/processed/model_input/sac_panel_model.csv"


# ---------------------------
# Legacy code mapping + suggested data sources
# ---------------------------
# This script intentionally uses cleaner ISO3/year inputs. If you are migrating from
# the former monolithic script, use this mapping:
#
#   Former object/file                 -> This script input
#   ---------------------------------------------------------------------------
#   FoodPriceLMIFull2.csv             -> COST_PER_KCAL_IN
#     - Former columns used: food_comp, PriceCal, country_title
#     - Former derivation: CostCal = sum(food_comp*PriceCal) / sum(food_comp)
#     - Here: precompute and provide iso3, year, usd_per_kcal
#
#   Parameters2.xlsx: "Parameters2" + "drug_efficacy"
#                                    -> MORBIDITY_PARAM_IN
#     - Former fields included weight diff / blood loss + efficacy by parasite
#     - Here: provide parasite-age burden threshold and below/above morbidity losses
#
#   annual_pop_by_age.csv / age splits -> AGE_POP_IN (optional)
#     - If absent, this script uses a single SAC age bucket with pop_req_pc
#
#   Chan et al. parasite k constants   -> K_PARAM_IN (optional; otherwise defaults)
#
# Suggested data sources if you need to reconstruct these tables:
#   - FAOSTAT Food Balance Sheets (kcal supply composition):
#       https://www.fao.org/faostat/en/#data/FBS
#   - International Comparison Program / World Bank PPP and CPI deflators
#     (for normalizing food price series to constant USD):
#       https://www.worldbank.org/en/programs/icp
#       https://data.worldbank.org/indicator/FP.CPI.TOTL
#   - WFP VAM / national market price series (commodity-level prices):
#       https://dataviz.vam.wfp.org/economic_explorer/prices
#   - WHO deworming efficacy and technical guidance:
#       https://www.who.int/teams/control-of-neglected-tropical-diseases
#   - Age-structured populations (UN WPP):
#       https://population.un.org/wpp/
#
# Practical note:
#   The most reproducible path is to build a frozen parameter folder under
#   data/parameters/food/ with explicit versioning metadata in filenames.

# ---------------------------
# Placeholder inputs you provide
# ---------------------------
# 1) Country-year USD per kcal (required)
#    required columns: iso3, year, usd_per_kcal
COST_PER_KCAL_IN <- Sys.getenv(
  "COST_PER_KCAL_IN",
  unset = "data/parameters/food/cost_per_kcal_iso3_year.csv"  # legacy-derived from FoodPriceLMIFull2.csv
)

# 2) Morbidity/calorie parameter table (required)
#    required columns:
#      parasite, age_class, burden_threshold,
#      blood_loss_l_below, tissue_loss_kg_below,
#      blood_loss_l_above, tissue_loss_kg_above,
#      drug_efficacy
MORBIDITY_PARAM_IN <- Sys.getenv(
  "MORBIDITY_PARAM_IN",
  unset = "data/parameters/food/morbidity_kcal_params_by_parasite_age.csv"  # legacy-derived from Parameters2.xlsx sheets
)

# 3) Optional SAC age population split
#    If present, required columns: iso3, year, age_class, pop_age
#    If missing, script uses single age_class = "SAC" with pop_req_pc.
AGE_POP_IN <- Sys.getenv(
  "AGE_POP_IN",
  unset = "data/processed/demography/sac_pop_by_iso3_year_age.csv"
)

# 4) Optional parasite-specific aggregation k priors
#    If present, required columns: parasite, k_mean
#    If missing, fallback defaults are used.
K_PARAM_IN <- Sys.getenv(
  "K_PARAM_IN",
  unset = "data/parameters/food/k_aggregation_by_parasite.csv"
)

# Legacy fallback files (optional)
LEGACY_COST_CAL_IN <- Sys.getenv("LEGACY_COST_CAL_IN", unset = "FoodPriceLMIFull2.csv")
LEGACY_PARAMETERS_XLSX <- Sys.getenv("LEGACY_PARAMETERS_XLSX", unset = "Parameters2.xlsx")

# ---------------------------
# Config
# ---------------------------
YEAR_MIN <- 2010L
YEAR_MAX <- 2023L
OUT_DIR <- Sys.getenv("OUT_DIR", unset = "outputs/cba_food")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Calorie conversion constants (override if needed)
KCAL_PER_L_BLOOD <- as.numeric(Sys.getenv("KCAL_PER_L_BLOOD", "1300"))
KCAL_PER_KG_TISSUE <- as.numeric(Sys.getenv("KCAL_PER_KG_TISSUE", "5000"))

# Fallback parasite-specific k values if K_PARAM_IN is absent
# (Can be replaced by your preferred table.)
DEFAULT_K <- tibble::tribble(
  ~parasite,              ~k_mean,
  "Ascariasis",           0.54,
  "Trichuriasis",         0.23,
  "Hookworm disease",     0.34,
  "Schistosomiasis",      0.23
)

clip01 <- function(x) pmin(pmax(as.numeric(x), 0), 1)
need_cols <- function(df, cols, nm) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(nm, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

# Convert prevalence (0..1) to mean burden M using NB relation:
# p = 1 - (k/(k+M))^k  =>  M = k * ((1-p)^(-1/k) - 1)
mean_burden_from_prev <- function(prev, k) {
  prev <- clip01(prev)
  k <- as.numeric(k)
  k * (((1 - prev)^(-1 / k)) - 1)
}

# Probability W > threshold under NB(size=k, mu=M)
# Returns conditional probability among infected persons: P(W>thr | W>0)
p_above_threshold_infected <- function(mu, k, threshold) {
  mu <- pmax(as.numeric(mu), 0)
  k <- as.numeric(k)
  threshold <- pmax(as.integer(round(threshold)), 0L)
  
  p0 <- dnbinom(0, size = k, mu = mu)
  p_pos <- pmax(1 - p0, 1e-12)
  p_above_all <- 1 - pnbinom(threshold, size = k, mu = mu)
  out <- p_above_all / p_pos
  pmin(pmax(out, 0), 1)
}

# ---------------------------
# Validate core files
# ---------------------------
stopifnot(file.exists(CF_IN), file.exists(ELIGIBLE_CY_IN), file.exists(WHO_PC_MODEL_IN))
if (!file.exists(PANEL_IN)) {
  warning("PANEL_IN not found; continuing because pop_req_pc comes from WHO_PC_MODEL_IN.")
}

build_cost_kcal_from_legacy <- function(path) {
  legacy <- read_csv(path, show_col_types = FALSE)
  need_cols(legacy, c("country_title", "food_comp", "PriceCal"), "LEGACY_COST_CAL_IN")
  legacy %>%
    mutate(
      iso3 = countrycode(country_title, "country.name", "iso3c"),
      prop = as.numeric(food_comp) * as.numeric(PriceCal)
    ) %>%
    filter(!is.na(iso3)) %>%
    group_by(iso3) %>%
    summarise(
      usd_per_kcal = sum(prop, na.rm = TRUE) / sum(as.numeric(food_comp), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::crossing(year = YEAR_MIN:YEAR_MAX) %>%
    select(iso3, year, usd_per_kcal)
}

build_morbidity_from_legacy <- function(path) {
  p2 <- read_excel(path, sheet = "Parameters2")
  eff <- read_excel(path, sheet = "drug_efficacy")
  need_cols(p2, c("Parasite", "WeightDiff", "BloodLoss"), "Parameters2.xlsx::Parameters2")
  need_cols(eff, c("Parasite", "Drug_Efficacy"), "Parameters2.xlsx::drug_efficacy")
  
  # conservative single SAC class defaults with parasite-specific thresholds
  thr <- tibble::tribble(
    ~parasite,            ~burden_threshold,
    "Ascariasis",         15,
    "Trichuriasis",       130,
    "Hookworm disease",   30,
    "Schistosomiasis",    171
  )
  
  p2 %>%
    transmute(
      parasite = as.character(Parasite),
      blood_loss_l_above = pmax(as.numeric(BloodLoss), 0),
      tissue_loss_kg_above = pmax(as.numeric(WeightDiff), 0)
    ) %>%
    group_by(parasite) %>%
    summarise(
      blood_loss_l_above = mean(blood_loss_l_above, na.rm = TRUE),
      tissue_loss_kg_above = mean(tissue_loss_kg_above, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      eff %>% transmute(parasite = as.character(Parasite), drug_efficacy = clip01(as.numeric(Drug_Efficacy))),
      by = "parasite"
    ) %>%
    left_join(thr, by = "parasite") %>%
    mutate(
      age_class = "SAC",
      blood_loss_l_below = 0,
      tissue_loss_kg_below = 0,
      burden_threshold = coalesce(burden_threshold, 50),
      drug_efficacy = coalesce(drug_efficacy, 0.75)
    ) %>%
    select(parasite, age_class, burden_threshold,
           blood_loss_l_below, tissue_loss_kg_below,
           blood_loss_l_above, tissue_loss_kg_above,
           drug_efficacy)
}

# ---------------------------
# Load modeled prevalence scenarios
# ---------------------------
cf <- read_csv(CF_IN, show_col_types = FALSE) %>%
  transmute(
    iso3 = toupper(iso3),
    year = as.integer(year),
    parasite = as.character(parasite),
    scenario_id = as.character(scenario),
    prev = clip01(prev_q50)
  ) %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX)

eligible_cy <- read_csv(ELIGIBLE_CY_IN, show_col_types = FALSE) %>%
  transmute(iso3 = toupper(iso3), year = as.integer(year)) %>%
  distinct() %>%
  filter(!is.na(iso3), nchar(iso3) == 3, !is.na(year), year >= YEAR_MIN, year <= YEAR_MAX)

if (nrow(eligible_cy) == 0) {
  stop("WHO eligible SAC country-year list is empty.", call. = FALSE)
}

cf <- cf %>% semi_join(eligible_cy, by = c("iso3", "year"))

scen <- sort(unique(cf$scenario_id))
if (!("NO_PC" %in% scen)) {
  stop("Counterfactual file must include scenario NO_PC.", call. = FALSE)
}

# ---------------------------
# Load population requiring PC (SAC)
# ---------------------------
pc_req <- read_csv(WHO_PC_MODEL_IN, show_col_types = FALSE) %>%
  transmute(
    iso3 = toupper(iso3),
    year = as.integer(year),
    parasite = case_when(
      parasite_pc == "SCH" ~ "Schistosomiasis",
      parasite_pc == "STH" ~ "STH",
      TRUE ~ as.character(parasite_pc)
    ),
    pop_req_pc = as.numeric(if ("pop_req_pc" %in% names(.)) pop_req_pc else requiring_pc)
  ) %>%
  # Split STH umbrella into species-level to match model parasite strings
  tidyr::uncount(if_else(parasite == "STH", 3L, 1L), .remove = FALSE, .id = "sth_id") %>%
  mutate(
    parasite = case_when(
      parasite == "STH" & sth_id == 1L ~ "Hookworm disease",
      parasite == "STH" & sth_id == 2L ~ "Ascariasis",
      parasite == "STH" & sth_id == 3L ~ "Trichuriasis",
      TRUE ~ parasite
    )
  ) %>%
  select(-sth_id) %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
  semi_join(eligible_cy, by = c("iso3", "year")) %>%
  distinct()

# ---------------------------
# Optional age split (if provided)
# ---------------------------
if (file.exists(AGE_POP_IN)) {
  age_pop <- read_csv(AGE_POP_IN, show_col_types = FALSE) %>%
    transmute(
      iso3 = toupper(iso3),
      year = as.integer(year),
      age_class = as.character(age_class),
      pop_age = as.numeric(pop_age)
    )
  need_cols(age_pop, c("iso3", "year", "age_class", "pop_age"), "AGE_POP_IN")
  
  age_weights <- age_pop %>%
    group_by(iso3, year) %>%
    mutate(pop_total_age = sum(pop_age, na.rm = TRUE),
           age_weight = if_else(pop_total_age > 0, pop_age / pop_total_age, NA_real_)) %>%
    ungroup() %>%
    select(iso3, year, age_class, age_weight)
} else {
  message("AGE_POP_IN not found. Using single age_class = SAC with weight 1.")
  age_weights <- pc_req %>%
    distinct(iso3, year) %>%
    mutate(age_class = "SAC", age_weight = 1)
}

# ---------------------------
# Load food cost/kcal and morbidity params
# ---------------------------
cost_kcal <- if (file.exists(COST_PER_KCAL_IN)) {
  read_csv(COST_PER_KCAL_IN, show_col_types = FALSE) %>%
    transmute(
      iso3 = toupper(iso3),
      year = as.integer(year),
      usd_per_kcal = as.numeric(usd_per_kcal)
    )
} else if (file.exists(LEGACY_COST_CAL_IN)) {
  message("COST_PER_KCAL_IN not found; deriving usd_per_kcal from legacy file: ", LEGACY_COST_CAL_IN)
  build_cost_kcal_from_legacy(LEGACY_COST_CAL_IN)
} else {
  stop(
    "Missing COST_PER_KCAL_IN and no legacy fallback found. Provide either:\n",
    "  - ", COST_PER_KCAL_IN, " (iso3, year, usd_per_kcal), or\n",
    "  - ", LEGACY_COST_CAL_IN, " (country_title, food_comp, PriceCal).",
    call. = FALSE
  )
}
need_cols(cost_kcal, c("iso3", "year", "usd_per_kcal"), "cost_kcal")

# If you only have the legacy country-level file (e.g., FoodPriceLMIFull2.csv),
# create COST_PER_KCAL_IN by:
#   1) mapping country_title -> iso3,
#   2) assigning/expanding to years (or using year-specific series),
#   3) computing usd_per_kcal = sum(food_comp*PriceCal)/sum(food_comp).

morb <- if (file.exists(MORBIDITY_PARAM_IN)) {
  read_csv(MORBIDITY_PARAM_IN, show_col_types = FALSE) %>%
    transmute(
      parasite = as.character(parasite),
      age_class = as.character(age_class),
      burden_threshold = as.numeric(burden_threshold),
      blood_loss_l_below = as.numeric(blood_loss_l_below),
      tissue_loss_kg_below = as.numeric(tissue_loss_kg_below),
      blood_loss_l_above = as.numeric(blood_loss_l_above),
      tissue_loss_kg_above = as.numeric(tissue_loss_kg_above),
      drug_efficacy = as.numeric(drug_efficacy)
    ) %>%
    mutate(drug_efficacy = clip01(drug_efficacy))
} else if (file.exists(LEGACY_PARAMETERS_XLSX)) {
  message("MORBIDITY_PARAM_IN not found; deriving morbidity parameters from legacy file: ", LEGACY_PARAMETERS_XLSX)
  build_morbidity_from_legacy(LEGACY_PARAMETERS_XLSX)
} else {
  stop(
    "Missing MORBIDITY_PARAM_IN and no legacy fallback found. Provide either:\n",
    "  - ", MORBIDITY_PARAM_IN, " (parasite-age morbidity params), or\n",
    "  - ", LEGACY_PARAMETERS_XLSX, " with sheets Parameters2 and drug_efficacy.",
    call. = FALSE
  )
}
need_cols(
  morb,
  c("parasite", "age_class", "burden_threshold", "blood_loss_l_below", "tissue_loss_kg_below",
    "blood_loss_l_above", "tissue_loss_kg_above", "drug_efficacy"),
  "MORBIDITY_PARAM_IN"
)

# Optional k input with fallback defaults
if (file.exists(K_PARAM_IN)) {
  k_tbl <- read_csv(K_PARAM_IN, show_col_types = FALSE) %>%
    transmute(parasite = as.character(parasite), k_mean = as.numeric(k_mean))
  need_cols(k_tbl, c("parasite", "k_mean"), "K_PARAM_IN")
} else {
  message("K_PARAM_IN not found. Using built-in fallback k values.")
  k_tbl <- DEFAULT_K
}

# ---------------------------
# Build analytic base table
# ---------------------------
base <- cf %>%
  inner_join(pc_req, by = c("iso3", "year", "parasite")) %>%
  filter(!is.na(pop_req_pc), pop_req_pc > 0) %>%
  inner_join(age_weights, by = c("iso3", "year")) %>%
  mutate(pop_age = pop_req_pc * age_weight) %>%
  inner_join(k_tbl, by = "parasite") %>%
  inner_join(morb, by = c("parasite", "age_class")) %>%
  inner_join(cost_kcal, by = c("iso3", "year")) %>%
  mutate(
    prev = clip01(prev),
    mean_burden = mean_burden_from_prev(prev, k_mean),
    infected_n = pop_age * prev,
    p_above_infected = p_above_threshold_infected(mean_burden, k_mean, burden_threshold),
    infected_above_n = infected_n * p_above_infected,
    infected_below_n = pmax(infected_n - infected_above_n, 0),
    kcal_per_person_below = blood_loss_l_below * KCAL_PER_L_BLOOD + tissue_loss_kg_below * KCAL_PER_KG_TISSUE,
    kcal_per_person_above = blood_loss_l_above * KCAL_PER_L_BLOOD + tissue_loss_kg_above * KCAL_PER_KG_TISSUE,
    # Gross calorie loss under scenario (before efficacy attribution)
    kcal_loss_scn = infected_below_n * kcal_per_person_below + infected_above_n * kcal_per_person_above,
    # Effectively avertable kcal under treatment effect assumption
    kcal_loss_scn_effective = kcal_loss_scn * drug_efficacy,
    food_loss_usd_scn = kcal_loss_scn_effective * usd_per_kcal
  ) %>%
  mutate(
    country = countrycode(iso3, "iso3c", "country.name"),
    region = countrycode(iso3, "iso3c", "region")
  )

if (nrow(base) == 0) {
  stop(
    "No rows after joins. Check key consistency across: cf prevalence, WHO PC, age weights, k params, morbidity params, cost/kcal.",
    call. = FALSE
  )
}

# ---------------------------
# Isolate attributable food savings vs NO_PC
# ---------------------------
anchor_nopc <- base %>%
  filter(scenario_id == "NO_PC") %>%
  select(iso3, year, parasite, age_class, food_loss_usd_nopc = food_loss_usd_scn, kcal_loss_nopc = kcal_loss_scn_effective)

food_by <- base %>%
  left_join(anchor_nopc, by = c("iso3", "year", "parasite", "age_class")) %>%
  mutate(
    food_savings_usd = food_loss_usd_nopc - food_loss_usd_scn,
    kcal_savings = kcal_loss_nopc - kcal_loss_scn_effective,
    food_savings_usd = if_else(scenario_id == "NO_PC", 0, food_savings_usd),
    kcal_savings = if_else(scenario_id == "NO_PC", 0, kcal_savings)
  )

# QC checks
if (any(food_by$year < YEAR_MIN | food_by$year > YEAR_MAX, na.rm = TRUE)) {
  stop("Years outside 2010-2023 in output.", call. = FALSE)
}
if (any(is.na(food_by$usd_per_kcal) | food_by$usd_per_kcal <= 0, na.rm = TRUE)) {
  stop("Non-positive or missing usd_per_kcal in analytic rows.", call. = FALSE)
}
if (any(is.na(food_by$k_mean) | food_by$k_mean <= 0, na.rm = TRUE)) {
  stop("Non-positive or missing k_mean in analytic rows.", call. = FALSE)
}
if (any(is.na(food_by$drug_efficacy), na.rm = TRUE)) {
  stop("Missing drug_efficacy in analytic rows.", call. = FALSE)
}

dup <- food_by %>% count(iso3, year, parasite, age_class, scenario_id) %>% filter(n > 1)
if (nrow(dup) > 0) {
  stop("Duplicate keys in food output: iso3-year-parasite-age_class-scenario_id.", call. = FALSE)
}

# ---------------------------
# Write outputs
# ---------------------------
food_country_year_parasite <- food_by %>%
  group_by(iso3, country, region, year, parasite, scenario_id) %>%
  summarise(
    pop_req_pc = sum(pop_req_pc, na.rm = TRUE),
    infected_n = sum(infected_n, na.rm = TRUE),
    infected_above_n = sum(infected_above_n, na.rm = TRUE),
    infected_below_n = sum(infected_below_n, na.rm = TRUE),
    kcal_loss_scn = sum(kcal_loss_scn_effective, na.rm = TRUE),
    food_loss_usd_scn = sum(food_loss_usd_scn, na.rm = TRUE),
    kcal_savings = sum(kcal_savings, na.rm = TRUE),
    food_savings_usd = sum(food_savings_usd, na.rm = TRUE),
    .groups = "drop"
  )

food_region_year <- food_country_year_parasite %>%
  group_by(region, year, scenario_id) %>%
  summarise(
    food_savings_usd = sum(food_savings_usd, na.rm = TRUE),
    kcal_savings = sum(kcal_savings, na.rm = TRUE),
    .groups = "drop"
  )

food_global_year <- food_country_year_parasite %>%
  group_by(year, scenario_id) %>%
  summarise(
    food_savings_usd = sum(food_savings_usd, na.rm = TRUE),
    kcal_savings = sum(kcal_savings, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(region = "Global")

write_csv(food_by, file.path(OUT_DIR, "food_chain_by_country_year_parasite_age_scenario.csv"))
write_csv(food_country_year_parasite, file.path(OUT_DIR, "food_savings_by_country_year_parasite_scenario.csv"))
write_csv(food_region_year, file.path(OUT_DIR, "food_savings_by_region_year_scenario.csv"))
write_csv(food_global_year, file.path(OUT_DIR, "food_savings_global_year_scenario.csv"))

message("Done. Outputs written to: ", OUT_DIR)