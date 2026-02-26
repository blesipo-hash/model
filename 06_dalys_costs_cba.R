# =============================================================================
# 06_dalys_costs_cba.R
# Q1-only DALYs + costs for WHO-eligible SAC country-years (2010-2023)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(countrycode)
})

CF_IN <- "outputs/counterfactuals/cf_prevalence_summary.csv"
PANEL_IN <- "data/processed/model_input/sac_panel_model.csv"
ELIGIBLE_CY_IN <- "data/processed/who_pc/country_year_sac_pc_required_2010_2023.csv"
WHO_PC_MODEL_IN <- "data/processed/who_pc/pc_sac_for_model.csv"
DALY_IN <- "data/processed/gbd/daly_sac_by_iso3_year_parasite.csv"
UNIT_COST_IN <- "data/parameters/unit_cost_by_iso3_parasite.csv"
WDI_GDP_IN <- "data/raw/econ/wdi_gdp_pcap_ppp_const2021.csv"

DALY_ONLY <- as.integer(Sys.getenv("DALY_ONLY", "0")) == 1
# Preferred output env var for Q1 runner; keep CBA_OUT_DIR fallback for compatibility.
OUT_DIR <- Sys.getenv("OUT_DIR", unset = Sys.getenv("CBA_OUT_DIR", "outputs/cba"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

YEAR_MIN <- 2010L
YEAR_MAX <- 2023L
# Monetization assumption (USD per DALY averted); override via env if needed.
DALY_VALUE <- as.numeric(Sys.getenv("DALY_VALUE_USD", "1000"))
DEFAULT_UNIT_COST <- 0.30

if (!is.finite(DALY_VALUE) || DALY_VALUE <= 0) {
  stop("DALY_VALUE_USD must be a positive finite number.", call. = FALSE)
}

if (DALY_ONLY) {
  message("Running 06_dalys_costs_cba.R in DALY-only mode.")
}

clip01 <- function(x) pmin(pmax(as.numeric(x), 0), 1)
need_cols <- function(df, cols, nm) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(nm, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

stopifnot(file.exists(CF_IN), file.exists(PANEL_IN), file.exists(ELIGIBLE_CY_IN), file.exists(DALY_IN))

cf <- read_csv(CF_IN, show_col_types = FALSE) %>%
  mutate(
    iso3 = toupper(iso3),
    year = as.integer(year),
    parasite = as.character(parasite),
    scenario = as.character(scenario),
    prev = clip01(prev_q50)
  ) %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX)

eligible_cy <- read_csv(ELIGIBLE_CY_IN, show_col_types = FALSE) %>%
  transmute(iso3 = toupper(iso3), year = as.integer(year)) %>%
  distinct() %>%
  filter(!is.na(iso3), nchar(iso3) == 3, !is.na(year), year >= YEAR_MIN, year <= YEAR_MAX)

if (nrow(eligible_cy) == 0) {
  stop("WHO eligible country-year list is empty: ", ELIGIBLE_CY_IN, call. = FALSE)
}

cf <- cf %>% semi_join(eligible_cy, by = c("iso3", "year"))

panel <- read_csv(PANEL_IN, show_col_types = FALSE) %>%
  mutate(iso3 = toupper(iso3), year = as.integer(year), parasite = as.character(parasite)) %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
  semi_join(eligible_cy, by = c("iso3", "year"))

daly <- read_csv(DALY_IN, show_col_types = FALSE) %>%
  mutate(iso3 = toupper(iso3), year = as.integer(year), parasite = as.character(parasite))
need_cols(daly, c("iso3", "year", "parasite", "daly_gbd"), "DALY_IN")

pc_req <- if (file.exists(WHO_PC_MODEL_IN)) {
  read_csv(WHO_PC_MODEL_IN, show_col_types = FALSE) %>%
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
    distinct()
} else {
  stop("Missing WHO PC model file: ", WHO_PC_MODEL_IN, call. = FALSE)
}

if (file.exists(UNIT_COST_IN)) {
  unit_cost <- read_csv(UNIT_COST_IN, show_col_types = FALSE) %>%
    mutate(iso3 = toupper(iso3), parasite = as.character(parasite), unit_cost = as.numeric(unit_cost)) %>%
    select(iso3, parasite, unit_cost)
} else {
  unit_cost <- tibble(iso3 = character(), parasite = character(), unit_cost = numeric())
}

costing <- panel %>%
  distinct(iso3, year, parasite, pop_sac, pc_cov_filled) %>%
  left_join(pc_req, by = c("iso3", "year", "parasite")) %>%
  mutate(
    pop_req_pc = as.numeric(pop_req_pc),
    cov_realized = clip01(pc_cov_filled),
    cov_scaleup_set_75 = 0.75,
    cov_scaleup_floor_75 = pmax(cov_realized, 0.75),
    cov_scaleup_100 = 1,
    cov_NO_PC = 0
  ) %>%
  pivot_longer(starts_with("cov_"), names_to = "scenario", values_to = "coverage") %>%
  mutate(scenario = str_remove(scenario, "^cov_"), coverage = clip01(coverage)) %>%
  left_join(unit_cost, by = c("iso3", "parasite")) %>%
  mutate(
    unit_cost = coalesce(unit_cost, DEFAULT_UNIT_COST),
    pop_sac = as.numeric(pop_sac),
    treated_n = coverage * pop_req_pc,
    cost = treated_n * unit_cost
  ) %>%
  select(iso3, year, parasite, scenario, pop_sac, pop_req_pc, coverage, treated_n, unit_cost, cost)

prev_anchor <- cf %>%
  filter(scenario %in% c("realized", "NO_PC")) %>%
  select(iso3, year, parasite, scenario, prev) %>%
  pivot_wider(names_from = scenario, values_from = prev, names_prefix = "prev_")

if (!all(c("prev_realized", "prev_NO_PC") %in% names(prev_anchor))) {
  stop("Counterfactual summary must include scenarios: realized and NO_PC.", call. = FALSE)
}

core <- cf %>%
  left_join(prev_anchor, by = c("iso3", "year", "parasite")) %>%
  left_join(daly, by = c("iso3", "year", "parasite")) %>%
  left_join(costing, by = c("iso3", "year", "parasite", "scenario")) %>%
  filter(!is.na(pop_req_pc), pop_req_pc > 0) %>%
  mutate(
    prev_realized = clip01(prev_realized),
    prev_NO_PC = clip01(prev_NO_PC),
    burden_ratio = if_else(prev_realized > 0, prev / prev_realized, NA_real_),
    dalys_pc = daly_gbd * burden_ratio,
    country = countrycode(iso3, "iso3c", "country.name"),
    region = countrycode(iso3, "iso3c", "region")
  ) %>%
  group_by(iso3, year, parasite) %>%
  mutate(
    dalys_nopc = dalys_pc[scenario == "NO_PC"][1],
    cost_nopc = cost[scenario == "NO_PC"][1],
    dalys_averted = dalys_nopc - dalys_pc,
    dalys_averted = if_else(scenario == "NO_PC", 0, dalys_averted),
    cost_usd = cost - cost_nopc,
    icer_cost_per_daly_averted = if_else(dalys_averted > 0, cost_usd / dalys_averted, NA_real_),
    usd_per_daly_assumption = DALY_VALUE,
    dalys_value_usd = dalys_averted * DALY_VALUE,
    nmb_dalys = dalys_value_usd - cost_usd,
    net_benefit_usd = nmb_dalys,
    bcr_dalys = if_else(cost_usd > 0, dalys_value_usd / cost_usd, NA_real_),
    scenario_id = scenario,
    coverage_scn = coverage
  ) %>%
  ungroup()

# QC assertions
if (any(core$year < YEAR_MIN | core$year > YEAR_MAX, na.rm = TRUE)) stop("Years outside 2010-2023 in output.", call. = FALSE)
if (any(is.na(core$pop_req_pc) | core$pop_req_pc <= 0, na.rm = TRUE)) stop("Non-eligible rows present.", call. = FALSE)
if (any(core$coverage_scn < 0 | core$coverage_scn > 1, na.rm = TRUE)) stop("coverage_scn outside [0,1].", call. = FALSE)
if (any(abs(core$treated_n - (core$coverage_scn * core$pop_req_pc)) > 1e-6, na.rm = TRUE)) stop("treated_n != coverage_scn * pop_req_pc", call. = FALSE)

dup <- core %>% count(iso3, year, parasite, scenario_id) %>% filter(n > 1)
if (nrow(dup) > 0) stop("Duplicate iso3-year-parasite-scenario_id keys.", call. = FALSE)

scen <- unique(core$scenario_id)
if (!("NO_PC" %in% scen && "realized" %in% scen && any(grepl("scaleup_", scen)))) {
  stop("Required scenarios missing in output.", call. = FALSE)
}

dalys_by <- core %>%
  select(
    iso3, country, region, year, parasite, scenario_id,
    pop_sac, pop_req_pc,
    coverage_scn, treated_n,
    cost_usd,
    dalys_pc, dalys_nopc, dalys_averted,
    icer_cost_per_daly_averted,
    usd_per_daly_assumption,
    dalys_value_usd, nmb_dalys, net_benefit_usd, bcr_dalys
  )

write_csv(dalys_by, file.path(OUT_DIR, "dalys_by_country_year_parasite_scenario.csv"))

if (DALY_ONLY) {
  if (file.exists(WDI_GDP_IN)) {
    analysis_iso3 <- dalys_by %>% distinct(iso3)

    gdp_thresholds <- read_csv(WDI_GDP_IN, show_col_types = FALSE) %>%
      filter(`Indicator Code` == "NY.GDP.PCAP.PP.KD") %>%
      select(`Country Code`, any_of(as.character(YEAR_MIN:YEAR_MAX))) %>%
      pivot_longer(
        cols = all_of(as.character(YEAR_MIN:YEAR_MAX)),
        names_to = "year",
        values_to = "gdp_pcap_ppp_const2021"
      ) %>%
      transmute(
        iso3 = toupper(`Country Code`),
        year = as.integer(year),
        gdp_pcap_ppp_const2021 = as.numeric(gdp_pcap_ppp_const2021)
      ) %>%
      filter(!is.na(iso3), nchar(iso3) == 3) %>%
      semi_join(analysis_iso3, by = "iso3")

    q1_nmb_gdp <- dalys_by %>%
      select(iso3, year, parasite, scenario_id, dalys_averted, cost_usd) %>%
      left_join(gdp_thresholds, by = c("iso3", "year")) %>%
      tidyr::crossing(lambda_mult = c(0.5, 1.0)) %>%
      mutate(
        lambda = lambda_mult * gdp_pcap_ppp_const2021,
        daly_benefit_usd = lambda * dalys_averted,
        nmb_usd = daly_benefit_usd - cost_usd
      ) %>%
      select(
        iso3, year, parasite, scenario_id,
        dalys_averted, cost_usd,
        gdp_pcap_ppp_const2021,
        lambda_mult, lambda,
        daly_benefit_usd, nmb_usd
      )

    write_csv(q1_nmb_gdp, file.path(OUT_DIR, "q1_nmb_gdp_thresholds.csv"))
  } else {
    warning(
      "Skipping GDP-threshold NMB output because WDI GDP file is missing: ",
      WDI_GDP_IN,
      call. = FALSE
    )
  }
}

global_summary <- dalys_by %>%
  group_by(scenario_id, year, parasite) %>%
  summarise(
    global_cost_usd = sum(cost_usd, na.rm = TRUE),
    global_dalys_pc = sum(dalys_pc, na.rm = TRUE),
    global_dalys_nopc = sum(dalys_nopc, na.rm = TRUE),
    global_dalys_averted = sum(dalys_averted, na.rm = TRUE),
    global_icer = if_else(sum(dalys_averted, na.rm = TRUE) > 0, sum(cost_usd, na.rm = TRUE) / sum(dalys_averted, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )
write_csv(global_summary, file.path(OUT_DIR, "global_summary_dalys.csv"))

impl_gap_country <- dalys_by %>%
  filter(scenario_id %in% c("scaleup_set_75", "scaleup_floor_75", "scaleup_100")) %>%
  left_join(
    dalys_by %>%
      filter(scenario_id == "realized") %>%
      select(iso3, year, parasite, cost_realized = cost_usd, dalys_realized = dalys_averted),
    by = c("iso3", "year", "parasite")
  ) %>%
  mutate(
    inc_cost_usd = cost_usd - cost_realized,
    inc_dalys_averted = dalys_averted - dalys_realized,
    inc_icer = if_else(inc_dalys_averted > 0, inc_cost_usd / inc_dalys_averted, NA_real_)
  ) %>%
  select(iso3, year, parasite, scenario_id, inc_cost_usd, inc_dalys_averted, inc_icer)

impl_gap <- impl_gap_country %>%
  group_by(year, parasite, scenario_id) %>%
  summarise(
    inc_cost_usd = sum(inc_cost_usd, na.rm = TRUE),
    inc_dalys_averted = sum(inc_dalys_averted, na.rm = TRUE),
    inc_icer = if_else(sum(inc_dalys_averted, na.rm = TRUE) > 0, sum(inc_cost_usd, na.rm = TRUE) / sum(inc_dalys_averted, na.rm = TRUE), NA_real_),
    .groups = "drop"
  )
write_csv(impl_gap, file.path(OUT_DIR, "implementation_gap_incremental_vs_realized.csv"))
write_csv(impl_gap_country, file.path(OUT_DIR, "implementation_gap_incremental_vs_realized_country.csv"))

top_countries <- dalys_by %>%
  filter(scenario_id != "NO_PC") %>%
  group_by(scenario_id, iso3, country, region) %>%
  summarise(dalys_averted = sum(dalys_averted, na.rm = TRUE), .groups = "drop") %>%
  group_by(scenario_id) %>%
  arrange(desc(dalys_averted), .by_group = TRUE) %>%
  slice_head(n = 25) %>%
  ungroup()
write_csv(top_countries, file.path(OUT_DIR, "top_countries_by_dalys_averted.csv"))

realized_total <- global_summary %>%
  filter(scenario_id == "realized") %>%
  summarise(v = sum(global_dalys_averted, na.rm = TRUE)) %>%
  pull(v)
if (length(realized_total) == 0 || is.na(realized_total) || realized_total == 0) {
  write_csv(global_summary, file.path(OUT_DIR, "diagnostic_global_summary_dalys_zero.csv"))
  stop("Global DALYs averted for realized are zero; diagnostics written.", call. = FALSE)
}

# Backward-compatibility: keep compact summary and main table for Q1 consumers
write_csv(
  global_summary %>%
    group_by(scenario_id) %>%
    summarise(
      daly_averted = sum(global_dalys_averted, na.rm = TRUE),
      costs = sum(global_cost_usd, na.rm = TRUE),
      icer = if_else(sum(global_dalys_averted, na.rm = TRUE) > 0, sum(global_cost_usd, na.rm = TRUE) / sum(global_dalys_averted, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    rename(scenario = scenario_id),
  file.path(OUT_DIR, "master_summary_by_scenario.csv")
)

# Canonical Q1 main table
write_csv(dalys_by, file.path(OUT_DIR, "q1_main_results.csv"))
# Backward-compatibility alias for prior downstream readers
write_csv(dalys_by, file.path(OUT_DIR, "q1_q2_q4_main_results.csv"))

message("Wrote Q1 DALY outputs to ", OUT_DIR)
