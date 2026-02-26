# =============================================================================
# 01b_build_who_pc_inputs.R
# Build WHO PC SAC inputs from WHO PCT exports (SCH + STH), 2010-2023.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(countrycode)
})

# ---------------------------
# Config
# ---------------------------
SCH_IN <- "data/raw/who_pc/who_pc_sch.xlsx"
STH_IN <- "data/raw/who_pc/who_pc_sth.xlsx"

OUT_DIR <- "data/processed/who_pc"
QC_DIR <- "data/processed/qc/01b_build_who_pc_inputs"
OUT_MODEL <- file.path(OUT_DIR, "pc_sac_for_model.csv")
OUT_MERGE <- file.path(OUT_DIR, "pc_sac_for_merge_gbd_causes.csv")

YEAR_MIN <- 2010L
YEAR_MAX <- 2023L

# ---------------------------
# Helpers
# ---------------------------
norm_txt <- function(x) {
  str_squish(str_replace_all(str_replace_all(as.character(x), "–", "-"), "’", "'"))
}

clamp01 <- function(x) if_else(is.na(x), NA_real_, pmin(pmax(as.numeric(x), 0), 1))

parse_pop_req <- function(x) {
  s <- str_to_lower(str_squish(as.character(x)))
  out <- rep(NA_real_, length(s))
  is_no_pc <- str_detect(s, "\\bno\\s*pc\\s*required\\b")
  is_missing_code <- str_detect(s, "\\b(to\\s*be\\s*defined|surveillance|no\\s*data\\s*available)\\b")
  is_na_token <- is.na(s) | s %in% str_to_lower(c("", "NA", "-", "—", "–"))
  out[is_no_pc] <- 0
  out[is_missing_code | is_na_token] <- NA_real_
  idx_num <- !(is_no_pc | is_missing_code | is_na_token) & str_detect(s, "\\d")
  out[idx_num] <- parse_number(s[idx_num], na = c("", "NA", "-", "—", "–"))
  out
}

map_iso3 <- function(country, country_code) {
  cname <- recode(country, "Malasya" = "Malaysia")
  code_clean <- str_to_upper(str_squish(as.character(country_code)))
  valid_iso3 <- unique(countrycode::codelist$iso3c)
  use_code <- !is.na(code_clean) & str_length(code_clean) == 3 & code_clean %in% valid_iso3
  out <- if_else(use_code, code_clean, countrycode(cname, "country.name", "iso3c"))
  str_to_upper(out)
}

unique_coverage_or_na <- function(x) {
  u <- unique(x[!is.na(x)])
  if (length(u) == 1) u[[1]] else NA_real_
}

# ---------------------------
# File checks
# ---------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SCH_IN)) {
  stop(
    "Missing SCH WHO export. Please place WHO SCH export at data/raw/who_pc/who_pc_sch.xlsx.",
    call. = FALSE
  )
}
if (!file.exists(STH_IN)) stop("Missing STH input: ", STH_IN, call. = FALSE)

# =============================================================================
# SCH: WHO export with exact expected columns
# =============================================================================
raw_sch <- read_excel(SCH_IN)
required_sch_cols <- c(
  "region", "country", "year",
  "SAC population requiring PC for SCH annually",
  "National coverage (%)",
  "country_code"
)
miss_sch <- setdiff(required_sch_cols, names(raw_sch))
if (length(miss_sch) > 0) {
  stop("SCH WHO export missing columns: ", paste(miss_sch, collapse = ", "), call. = FALSE)
}

sch_line <- raw_sch %>%
  transmute(
    year = as.integer(year),
    country = norm_txt(country),
    iso3 = map_iso3(country, country_code),
    region = norm_txt(region),
    pop_req_pc = parse_pop_req(`SAC population requiring PC for SCH annually`),
    coverage_line = parse_number(norm_txt(`National coverage (%)`), na = c("", "NA", "-", "—", "–")) / 100
  ) %>%
  filter(year %in% YEAR_MIN:YEAR_MAX, !is.na(iso3), nchar(iso3) == 3)

sch_cov_conflicts <- sch_line %>%
  group_by(iso3, country, year) %>%
  summarise(
    n_rows = n(),
    n_unique_coverage = n_distinct(coverage_line[!is.na(coverage_line)]),
    coverage_values = paste(sort(unique(coverage_line[!is.na(coverage_line)])), collapse = "|"),
    .groups = "drop"
  ) %>%
  filter(n_unique_coverage > 1)
write_csv(sch_cov_conflicts, file.path(QC_DIR, "sch_coverage_conflicts.csv"))

pc_sch <- sch_line %>%
  group_by(iso3, country, year, region) %>%
  summarise(
    pop_req_pc = if (all(is.na(pop_req_pc))) NA_real_ else max(pop_req_pc, na.rm = TRUE),
    coverage_obs = unique_coverage_or_na(coverage_line),
    .groups = "drop"
  ) %>%
  mutate(
    parasite_pc = "SCH",
    pop_treated = NA_real_,
    coverage_obs = clamp01(coverage_obs),
    coverage_obs_source = "sch_national_coverage_percent_who",
    eligible = !is.na(pop_req_pc) & pop_req_pc > 0,
    delivered = !is.na(coverage_obs) & coverage_obs > 0
  ) %>%
  select(iso3, region, country, year, parasite_pc, pop_req_pc, pop_treated,
         coverage_obs, coverage_obs_source, eligible, delivered)

# =============================================================================
# STH: keep existing logic pattern, output to same schema
# =============================================================================
raw_sth <- read_excel(STH_IN)
required_sth_cols <- c(
  "year", "country", "country_code",
  "Population requiring PC for STH, SAC",
  "National coverage, SAC (%)"
)
miss_sth <- setdiff(required_sth_cols, names(raw_sth))
if (length(miss_sth) > 0) {
  stop("WHO STH missing columns: ", paste(miss_sth, collapse = ", "), call. = FALSE)
}

sth_line <- raw_sth %>%
  transmute(
    year = as.integer(year),
    country = norm_txt(country),
    iso3 = map_iso3(country, country_code),
    region = if ("region" %in% names(raw_sth)) norm_txt(raw_sth$region) else countrycode(iso3, "iso3c", "region"),
    pop_req_pc = parse_pop_req(`Population requiring PC for STH, SAC`),
    coverage_line = parse_number(norm_txt(`National coverage, SAC (%)`), na = c("", "NA", "-", "—", "–")) / 100,
    pop_treated_line = if ("Reported number of SAC treated" %in% names(raw_sth)) {
      parse_number(norm_txt(`Reported number of SAC treated`), na = c("", "NA", "-", "—", "–"))
    } else {
      NA_real_
    }
  ) %>%
  filter(year %in% YEAR_MIN:YEAR_MAX, !is.na(iso3), nchar(iso3) == 3)

sth_cov_conflicts <- sth_line %>%
  group_by(iso3, country, year) %>%
  summarise(
    n_rows = n(),
    n_unique_coverage = n_distinct(coverage_line[!is.na(coverage_line)]),
    coverage_values = paste(sort(unique(coverage_line[!is.na(coverage_line)])), collapse = "|"),
    .groups = "drop"
  ) %>%
  filter(n_unique_coverage > 1)
write_csv(sth_cov_conflicts, file.path(QC_DIR, "sth_coverage_conflicts.csv"))

pc_sth <- sth_line %>%
  group_by(iso3, country, year, region) %>%
  summarise(
    pop_req_pc = if (all(is.na(pop_req_pc))) NA_real_ else max(pop_req_pc, na.rm = TRUE),
    pop_treated = if (all(is.na(pop_treated_line))) NA_real_ else sum(pop_treated_line, na.rm = TRUE),
    coverage_obs = unique_coverage_or_na(coverage_line),
    .groups = "drop"
  ) %>%
  mutate(
    parasite_pc = "STH",
    coverage_obs = clamp01(coverage_obs),
    coverage_obs_source = "sth_national_coverage_percent_who",
    eligible = !is.na(pop_req_pc) & pop_req_pc > 0,
    delivered = !is.na(coverage_obs) & coverage_obs > 0
  ) %>%
  select(iso3, region, country, year, parasite_pc, pop_req_pc, pop_treated,
         coverage_obs, coverage_obs_source, eligible, delivered)

# =============================================================================
# Combine + output model table
# =============================================================================
pc_sac_for_model <- bind_rows(pc_sch, pc_sth) %>%
  mutate(
    pop_req_pc = as.numeric(pop_req_pc),
    pop_treated = as.numeric(pop_treated),
    coverage_obs = clamp01(coverage_obs),
    pop_treated = if_else(is.na(pop_treated), 0, pop_treated),
    eligible = !is.na(pop_req_pc) & pop_req_pc > 0,
    delivered = !is.na(coverage_obs) & coverage_obs > 0
  ) %>%
  arrange(iso3, year, parasite_pc)

if (any(!is.na(pc_sac_for_model$pop_req_pc) & pc_sac_for_model$pop_req_pc < 0)) {
  stop("pop_req_pc must be >= 0.", call. = FALSE)
}

dup_model <- pc_sac_for_model %>% count(iso3, year, parasite_pc, name = "n") %>% filter(n > 1)
if (nrow(dup_model) > 0) {
  print(dup_model, n = min(50, nrow(dup_model)))
  stop("pc_sac_for_model has duplicate keys: iso3, year, parasite_pc", call. = FALSE)
}

write_csv(pc_sac_for_model, OUT_MODEL)
message("Wrote: ", OUT_MODEL)

# =============================================================================
# Build merge table keyed to GBD parasites (backward compatible fields)
# =============================================================================
pc_sac_for_merge <- bind_rows(
  pc_sac_for_model %>%
    filter(parasite_pc == "SCH") %>%
    transmute(
      iso3,
      year,
      parasite = "Schistosomiasis",
      coverage_use = coverage_obs,
      coverage_use_source = coverage_obs_source,
      coverage_use_raw = coverage_obs,
      rounds_use = if_else(delivered, 1L, 0L),
      rounds_use_source = "derived_from_annual_coverage",
      requiring_pc = pop_req_pc,
      targeted = NA_real_,
      treated = pop_treated,
      flag_cov_over_100 = FALSE,
      flag_cov_negative = FALSE,
      flag_treated_gt_requiring = FALSE,
      flag_treated_when_no_pc_required = FALSE,
      flag_requiring_conflict = FALSE
    ),
  pc_sac_for_model %>%
    filter(parasite_pc == "STH") %>%
    transmute(
      iso3,
      year,
      coverage_use = coverage_obs,
      coverage_use_source = coverage_obs_source,
      coverage_use_raw = coverage_obs,
      rounds_use = if_else(delivered, 1L, 0L),
      rounds_use_source = "derived_from_annual_coverage",
      requiring_pc = pop_req_pc,
      targeted = NA_real_,
      treated = pop_treated,
      flag_cov_over_100 = FALSE,
      flag_cov_negative = FALSE,
      flag_treated_gt_requiring = FALSE,
      flag_treated_when_no_pc_required = FALSE,
      flag_requiring_conflict = FALSE
    ) %>%
    crossing(parasite = c("Hookworm disease", "Ascariasis", "Trichuriasis"))
) %>%
  arrange(iso3, year, parasite)

dup_merge <- pc_sac_for_merge %>% count(iso3, year, parasite, name = "n") %>% filter(n > 1)
if (nrow(dup_merge) > 0) {
  print(dup_merge, n = min(50, nrow(dup_merge)))
  stop("pc_sac_for_merge has duplicate keys: iso3, year, parasite", call. = FALSE)
}

write_csv(pc_sac_for_merge, OUT_MERGE)
message("Wrote: ", OUT_MERGE)

# Short run summary for reviewers
elig_summary <- pc_sac_for_model %>%
  filter(year %in% YEAR_MIN:YEAR_MAX) %>%
  count(parasite_pc, wt = as.integer(eligible), name = "eligible_country_years") %>%
  left_join(
    pc_sac_for_model %>%
      filter(year %in% YEAR_MIN:YEAR_MAX, eligible) %>%
      distinct(parasite_pc, iso3) %>%
      count(parasite_pc, name = "eligible_countries"),
    by = "parasite_pc"
  )
message("Eligibility summary (2010-2023):")
print(elig_summary)
