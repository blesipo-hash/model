# =============================================================================
# 01_build_inputs.R
# Purpose:
#   Build prevalence inputs from GBD draws using draw-by-draw WPP population weighting.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(countrycode)
})

# ---------------------------
# Helpers
# ---------------------------
norm_text <- function(x) {
  x %>% as.character() %>%
    str_replace_all("–", "-") %>%
    str_replace_all("’", "'") %>%
    str_squish()
}

# Parse ONLY the age labels that actually exist in raw GBD file
# NOTE: "12 to 23 months" is treated as age 1 (approximation for WPP single-age mapping).
parse_age_bounds <- function(age_name) {
  a <- norm_text(age_name)
  if (a == "12 to 23 months") return(tibble(age_start = 1L, age_end = 1L))
  if (a == "2 to 4") return(tibble(age_start = 2L, age_end = 4L))
  if (a == "95 plus") return(tibble(age_start = 95L, age_end = 100L))
  if (str_detect(a, "^\\d+ to \\d+$")) {
    m <- str_match(a, "^(\\d+) to (\\d+)$")
    return(tibble(age_start = as.integer(m[, 2]), age_end = as.integer(m[, 3])))
  }
  tibble(age_start = NA_integer_, age_end = NA_integer_)
}

clip01 <- function(x, eps_local = 1e-6) pmin(pmax(x, eps_local), 1 - eps_local)
logit  <- function(p) log(p / (1 - p))

# ---------------------------
# Config
# ---------------------------
EPS <- 1e-6
TARGET_AGES <- c("5 to 9", "10 to 14")  # SAC construction bins
DROP_AGES   <- c("Early Neonatal", "Late Neonatal", "1-5 months", "6-11 months")

# WHO PC country subselect for faster, MDA-focused runs.
# Examples:
#   SUBSELECT_TO_PC_COUNTRIES=1 Rscript R/01_build_inputs.R
#   PC_COUNTRY_LIST_PATH=data/processed/who_pc/countries_sac_pc_required_2010_2023.csv SUBSELECT_TO_PC_COUNTRIES=1 Rscript R/01_build_inputs.R
#   SUBSELECT_TO_PC_COUNTRIES=0 Rscript R/01_build_inputs.R
PC_COUNTRY_LIST_PATH <- Sys.getenv("PC_COUNTRY_LIST_PATH", unset = "data/processed/who_pc/countries_sac_pc_required_2010_2023.csv")
SUBSELECT_TO_PC_COUNTRIES <- as.integer(Sys.getenv("SUBSELECT_TO_PC_COUNTRIES", unset = "1")) == 1

pc_iso3 <- character(0)
if (SUBSELECT_TO_PC_COUNTRIES) {
  if (!file.exists(PC_COUNTRY_LIST_PATH)) {
    stop(
      "SUBSELECT_TO_PC_COUNTRIES=1 but PC country list file does not exist: ",
      PC_COUNTRY_LIST_PATH,
      call. = FALSE
    )
  }

  pc_iso3 <- read_csv(PC_COUNTRY_LIST_PATH, show_col_types = FALSE) %>%
    transmute(iso3 = toupper(str_squish(iso3))) %>%
    filter(!is.na(iso3), nchar(iso3) == 3) %>%
    distinct() %>%
    pull(iso3)

  if (length(pc_iso3) == 0) {
    stop("PC country ISO3 list is empty after cleaning; check country list file.", call. = FALSE)
  }
}

QC_DIR <- "data/processed/qc/01_build_inputs"
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Minimal QC helpers (delete later)
# ---------------------------
qc_check_unique <- function(df, keys, label, n_show = 20) {
  dups <- df %>% count(across(all_of(keys)), name = "n") %>% filter(n > 1) %>% arrange(desc(n))
  if (nrow(dups) > 0) {
    message("=== QC FAIL: duplicated keys at ", label, " ===")
    print(head(dups, n_show))
    stop("Duplicated keys detected at: ", label, call. = FALSE)
  } else {
    message("QC OK (unique keys): ", label)
  }
}

qc_draw_summary <- function(df, draw_cols, label) {
  m <- as.matrix(df[, draw_cols, drop = FALSE])
  message(
    "=== QC: ", label, " ===\n",
    "rows=", nrow(df),
    " | draw_cells=", length(m),
    " | NA_cells=", sum(is.na(m)),
    " | min=", signif(suppressWarnings(min(m, na.rm = TRUE)), 6),
    " | max=", signif(suppressWarnings(max(m, na.rm = TRUE)), 6),
    " | cells<0=", sum(m < 0, na.rm = TRUE),
    " | cells>1=", sum(m > 1, na.rm = TRUE)
  )
}

# =============================================================================
# 1) WPP single-age populations for 2010–2023 (ages 1+)
# =============================================================================
wpp_estimates <- read_excel(
  "data/raw/pop.xlsx",
  sheet = "Estimates",
  skip = 16,
  col_types = "text"
)

age_cols <- intersect(names(wpp_estimates), c(as.character(0:99), "100+"))

wpp_long <- wpp_estimates %>%
  rename(
    iso3_raw = `ISO3 Alpha-code`,
    country  = `Region, subregion, country or area *`,
    year     = Year,
    type     = Type
  ) %>%
  filter(norm_text(type) == "Country/Area") %>%
  mutate(country = norm_text(country)) %>%
  pivot_longer(cols = all_of(age_cols), names_to = "age_col", values_to = "pop_thousands") %>%
  mutate(
    iso3 = toupper(coalesce(
      na_if(str_trim(iso3_raw), ""),
      countrycode(country, "country.name", "iso3c")
    )),
    year = as.integer(year),
    age  = if_else(str_trim(age_col) == "100+", 100L, as.integer(str_trim(age_col))),
    pop  = parse_number(pop_thousands) * 1000
  ) %>%
  filter(
    year >= 2010, year <= 2023,
    !is.na(iso3), nchar(iso3) == 3,
    !is.na(age), age >= 1
  )

wpp_single <- wpp_long %>%
  filter(!is.na(pop), pop >= 0) %>%
  select(iso3, country, year, age, pop)

if (SUBSELECT_TO_PC_COUNTRIES) {
  wpp_single <- wpp_single %>% filter(iso3 %in% pc_iso3)
}

# QC_START wpp_single
qc_check_unique(wpp_single, c("iso3","year","age"), "wpp_single (iso3-year-age)")
if (any(is.na(wpp_single$pop)) || any(wpp_single$pop < 0, na.rm = TRUE)) {
  stop("WPP population has NA or negative values after parsing.", call. = FALSE)
}
# QC_END wpp_single

write_csv(wpp_single, "data/processed/un_wpp/wpp_single_age_2010_2023.csv")

# =============================================================================
# 2) GBD prevalence draws: drop unwanted ages + parse age bounds + rate->prevalence transform
# =============================================================================
gbd_draws_raw <- read_csv("data/raw/gbd/prev_proportion_draw_both.csv", show_col_types = FALSE)

draw_cols <- names(gbd_draws_raw)[str_detect(names(gbd_draws_raw), "^draw_\\d+$")]
if (length(draw_cols) == 0) {
  stop("No draw columns detected in raw GBD input (expected columns like draw_0, draw_1, ...).", call. = FALSE)
}

gbd_draws_raw <- gbd_draws_raw %>% mutate(across(all_of(draw_cols), as.numeric))

# QC_START gbd_raw_unique
qc_check_unique(gbd_draws_raw,
                c("location_id","year_id","cause_id","sex","age_group_id"),
                "gbd_draws_raw (location_id-year_id-cause_id-sex-age_group_id)")
# QC_END gbd_raw_unique

drop_ages <- DROP_AGES

gbd_draws_mapped <- gbd_draws_raw %>%
  mutate(
    location_name = norm_text(location_name),
    iso3 = toupper(countrycode(location_name, "country.name", "iso3c")),
    age_name = norm_text(age_group_name),
    year = as.integer(year_id)
  )

gbd_draws <- gbd_draws_mapped %>%
  filter(
    year >= 2010, year <= 2023,
    !is.na(iso3), nchar(iso3) == 3,
    !age_name %in% drop_ages
  )

if (SUBSELECT_TO_PC_COUNTRIES) {
  gbd_draws <- gbd_draws %>% filter(iso3 %in% pc_iso3)

  if (nrow(gbd_draws) == 0) {
    stop(
      "GBD rows after PC country subselect is zero; check ISO3 mapping or the country list file.",
      call. = FALSE
    )
  }

  gbd_iso3_after <- gbd_draws %>% distinct(iso3) %>% pull(iso3)
  missing_pc_iso3 <- setdiff(pc_iso3, gbd_iso3_after)

  message("=== PC country subselect summary ===")
  message("n_pc_countries = ", length(pc_iso3))
  message("n_gbd_countries_after_filter = ", length(gbd_iso3_after))
  if (length(missing_pc_iso3) > 0) {
    message("top_10_missing_pc_iso3_in_gbd = ", paste(head(sort(missing_pc_iso3), 10), collapse = ", "))
  } else {
    message("top_10_missing_pc_iso3_in_gbd = none")
  }
}

age_bounds <- gbd_draws %>%
  distinct(age_name) %>%
  mutate(tmp = map(age_name, parse_age_bounds)) %>%
  unnest(tmp)

# QC_START age_bounds_parse
bad_age <- age_bounds %>% filter(is.na(age_start) | is.na(age_end))
if (nrow(bad_age) > 0) {
  write_csv(bad_age, file.path(QC_DIR, "age_bounds_unparsed.csv"))
  stop("Some age_group_name values could not be parsed to age bounds. See QC output.", call. = FALSE)
}
# QC_END age_bounds_parse

# Scale detection (prefer metric column)
metric_vals <- unique(norm_text(gbd_draws$metric))
is_rate <- any(str_detect(tolower(metric_vals), "rate"))
if (!is_rate) {
  max_draw_val <- suppressWarnings(max(as.matrix(gbd_draws[, draw_cols, drop = FALSE]), na.rm = TRUE))
  is_rate <- is.finite(max_draw_val) && max_draw_val > 1
}

# QC_START scale_detection
message("=== QC: scale_detection ===")
message("metric values: ", paste(metric_vals, collapse = " | "))
message("is_rate = ", is_rate)
# QC_END scale_detection

# Convert rates -> prevalence proportions (NO capping here)
if (is_rate) {
  gbd_draws <- gbd_draws %>%
    mutate(across(all_of(draw_cols), ~ .x / 100000))
} else {
  # already proportions; do not clip yet (we want to *detect* violations first)
  gbd_draws <- gbd_draws
}

# Attach age bounds
gbd_draws <- gbd_draws %>% left_join(age_bounds, by = "age_name")

# QC_START post_conversion_summary
qc_draw_summary(gbd_draws, draw_cols, "post_conversion (before exclude)")
# QC_END post_conversion_summary

# ---------------------------
# EXCLUDE POLICY (SAC only):
# drop SAC keys where ANY draw is outside [0,1] in either bin
# ---------------------------
gbd_sac <- gbd_draws %>% filter(age_name %in% TARGET_AGES)

# QC_START sac_presence
if (nrow(gbd_sac) == 0) stop("No SAC rows found for TARGET_AGES; check labels.", call. = FALSE)
# QC_END sac_presence

m_sac <- as.matrix(gbd_sac[, draw_cols, drop = FALSE])
row_min <- apply(m_sac, 1, min, na.rm = TRUE)
row_max <- apply(m_sac, 1, max, na.rm = TRUE)
viol_row <- (row_min < 0) | (row_max > 1)

viol_keys <- gbd_sac[viol_row, c("iso3","year","cause_id","cause_name","sex")] %>%
  distinct()

# QC_START exclude_impact
message("=== QC: exclude_impact (SAC) ===")
message("SAC viol rows = ", sum(viol_row), " / ", nrow(gbd_sac))
message("SAC viol keys = ", nrow(viol_keys), " / ",
        nrow(gbd_sac %>% distinct(iso3, year, cause_id, sex)))
write_csv(viol_keys, file.path(QC_DIR, "viol_keys_sac_exclude.csv"))
# QC_END exclude_impact

# Drop only SAC rows for those keys (both bins), keep other ages intact
gbd_draws <- gbd_draws %>%
  filter(!age_name %in% TARGET_AGES) %>%
  bind_rows(
    gbd_draws %>%
      filter(age_name %in% TARGET_AGES) %>%
      anti_join(viol_keys %>% select(iso3, year, cause_id, sex), by = c("iso3","year","cause_id","sex"))
  ) %>%
  arrange(iso3, year, cause_id, sex, age_name)

# Now stability clamp for modeling/logit
gbd_draws <- gbd_draws %>%
  mutate(across(all_of(draw_cols), ~ clip01(.x, EPS)))

# QC_START post_exclude_bounds
qc_draw_summary(gbd_draws %>% filter(age_name %in% TARGET_AGES), draw_cols, "post_exclude (SAC after clip)")
# QC_END post_exclude_bounds

# =============================================================================
# 3) Map WPP single-age population into each kept GBD age bin
# =============================================================================
age_bounds_clean <- age_bounds %>%
  filter(!is.na(age_start), !is.na(age_end)) %>%
  distinct(age_name, age_start, age_end)

age_map <- age_bounds_clean %>%
  mutate(age = map2(age_start, age_end, ~ seq(.x, .y))) %>%
  unnest(age)

wpp_in_gbd_bins <- wpp_single %>%
  inner_join(age_map, by = "age") %>%
  group_by(iso3, year, age_name, age_start, age_end) %>%
  summarise(population = sum(pop, na.rm = TRUE), .groups = "drop")

# QC_START wpp_bins
qc_check_unique(wpp_in_gbd_bins, c("iso3","year","age_name"), "wpp_in_gbd_bins (iso3-year-age_name)")
bad_target <- wpp_in_gbd_bins %>% filter(age_name %in% TARGET_AGES, population <= 0)
if (nrow(bad_target) > 0) {
  write_csv(bad_target, file.path(QC_DIR, "wpp_bins_bad_target_ages.csv"))
  stop("WPP population <= 0 for TARGET_AGES. See QC output.", call. = FALSE)
}
# QC_END wpp_bins

write_csv(wpp_in_gbd_bins, "data/processed/un_wpp/wpp_population_gbd_bins.csv")

# =============================================================================
# 4) Draw-by-draw weighted aggregation into preSAC / SAC / adult
# =============================================================================
gbd_joined_pop <- gbd_draws %>%
  left_join(
    wpp_in_gbd_bins %>% select(iso3, year, age_name, population),
    by = c("iso3", "year", "age_name")
  )

# QC_START join_pop
if (any(is.na(gbd_joined_pop$population))) {
  miss <- gbd_joined_pop %>% filter(is.na(population)) %>%
    select(iso3, year, age_name, location_name) %>% distinct()
  write_csv(miss, file.path(QC_DIR, "missing_population_after_join.csv"))
  stop("Missing population after join. See QC output.", call. = FALSE)
}
# QC_END join_pop

gbd_weighted_base <- gbd_joined_pop %>%
  filter(population > 0) %>%
  mutate(
    target_age_group = case_when(
      age_start >= 1L  & age_end <= 4L  ~ "preSAC",
      age_start >= 5L  & age_end <= 14L ~ "SAC",
      age_start >= 15L                 ~ "adult",
      TRUE                             ~ NA_character_
    )
  ) %>%
  filter(!is.na(target_age_group))

# Population-weighted within-draw aggregation:
# cases_draw = p_draw * pop; p_agg_draw = sum(cases_draw)/sum(pop)
weighted_draws <- gbd_weighted_base %>%
  group_by(iso3, location_name, year, cause_name, cause_id, sex, target_age_group) %>%
  group_modify(~{
    P   <- as.matrix(.x[, draw_cols, drop = FALSE])
    pop <- .x$population
    num <- colSums(P * pop, na.rm = TRUE)
    den <- rep(sum(pop), length(num))
    p_agg <- num / den
    tibble(total_population = sum(pop), !!!setNames(as.list(p_agg), draw_cols))
  }) %>%
  ungroup() %>%
  rename(parasite = cause_name, age_group = target_age_group)

# QC_START weighted_bounds
qc_check_unique(weighted_draws, c("iso3","year","cause_id","sex","age_group"), "weighted_draws (iso3-year-cause_id-sex-age_group)")
qc_draw_summary(weighted_draws, draw_cols, "weighted_draws (post-aggregation)")
m_w <- as.matrix(weighted_draws[, draw_cols, drop = FALSE])
if (sum(m_w < 0 | m_w > 1, na.rm = TRUE) > 0) stop("Weighted draws out of bounds.", call. = FALSE)
# QC_END weighted_bounds

write_csv(weighted_draws, "data/processed/gbd/gbd_country_age_draws_weighted.csv")

# =============================================================================
# 5) Summaries from draw-level outputs
# =============================================================================
draw_mat <- as.matrix(weighted_draws[, draw_cols, drop = FALSE])

gbd_country_age <- weighted_draws %>%
  transmute(
    iso3, location_name, year, parasite, age_group, sex, total_population,
    weighted_prev    = rowMeans(draw_mat, na.rm = TRUE),
    weighted_prev_lo = apply(draw_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    weighted_prev_hi = apply(draw_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    weighted_se_p    = apply(draw_mat, 1, sd, na.rm = TRUE),
    sigma_logit      = apply(draw_mat, 1, function(v) sd(logit(clip01(v, EPS)), na.rm = TRUE))
  )

write_csv(gbd_country_age, "data/processed/gbd/gbd_country_age_with_uncertainty.csv")

# =============================================================================
# 6) SAC-only panel + contract file for downstream Stan prep
# =============================================================================
gbd_sac_panel <- gbd_country_age %>%
  filter(age_group == "SAC") %>%
  transmute(
    iso3, location_name, year, parasite, sex,
    prev_obs    = clip01(weighted_prev, EPS),
    prev_lo     = clip01(weighted_prev_lo, EPS),
    prev_hi     = clip01(weighted_prev_hi, EPS),
    sigma_logit = as.numeric(sigma_logit),
    pop_sac     = as.numeric(total_population)
  )

# QC_START final_panel
if (any(is.na(gbd_sac_panel$prev_obs))) stop("prev_obs has NA in final SAC panel.", call. = FALSE)
# QC_END final_panel

write_csv(gbd_sac_panel, "data/processed/gbd/gbd_sac_panel_for_inference.csv")

who_pc_path <- "data/processed/who_pc/pc_sac_for_merge_gbd_causes.csv"
if (!file.exists(who_pc_path)) {
  stop("WHO PC merge file not found: data/processed/who_pc/pc_sac_for_merge_gbd_causes.csv. Run 01b_build_who_pc_inputs.R first.", call. = FALSE)
}

who_pc <- read_csv(who_pc_path, show_col_types = FALSE) %>%
  transmute(
    iso3 = toupper(str_squish(iso3)),
    year = as.integer(year),
    parasite = norm_text(parasite),
    pc_cov = as.numeric(coverage_use),
    pc_status = as.character(coverage_use_source)
  )

qc_check_unique(who_pc, c("iso3","year","parasite"), "who_pc (iso3-year-parasite)")

panel_keys <- gbd_sac_panel %>% distinct(iso3, year, parasite)
overlap_n <- panel_keys %>%
  inner_join(who_pc %>% distinct(iso3, year, parasite), by = c("iso3", "year", "parasite")) %>%
  nrow()
overlap_rate <- overlap_n / nrow(panel_keys)

if (overlap_n == 0) stop("WHO PC merge overlap with GBD SAC panel is zero.", call. = FALSE)
if (overlap_rate < 0.05) stop(sprintf("WHO PC merge overlap is %.1f%% (<5%%).", 100 * overlap_rate), call. = FALSE)
if (overlap_rate < 0.50) warning(sprintf("WHO PC merge overlap is %.1f%% (<50%%).", 100 * overlap_rate), call. = FALSE)

panel <- gbd_sac_panel %>%
  mutate(y_obs = logit(prev_obs)) %>%
  left_join(who_pc, by = c("iso3", "year", "parasite"))

if (nrow(panel) != nrow(gbd_sac_panel)) {
  stop("WHO PC join multiplied or dropped rows in final panel.", call. = FALSE)
}

write_csv(panel, "data/processed/model_input/sac_panel_merged.csv")
