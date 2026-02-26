# =============================================================================
# R/01c_make_pc_required_country_lists.R
# Build WHO SAC-PC-required country lists for 2010-2023.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(countrycode)
})

YEAR_MIN <- 2010L
YEAR_MAX <- 2023L
YEARS <- YEAR_MIN:YEAR_MAX

IN_PROCESSED <- "data/processed/who_pc/pc_sac_for_model.csv"
SCH_RAW <- "data/raw/who_pc/schistosomiasis-treatment-coverage.csv"
STH_RAW <- "data/raw/who_pc/who_pc_sth.xlsx"

OUT_DIR <- "data/processed/who_pc"
OUT_COUNTRY_LIST <- file.path(OUT_DIR, "countries_sac_pc_required_2010_2023.csv")
OUT_COUNTRY_YEAR <- file.path(OUT_DIR, "country_year_sac_pc_required_2010_2023.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

norm_txt <- function(x) {
  str_squish(str_replace_all(str_replace_all(as.character(x), "–", "-"), "’", "'"))
}

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

load_from_processed <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)

  req_cols <- c("iso3", "year", "parasite_pc")
  if (!any(c("pop_req_pc", "requiring_pc") %in% names(df))) {
    stop("Processed WHO PC file missing pop requirement column (pop_req_pc or requiring_pc).", call. = FALSE)
  }
  miss <- setdiff(req_cols, names(df))
  if (length(miss) > 0) {
    stop("Processed WHO PC file missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df %>%
    transmute(
      iso3 = str_to_upper(str_squish(as.character(iso3))),
      country = if ("country" %in% names(df)) norm_txt(country) else NA_character_,
      region = if ("region" %in% names(df)) norm_txt(region) else countrycode(iso3, "iso3c", "region"),
      year = as.integer(year),
      parasite_pc = str_to_upper(str_squish(as.character(parasite_pc))),
      pop_req_pc = as.numeric(if ("pop_req_pc" %in% names(df)) pop_req_pc else requiring_pc)
    )
}

load_from_raw <- function(sch_path, sth_path) {
  if (!file.exists(sch_path)) stop("Missing raw SCH file: ", sch_path, call. = FALSE)
  if (!file.exists(sth_path)) stop("Missing raw STH file: ", sth_path, call. = FALSE)

  raw_sch <- read_csv(sch_path, show_col_types = FALSE)
  sch_col <- names(raw_sch)[
    str_detect(str_to_lower(names(raw_sch)), "population\\s*requiring") &
      str_detect(str_to_lower(names(raw_sch)), "pc")
  ]
  if (length(sch_col) == 0) {
    stop("No SCH 'population requiring PC' column found in ", sch_path, call. = FALSE)
  }

  sch <- raw_sch %>%
    transmute(
      iso3 = str_to_upper(str_squish(as.character(Code))),
      country = norm_txt(Entity),
      region = countrycode(iso3, "iso3c", "region"),
      year = as.integer(Year),
      parasite_pc = "SCH",
      pop_req_pc = parse_pop_req(.data[[sch_col[[1]]]])
    )

  raw_sth <- read_excel(sth_path)
  req_sth_col <- "Population requiring PC for STH, SAC"
  if (!req_sth_col %in% names(raw_sth)) {
    stop("STH file missing required column: ", req_sth_col, call. = FALSE)
  }

  sth <- raw_sth %>%
    transmute(
      country = norm_txt(country),
      iso3 = {
        cname <- recode(country, "Malasya" = "Malaysia")
        code_clean <- str_to_upper(str_squish(as.character(country_code)))
        valid_iso3 <- unique(countrycode::codelist$iso3c)
        use_code <- !is.na(code_clean) & str_length(code_clean) == 3 & code_clean %in% valid_iso3
        out <- if_else(use_code, code_clean, countrycode(cname, "country.name", "iso3c"))
        str_to_upper(out)
      },
      region = if ("region" %in% names(raw_sth)) norm_txt(raw_sth$region) else countrycode(iso3, "iso3c", "region"),
      year = as.integer(year),
      parasite_pc = "STH",
      pop_req_pc = parse_pop_req(.data[[req_sth_col]])
    )

  bind_rows(sch, sth)
}

pc <- if (file.exists(IN_PROCESSED)) {
  message("Using processed input: ", IN_PROCESSED)
  load_from_processed(IN_PROCESSED)
} else {
  message("Processed input not found; using fallback raw WHO files.")
  load_from_raw(SCH_RAW, STH_RAW)
}

pc <- pc %>%
  filter(year %in% YEARS) %>%
  mutate(
    parasite_pc = str_to_upper(parasite_pc),
    iso3 = str_to_upper(iso3)
  )

# Assertions
if (any(is.na(pc$iso3) | nchar(pc$iso3) != 3L)) {
  bad <- pc %>% filter(is.na(iso3) | nchar(iso3) != 3L) %>% head(20)
  print(bad)
  stop("Assertion failed: iso3 must be non-missing and length 3", call. = FALSE)
}

if (any(is.na(pc$year) | pc$year < YEAR_MIN | pc$year > YEAR_MAX)) {
  bad <- pc %>% filter(is.na(year) | year < YEAR_MIN | year > YEAR_MAX) %>% head(20)
  print(bad)
  stop("Assertion failed: year must be integer and within 2010-2023", call. = FALSE)
}

if (any(!is.na(pc$pop_req_pc) & !is.finite(pc$pop_req_pc))) {
  stop("Assertion failed: pop_req_pc must be numeric", call. = FALSE)
}
if (any(!is.na(pc$pop_req_pc) & pc$pop_req_pc < 0)) {
  bad <- pc %>% filter(!is.na(pop_req_pc) & pop_req_pc < 0) %>% head(20)
  print(bad)
  stop("Assertion failed: pop_req_pc must be >= 0", call. = FALSE)
}

country_year <- pc %>%
  filter(!is.na(pop_req_pc), pop_req_pc > 0) %>%
  select(iso3, year, parasite_pc, pop_req_pc) %>%
  arrange(parasite_pc, iso3, year)

# Uniqueness assertion for country-year panel
cy_dups <- country_year %>% count(iso3, year, parasite_pc, name = "n") %>% filter(n > 1)
if (nrow(cy_dups) > 0) {
  print(cy_dups, n = min(50, nrow(cy_dups)))
  stop("Assertion failed: country-year panel is not unique on iso3-year-parasite_pc", call. = FALSE)
}

pc_missing_summary <- pc %>%
  group_by(iso3, parasite_pc) %>%
  summarise(any_years_missing_pop_req = any(is.na(pop_req_pc)), .groups = "drop")

country_list <- country_year %>%
  group_by(iso3, parasite_pc) %>%
  summarise(
    first_year_required = min(year),
    last_year_required = max(year),
    n_years_required = n_distinct(year),
    .groups = "drop"
  ) %>%
  left_join(
    pc %>%
      group_by(iso3, parasite_pc) %>%
      summarise(
        country = first(na.omit(country)),
        region = first(na.omit(region)),
        .groups = "drop"
      ),
    by = c("iso3", "parasite_pc")
  ) %>%
  mutate(
    country = if_else(is.na(country), countrycode(iso3, "iso3c", "country.name"), country),
    region = if_else(is.na(region), countrycode(iso3, "iso3c", "region"), region)
  ) %>%
  left_join(pc_missing_summary, by = c("iso3", "parasite_pc")) %>%
  select(iso3, country, region, parasite_pc, first_year_required, last_year_required,
         n_years_required, any_years_missing_pop_req) %>%
  arrange(parasite_pc, iso3)

# Uniqueness assertion for country list
cl_dups <- country_list %>% count(iso3, parasite_pc, name = "n") %>% filter(n > 1)
if (nrow(cl_dups) > 0) {
  print(cl_dups, n = min(50, nrow(cl_dups)))
  stop("Assertion failed: country list is not unique on iso3-parasite_pc", call. = FALSE)
}

write_csv(country_list, OUT_COUNTRY_LIST)
write_csv(country_year, OUT_COUNTRY_YEAR)

message("Wrote: ", OUT_COUNTRY_LIST)
message("Wrote: ", OUT_COUNTRY_YEAR)

message("Eligible country-years by parasite_pc")
print(country_year %>% count(parasite_pc, name = "n_country_years") %>% arrange(parasite_pc))

message("Eligible countries by parasite_pc")
print(country_list %>% count(parasite_pc, name = "n_countries") %>% arrange(parasite_pc))

message("Top 10 countries by n_years_required per parasite")
print(
  country_list %>%
    group_by(parasite_pc) %>%
    arrange(desc(n_years_required), iso3, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    select(parasite_pc, iso3, country, n_years_required, first_year_required, last_year_required)
)
