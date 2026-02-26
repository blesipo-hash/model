# =============================================================================
# 01c_build_daly_in.R
# Build DALY_IN for SAC from raw GBD burden extracts.
# Output: iso3, year, parasite, daly_gbd (unique by iso3-year-parasite)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(janitor)
  library(countrycode)
})

DALY_DIR <- "data/raw/daly"
PATTERN <- "^gbd_daly_.*\\.csv$"
DALY_IN <- "data/processed/gbd/daly_sac_by_iso3_year_parasite.csv"
QC_DIR <- "data/processed/model_input/qc_daly_build"
dir.create(dirname(DALY_IN), recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

CAUSES <- c("Hookworm disease", "Ascariasis", "Trichuriasis", "Schistosomiasis")

need_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(name, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

assert_unique <- function(df, keys, nm) {
  dup <- df %>% count(across(all_of(keys)), name = "n") %>% filter(n > 1)
  if (nrow(dup) > 0) {
    readr::write_csv(dup, file.path(QC_DIR, paste0("dup_", nm, ".csv")))
    stop("Duplicate keys in ", nm, ". See QC output in ", QC_DIR, call. = FALSE)
  }
}

daly_files <- list.files(DALY_DIR, full.names = TRUE, pattern = PATTERN)
if (length(daly_files) == 0) stop("No DALY files found in: ", DALY_DIR, call. = FALSE)

gbd_raw <- map_dfr(daly_files, ~ read_csv(.x, show_col_types = FALSE)) %>% clean_names()

need_cols(
  gbd_raw,
  c("measure_name", "metric_name", "sex_name", "age_id", "cause_name", "location_name", "val"),
  name = "GBD burden extract(s)"
)

if (!any(c("year", "year_id") %in% names(gbd_raw))) {
  stop("GBD burden extract(s) must include either year or year_id.", call. = FALSE)
}

gbd_burden <- gbd_raw %>%
  mutate(
    year = as.integer(coalesce(
      if ("year" %in% names(gbd_raw)) as.character(year) else NA_character_,
      if ("year_id" %in% names(gbd_raw)) as.character(year_id) else NA_character_
    )),
    measure_name_l = str_to_lower(str_squish(measure_name)),
    metric_name_l = str_to_lower(str_squish(metric_name)),
    sex_name_l = str_to_lower(str_squish(sex_name)),
    age_id = as.integer(age_id),
    cause_name = as.character(cause_name),
    location_name = as.character(location_name),
    val = as.numeric(val),
    lower = as.numeric(if ("lower" %in% names(gbd_raw)) lower else NA_real_),
    upper = as.numeric(if ("upper" %in% names(gbd_raw)) upper else NA_real_)
  )

# DALY only (exclude YLD/YLL and any other burden measures), number metric, both sexes, SAC age_id 23, 2010-2023.
sac_daly_base <- gbd_burden %>%
  filter(
    year >= 2010, year <= 2023,
    str_detect(measure_name_l, "\\bdaly"),
    !str_detect(measure_name_l, "\\byld"),
    metric_name_l == "number",
    sex_name_l %in% c("both", "both sexes"),
    age_id == 23,
    cause_name %in% CAUSES
  )

iso3_fail <- sac_daly_base %>%
  mutate(iso3_chk = countrycode(location_name, "country.name", "iso3c")) %>%
  filter(is.na(iso3_chk)) %>%
  distinct(location_name) %>%
  arrange(location_name)

sac_daly_raw <- sac_daly_base %>%
  mutate(
    iso3 = countrycode(location_name, "country.name", "iso3c"),
    parasite = cause_name
  ) %>%
  filter(!is.na(iso3)) %>%
  select(iso3, year, parasite, val, lower, upper)

if (nrow(sac_daly_raw) == 0) {
  stop("No DALY rows remained after filters (2010-2023, DALY measure, number metric, both sexes, age_id 23, target causes).", call. = FALSE)
}

dup_diag <- sac_daly_raw %>%
  group_by(iso3, year, parasite) %>%
  summarise(
    n = n(),
    nd_val = n_distinct(val),
    nd_lower = n_distinct(lower),
    nd_upper = n_distinct(upper),
    .groups = "drop"
  ) %>%
  filter(n > 1)

if (nrow(dup_diag) > 0) {
  readr::write_csv(dup_diag, file.path(QC_DIR, "dup_iso3_year_parasite.csv"))
  bad <- dup_diag %>% filter(nd_val > 1 | nd_lower > 1 | nd_upper > 1)
  if (nrow(bad) > 0) {
    readr::write_csv(bad, file.path(QC_DIR, "dup_conflicting_values.csv"))
    stop("Duplicate DALY rows with conflicting values detected. See QC output in ", QC_DIR, call. = FALSE)
  }
  message("Duplicate DALY keys with identical values found; deduplicating safely.")
}

daly_in <- sac_daly_raw %>%
  distinct(iso3, year, parasite, val) %>%
  group_by(iso3, year, parasite) %>%
  summarise(daly_gbd = first(val), .groups = "drop") %>%
  arrange(iso3, year, parasite)

assert_unique(daly_in, c("iso3", "year", "parasite"), "daly_in_iso3_year_parasite")

coverage_qc <- daly_in %>%
  summarise(
    rows = n(),
    years_min = min(year, na.rm = TRUE),
    years_max = max(year, na.rm = TRUE),
    n_iso3 = n_distinct(iso3),
    n_parasite = n_distinct(parasite),
    n_missing_daly = sum(is.na(daly_gbd))
  )
print(coverage_qc)

readr::write_csv(daly_in, DALY_IN)
message("Wrote DALY_IN: ", DALY_IN)

if (nrow(iso3_fail) > 0) {
  write_csv(iso3_fail, file.path(QC_DIR, "location_name_iso3_mapping_failed.csv"))
  message("NOTE: Some location_name could not be mapped to iso3. QC: ",
          file.path(QC_DIR, "location_name_iso3_mapping_failed.csv"))
}

# ---------------------------
# QC: require full coverage of panel keys (prevents 02 from failing later)
# ---------------------------
PANEL_IN <- "data/processed/model_input/sac_panel_merged.csv"
if (!file.exists(PANEL_IN)) stop("Missing panel for coverage QC: ", PANEL_IN, call. = FALSE)

panel_keys <- read_csv(PANEL_IN, show_col_types = FALSE) %>%
  transmute(
    iso3 = toupper(str_squish(as.character(iso3))),
    year = as.integer(year),
    parasite = as.character(parasite)
  ) %>%
  distinct()

missing_keys <- panel_keys %>%
  anti_join(daly_in, by = c("iso3","year","parasite")) %>%
  arrange(iso3, year, parasite)

if (nrow(missing_keys) > 0) {
  write_csv(missing_keys, file.path(QC_DIR, "panel_keys_missing_daly.csv"))
  stop("DALY_IN does not cover all panel keys. QC: ",
       file.path(QC_DIR, "panel_keys_missing_daly.csv"),
       call. = FALSE)
}
