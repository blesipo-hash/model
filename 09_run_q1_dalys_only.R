# =============================================================================
# R/09_run_q1_dalys_only.R
# Q1-only end-to-end runner
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

YEAR_MIN <- 2010L
YEAR_MAX <- 2023L
SCENARIOS <- c("NO_PC", "realized", "scaleup_set_75", "scaleup_floor_75", "scaleup_100")
RUN_PPC <- FALSE
DALY_ONLY <- TRUE
DALY_VALUE_USD <- Sys.getenv("DALY_VALUE_USD", "1000")
OUT_DIR <- "outputs/q1_dalys"
LOG_PATH <- file.path(OUT_DIR, "run_log_q1.txt")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(character(), LOG_PATH)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  write(msg, file = LOG_PATH, append = TRUE)
}

run_step <- function(name, script, required_outputs = character(), skip_if_outputs_exist = TRUE, env = character()) {
  already <- length(required_outputs) > 0 && all(file.exists(required_outputs))
  if (skip_if_outputs_exist && already) {
    log_msg("SKIP ", name, " (outputs already exist)")
    return(invisible(TRUE))
  }

  log_msg("RUN  ", name, " -> ", script)
  status <- system2("Rscript", c(script), stdout = TRUE, stderr = TRUE, env = env)
  attr_status <- attr(status, "status")
  if (!is.null(attr_status) && attr_status != 0) {
    writeLines(status)
    stop("Step failed: ", name, " (", script, ")", call. = FALSE)
  }

  if (length(required_outputs) > 0 && !all(file.exists(required_outputs))) {
    stop("Step completed but required outputs missing for ", name, ": ",
         paste(required_outputs[!file.exists(required_outputs)], collapse = ", "), call. = FALSE)
  }

  log_msg("DONE ", name)
  invisible(TRUE)
}

log_msg("Starting Q1 DALY-only run")

run_step("00_make_err_priors", "00_make_err_priors.R", c("data/parameters/drug_err_default.csv"))
run_step("00_make_k_priors", "00_make_k_priors.R", c("data/parameters/k_priors_default_minimal.csv"))
run_step("01b_build_who_pc_inputs", "01b_build_who_pc_inputs.R", c("data/processed/who_pc/pc_sac_for_model.csv", "data/processed/who_pc/pc_sac_for_merge_gbd_causes.csv"))
run_step("01c_make_pc_required_country_lists", "R/01c_make_pc_required_country_lists.R", c("data/processed/who_pc/countries_sac_pc_required_2010_2023.csv", "data/processed/who_pc/country_year_sac_pc_required_2010_2023.csv"))
run_step("01_build_inputs", "01_build_inputs.R", c("data/processed/gbd/gbd_sac_panel_for_inference.csv"))
run_step("01c_build_adjacency_edges", "01c_build_adjacency_edges.R", c("data/input/adjacency_edges.csv"))
run_step("01c_build_daly_in", "01c_build_daly_in.R", c("data/processed/gbd/daly_sac_by_iso3_year_parasite.csv"))
run_step("02_prepare_stan_inputs", "02_prepare_stan_inputs.R", c("data/processed/stan_data_sac.rds", "data/processed/model_input/sac_panel_model.csv"))
run_step("03_fit_model", "03_fit_model.R", c("data/processed/fit_sac.rds", "data/processed/fit_sac_summary.csv"))

if (RUN_PPC) {
  run_step("04_ppc", "04_ppc.R", c("outputs/ppc/ppc_pred_summary_selected.csv"), skip_if_outputs_exist = FALSE)
} else {
  log_msg("SKIP 04_ppc (RUN_PPC=FALSE)")
}

run_step("05_counterfactuals", "05_counterfactuals.R", c("outputs/counterfactuals/cf_prevalence_summary.csv"), skip_if_outputs_exist = FALSE)
run_step(
  "06_dalys_costs_cba DALY_ONLY",
  "06_dalys_costs_cba.R",
  c(
    file.path(OUT_DIR, "dalys_by_country_year_parasite_scenario.csv"),
    file.path(OUT_DIR, "global_summary_dalys.csv"),
    file.path(OUT_DIR, "implementation_gap_incremental_vs_realized.csv"),
    file.path(OUT_DIR, "top_countries_by_dalys_averted.csv")
  ),
  skip_if_outputs_exist = FALSE,
  env = c("DALY_ONLY=1", paste0("OUT_DIR=", OUT_DIR), paste0("DALY_VALUE_USD=", DALY_VALUE_USD))
)

# QC checks
cf <- read_csv("outputs/counterfactuals/cf_prevalence_summary.csv", show_col_types = FALSE)
scn_missing <- setdiff(SCENARIOS, unique(cf$scenario))
if (length(scn_missing) > 0) {
  stop("Counterfactual scenarios missing: ", paste(scn_missing, collapse = ", "), call. = FALSE)
}

q1 <- read_csv(file.path(OUT_DIR, "dalys_by_country_year_parasite_scenario.csv"), show_col_types = FALSE)
if (nrow(q1) == 0) stop("Q1 output table is empty.", call. = FALSE)
if (any(q1$year < YEAR_MIN | q1$year > YEAR_MAX, na.rm = TRUE)) stop("Q1 contains years outside 2010-2023.", call. = FALSE)
if (any(is.na(q1$pop_req_pc) | q1$pop_req_pc <= 0, na.rm = TRUE)) stop("Q1 contains non-eligible rows.", call. = FALSE)
if (any(q1$coverage_scn < 0 | q1$coverage_scn > 1, na.rm = TRUE)) stop("coverage_scn outside [0,1].", call. = FALSE)
if (any(abs(q1$treated_n - (q1$coverage_scn * q1$pop_req_pc)) > 1e-6, na.rm = TRUE)) stop("treated_n != coverage_scn * pop_req_pc.", call. = FALSE)

dup <- q1 %>% count(iso3, year, parasite, scenario_id) %>% filter(n > 1)
if (nrow(dup) > 0) stop("Duplicate keys in q1 output.", call. = FALSE)

glob <- read_csv(file.path(OUT_DIR, "global_summary_dalys.csv"), show_col_types = FALSE)
realized_total <- glob %>% filter(scenario_id == "realized") %>% summarize(v = sum(global_dalys_averted, na.rm = TRUE)) %>% pull(v)
if (length(realized_total) == 0 || is.na(realized_total) || realized_total == 0) {
  stop("Global realized DALYs averted is zero.", call. = FALSE)
}

log_msg("QC passed for Q1 DALY-only outputs")
log_msg("Completed Q1 DALY-only run")
