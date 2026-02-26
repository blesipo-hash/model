# =============================================================================
# 03_fit_model.R
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(readr)
  library(dplyr)
  library(stringr)
})

# ---------------------------
# helpers
# ---------------------------
get_script_path <- function() {
  args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(args_file) > 0) return(sub("^--file=", "", args_file[[1]]))
  if (!is.null(sys.frames()[[1]]$ofile)) return(sys.frames()[[1]]$ofile)
  return(NA_character_)
}

pick_root_dir <- function(script_path) {
  candidates <- character(0)
  
  if (is.finite(nchar(script_path)) && !is.na(script_path) && nzchar(script_path)) {
    script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
    candidates <- c(
      script_dir,                               # e.g., repo/R
      normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE) # e.g., repo
    )
  }
  
  candidates <- c(candidates, normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  candidates <- unique(candidates)
  
  # choose first candidate that contains the expected stan_data
  for (root in candidates) {
    if (file.exists(file.path(root, "data/processed/stan_data_sac.rds"))) return(root)
  }
  
  # fallback: first candidate
  candidates[[1]]
}

# ---------------------------
# paths
# ---------------------------
script_path <- get_script_path()
ROOT_DIR <- pick_root_dir(script_path)

DATA_RDS <- file.path(ROOT_DIR, "data/processed/stan_data_sac.rds")

MODEL_STAN_CANDIDATES <- c(
  file.path(ROOT_DIR, "stan/03_model_sac.stan"),
  file.path(ROOT_DIR, "03_model_sac.stan")
)
MODEL_STAN <- MODEL_STAN_CANDIDATES[file.exists(MODEL_STAN_CANDIDATES)][1]

OUT_FIT_RDS <- file.path(ROOT_DIR, "data/processed/fit_sac.rds")
OUT_SUMMARY <- file.path(ROOT_DIR, "data/processed/fit_sac_summary.csv")
OUT_DIAG    <- file.path(ROOT_DIR, "data/processed/fit_sac_cmdstan_diagnose.txt")

dir.create(file.path(ROOT_DIR, "data/processed"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# checks
# ---------------------------
if (!file.exists(DATA_RDS)) {
  stop("Missing Stan data. Expected at: ", DATA_RDS,
       "\nRun 02_prepare_stan_inputs.R first (and ensure you are in the repo root).",
       call. = FALSE)
}

if (length(MODEL_STAN) == 0 || is.na(MODEL_STAN)) {
  stop("Missing Stan model file. Tried:\n- ",
       paste(MODEL_STAN_CANDIDATES, collapse = "\n- "),
       call. = FALSE)
}

stan_data <- readRDS(DATA_RDS)

# ---------------------------
# compile
# ---------------------------
mod <- cmdstan_model(
  MODEL_STAN,
  cpp_options = list(stan_threads = TRUE)
)

# ---------------------------
# sample settings
# ---------------------------
SEED <- 1
CHAINS <- 4
WARMUP <- 1000
SAMPLES <- 1000
THREADS_PER_CHAIN <- 1

fit <- mod$sample(
  data = stan_data,
  seed = SEED,
  chains = CHAINS,
  parallel_chains = CHAINS,
  iter_warmup = WARMUP,
  iter_sampling = SAMPLES,
  adapt_delta = 0.95,
  max_treedepth = 12,
  threads_per_chain = THREADS_PER_CHAIN,
  refresh = 200
)

# ---------------------------
# save fit + diagnostics
# ---------------------------
fit$save_object(file = OUT_FIT_RDS)

sum_df <- fit$summary()
write_csv(sum_df, OUT_SUMMARY)

diag_txt <- fit$cmdstan_diagnose()
writeLines(diag_txt, OUT_DIAG)

message("ROOT_DIR: ", ROOT_DIR)
message("Wrote: ", OUT_FIT_RDS)
message("Wrote: ", OUT_SUMMARY)
message("Wrote: ", OUT_DIAG)

# quick sanity subset
print(
  sum_df %>%
    filter(str_detect(variable, "sigma_proc|sigma_v|sigma_u|sigma_init|sigma_obs_param|ERR\\[|k\\[")) %>%
    slice_head(n = 30)
)


