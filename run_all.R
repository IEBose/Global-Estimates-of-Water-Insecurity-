################################################################################
#  run_all.R — RUN THE FULL PIPELINE IN ORDER
#
#  Each script reads what it needs from OUTPUT/ (written by earlier scripts),
#  so any script can also be run on its own once its inputs exist.
#  Run times for every script are appended to OUTPUT/logs/run_times.csv.
#
#  Run times from the October 2026 run on the full data (logs/run_times.csv):
#    A0  World Bank pull (API)             ~1-3 min   (only when refreshing data)
#    A1  Data cleaning                     ~2-5 min
#    A4  MICE imputation                   ~4 h (individual-level MICE dominates)
#    A5  Descriptive statistics            ~2-5 min (mostly the survey-package check)
#    B1  Individual-level training         ~45 min (RF + mixed model in 5-fold LCO CV; SHAP);
#                                          longer from Nov 2026 (RF and XGBoost single-threaded
#                                          so results are identical on every computer)
#    B2  Individual ENET predictions       ~1-2 min
#    B2b Imputation sensitivity (5 MICE imputations; no change to results) ~few min (new, Nov 2026)
#    B3  Individual XGB predictions        <1 min
#    B4  Country-level training (5-fold)   ~1-3 min
#    B5  Country ENET predictions          <1 min
#    B6  Country XGB predictions           <1 min
#    B7  Validation graphs                 ~1 min
#    B8  Global estimates + CIs            ~5-15 min (bootstraps; B8 intervals = sensitivity)
#    B9  Region sensitivity                <1 min
#    B10 Intervals, countries fixed        <1 min    (sensitivity)
#    B11 Model bootstrap UIs               ~20 min (MAIN group UIs; case + residual)
#    C1  Figure 2, S1-S3 (+ alternative maps) ~2-5 min  (first run downloads UN boundaries)
#    C2  Bivariate maps (alternative)      ~1-2 min
#    C3  Figure 4, S4-S6 (correlations)    <1 min
#    C4  Maps on UN boundaries (alternative) ~2-5 min
#    C5  Figures 1, 3, S7 (WHO map SOPs)   ~1-2 min  (first run downloads Natural Earth layers)
#    TOTAL (main pipeline, excl. A0)       ~5-5.5 hours (logged run, Oct 2026)
#
#  Optional, for interest only (NOT needed for the paper, NOT in the total):
#    A2  Original ML training (single imputation)   ~40 min
#    A3  Imputation sensitivity                     ~35 min (needs A2)
#    A4b MICE for all national variables (sensitivity; run on its own, writes to
#        OUTPUT_sensitivity_MICE_all_national/)     ~1-2 h
################################################################################

# Run from the scripts folder (all scripts call source("00_setup.R")).
# In RStudio: Session > Set Working Directory > To Source File Location.
source("00_setup.R")

run_A0       <- FALSE   # TRUE only when deliberately pulling a NEW World Bank snapshot
run_optional <- FALSE   # TRUE to also run the optional A2 and A3 scripts

pipeline <- c(
  "A1_Data_cleaning.R",
  "A4_MICE_Imputation.R",
  "A5_Descriptive_statistics.R",
  "B1_Individual_level_Training.R",
  "B2_Individual_ENET_predictions.R",
  "B2b_Imputation_sensitivity.R",
  "B3_Individual_XGB_predictions.R",
  "B4_Country_level_Training.R",
  "B5_Country_ENET_predictions_comparison.R",
  "B6_Country_XGB_predictions.R",
  "B7_Validation_Graphs.R",
  "B8_Global_Estimates_CIs.R",
  "B9_Sensitivity_region_estimates.R",
  "B10_Sensitivity_CI_without_replacement.R",
  "B11_Bootstrap_model_CIs.R",
  "C1_Global_Estimates_Figures_Maps.R",
  "C2_Bivariate_maps.R",
  "C3_Key_Indicator_correlations.R",
  "C4_UN_boundary_maps.R",
  "C5_WHO_style_maps.R"
)
optional_scripts <- c("A2_Original_ML_Training.R", "A3_Imputation_sensitivity.R")

if (run_A0) pipeline <- c("A0_External_data_merge.R", pipeline)

# Each script runs in a fresh environment so that no script silently relies on
# objects left in memory by another one.
run_script <- function(f) {
  message("\n================ ", f, " ================")
  source(f, local = new.env(parent = globalenv()), echo = FALSE)
}

pipeline_start <- Sys.time()
invisible(lapply(pipeline, run_script))
pipeline_mins <- as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))
message("\nMain pipeline finished in ", round(pipeline_mins, 1), " minutes")

if (run_optional) invisible(lapply(optional_scripts, run_script))

# Record package versions for reproducibility (each script also writes
# logs/sessionInfo_<script>.txt)
writeLines(capture.output(sessionInfo()), file.path(dir_logs, "sessionInfo.txt"))
