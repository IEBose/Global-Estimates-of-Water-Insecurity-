################################################################################
#  A4b_Sensitivity_MICE_all_national_vars.R — SENSITIVITY: MICE FOR ALL
#                                             INCOMPLETE NATIONAL VARIABLES
#
#  OPTIONAL — not part of run_all.R. Main outputs (OUTPUT/) are not changed.
#
#  Why    : In A4, MICE imputes only 4 national variables (corruption, WASH
#           mortality, climate vulnerability, log GNI) and uses the other 7 as
#           predictors. Four of those 7 have gaps (water stress 17 countries,
#           precipitation 13, basic water 1, under-5 mortality 1). MICE cannot
#           impute a country with a missing predictor, so 12 countries (mainly
#           Pacific islands and European microstates) got the model-region median
#           for some of the 4 variables, and the gaps in the 7 predictors were
#           also filled with the regional median.
#  Change : Here ALL incomplete national variables are imputed together by MICE
#           (pmm; all national variables + model region as predictors; same seed,
#           5 imputations x 10 iterations; mean of the 5 used, as in A4). The
#           regional median is kept only as a safety net (should not be needed).
#           Individual-level imputation is unchanged, so it is NOT re-run: the
#           individual variables are taken from the main A4 output (identical,
#           as A4 re-sets the seed before the individual-level MICE).
#  Then   : re-runs B1-B6 and B8 (optionally B7, B11) into a SEPARATE folder
#           OUTPUT_sensitivity_MICE_all_national/ and compares the estimates
#           with the main run.
#
#  Depends : main pipeline already run (A1, A4, B1-B8) — reads OUTPUT/
#  Outputs : OUTPUT_sensitivity_MICE_all_national/ (full set of B outputs) and,
#            in its tables/ folder:
#              A4b_national_values_comparison.csv  imputed values, main vs this
#              A4b_country_comparison.csv          country estimates, main vs this
#              A4b_group_comparison.csv            region / income / global
#              A4b_validation_comparison.csv       LCO CV metrics, main vs this
#            and, from A4c (run automatically at the end), the supplement table
#            and figure: A4b_table_group_comparison.csv,
#            A4b_table_countries_changed.csv, figures/A4b_fig_country_comparison.png/.pdf
#  Run time: ~1-2 h (mostly B1); add ~1-2 h if run_B11 = TRUE
#
#  Run from the scripts folder, in a fresh R session.
################################################################################

source("00_setup.R")    # main folders first
main_output_dir <- output_dir
main_dir_data   <- dir_data
main_dir_tables <- dir_tables

run_B7  <- FALSE   # validation figures (not needed for the comparison)
run_B11 <- FALSE   # bootstrap UIs (slow; set TRUE to compare intervals too)

sens_output_dir <- file.path(root_dir, "OUTPUT_sensitivity_MICE_all_national")

# Everything below writes to the sensitivity folder; the override is always
# removed at the end (even after an error) so later runs write to OUTPUT/ again.
assign("output_dir_override", sens_output_dir, envir = globalenv())
tryCatch({

source("00_setup.R")    # now points to the sensitivity folder
load_pkgs("mice")
select <- dplyr::select
filter <- dplyr::filter
t_start <- Sys.time()

# ══════════════════════════════════════════════════════════════════════════════
# 1. COPY THE UNCHANGED A1 OUTPUTS INTO THE SENSITIVITY FOLDER
# ══════════════════════════════════════════════════════════════════════════════
for (f in c("df_clean.rds", "new_data_clean.rds", "wb_clean.rds")) {
  file.copy(file.path(main_dir_data, f), file.path(dir_data, f), overwrite = TRUE)
}
a1_tables <- list.files(main_dir_tables, pattern = "^A1_", full.names = TRUE)
file.copy(a1_tables, dir_tables, overwrite = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 2. COUNTRY-LEVEL MICE FOR ALL INCOMPLETE NATIONAL VARIABLES
# ══════════════════════════════════════════════════════════════════════════════
wb       <- load_data("wb_clean")
df       <- readRDS(file.path(main_dir_data, "df_imputed.rds"))        # individual vars from main A4
new_data <- readRDS(file.path(main_dir_data, "new_data_imputed.rds"))

m_imp <- 5; maxit <- 10
nat_vars <- c(ctry_vars_wb, ctry_aux_wb)                           # all 11 national predictors
imp_vars <- nat_vars[sapply(nat_vars, function(v) any(is.na(wb[[v]])))]
message("National variables imputed by MICE here: ", paste(imp_vars, collapse = ", "))

ctry_mice_input <- wb %>% select(country_code = iso3, region_final, all_of(nat_vars))
ctry_mice_dat   <- ctry_mice_input %>% select(-country_code)

meth_c <- make.method(ctry_mice_dat)                 # pmm for numeric variables with gaps
meth_c[!names(meth_c) %in% imp_vars] <- ""
pred_c <- make.predictorMatrix(ctry_mice_dat)         # all other variables (incl. region) predict
print(meth_c)

set.seed(seed)
imp_ctry <- mice(ctry_mice_dat, m = m_imp, maxit = maxit,
                 method = meth_c, predictorMatrix = pred_c, printFlag = FALSE)
message("Logged events: ", if (is.null(imp_ctry$loggedEvents)) "none" else nrow(imp_ctry$loggedEvents))
if (!is.null(imp_ctry$loggedEvents)) print(imp_ctry$loggedEvents)

complete_ctry <- function(k) {
  out <- complete(imp_ctry, k) %>%
    mutate(country_code = ctry_mice_input$country_code, region_final = ctry_mice_input$region_final)
  n_left <- sum(is.na(out[imp_vars]))
  if (n_left > 0) message("Imputation ", k, ": ", n_left, " values not imputed by MICE - regional median used")
  out %>%
    group_by(region_final) %>%
    mutate(across(all_of(imp_vars), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup() %>%
    select(country_code, all_of(imp_vars))
}
ctry_completed <- lapply(seq_len(m_imp), complete_ctry) %>%
  bind_rows() %>%
  group_by(country_code) %>%
  summarise(across(all_of(imp_vars), mean), .groups = "drop")

# Put the imputed values into wb, df (model names) and new_data (World Bank names)
for (v in imp_vars) {
  wb[[v]]                     <- ctry_completed[[v]][match(wb$iso3, ctry_completed$country_code)]
  df[[ctry_name_map[[v]]]]    <- ctry_completed[[v]][match(as.character(df$country_code), ctry_completed$country_code)]
  new_data[[v]]               <- ctry_completed[[v]][match(as.character(new_data$country_code), ctry_completed$country_code)]
}

# CHECK: nothing missing (expect 0)
message("Remaining NAs - wb: ", sum(is.na(wb[nat_vars])),
        " | df: ", sum(is.na(df[unname(ctry_name_map)])),
        " | new_data: ", sum(is.na(new_data[nat_vars])))

save_data(df,       "df_imputed")
save_data(new_data, "new_data_imputed")
save_data(wb,       "wb_imputed")

# Imputed national values: main run vs this run (countries with any gap)
wb_main  <- readRDS(file.path(main_dir_data, "wb_imputed.rds"))
gap_iso3 <- wb$iso3[!complete.cases(load_data("wb_clean")[nat_vars])]
values_cmp <- bind_rows(
  wb_main %>% filter(iso3 %in% gap_iso3) %>% select(iso3, region_final, all_of(nat_vars)) %>% mutate(run = "Main (A4)"),
  wb      %>% filter(iso3 %in% gap_iso3) %>% select(iso3, region_final, all_of(nat_vars)) %>% mutate(run = "MICE all national (A4b)")
) %>% arrange(iso3, run)
save_table(values_cmp, "A4b_national_values_comparison")


# ══════════════════════════════════════════════════════════════════════════════
# 3. RE-RUN THE MODELS AND ESTIMATES INTO THE SENSITIVITY FOLDER
# ══════════════════════════════════════════════════════════════════════════════
# Clear this script's working objects so the re-run scripts only use what they load
rm(df, new_data, wb, wb_main, imp_ctry, ctry_completed, ctry_mice_input, ctry_mice_dat)

downstream <- c("B1_Individual_level_Training.R", "B2_Individual_ENET_predictions.R",
                "B3_Individual_XGB_predictions.R", "B4_Country_level_Training.R",
                "B5_Country_ENET_predictions_comparison.R", "B6_Country_XGB_predictions.R",
                if (run_B7) "B7_Validation_Graphs.R",
                "B8_Global_Estimates_CIs.R",
                if (run_B11) "B11_Bootstrap_model_CIs.R")
for (f in downstream) {
  message("\n================ ", f, " (sensitivity) ================")
  source(f, local = new.env(parent = globalenv()), echo = FALSE)
}
# B scripts re-assign select/filter; re-assert dplyr versions
select <- dplyr::select; filter <- dplyr::filter


# ══════════════════════════════════════════════════════════════════════════════
# 4. COMPARE WITH THE MAIN RUN
# ══════════════════════════════════════════════════════════════════════════════
main_tiers <- readRDS(file.path(main_dir_data, "all_tiers.rds"))
sens_tiers <- load_data("all_tiers")

country_cmp <- main_tiers %>%
  select(country_code, country_name, tier, m49_region, income_group, pop_15plus, prev_main = prev_best) %>%
  left_join(sens_tiers %>% select(country_code, prev_sens = prev_best), by = "country_code") %>%
  mutate(diff_pp = 100 * (prev_sens - prev_main),
         had_national_gap = country_code %in% gap_iso3) %>%
  arrange(desc(abs(diff_pp)))
save_table(country_cmp, "A4b_country_comparison")

group_prev <- function(d, prev_col, g) {
  bind_rows(
    d %>% group_by(group = as.character(.data[[g]])) %>%
      summarise(prev = weighted.mean(.data[[prev_col]], pop_15plus, na.rm = TRUE), .groups = "drop"),
    d %>% summarise(group = "Global", prev = weighted.mean(.data[[prev_col]], pop_15plus, na.rm = TRUE))
  ) %>% mutate(grouping = g)
}
group_cmp <- bind_rows(lapply(c("m49_region", "income_group"), function(g)
  group_prev(country_cmp, "prev_main", g) %>% rename(prev_main = prev) %>%
    left_join(group_prev(country_cmp, "prev_sens", g) %>% rename(prev_sens = prev), by = c("group", "grouping")))) %>%
  distinct(grouping, group, .keep_all = TRUE) %>%
  mutate(diff_pp = 100 * (prev_sens - prev_main)) %>%
  select(grouping, group, prev_main, prev_sens, diff_pp)
save_table(group_cmp, "A4b_group_comparison")

# LCO CV validation metrics (country-level models and individual model, country scale)
read_metrics <- function(dir, file) {
  f <- file.path(dir, file)
  if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE) else NULL
}
validation_cmp <- bind_rows(
  read_metrics(main_dir_tables, "B4_model_comparison_all_sets.csv") %>% mutate(run = "Main (A4)"),
  read_metrics(dir_tables,      "B4_model_comparison_all_sets.csv") %>% mutate(run = "MICE all national (A4b)")
)
if (nrow(validation_cmp) > 0) {
  validation_cmp <- validation_cmp %>%
    filter(validation == "LCO CV (70)", model_level != "Baseline") %>%
    select(run, model_level, model, R2, MAE, RMSE, calibration_slope, calibration_intercept) %>%
    arrange(model_level, model, run)
  save_table(validation_cmp, "A4b_validation_comparison")
}

message("\nSummary - countries whose estimate changed by more than 1 percentage point: ",
        sum(abs(country_cmp$diff_pp) > 1, na.rm = TRUE),
        " (of which with national data gaps: ",
        sum(abs(country_cmp$diff_pp) > 1 & country_cmp$had_national_gap, na.rm = TRUE), ")")
print(group_cmp %>% mutate(across(where(is.numeric), ~ round(.x, 4))), n = Inf)
print(country_cmp %>% filter(had_national_gap) %>%
        select(country_code, tier, prev_main, prev_sens, diff_pp) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf)
message("A4b finished in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
        " minutes. Outputs in: ", sens_output_dir)

}, finally = {
  if (exists("output_dir_override", envir = globalenv(), inherits = FALSE))
    rm("output_dir_override", envir = globalenv())
  message("Output folder reset to the main OUTPUT/ for any later scripts in this session.")
})

# Supplement table and figure (also runs on its own: source("A4c_MICE_sensitivity_table_figure.R"))
source("A4c_MICE_sensitivity_table_figure.R", local = new.env(parent = globalenv()))
