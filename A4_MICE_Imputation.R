################################################################################
#  A4_MICE_Imputation.R — MULTIPLE IMPUTATION (MICE) OF MISSING PREDICTORS
#
#  Purpose : Impute missing predictor values with MICE, at two levels:
#              1. Country level — all 194 countries (more data points), imputing
#                 Corruption (CPI), WASH mortality, ND-GAIN (cvi) and log GNI
#                 from 7 other national indicators + model region. These 7 are
#                 NOT imputed by MICE; 4 have gaps (water stress 17 countries,
#                 precipitation 13, basic water 1, under-5 mortality 1).
#              2. Individual level — df + new_data combined (161 countries),
#                 imputing 9 individual variables. The OUTCOME IS EXCLUDED to
#                 prevent leakage. Model region is used instead of country
#                 (161 country dummies caused convergence failures).
#            Countries MICE cannot resolve get regional medians/modes:
#              national: a country missing one of the 7 predictors cannot be
#                imputed by MICE (mice leaves it missing, with no logged event),
#                so its missing target values get the model-region median
#                (12 countries, mainly Pacific islands and European microstates);
#                gaps in the 7 predictors themselves are then also filled with
#                the model-region median (Section 3).
#              individual: Cuba, Djibouti, Guyana (no income data).
#            Sensitivity analysis with ALL incomplete national variables imputed
#            by MICE: A4b_Sensitivity_MICE_all_national_vars.R (optional).
#
#  Depends : A1 only (not A2 or A3)
#  Inputs  : OUTPUT/data/df_clean.rds, new_data_clean.rds, wb_clean.rds
#  Outputs : OUTPUT/data/df_imputed.rds, new_data_imputed.rds, wb_imputed.rds
#            OUTPUT/data/mice_completed_all.rds (all 5 completed datasets)
#            OUTPUT/tables/A4_mice_logged_events.csv, wb_country_vars_MICE_imputed.csv
#
#  Run time: the individual-level MICE dominates (~95 s per iteration;
#            5 imputations x 10 iterations ~ 80 min).
#
#  What to look out for:
#    * "Logged events" — MICE lists variables it dropped or could not impute
#      (e.g. collinear predictors). None, or only expected ones, is good.
#    * "Remaining NAs" after imputation must be 0 for every model variable.
#
#  Which imputation is used:
#    * National variables: the MEAN of the 5 imputations (after regional
#      fallback) — a single best estimate, less dependent on one random draw.
#    * Individual variables: imputation 1. Most are categorical, so they cannot
#      be averaged; robustness across all 5 imputations is checked in A3.
################################################################################

source("00_setup.R")
start_timer("A4")
load_pkgs("mice")
select <- dplyr::select   # mice is loaded after setup, so re-assert dplyr versions
filter <- dplyr::filter

df       <- load_data("df_clean")
new_data <- load_data("new_data_clean")
wb       <- load_data("wb_clean")

m_imp <- 5    # number of imputations
maxit <- 10   # iterations per imputation

# Individual-level variables to impute (same names in df and new_data)
ind_vars <- c("age", "INCOME_5", "inadincome", "education", "waterdissatisfied",
              "maritalstatus", "urban_final", "employment", "hhsize")
aux_vars <- c("female")   # complete, used only as a predictor

# CHECK: df must have the columns that receive the imputed national values
missing_in_df <- setdiff(unname(ctry_name_map), names(df))
if (length(missing_in_df) > 0) stop("Columns not found in df: ", paste(missing_in_df, collapse = ", "))

df       <- df %>% mutate(row_id = row_number())
new_data <- new_data %>% mutate(row_id = max(df$row_id) + row_number())   # no overlap with df


# ══════════════════════════════════════════════════════════════════════════════
# 1. COUNTRY-LEVEL MICE — ALL 194 COUNTRIES
# ══════════════════════════════════════════════════════════════════════════════

ctry_mice_input <- wb %>%
  select(country_code = iso3, region_final, all_of(c(ctry_vars_wb, ctry_aux_wb)))

ctry_mice_dat <- ctry_mice_input %>% select(-country_code)

meth_c <- make.method(ctry_mice_dat)
meth_c[!names(meth_c) %in% ctry_vars_wb] <- ""      # only impute the 4 target variables
pred_c <- make.predictorMatrix(ctry_mice_dat)        # all other columns (incl. region) predict

set.seed(seed)
imp_ctry <- mice(ctry_mice_dat, m = m_imp, maxit = maxit,
                 method = meth_c, predictorMatrix = pred_c, printFlag = FALSE)

if (!is.null(imp_ctry$loggedEvents)) {
  message("Country-level MICE logged events:"); print(imp_ctry$loggedEvents)
} else {
  message("Country-level MICE complete. No issues.")
}

# Completed national data for imputation k, with regional-median fallback for
# countries MICE could not impute
complete_ctry <- function(k) {
  complete(imp_ctry, k) %>%
    mutate(country_code = ctry_mice_input$country_code,
           region_final = ctry_mice_input$region_final) %>%
    group_by(region_final) %>%
    mutate(across(all_of(ctry_vars_wb), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup() %>%
    select(country_code, all_of(ctry_vars_wb))
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. INDIVIDUAL-LEVEL MICE — df + new_data COMBINED (161 COUNTRIES)
# ══════════════════════════════════════════════════════════════════════════════

ind_vars_both <- intersect(intersect(ind_vars, names(df)), names(new_data))
if (length(ind_vars_both) < length(ind_vars)) {
  warning("Not all ind_vars are in both datasets: missing ",
          paste(setdiff(ind_vars, ind_vars_both), collapse = ", "))
}

combined <- bind_rows(
  df       %>% select(all_of(c("row_id", "country_code", "region_final", ind_vars_both, aux_vars))),
  new_data %>% select(all_of(c("row_id", "country_code", "region_final", ind_vars_both, aux_vars)))
) %>%
  mutate(
    country_code      = as.factor(country_code),
    region_final      = factor(region_final, levels = region_final_levels),
    INCOME_5          = as.numeric(as.character(INCOME_5)),   # imputed as numeric (pmm)
    inadincome        = factor(inadincome),
    waterdissatisfied = factor(waterdissatisfied),
    education         = factor(education,     levels = lv_edu),
    maritalstatus     = factor(maritalstatus, levels = lv_mar),
    urban_final       = factor(urban_final,   levels = lv_urban),
    employment        = factor(employment,    levels = lv_emp),
    hhsize            = as.numeric(hhsize),
    female            = factor(female)
  )

# Variables that are completely missing in at least one country (these rely
# on the regional fallback in Section 3)
vars_regional <- combined %>%
  group_by(country_code) %>%
  summarise(across(all_of(ind_vars_both), ~ all(is.na(.)))) %>%
  pivot_longer(-country_code, names_to = "variable", values_to = "all_na") %>%
  filter(all_na)
save_table(vars_regional, "A4_vars_fully_missing_by_country")

ind_mice_dat <- combined %>%
  select(all_of(c("row_id", "country_code", "region_final", ind_vars_both, aux_vars)))

# Imputation methods by variable type:
#   pmm     — continuous (age, hhsize) and INCOME_5 (treated as numeric)
#   logreg  — binary (waterdissatisfied); inadincome is also binary but uses
#             polyreg, which for two categories is the same logistic model
#   polyreg — other categorical variables
meth_i <- make.method(ind_mice_dat)
meth_i[!names(meth_i) %in% ind_vars_both] <- ""
meth_i[c("inadincome", "education", "maritalstatus", "urban_final", "employment")] <- "polyreg"
meth_i["waterdissatisfied"] <- "logreg"
print(meth_i[ind_vars_both])   # CHECK: methods as listed above

pred_i <- make.predictorMatrix(ind_mice_dat)
pred_i[, "row_id"]       <- 0
pred_i[, "country_code"] <- 0   # 161 levels caused convergence failures
pred_i[, "region_final"] <- 1   # use the 9 model regions instead

set.seed(seed)
imp_ind <- mice(ind_mice_dat, m = m_imp, maxit = maxit,
                method = meth_i, predictorMatrix = pred_i, printFlag = TRUE)

if (!is.null(imp_ind$loggedEvents)) {
  message("Individual-level MICE logged events:"); print(imp_ind$loggedEvents)
} else {
  message("Individual-level MICE complete. No issues.")
}

save_table(bind_rows(
  if (!is.null(imp_ctry$loggedEvents)) mutate(as_tibble(imp_ctry$loggedEvents), level = "country"),
  if (!is.null(imp_ind$loggedEvents))  mutate(as_tibble(imp_ind$loggedEvents),  level = "individual")
), "A4_mice_logged_events")


# ══════════════════════════════════════════════════════════════════════════════
# 3. COMPLETED DATA, WITH REGIONAL FALLBACK
# ══════════════════════════════════════════════════════════════════════════════

get_mode <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

# Completed individual data for imputation k. Anything MICE could not impute
# (e.g. Cuba, Djibouti, Guyana — no income data) gets the regional median/mode.
complete_full <- function(k) {
  out <- complete(imp_ind, k) %>%
    group_by(region_final) %>%
    mutate(
      age      = ifelse(is.na(age), median(age, na.rm = TRUE), age),
      # round(): a median of quintiles can be x.5, which is not a valid level
      INCOME_5 = ifelse(is.na(INCOME_5), round(median(INCOME_5, na.rm = TRUE)), INCOME_5)
    ) %>%
    ungroup()

  for (v in c("inadincome", "education", "waterdissatisfied", "maritalstatus", "urban_final")) {
    if (any(is.na(out[[v]]))) {
      out <- out %>%
        group_by(region_final) %>%
        mutate(!!v := ifelse(is.na(.data[[v]]), get_mode(.data[[v]]), as.character(.data[[v]]))) %>%
        ungroup()
    }
  }
  out
}

completed_1 <- complete_full(1)
df_rows <- completed_1 %>% filter(row_id %in% df$row_id)
nd_rows <- completed_1 %>% filter(row_id %in% new_data$row_id)

# ── Individual variables ──
# df: imputed values go into the _imp columns (what the model uses)
# new_data: imputed values go into the original columns (B2 copies them to _imp)
for (v in ind_vars_both) {
  v_imp <- paste0(v, "_imp")
  if (v_imp %in% names(df)) {
    df[[v_imp]] <- df_rows[[v]][match(df$row_id, df_rows$row_id)]
  } else {
    df[[v]] <- df_rows[[v]][match(df$row_id, df_rows$row_id)]
  }
  new_data[[v]] <- nd_rows[[v]][match(new_data$row_id, nd_rows$row_id)]
}

# ── National variables: mean of the 5 imputations ──
ctry_completed <- lapply(seq_len(m_imp), complete_ctry) %>%
  bind_rows() %>%
  group_by(country_code) %>%
  summarise(across(all_of(ctry_vars_wb), mean), .groups = "drop")

# CHECK: spread of the 5 imputed values for imputed countries. A large range
# (relative to the variable's scale) means the value is poorly determined.
ctry_imp_spread <- lapply(seq_len(m_imp), complete_ctry) %>%
  bind_rows() %>%
  filter(country_code %in% ctry_mice_input$country_code[!complete.cases(ctry_mice_input[ctry_vars_wb])]) %>%
  pivot_longer(all_of(ctry_vars_wb), names_to = "variable") %>%
  group_by(country_code, variable) %>%
  summarise(mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>%
  filter(max > min)
save_table(ctry_imp_spread, "A4_country_imputation_spread")

for (wb_name in ctry_vars_wb) {
  idx_df <- match(as.character(df$country_code),       ctry_completed$country_code)
  idx_nd <- match(as.character(new_data$country_code), ctry_completed$country_code)
  idx_wb <- match(wb$iso3,                             ctry_completed$country_code)
  df[[ctry_name_map[[wb_name]]]] <- ctry_completed[[wb_name]][idx_df]
  new_data[[wb_name]]            <- ctry_completed[[wb_name]][idx_nd]
  wb[[wb_name]]                  <- ctry_completed[[wb_name]][idx_wb]
}

# Auxiliary national variables (not imputed by MICE): fill their gaps with the
# model-region median (water stress 17, precipitation 13, basic water 1,
# under-5 mortality 1 country in the current snapshot)
wb <- wb %>%
  group_by(region_final) %>%
  mutate(across(all_of(ctry_aux_wb), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  ungroup()

for (wb_name in ctry_aux_wb) {
  new_data[[wb_name]] <- wb[[wb_name]][match(as.character(new_data$country_code), wb$iso3)]
}

# ── Factor levels ──
df <- df %>% mutate(
  INCOME_5_imp          = factor(INCOME_5_imp, levels = 1:5),
  inadincome_imp        = factor(inadincome_imp),
  waterdissatisfied_imp = factor(waterdissatisfied_imp),
  education_imp         = factor(education_imp,     levels = lv_edu),
  maritalstatus_imp     = factor(maritalstatus_imp, levels = lv_mar),
  urban_final_imp       = factor(urban_final_imp,   levels = lv_urban),
  employment            = factor(employment,        levels = lv_emp)
)
new_data <- new_data %>% mutate(
  INCOME_5          = factor(INCOME_5, levels = 1:5),
  inadincome        = factor(inadincome),
  waterdissatisfied = factor(waterdissatisfied),
  education         = factor(education,     levels = lv_edu),
  maritalstatus     = factor(maritalstatus, levels = lv_mar),
  urban_final       = factor(urban_final,   levels = lv_urban),
  employment        = factor(employment,    levels = lv_emp)
)


# ══════════════════════════════════════════════════════════════════════════════
# 4. CHECK: NO MISSING VALUES LEFT (all counts should be 0)
# ══════════════════════════════════════════════════════════════════════════════

count_na <- function(d, vars, label) {
  d %>% summarise(across(any_of(vars), ~ sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_na") %>%
    mutate(dataset = label)
}
remaining_na <- bind_rows(
  count_na(df,       c(paste0(ind_vars, "_imp"), "employment", "hhsize", unname(ctry_name_map)), "df"),
  count_na(new_data, c(ind_vars, ctry_vars_wb, ctry_aux_wb), "new_data"),
  count_na(wb,       c(ctry_vars_wb, ctry_aux_wb), "wb")
) %>% filter(n_na > 0)

if (nrow(remaining_na) == 0) {
  message("Remaining NAs after imputation: none")
} else {
  warning("Remaining NAs after imputation — check before modelling:")
  print(remaining_na)
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. SAVE
# ══════════════════════════════════════════════════════════════════════════════

save_data(df,       "df_imputed")
save_data(new_data, "new_data_imputed")
save_data(wb,       "wb_imputed")

# All 5 completed datasets (for checking robustness across imputations)
save_data(list(individual = lapply(seq_len(m_imp), complete_full),
               country    = lapply(seq_len(m_imp), complete_ctry)),
          "mice_completed_all")

save_table(wb %>% select(iso3, region_final, all_of(c(ctry_vars_wb, ctry_aux_wb))),
           "wb_country_vars_MICE_imputed")

stop_timer()
