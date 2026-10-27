################################################################################
#  B2_Individual_ENET_predictions.R — TIER 2 PREDICTIONS (85 COUNTRIES, ENET)
#
#  Purpose : Apply the final individual-level Elastic Net (trained in B1) to the
#            85 countries with GWP data but no IWISE data. Each respondent gets
#            a predicted probability of water insecurity; these are averaged
#            within country using the survey weights to give the national
#            predicted prevalence (Tier 2 in the final estimates).
#            Also checks extrapolation: predictor values outside the range of
#            the training data (used for Table S9).
#
#  Depends : A4 (imputed new_data), B1 (recipe + model)
#  Inputs  : data/new_data_imputed.rds, models/ind_recipe.rds,
#            models/ind_enet_final.rds, models/ind_settings.rds
#  Outputs : data/new_data_pred.rds (respondent-level, with pred_prob)
#            tables/B2_tier2_country_prevalence_enet.csv
#            tables/B2_tier2_extrapolation_by_variable.csv
#            tables/B2_tier2_extrapolation_by_country.csv   (Table S9 input)
#
#  What to look out for:
#    * "Columns match training data: TRUE" — the prediction matrix must have
#      exactly the same columns, in the same order, as the training matrix.
#    * 85 countries, none of them in the training data.
#    * Extrapolation: countries with predictors outside the training range
#      (e.g. very high water stress in KWT/ARE/SAU, very high population
#      density in SGP/BHR/MDV/MLT) have less certain predictions.
################################################################################

source("00_setup.R")
start_timer("B2")

new_data   <- load_data("new_data_imputed")
rec        <- load_model("ind_recipe")
enet_final <- load_model("ind_enet_final")
settings   <- load_model("ind_settings")


# ══════════════════════════════════════════════════════════════════════════════
# 1. PUT new_data INTO THE SAME FORMAT AS THE TRAINING DATA
# ══════════════════════════════════════════════════════════════════════════════
# National variables: MICE-imputed World Bank names -> individual-model names.
# Individual variables: MICE-imputed originals -> _imp columns.
# The missingness flags (income_missing, urban_missing) were created in A1 from
# the ORIGINAL values and are NOT recalculated here (after imputation they would
# all be 0).

new_data <- new_data %>%
  mutate(
    log_GNI_imp    = log_gni_pc_ppp,
    Corruption_imp = cpi_score,
    Mort_WASH_imp  = wash_mort,
    CVI_imp        = cvi,
    Drinking_basic = basic_water,
    age_imp               = age,
    INCOME_5_imp          = factor(INCOME_5, levels = 1:5),
    inadincome_imp        = inadincome,
    education_imp         = education,
    waterdissatisfied_imp = waterdissatisfied,
    maritalstatus_imp     = maritalstatus,
    urban_final_imp       = urban_final,
    iwise12_imp           = NA_real_          # placeholder required by the recipe
  )

recipe_preds <- rec$var_info %>% filter(role == "predictor") %>% pull(variable)
missing_preds <- setdiff(recipe_preds, names(new_data))
if (length(missing_preds) > 0) stop("Predictors missing in new_data: ", paste(missing_preds, collapse = ", "))


# ══════════════════════════════════════════════════════════════════════════════
# 2. PREDICT
# ══════════════════════════════════════════════════════════════════════════════

X_new <- bake(rec, new_data) %>% select(-iwise12_imp) %>% as.matrix()

# CHECK: identical columns to the training matrix (should print TRUE)
message("Columns match training data: ", identical(colnames(X_new), settings$train_colnames))
if (!identical(colnames(X_new), settings$train_colnames)) stop("Column mismatch — check factor levels.")

new_data$pred_prob <- predict(enet_final, X_new, type = "response")[, 1]

# CHECK: no missing predictions; probabilities between 0 and 1
message("Predicted probabilities: n = ", nrow(new_data), ", NAs = ", sum(is.na(new_data$pred_prob)),
        ", range ", paste(round(range(new_data$pred_prob, na.rm = TRUE), 3), collapse = " to "))

# CHECK: 85 countries, none of which were training/observed countries
new_countries <- unique(as.character(new_data$country_code))
message("Tier 2 countries: ", length(new_countries), " | also in training data: ",
        sum(new_countries %in% settings$train_countries))


# ══════════════════════════════════════════════════════════════════════════════
# 3. EXTRAPOLATION CHECK (Table S9 input)
# ══════════════════════════════════════════════════════════════════════════════
# Compares each numeric predictor with its range in the 70 training countries.

ranges <- settings$numeric_ranges %>% rename(train_min = min, train_max = max)

extrap_by_variable <- ranges %>%
  rowwise() %>%
  mutate(new_min = min(new_data[[variable]], na.rm = TRUE),
         new_max = max(new_data[[variable]], na.rm = TRUE),
         extrapolates = new_min < train_min | new_max > train_max) %>%
  ungroup()
save_table(extrap_by_variable, "B2_tier2_extrapolation_by_variable")

# Per country: national predictors only (one value per country), listing which
# variables fall outside the training range
ctry_vars_check <- intersect(ctry_predictors_ind, ranges$variable)
extrap_by_country <- new_data %>%
  select(country_code, all_of(ctry_vars_check)) %>%
  distinct() %>%
  pivot_longer(-country_code, names_to = "variable", values_to = "value") %>%
  left_join(ranges, by = "variable") %>%
  mutate(outside = value < train_min | value > train_max) %>%
  group_by(country_code) %>%
  summarise(n_extrapolating = sum(outside, na.rm = TRUE),
            variables       = paste(variable[outside %in% TRUE], collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(n_extrapolating))
save_table(extrap_by_country, "B2_tier2_extrapolation_by_country")


# ══════════════════════════════════════════════════════════════════════════════
# 4. COUNTRY PREVALENCE (survey-weighted mean of predicted probabilities)
# ══════════════════════════════════════════════════════════════════════════════

country_prevalence <- new_data %>%
  group_by(country_code, region_final, m49_region, m49_sub_region, YEAR_WAVE, pre_2015) %>%
  summarise(n = n(),
            prev_enet_individual = weighted.mean(pred_prob, .data[[weight_var]]),
            .groups = "drop") %>%
  arrange(desc(prev_enet_individual))
save_table(country_prevalence, "B2_tier2_country_prevalence_enet")

# Diagnostic plot: predicted prevalence by country (pre-2015 data faded)
p_bar <- country_prevalence %>%
  mutate(country_code = fct_reorder(as.character(country_code), prev_enet_individual),
         caution = factor(ifelse(pre_2015 == 1, "Pre-2015 (caution)", "2015+")),
         m49_region = factor(m49_region, levels = m49_region_levels)) %>%
  ggplot(aes(x = prev_enet_individual, y = country_code, fill = m49_region, alpha = caution)) +
  geom_col() +
  scale_fill_manual(values = m49_col) +
  scale_alpha_manual(values = c("2015+" = 1, "Pre-2015 (caution)" = 0.5)) +
  labs(title = "Predicted Water Insecurity Prevalence - Tier 2 countries (Elastic Net)",
       x = "Predicted Prevalence", y = NULL, fill = "M49 region", alpha = "Data era") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6), legend.position = "bottom", legend.box = "vertical")
save_fig(p_bar, "B2_tier2_prevalence_barchart.png", width = 10, height = 14, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 5. SAVE
# ══════════════════════════════════════════════════════════════════════════════

save_data(new_data, "new_data_pred")

stop_timer()
