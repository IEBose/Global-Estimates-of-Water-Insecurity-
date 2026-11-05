################################################################################
#  B4_Country_level_Training.R — TRAIN AND VALIDATE COUNTRY-LEVEL MODELS
#
#  Purpose : Model country prevalence of moderate-to-high water insecurity in the
#            76 observed countries from national predictors only (11 national
#            indicators + model region). Used for the 33 countries with no GWP
#            individual data (Tier 3). Four models are compared:
#              Beta regression and quasi-binomial ("fractional logit") — no penalty
#              Elastic Net on logit(prevalence)                       <- final model
#              XGBoost on logit(prevalence)                           <- sensitivity
#
#  Validation — the SAME design as the individual-level model (B1):
#    * The same 6 hold-out countries (India + 5) are set aside as a test set.
#    * 5-fold leave-countries-out CV on the 70 training countries, using the
#      SAME country folds as B1, so Tables S7 and S8 compare like with like.
#    * Hyperparameters (ENET alpha/lambda, XGB rounds) are tuned on the 70
#      training countries using these folds; standardisation is re-estimated
#      within each fold.
#    * Out-of-sample predictions therefore exist for all 76 countries
#      (70 out-of-fold + 6 hold-out). Their residuals are used for the Tier 3
#      confidence intervals in B8.
#  Final model for prediction (B5/B6): refitted on ALL 76 countries with the
#  tuned hyperparameters (more data; validation figures are slightly conservative).
#
#  Depends : A1 (split), A4, B1 (country folds; out-of-sample individual predictions and fit for Tables S7/S8)
#  Inputs  : data/wb_imputed.rds, data/df_clean.rds, data/train_df_oof.rds
#  Outputs : models/country_*.rds, country_xgb_model.json, country_settings.rds
#            data/country_df_oos.rds (76 countries with out-of-sample predictions)
#            tables/B4_tableS07_model_comparison_LCO_CV.csv (Table S7: individual
#              and country models + model-region mean, LCO CV 70: R², MAE, RMSE,
#              calibration slope / intercept)
#            tables/B4_tableS08_holdout_metrics.csv (Table S8: the 6 hold-out
#              countries: MAE, RMSE, mean bias)
#            Full versions (all metrics and baselines; hold-out and all 76):
#              B4_model_comparison_all_sets.csv, B4_model_comparison_LCO_CV_full.csv,
#              B4_model_comparison_holdout_all76_full.csv, B4_country_cv_metrics.csv,
#              B4_baseline_oos_predictions.csv, B4_prediction_interval_coverage.csv
#            figures: B4_figS14_country_cv_all_models.png,
#                     B4_figS18_country_feature_importance.png,
#                     B4_figS19_shap_country_xgb.png
#
#  HOW TO READ THE METRICS (Tables S7, S8)
#    R²         squared correlation of observed vs out-of-sample predicted
#               prevalence (higher = more between-country variation captured;
#               does not detect systematic bias)
#    MAE        mean absolute error (proportion scale; 0.08 = 8 percentage points)
#    RMSE       like MAE but penalises large errors more
#    mean_bias  mean(predicted - observed): > 0 over-predicts on average
#    slope/intercept  from regressing observed on predicted; ideal = 1 and 0.
#               slope < 1: predictions too spread out; slope > 1: too compressed.
#    Best model = highest R², lowest MAE/RMSE, slope closest to 1.
#    The individual-level rows use the same 70 LCO CV / 6 hold-out countries.
################################################################################

source("00_setup.R")
start_timer("B4")
load_pkgs(c("betareg", "shapviz", "survey"))
select <- dplyr::select
filter <- dplyr::filter

wb <- load_data("wb_imputed")


# ══════════════════════════════════════════════════════════════════════════════
# 1. DATA: OBSERVED COUNTRY PREVALENCE + IMPUTED NATIONAL PREDICTORS
# ══════════════════════════════════════════════════════════════════════════════
# Observed prevalence (outcome) is calculated here directly from the cleaned
# survey data: survey-weighted proportion with IWISE >= 12 in each country,
# with PSU and strata (observed_country_prevalence() in 00_setup.R). This
# replaces the Stata-produced spreadsheet and is identical to Tier 1 in B8.

country_prev <- observed_country_prevalence(load_data("df_clean")) %>%
  transmute(Country_Code = country_code, Mod_high_WI = prevalence)

# merge() sorts by country code
country_df_merge <- merge(
  country_prev,
  wb %>% select(iso3, un_country_name, m49_region, m49_sub_region, m49_intermediate_region,
                region_final, all_of(c(ctry_vars_wb, ctry_aux_wb))),
  by.x = "Country_Code", by.y = "iso3", all.x = TRUE
) %>%
  filter(!is.na(Mod_high_WI)) %>%
  arrange(desc(Mod_high_WI))

# CHECK: 76 countries, no missing predictors after MICE (both should be as stated)
message("Observed countries: ", nrow(country_df_merge), " (expect 76) | missing predictor values: ",
        sum(is.na(country_df_merge[country_vars])))

# Beta regression needs values strictly between 0 and 1
country_df_merge <- country_df_merge %>%
  mutate(prevalence_beta = pmin(pmax(Mod_high_WI, 0.001), 0.999),
         logit_prev      = log(prevalence_beta / (1 - prevalence_beta)),
         region_final    = droplevels(factor(region_final, levels = region_final_levels)))

save_table(as_tibble(cor(country_df_merge[country_vars], use = "pairwise.complete.obs"), rownames = "variable"),
           "B4_country_predictor_correlations")

# ── Same split and folds as the individual-level model ──
split_lookup <- load_data("df_clean") %>% distinct(country_code, split) %>%
  mutate(country_code = as.character(country_code))
fold_lookup  <- load_data("train_df_oof") %>% distinct(country_code, fold) %>%
  mutate(country_code = as.character(country_code))

country_df_merge <- country_df_merge %>%
  left_join(split_lookup, by = c("Country_Code" = "country_code")) %>%
  left_join(fold_lookup,  by = c("Country_Code" = "country_code"))

# CHECK: 70 training countries in 5 folds + 6 hold-outs, nothing missing
message("Country split — train: ", sum(country_df_merge$split == "train", na.rm = TRUE),
        " | hold-out: ", sum(country_df_merge$split != "train", na.rm = TRUE),
        " | missing split: ", sum(is.na(country_df_merge$split)),
        " | training countries without a fold: ",
        sum(country_df_merge$split == "train" & is.na(country_df_merge$fold), na.rm = TRUE))

train_c   <- country_df_merge %>% filter(split == "train")
holdout_c <- country_df_merge %>% filter(split != "train")


# ══════════════════════════════════════════════════════════════════════════════
# 2. MODEL SPECIFICATIONS AND HELPERS
# ══════════════════════════════════════════════════════════════════════════════

beta_formula <- as.formula(paste("prevalence_beta ~", paste(c(country_vars, "region_final"), collapse = " + ")))
frac_formula <- as.formula(paste("Mod_high_WI ~",     paste(c(country_vars, "region_final"), collapse = " + ")))

# Recipe: dummies for region, drop zero-variance columns, standardise
make_country_recipe <- function(d) {
  recipe(logit_prev ~ ., data = d %>% select(logit_prev, all_of(country_vars), region_final)) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep()
}
bake_Xc <- function(rec, d) bake(rec, d) %>% select(-logit_prev) %>% as.matrix()

# XGBoost: shallow trees and strong regularisation because n is small
xgb_country_params <- list(objective = "reg:squarederror", eval_metric = "rmse",
                           eta = 0.05, max_depth = 3, subsample = 0.8,
                           colsample_bytree = 0.7, min_child_weight = 3, gamma = 2,
                           tree_method = "hist", nthread = 1,   # identical on every computer
                           seed = seed)

# Fit all 4 models on `train`, return predicted prevalence for `test`
fit_predict_country <- function(train, test, alpha, lambda, nrounds) {
  rec  <- make_country_recipe(train)
  X_tr <- bake_Xc(rec, train); X_te <- bake_Xc(rec, test)
  tibble(
    Country_Code = test$Country_Code,
    oos_beta = tryCatch(predict(betareg(beta_formula, data = train, link = "logit"),
                                newdata = test, type = "response"), error = function(e) NA_real_),
    oos_frac = tryCatch(predict(glm(frac_formula, data = train, family = quasibinomial(link = "logit")),
                                newdata = test, type = "response"), error = function(e) NA_real_),
    oos_enet = inv_logit(predict(glmnet(X_tr, train$logit_prev, alpha = alpha, lambda = lambda),
                                 X_te, type = "response")[, 1]),
    oos_xgb  = inv_logit(predict(xgb.train(params = xgb_country_params,
                                           data = xgb.DMatrix(X_tr, label = train$logit_prev),
                                           nrounds = nrounds, verbose = 0),
                                 xgb.DMatrix(X_te)))
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. TUNE HYPERPARAMETERS ON THE 70 TRAINING COUNTRIES (same 5 folds as B1)
# ══════════════════════════════════════════════════════════════════════════════

rec_train <- make_country_recipe(train_c)
X_train_c <- bake_Xc(rec_train, train_c)
y_train_c <- train_c$logit_prev

# Elastic Net: alpha grid x lambda path, minimum CV mean squared error
enet_country_cv <- lapply(c(0, 0.25, 0.5, 0.75, 1), function(a) {
  cv_fit <- cv.glmnet(X_train_c, y_train_c, alpha = a, foldid = train_c$fold, type.measure = "mse")
  tibble(alpha = a, lambda_min = cv_fit$lambda.min, mse_min = min(cv_fit$cvm))
}) %>% bind_rows()
save_table(enet_country_cv, "B4_enet_country_alpha_search")

best_alpha_c  <- enet_country_cv %>% slice_min(mse_min, with_ties = FALSE) %>% pull(alpha)
best_lambda_c <- enet_country_cv %>% slice_min(mse_min, with_ties = FALSE) %>% pull(lambda_min)

# XGBoost: number of rounds by the same folds, early stopping on RMSE
set.seed(seed)
xgb_country_cv <- xgb.cv(params = xgb_country_params,
                         data = xgb.DMatrix(X_train_c, label = y_train_c), nrounds = 500,
                         folds = split(seq_len(nrow(train_c)), train_c$fold),
                         early_stopping_rounds = 20, verbose = 0)
best_nrounds_c <- xgb_country_cv$evaluation_log$iter[which.min(xgb_country_cv$evaluation_log$test_rmse_mean)]

save_table(tibble(best_alpha_c, best_lambda_c, best_nrounds_c), "B4_country_hyperparameters")


# ══════════════════════════════════════════════════════════════════════════════
# 4. OUT-OF-SAMPLE PREDICTIONS FOR ALL 76 COUNTRIES
# ══════════════════════════════════════════════════════════════════════════════
# (a) 5-fold LCO CV: each training fold predicted by models fitted on the other 4
# (b) hold-out countries predicted by models fitted on all 70 training countries

set.seed(seed)
oof_preds <- map_dfr(sort(unique(train_c$fold)), function(f) {
  fit_predict_country(train_c %>% filter(fold != f), train_c %>% filter(fold == f),
                      best_alpha_c, best_lambda_c, best_nrounds_c)
}) %>% mutate(oos_source = "LCO CV (70)")

holdout_preds <- fit_predict_country(train_c, holdout_c, best_alpha_c, best_lambda_c, best_nrounds_c) %>%
  mutate(oos_source = "Hold-out (6)")

country_df_oos <- country_df_merge %>%
  left_join(bind_rows(oof_preds, holdout_preds), by = "Country_Code")

# CHECK: every country has an out-of-sample prediction (should be 0 for each model)
message("Missing out-of-sample predictions - beta: ", sum(is.na(country_df_oos$oos_beta)),
        " | frac: ", sum(is.na(country_df_oos$oos_frac)),
        " | enet: ", sum(is.na(country_df_oos$oos_enet)),
        " | xgb: ", sum(is.na(country_df_oos$oos_xgb)))


# ══════════════════════════════════════════════════════════════════════════════
# 5. VALIDATION METRICS (country-level models)
# ══════════════════════════════════════════════════════════════════════════════

evaluate_country <- function(observed, predicted, model_name, set_name) {
  ok  <- !is.na(predicted)
  obs <- observed[ok]; pred <- predicted[ok]
  cal <- coef(lm(obs ~ pred))
  tibble(model = model_name, validation = set_name, n_countries = sum(ok),
         R2 = cor(obs, pred)^2, MAE = mean(abs(pred - obs)),
         RMSE = sqrt(mean((pred - obs)^2)), mean_bias = mean(pred - obs),
         calibration_intercept = unname(cal[1]), calibration_slope = unname(cal[2]))
}

oos_cols <- c("Beta Regression" = "oos_beta", "Fractional Logit" = "oos_frac",
              "Elastic Net (logit)" = "oos_enet", "XGBoost" = "oos_xgb")
validation_sets <- list("LCO CV (70)" = country_df_oos$oos_source == "LCO CV (70)",
                        "Hold-out (6)" = country_df_oos$oos_source == "Hold-out (6)",
                        "All out-of-sample (76)" = rep(TRUE, nrow(country_df_oos)))

cv_results <- map_dfr(names(validation_sets), function(v) {
  d <- country_df_oos[validation_sets[[v]], ]
  map_dfr(names(oos_cols), ~ evaluate_country(d$Mod_high_WI, d[[oos_cols[[.x]]]], .x, v))
})
save_table(cv_results, "B4_country_cv_metrics")


# ══════════════════════════════════════════════════════════════════════════════
# 6. FINAL MODELS ON ALL 76 COUNTRIES (used for prediction in B5/B6)
# ══════════════════════════════════════════════════════════════════════════════

enet_country_rec <- make_country_recipe(country_df_merge)
X_country <- bake_Xc(enet_country_rec, country_df_merge)
y_country <- country_df_merge$logit_prev

beta_model <- betareg(beta_formula, data = country_df_merge, link = "logit")
frac_model <- glm(frac_formula, data = country_df_merge, family = quasibinomial(link = "logit"))
save_table(tibble(model = "Beta regression", pseudo_R2 = summary(beta_model)$pseudo.r.squared),
           "B4_beta_pseudo_R2")

enet_country_model <- glmnet(X_country, y_country, alpha = best_alpha_c, lambda = best_lambda_c)
# All coefficients, including those shrunk to 0 (not selected). The outcome is
# the logit of country prevalence, so odds_ratio = exp(coefficient) is the
# multiplicative change in the odds of country prevalence per 1 SD of the
# predictor (all predictors, including region dummies, are standardised).
enet_country_coefs_full <- coef(enet_country_model) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>% rename(coefficient = s0) %>%
  filter(feature != "(Intercept)") %>%
  mutate(feature_label = relabel(feature),
         odds_ratio    = exp(coefficient),
         selected      = coefficient != 0) %>%
  arrange(desc(abs(coefficient)))
save_table(enet_country_coefs_full, "B4_enet_country_coefficients")

enet_country_coefs <- enet_country_coefs_full %>% filter(selected)   # selected only (Figure S18)

set.seed(seed)
xgb_country_model <- xgb.train(params = xgb_country_params, data = xgb.DMatrix(X_country, label = y_country),
                               nrounds = best_nrounds_c, verbose = 0)
xgb_country_imp <- xgb.importance(model = xgb_country_model) %>% as_tibble()
save_table(xgb_country_imp, "B4_xgb_country_importance")


# ══════════════════════════════════════════════════════════════════════════════
# 7. FEATURE IMPORTANCE (Figures S18, S19)
# ══════════════════════════════════════════════════════════════════════════════

p_coef_country <- enet_country_coefs %>%
  mutate(feature_label = fct_reorder(feature_label, abs(coefficient)),
         direction = ifelse(coefficient > 0, "Higher prevalence", "Lower prevalence")) %>%
  ggplot(aes(x = coefficient, y = feature_label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Higher prevalence" = "darkred", "Lower prevalence" = "navy")) +
  labs(title = "a) Elastic Net coefficients (country-level)", x = "Coefficient (logit scale)",
       y = NULL, fill = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

p_xgb_country <- xgb_country_imp %>%
  head(15) %>%
  mutate(Feature = fct_reorder(relabel(Feature), Gain)) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(fill = "#4DAF4A") +
  labs(title = "b) XGBoost feature importance (country-level)", x = "Gain", y = NULL) +
  theme_minimal()

save_fig(p_coef_country | p_xgb_country, "B4_figS18_country_feature_importance.png",
         width = 14, height = 7, dpi = 300)

shap_country <- shapviz(xgb_country_model, X_pred = X_country)
colnames(shap_country$S) <- colnames(shap_country$X) <- relabel(colnames(shap_country$S))
save_fig(sv_importance(shap_country, kind = "beeswarm", max_display = 15) +
           labs(title = "SHAP Feature Importance - XGBoost (country-level)") + theme_minimal(),
         "B4_figS19_shap_country_xgb.png", width = 10, height = 7, dpi = 300)
save_fig(sv_importance(shap_country, kind = "bar", max_display = 15) +
           labs(title = "Mean |SHAP| - XGBoost (country-level)") + theme_minimal(),
         "B4_shap_bar_country_xgb.png", width = 8, height = 6, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 8. PLOTS: OUT-OF-SAMPLE PREDICTED vs OBSERVED (Figure S14) AND RESIDUALS
# ══════════════════════════════════════════════════════════════════════════════
# Circles = LCO CV countries; triangles = hold-out countries

oos_long <- country_df_oos %>%
  select(Country_Code, m49_region, Mod_high_WI, oos_source, all_of(unname(oos_cols))) %>%
  pivot_longer(all_of(unname(oos_cols)), names_to = "model", values_to = "predicted") %>%
  mutate(model = names(oos_cols)[match(model, oos_cols)],
         m49_region = factor(m49_region, levels = m49_region_levels))

ax_max <- max(country_df_oos$Mod_high_WI) * 1.1
p_oos <- ggplot(oos_long, aes(x = Mod_high_WI, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = m49_region, shape = oos_source), size = 3, alpha = 0.7) +
  geom_text(aes(label = Country_Code), size = 2.2, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~ model) +
  scale_colour_manual(values = m49_col) +
  coord_equal(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  labs(title = "Country-Level Predicted vs Observed Prevalence (out-of-sample)",
       x = "Observed Prevalence (weighted)", y = "Out-of-sample Predicted Prevalence",
       colour = "M49 region", shape = NULL) +
  theme_minimal() + theme(legend.position = "bottom", legend.box = "vertical")
save_fig(p_oos, "B4_figS14_country_cv_all_models.png", width = 14, height = 7, dpi = 300)

# Residuals for the best model (lowest MAE in LCO CV, comparable with Table S7)
best_model_name <- cv_results %>% filter(validation == "LCO CV (70)") %>%
  slice_min(MAE, with_ties = FALSE) %>% pull(model)
message("Best country-level model by LCO CV MAE: ", best_model_name)

p_resid <- country_df_oos %>%
  mutate(residual = .data[[oos_cols[[best_model_name]]]] - Mod_high_WI,
         m49_region = factor(m49_region, levels = m49_region_levels)) %>%
  ggplot(aes(x = Mod_high_WI, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = m49_region, shape = oos_source), size = 3, alpha = 0.7) +
  geom_text(aes(label = Country_Code), size = 2.5, vjust = -0.8, check_overlap = TRUE) +
  scale_colour_manual(values = m49_col) +
  labs(title = paste("Residuals -", best_model_name, "(country level, out-of-sample)"),
       subtitle = "Positive = overestimates, Negative = underestimates",
       x = "Observed Prevalence", y = "Residual (Predicted - Observed)",
       colour = "M49 region", shape = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
save_fig(p_resid, "B4_country_residuals_best_model.png", width = 10, height = 7, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 9. SAVE MODELS AND DATA
# ══════════════════════════════════════════════════════════════════════════════

save_table(country_df_oos %>%
             select(Country_Code, region_final, m49_region, split, fold, oos_source, Mod_high_WI,
                    all_of(unname(oos_cols))),
           "B4_country_out_of_sample_predictions")

save_model(enet_country_rec,   "country_recipe")
save_model(enet_country_model, "country_enet_model")
save_model(beta_model,         "country_beta_model")
save_model(frac_model,         "country_frac_model")
invisible(xgb.save(xgb_country_model, file.path(dir_models, "country_xgb_model.json")))
save_model(list(best_alpha_c = best_alpha_c, best_lambda_c = best_lambda_c,
                best_nrounds_c = best_nrounds_c, xgb_country_params = xgb_country_params,
                region_levels = levels(country_df_merge$region_final),
                train_colnames = colnames(X_country)),
           "country_settings")

save_data(country_df_oos, "country_df_oos")


# ══════════════════════════════════════════════════════════════════════════════
# 10. BASELINES, TABLES S7/S8 AND INTERVAL COVERAGE
# ══════════════════════════════════════════════════════════════════════════════

# ── (a) Baseline models: is the ML better than a simple average? ──
# Each test country is predicted by the mean observed prevalence of the
# training countries overall, in its model region, or in its income group,
# using exactly the same LCO CV folds and hold-outs as the models above.
income_lookup <- wb %>%
  transmute(Country_Code = iso3,
            income_group = coalesce(as.character(wb_income_fy27), "Unclassified"))
base_df <- country_df_oos %>%
  select(Country_Code, region_final, oos_source, fold, Mod_high_WI) %>%
  left_join(income_lookup, by = "Country_Code") %>%
  mutate(region_final = as.character(region_final))

baseline_predict <- function(train, test) {
  overall <- mean(train$Mod_high_WI)
  by_reg  <- train %>% group_by(region_final) %>% summarise(m = mean(Mod_high_WI), .groups = "drop")
  by_inc  <- train %>% group_by(income_group) %>% summarise(m = mean(Mod_high_WI), .groups = "drop")
  test %>%
    mutate(base_overall = overall,
           base_region  = coalesce(by_reg$m[match(region_final, by_reg$region_final)], overall),
           base_income  = coalesce(by_inc$m[match(income_group, by_inc$income_group)], overall))
}
base_train <- filter(base_df, oos_source == "LCO CV (70)")
baseline_oos <- bind_rows(
  map_dfr(sort(unique(base_train$fold)), function(f)
    baseline_predict(filter(base_train, fold != f), filter(base_train, fold == f))),
  baseline_predict(base_train, filter(base_df, oos_source == "Hold-out (6)"))
)
save_table(baseline_oos, "B4_baseline_oos_predictions")

base_cols <- c("Baseline: overall mean" = "base_overall",
               "Baseline: model-region mean" = "base_region",
               "Baseline: income-group mean" = "base_income")
base_sets <- list("LCO CV (70)" = baseline_oos$oos_source == "LCO CV (70)",
                  "Hold-out (6)" = baseline_oos$oos_source == "Hold-out (6)",
                  "All (76)" = rep(TRUE, nrow(baseline_oos)))
baseline_results <- map_dfr(names(base_sets), function(v) {
  d <- baseline_oos[base_sets[[v]], ]
  map_dfr(names(base_cols), function(m) {
    r <- suppressWarnings(evaluate_country(d$Mod_high_WI, d[[base_cols[[m]]]], m, v))
    # The overall mean is one value per fold, so R2 and calibration are not
    # meaningful for it (or for any baseline with a single predicted value)
    if (m == "Baseline: overall mean" || sd(d[[base_cols[[m]]]]) == 0) {
      r <- r %>% mutate(R2 = NA_real_, calibration_intercept = NA_real_, calibration_slope = NA_real_)
    }
    r
  })
}) %>% mutate(model_level = "Baseline")

# ── (b) Combined tables: LCO CV (70) -> Table S7; hold-out (6) -> Table S8 ──
ind_fit <- load_data("individual_country_fit") %>%
  rename(validation = set)
combined <- bind_rows(
  ind_fit,
  cv_results %>% mutate(model_level = "Country level",
                        validation = recode(validation, "All out-of-sample (76)" = "All (76)")),
  baseline_results
) %>%
  group_by(validation) %>%
  mutate(error_reduction_vs_region_mean =
           1 - MAE / MAE[model == "Baseline: model-region mean"][1]) %>%
  ungroup() %>%
  mutate(model_level = factor(model_level, levels = c("Individual level", "Country level", "Baseline"))) %>%
  arrange(validation, model_level) %>%
  select(validation, model_level, model, n_countries, n_respondents, R2, MAE, RMSE, mean_bias,
         calibration_slope, calibration_intercept, error_reduction_vs_region_mean)

fmt_tab <- function(d) d %>%
  mutate(across(c(R2, MAE, RMSE, mean_bias, calibration_slope, calibration_intercept,
                  error_reduction_vs_region_mean), ~ ifelse(is.na(.x), "", sprintf("%.3f", .x))),
         n_respondents = ifelse(is.na(n_respondents), "", format(n_respondents, big.mark = ","))) %>%
  rename(Validation = validation, `Model level` = model_level, Model = model,
         `N countries` = n_countries, `N respondents` = n_respondents, `Mean bias` = mean_bias,
         `Calibration slope` = calibration_slope, `Calibration intercept` = calibration_intercept,
         `Error reduction vs model-region mean` = error_reduction_vs_region_mean)
save_table(combined, "B4_model_comparison_all_sets")
save_table(fmt_tab(filter(combined, validation == "LCO CV (70)")), "B4_model_comparison_LCO_CV_full")
save_table(fmt_tab(filter(combined, validation != "LCO CV (70)")), "B4_model_comparison_holdout_all76_full")
print(filter(combined, validation == "LCO CV (70)") %>% select(model_level, model, R2, MAE), n = Inf)

# Tables S7 and S8 as in the paper: the 4 individual and 4 country models plus
# the model-region mean (reference for how much the models improve on a simple
# regional average). N countries / respondents go in the table footnote.
keep_baseline <- TRUE   # FALSE to drop the model-region mean row
f3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
paper_rows <- combined %>%
  filter(model_level != "Baseline" | (keep_baseline & model == "Baseline: model-region mean")) %>%
  mutate(model = recode(model, "Baseline: model-region mean" = "Model-region mean"))
tableS7 <- paper_rows %>% filter(validation == "LCO CV (70)") %>%
  transmute(`Model level` = model_level, Model = model, `R²` = f3(R2), MAE = f3(MAE), RMSE = f3(RMSE),
            `Calibration slope` = f3(calibration_slope), `Calibration intercept` = f3(calibration_intercept))
tableS8 <- paper_rows %>% filter(validation == "Hold-out (6)") %>%
  transmute(`Model level` = model_level, Model = model, MAE = f3(MAE), RMSE = f3(RMSE),
            `Mean bias` = f3(mean_bias))
readr::write_excel_csv(tableS7, file.path(dir_tables, "B4_tableS07_model_comparison_LCO_CV.csv"))
readr::write_excel_csv(tableS8, file.path(dir_tables, "B4_tableS08_holdout_metrics.csv"))
print(tableS7, width = Inf); print(tableS8, width = Inf)
message("Footnote numbers - LCO CV: ",
        paste(unique(paste0(paper_rows$model_level, " n_countries = ", paper_rows$n_countries,
                            ", n_respondents = ", paper_rows$n_respondents)[paper_rows$validation == "LCO CV (70)"]),
              collapse = " | "))

# ── (c) Coverage of 95% prediction intervals (observed countries) ──
# For each observed country: out-of-sample prediction + 2.5th / 97.5th
# percentiles of the out-of-sample residuals of the OTHER 75 countries (leave
# own residual out). Coverage = % of countries whose observed prevalence lies
# inside. About 95% (roughly 90-98% given 76 countries) means the
# residual-based model uncertainty used for predicted countries is about right.
# Two versions:
#   "logit"      residuals on the logit scale, added to logit(prediction) and
#                transformed back (as used for the B8 / B11 intervals)
#   "prevalence" residuals in percentage points, bounds clipped to [0, 1]
#                (for comparison)
# Checks the model-prediction component only.
coverage_one <- function(obs, pred, region, model_name) {
  res   <- obs - pred
  res_l <- logit_c(obs) - logit_c(pred)
  lo <- hi <- lo_l <- hi_l <- numeric(length(obs))
  for (i in seq_along(obs)) {
    q <- quantile(res[-i], c(0.025, 0.975))
    lo[i] <- max(pred[i] + q[1], 0); hi[i] <- min(pred[i] + q[2], 1)
    ql <- quantile(res_l[-i], c(0.025, 0.975))
    lo_l[i] <- plogis(logit_c(pred[i]) + ql[1]); hi_l[i] <- plogis(logit_c(pred[i]) + ql[2])
  }
  bind_rows(
    tibble(model = model_name, scale = "logit", region = region,
           inside = obs >= lo_l & obs <= hi_l, width = hi_l - lo_l, lower_at_0 = lo_l <= 0),
    tibble(model = model_name, scale = "prevalence", region = region,
           inside = obs >= lo & obs <= hi, width = hi - lo, lower_at_0 = lo <= 0))
}
ind_oos_cov <- load_data("country_pred_individual_oos")
cov_all <- bind_rows(
  coverage_one(ind_oos_cov$observed, ind_oos_cov$pred_enet, as.character(ind_oos_cov$m49_region),
               "Individual-level Elastic Net"),
  coverage_one(country_df_oos$Mod_high_WI, country_df_oos$oos_enet, as.character(country_df_oos$m49_region),
               "Country-level Elastic Net"))
coverage <- bind_rows(
  cov_all %>% group_by(model, scale) %>%
    summarise(region = "All", n_countries = n(), coverage_pct = 100 * mean(inside),
              mean_width_pp = 100 * mean(width), n_lower_at_0 = sum(lower_at_0), .groups = "drop"),
  cov_all %>% group_by(model, scale, region) %>%
    summarise(n_countries = n(), coverage_pct = 100 * mean(inside),
              mean_width_pp = 100 * mean(width), n_lower_at_0 = sum(lower_at_0), .groups = "drop")
) %>% arrange(model, desc(scale == "logit"), region != "All", region)
save_table(coverage, "B4_prediction_interval_coverage")
print(coverage, n = Inf)

stop_timer()
