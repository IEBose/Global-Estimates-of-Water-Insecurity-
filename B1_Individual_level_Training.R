################################################################################
#  B1_Individual_level_Training.R — TRAIN AND VALIDATE INDIVIDUAL-LEVEL MODELS
#
#  Purpose : Using the 76 observed countries (individual + national predictors),
#            fit and compare four models of moderate-to-high water insecurity:
#              Elastic Net logistic regression (ENET)   <- final model
#              Random Forest (RF)
#              XGBoost (XGB)                            <- sensitivity model
#              Mixed-effects logistic regression (ME)   <- checks between-country variation (ICC)
#            Validation:
#              (a) Country hold-out: India + 5 random countries never used in training
#              (b) 5-fold leave-countries-out cross-validation (LCO CV) on the 70
#                  training countries (all respondents of a country in one fold);
#                  out-of-fold (OOF) predictions are aggregated to country prevalence.
#            Country-level results (Figures S9, S10; Figure 5 uses the 70 LCO
#            CV countries) use all 76 observed countries: 70 from LCO CV + the
#            6 hold-outs, as for the country-level model (B4, Tables S7/S8).
#            Sensitivity: excluding water dissatisfaction as a predictor.
#
#  NOTE    : The final ENET and XGB models are fitted on the 70 TRAINING countries
#            (the 6 hold-out countries are not used), and these models are the
#            ones used for prediction in B2/B3.
#
#  Depends : A4
#  Inputs  : OUTPUT/data/df_imputed.rds
#  Outputs : models/ind_recipe.rds, ind_enet_final.rds, ind_xgb_model.json,
#            ind_settings.rds, ind_enet_nowd.rds
#            data/train_df_oof.rds (training data + fold + OOF predictions)
#            data/country_pred_individual_oos.rds (76 countries, out-of-sample; B7)
#            data/holdout_country_pred.rds (6 hold-out countries; B5, B8)
#            tables/B1_*.csv, figures/B1_*.png  (see "OUTPUT MAP" below)
#
#  OUTPUT MAP (paper; numbering as in the manuscript):
#    Figure 6   ENET coefficients (bracketed)       B1_fig6_enet_coefficients_bracketed.png
#    Figure S8  ROC curves, hold-out countries      B1_figS08_roc_holdouts.png
#    Figure S9  Pred vs obs (ENET & XGB, 76)        B1_figS09_pred_vs_obs_enet_xgb.png
#    Figure S10 Pred vs obs all 4 models (76)       B1_figS10_pred_vs_obs_all_models.png
#    Figure S15 SHAP importance (XGB)               B1_figS15_shap_beeswarm_xgb.png (+ grouped)
#    Figure S16 ENET coefficients without WD        B1_figS16_enet_coef_noWD.png
#    Figure S17 SHAP without WD                     B1_figS17_shap_noWD.png
#    Table S5   Test (hold-out) metrics             B1_tableS05_holdout_metrics.csv
#    Table S6   LCO CV metrics (individual level),  B1_tableS06_lco_cv_metrics.csv
#               mean (SD) across folds + pooled     (footnote: ICC, B1_mixed_effects_ICC.csv)
#    Check tables (unformatted / extra detail):
#               B1_holdout_metrics.csv, B1_holdout_metrics_extended.csv,
#               B1_lco_cv_by_fold.csv, B1_lco_cv_by_fold_extended.csv,
#               B1_lco_cv_fold_summary.csv, B1_lco_cv_pooled.csv,
#               B1_lco_country_R2_MAE.csv, B1_country_level_fit_extended.csv
#               (country-scale fit of the individual models; combined with the
#               country models in B4 -> Tables S7/S8)
#    ICC (Table S6 footnote)                        B1_mixed_effects_ICC.csv (with and without predictors)
#
#  HOW TO READ THE METRICS
#    AUC-ROC   Probability the model ranks a random water-insecure respondent above
#              a random secure one. 0.5 = chance, 1 = perfect; ~0.7-0.8 acceptable,
#              >0.8 good. Measures discrimination, not whether prevalence is right.
#    AUC-PR    Area under precision-recall curve; compare with the outcome
#              prevalence (the value for a random model). Higher is better.
#    Brier     Mean squared error of predicted probabilities; lower is better.
#    R²        Here the SQUARED CORRELATION between observed and predicted country
#              prevalence (0-1). Higher = the model captures more of the
#              between-country variation. It does NOT detect systematic bias
#              (predictions consistently too high/low can still give a high R²),
#              so always read it together with MAE and calibration.
#    MAE       Mean absolute error of country prevalence, on the proportion scale
#              (0.07 = 7 percentage points). Lower is better.
#    ICC       Share of the remaining (latent) variance that lies between
#              countries after the predictors. Lower = predictors capture more
#              of the between-country differences.
################################################################################

source("00_setup.R")
start_timer("B1")
load_pkgs(c("ranger", "lme4", "shapviz"))
select <- dplyr::select
filter <- dplyr::filter

run_alpha_sensitivity <- TRUE   # slow (~5 x cv.glmnet); set FALSE to skip

df <- load_data("df_imputed")

train_df      <- df %>% filter(split == "train")
test_focal_df <- df %>% filter(split == "holdout_focal")
test_other_df <- df %>% filter(split == "holdout_other")

# CHECK: 70 / 1 / 5 countries, and every predictor is present
message("Train: ", n_distinct(train_df$country_code), " countries | focal hold-out: ",
        paste(unique(test_focal_df$country_code), collapse = ", "), " | other hold-outs: ",
        paste(unique(test_other_df$country_code), collapse = ", "))
missing_pred <- setdiff(all_predictors, names(train_df))
if (length(missing_pred) > 0) stop("Missing predictors: ", paste(missing_pred, collapse = ", "))

num_outcome <- function(d) d %>% mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))


# ══════════════════════════════════════════════════════════════════════════════
# 1. PRE-PROCESSING (recipe fitted on the training countries only)
# ══════════════════════════════════════════════════════════════════════════════
# Categorical predictors -> dummy variables (first level = reference);
# zero-variance columns removed; numeric predictors standardised (mean 0, SD 1),
# so ENET coefficients are per 1 SD change.

make_recipe <- function(train_data, predictors) {
  recipe(iwise12_imp ~ ., data = train_data %>% select(all_of(c("iwise12_imp", predictors))) %>% num_outcome()) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = train_data %>% num_outcome())
}
bake_X <- function(rec, d) bake(rec, d %>% num_outcome()) %>% select(-iwise12_imp) %>% as.matrix()

rec <- make_recipe(train_df, all_predictors)

X_train      <- bake_X(rec, train_df)
X_test_focal <- bake_X(rec, test_focal_df)
X_test_other <- bake_X(rec, test_other_df)
y_train      <- as.numeric(as.character(train_df$iwise12_imp))
y_test_focal <- as.numeric(as.character(test_focal_df$iwise12_imp))
y_test_other <- as.numeric(as.character(test_other_df$iwise12_imp))
w_train      <- train_df[[weight_var]]


# ══════════════════════════════════════════════════════════════════════════════
# 2. EVALUATION FUNCTION
# ══════════════════════════════════════════════════════════════════════════════
# Classification metrics use the threshold that maximises Youden's J
# (sensitivity + specificity) on the evaluated data. Metrics are unweighted.

evaluate_model <- function(y_true, y_prob, model_name, test_name) {
  roc_obj <- roc(y_true, y_prob, quiet = TRUE)
  best    <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  y_pred  <- ifelse(y_prob >= best$threshold, 1, 0)
  cm <- confusionMatrix(factor(y_pred, levels = c("0", "1")),
                        factor(y_true, levels = c("0", "1")), positive = "1")
  pr <- pr.curve(scores.class0 = y_prob[y_true == 1],
                 scores.class1 = y_prob[y_true == 0], curve = FALSE)
  tibble(model       = model_name,
         test_set    = test_name,
         AUC_ROC     = round(as.numeric(auc(roc_obj)), 4),
         AUC_PR      = round(pr$auc.integral, 4),
         Accuracy    = round(cm$overall["Accuracy"], 4),
         Sensitivity = round(cm$byClass["Sensitivity"], 4),
         Specificity = round(cm$byClass["Specificity"], 4),
         F1          = round(cm$byClass["F1"], 4),
         Brier       = round(mean((y_prob - y_true)^2), 4),
         Threshold   = round(best$threshold, 4))
}

# Additional metrics for the supplementary tables (Tables S5/S6, extended):
#   observed prevalence and mean predicted probability (survey-weighted; close
#   values = good calibration-in-the-large), 95% CI of AUC-ROC (DeLong), and the
#   scaled Brier score = 1 - Brier / Brier of predicting the (unweighted)
#   observed prevalence for everyone (0 = no better than that; 1 = perfect).
extra_metrics <- function(y_true, y_prob, w, n_countries) {
  ci <- as.numeric(ci.auc(roc(y_true, y_prob, quiet = TRUE), method = "delong"))
  brier     <- mean((y_prob - y_true)^2)
  brier_ref <- mean((mean(y_true) - y_true)^2)
  tibble(n_respondents = length(y_true), n_countries = n_countries,
         observed_prev = weighted.mean(y_true, w), mean_predicted = weighted.mean(y_prob, w),
         AUC_ROC_lower = ci[1], AUC_ROC_upper = ci[3],
         Brier_scaled  = 1 - brier / brier_ref)
}
fmt3 <- function(x) sprintf("%.3f", x)


# ══════════════════════════════════════════════════════════════════════════════
# 3. MODEL 1: ELASTIC NET (final model)
# ══════════════════════════════════════════════════════════════════════════════
# alpha = 0.5 (equal lasso/ridge mix). lambda chosen by 5-fold CV (random folds
# of respondents, maximising AUC) on the training countries.

best_alpha <- 0.5

set.seed(seed)   # makes the cv.glmnet folds reproducible
cv_fit <- cv.glmnet(X_train, y_train, family = "binomial", weights = w_train,
                    alpha = best_alpha, nfolds = 5, type.measure = "auc")
best_lambda <- cv_fit$lambda.min

enet_final <- glmnet(X_train, y_train, family = "binomial", weights = w_train,
                     alpha = best_alpha, lambda = best_lambda)

enet_prob_focal <- predict(enet_final, X_test_focal, type = "response")[, 1]
enet_prob_other <- predict(enet_final, X_test_other, type = "response")[, 1]

# Coefficients (log-odds per 1 SD for numeric predictors) and odds ratios
enet_coefs_full <- coef(enet_final) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>%
  rename(coefficient = s0) %>%
  filter(feature != "(Intercept)") %>%
  mutate(feature_label = relabel(feature),
         odds_ratio    = exp(coefficient)) %>%
  arrange(desc(abs(coefficient)))
save_table(enet_coefs_full, "B1_enet_coefficients_full")

enet_coefs <- enet_coefs_full %>% filter(coefficient != 0)   # selected features
message("ENET: lambda = ", signif(best_lambda, 3), "; non-zero features: ", nrow(enet_coefs))

# Odds-ratio plot (top 15)
p_or <- enet_coefs %>%
  head(15) %>%
  mutate(feature_label = fct_reorder(feature_label, abs(coefficient))) %>%
  ggplot(aes(x = odds_ratio, y = feature_label)) +
  geom_point(size = 3, colour = "#377EB8") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "Elastic Net - Odds Ratios (Top 15 Predictors)",
       subtitle = "Per 1-SD change in standardised predictors",
       x = "Odds Ratio", y = NULL) +
  theme_minimal()
save_fig(p_or, "B1_enet_odds_ratios.png", width = 9, height = 7, dpi = 300)

# Alpha sensitivity: does the lasso/ridge mix change performance?
# Look for: similar CV and hold-out AUC across alpha -> choice of 0.5 not critical.
if (run_alpha_sensitivity) {
  alpha_sensitivity <- lapply(c(0, 0.25, 0.5, 0.75, 1), function(a) {
    set.seed(seed)
    cv_tmp   <- cv.glmnet(X_train, y_train, family = "binomial", weights = w_train,
                          alpha = a, nfolds = 5, type.measure = "auc")
    enet_tmp <- glmnet(X_train, y_train, family = "binomial", weights = w_train,
                       alpha = a, lambda = cv_tmp$lambda.min)
    tibble(alpha     = a,
           cv_auc    = round(max(cv_tmp$cvm), 4),
           auc_india = round(as.numeric(auc(roc(y_test_focal, predict(enet_tmp, X_test_focal, type = "response")[, 1], quiet = TRUE))), 4),
           auc_other = round(as.numeric(auc(roc(y_test_other, predict(enet_tmp, X_test_other, type = "response")[, 1], quiet = TRUE))), 4))
  }) %>% bind_rows()
  save_table(alpha_sensitivity, "B1_sensitivity_alpha")
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. MODEL 2: RANDOM FOREST (comparison only)
# ══════════════════════════════════════════════════════════════════════════════

rf_settings <- list(num.trees = 1000, mtry = floor(sqrt(length(all_predictors))),
                    min.node.size = 10)

fit_rf <- function(train_data, w, importance = "none") {
  ranger(iwise12_imp ~ ., data = train_data %>% select(all_of(c("iwise12_imp", all_predictors))),
         case.weights = w, num.trees = rf_settings$num.trees, mtry = rf_settings$mtry,
         min.node.size = rf_settings$min.node.size, probability = TRUE,
         importance = importance, seed = seed,
         num.threads = 1)   # single thread: with a seed, results are then the same on every computer
}

rf_model <- fit_rf(train_df, w_train, importance = "impurity")
rf_prob_focal <- predict(rf_model, test_focal_df)$predictions[, "1"]
rf_prob_other <- predict(rf_model, test_other_df)$predictions[, "1"]

rf_importance <- tibble(feature = names(rf_model$variable.importance),
                        importance = rf_model$variable.importance) %>%
  arrange(desc(importance))
save_table(rf_importance, "B1_rf_feature_importance")


# ══════════════════════════════════════════════════════════════════════════════
# 5. MODEL 3: XGBOOST (sensitivity model)
# ══════════════════════════════════════════════════════════════════════════════
# Weights rescaled to mean 1 (min_child_weight is on the weight scale).
# Number of rounds chosen by 5-fold CV with early stopping (max AUC).

w_train_scaled <- w_train / mean(w_train)
dtrain <- xgb.DMatrix(data = X_train, label = y_train, weight = w_train_scaled)

# nthread = 1 and an explicit tree_method make the trees identical on every
# computer (with several threads, results depend on the number of CPU cores)
xgb_params <- list(objective = "binary:logistic", eval_metric = "auc",
                   eta = 0.05, max_depth = 6, subsample = 0.8,
                   colsample_bytree = 0.8, min_child_weight = 5, gamma = 1,
                   tree_method = "hist", nthread = 1, seed = seed)

set.seed(seed)
xgb_cv <- xgb.cv(params = xgb_params, data = dtrain, nrounds = 1000, nfold = 5,
                 early_stopping_rounds = 30, verbose = 0)
best_nrounds <- xgb_cv$evaluation_log$iter[which.max(xgb_cv$evaluation_log$test_auc_mean)]

xgb_model <- xgb.train(params = xgb_params, data = dtrain, nrounds = best_nrounds, verbose = 0)

xgb_prob_focal <- predict(xgb_model, xgb.DMatrix(data = X_test_focal))
xgb_prob_other <- predict(xgb_model, xgb.DMatrix(data = X_test_other))

save_table(xgb.importance(model = xgb_model) %>% as_tibble(), "B1_xgb_feature_importance")


# ══════════════════════════════════════════════════════════════════════════════
# 6. MODEL 4: MIXED-EFFECTS LOGISTIC REGRESSION (random intercept for country)
# ══════════════════════════════════════════════════════════════════════════════
# Reduced predictor set for convergence. Predictions for new countries use the
# fixed effects only (random effect = 0). The ICC indicates how much
# between-country variation remains after the predictors.

me_formula <- outcome_num ~ waterdissatisfied_imp + inadincome_imp + age_imp +
  hhsize + female + education_imp + INCOME_5_imp + urban_final_imp +
  Drinking_basic + u5_mort + CVI_imp + log_GNI_imp + log_pop_density +
  region_final + (1 | country_code)

fit_me <- function(train_data) {
  d <- train_data %>%
    mutate(outcome_num = as.numeric(as.character(iwise12_imp)),
           wt_scaled   = .data[[weight_var]] / mean(.data[[weight_var]]))
  suppressWarnings(glmer(me_formula, data = d, family = binomial(link = "logit"),
                         weights = wt_scaled,   # column of d
                         control = glmerControl(optimizer = "nloptwrap",
                                                optCtrl = list(maxeval = 20000),
                                                calc.derivs = FALSE),
                         nAGQ = 0))
}

me_model <- fit_me(train_df)

country_var <- as.data.frame(VarCorr(me_model))$vcov[1]
icc <- country_var / (country_var + (pi^2 / 3))   # latent-scale ICC for logistic models

# Same model without predictors (random intercept only): how much of the
# between-country variance the predictors explain (Table S6 footnote).
null_me <- local({
  d <- train_df %>%
    mutate(outcome_num = as.numeric(as.character(iwise12_imp)),
           wt_scaled   = .data[[weight_var]] / mean(.data[[weight_var]]))
  suppressWarnings(glmer(outcome_num ~ 1 + (1 | country_code), data = d,
                         family = binomial(link = "logit"), weights = wt_scaled, nAGQ = 0))
})
country_var_null <- as.data.frame(VarCorr(null_me))$vcov[1]
icc_null <- country_var_null / (country_var_null + (pi^2 / 3))
save_table(tibble(country_variance = country_var, ICC = icc,
                  country_variance_null = country_var_null, ICC_null = icc_null,
                  between_country_variance_explained = 1 - country_var / country_var_null),
           "B1_mixed_effects_ICC")
message("ICC: ", round(100 * icc_null, 1), "% without predictors, ", round(100 * icc, 1),
        "% with predictors (", round(100 * (1 - country_var / country_var_null)), "% of between-country variance explained)")

me_prob_focal <- predict(me_model, newdata = test_focal_df, type = "response", re.form = NA, allow.new.levels = TRUE)
me_prob_other <- predict(me_model, newdata = test_other_df, type = "response", re.form = NA, allow.new.levels = TRUE)

me_sum <- summary(me_model)$coefficients
save_table(tibble(feature = rownames(me_sum), estimate = me_sum[, "Estimate"],
                  std_error = me_sum[, "Std. Error"], z_value = me_sum[, "z value"],
                  p_value = me_sum[, "Pr(>|z|)"], odds_ratio = exp(me_sum[, "Estimate"])),
           "B1_mixed_effects_coefficients")


# ══════════════════════════════════════════════════════════════════════════════
# 7. HOLD-OUT RESULTS (Table S5) AND ROC CURVES (Figure S8)
# ══════════════════════════════════════════════════════════════════════════════

holdout_results <- bind_rows(
  evaluate_model(y_test_focal, enet_prob_focal, "Elastic Net",      "Focal (India)"),
  evaluate_model(y_test_other, enet_prob_other, "Elastic Net",      "Other countries"),
  evaluate_model(y_test_focal, rf_prob_focal,   "Random Forest",    "Focal (India)"),
  evaluate_model(y_test_other, rf_prob_other,   "Random Forest",    "Other countries"),
  evaluate_model(y_test_focal, xgb_prob_focal,  "XGBoost",          "Focal (India)"),
  evaluate_model(y_test_other, xgb_prob_other,  "XGBoost",          "Other countries"),
  evaluate_model(y_test_focal, me_prob_focal,   "Mixed-Effects LR", "Focal (India)"),
  evaluate_model(y_test_other, me_prob_other,   "Mixed-Effects LR", "Other countries")
)
save_table(holdout_results, "B1_holdout_metrics")

# Extended hold-out table (all columns above kept, extra columns added) -> Table S5
w_focal <- test_focal_df[[weight_var]]; w_other <- test_other_df[[weight_var]]
nc_focal <- n_distinct(test_focal_df$country_code); nc_other <- n_distinct(test_other_df$country_code)
holdout_probs <- list(
  "Elastic Net"      = list(enet_prob_focal, enet_prob_other),
  "Random Forest"    = list(rf_prob_focal,   rf_prob_other),
  "XGBoost"          = list(xgb_prob_focal,  xgb_prob_other),
  "Mixed-Effects LR" = list(me_prob_focal,   me_prob_other))
holdout_extra <- bind_rows(lapply(names(holdout_probs), function(m) bind_rows(
  extra_metrics(y_test_focal, holdout_probs[[m]][[1]], w_focal, nc_focal) %>% mutate(model = m, test_set = "Focal (India)"),
  extra_metrics(y_test_other, holdout_probs[[m]][[2]], w_other, nc_other) %>% mutate(model = m, test_set = "Other countries"))))
tableS3_ext <- holdout_results %>%
  left_join(holdout_extra, by = c("model", "test_set")) %>%
  select(model, test_set, n_respondents, n_countries, observed_prev, mean_predicted,
         AUC_ROC, AUC_ROC_lower, AUC_ROC_upper, AUC_PR, Brier, Brier_scaled,
         Threshold, Accuracy, Sensitivity, Specificity, F1)
save_table(tableS3_ext, "B1_holdout_metrics_extended")
save_table(tableS3_ext %>%
             transmute(Model = model, `Test set` = test_set, `N respondents` = n_respondents,
                       `N countries` = n_countries, `Observed prevalence` = fmt3(observed_prev),
                       `Mean predicted` = fmt3(mean_predicted),
                       `AUC-ROC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", AUC_ROC, AUC_ROC_lower, AUC_ROC_upper),
                       `AUC-PR` = fmt3(AUC_PR), Brier = fmt3(Brier), `Scaled Brier` = fmt3(Brier_scaled),
                       Threshold = fmt3(Threshold), Accuracy = fmt3(Accuracy), Sensitivity = fmt3(Sensitivity),
                       Specificity = fmt3(Specificity), F1 = fmt3(F1)),
           "B1_tableS05_holdout_metrics")

# Country prevalence predicted for the 6 hold-out countries (used in Tables S7/S8,
# S10, Figures S9, S10, and for the residuals in B8/B11)
holdout_country_pred <- bind_rows(test_focal_df, test_other_df) %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp)),
         enet_prob = c(enet_prob_focal, enet_prob_other),
         xgb_prob  = c(xgb_prob_focal,  xgb_prob_other),
         rf_prob   = c(rf_prob_focal,   rf_prob_other),
         me_prob   = c(me_prob_focal,   me_prob_other)) %>%
  group_by(country_code, region_final, m49_region) %>%
  summarise(n                    = n(),
            observed             = weighted.mean(outcome_num, .data[[weight_var]]),
            prev_enet_individual = weighted.mean(enet_prob, .data[[weight_var]]),
            prev_xgb_individual  = weighted.mean(xgb_prob,  .data[[weight_var]]),
            prev_rf_individual   = weighted.mean(rf_prob,   .data[[weight_var]]),
            prev_me_individual   = weighted.mean(me_prob,   .data[[weight_var]]),
            .groups = "drop") %>%
  mutate(country_code = as.character(country_code), source = "hold-out prediction")
save_data(holdout_country_pred, "holdout_country_pred")

plot_roc_comparison <- function(y_true, prob_list, model_names, title) {
  roc_list <- lapply(prob_list, function(p) roc(y_true, p, quiet = TRUE))
  cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  plot(roc_list[[1]], col = cols[1], lwd = 2, main = title)
  for (i in 2:length(roc_list)) lines(roc_list[[i]], col = cols[i], lwd = 2)
  legend("bottomright",
         legend = paste0(model_names, " (AUC=", sapply(roc_list, function(r) round(auc(r), 3)), ")"),
         col = cols[seq_along(roc_list)], lwd = 2, cex = 0.9)
}
model_names <- c("Elastic Net", "Random Forest", "XGBoost", "Mixed-Effects LR")

png(file.path(dir_figures, "B1_figS08_roc_holdouts.png"), width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
plot_roc_comparison(y_test_focal, list(enet_prob_focal, rf_prob_focal, xgb_prob_focal, me_prob_focal),
                    model_names, "ROC - India Holdout")
plot_roc_comparison(y_test_other, list(enet_prob_other, rf_prob_other, xgb_prob_other, me_prob_other),
                    model_names, "ROC - Other Countries Holdout")
invisible(dev.off())


# ══════════════════════════════════════════════════════════════════════════════
# 8. LEAVE-COUNTRIES-OUT CROSS-VALIDATION (5 folds, 70 training countries)
# ══════════════════════════════════════════════════════════════════════════════
# Countries (not respondents) are assigned to folds, so each fold's countries
# are entirely unseen when predicted. Hyperparameters (alpha, lambda, nrounds)
# are those chosen above on all training data (not re-tuned within folds).
# One loop stores the OOF predictions; fold metrics are computed from them.

k_folds <- 5
set.seed(seed)
country_folds <- tibble(country_code = unique(as.character(train_df$country_code))) %>%
  mutate(fold = sample(rep(1:k_folds, length.out = n())))

train_df <- train_df %>%
  mutate(country_code = as.character(country_code)) %>%
  select(-any_of("fold")) %>%
  left_join(country_folds, by = "country_code") %>%
  mutate(country_code = factor(country_code),
         enet_oof_prob = NA_real_, xgb_oof_prob = NA_real_,
         rf_oof_prob   = NA_real_, me_oof_prob  = NA_real_)

for (f in 1:k_folds) {
  message("  LCO CV fold ", f, " of ", k_folds)
  idx_test  <- which(train_df$fold == f)
  idx_train <- which(train_df$fold != f)
  cv_train  <- train_df[idx_train, ]
  cv_test   <- train_df[idx_test, ]

  cv_rec  <- make_recipe(cv_train, all_predictors)
  X_cv_tr <- bake_X(cv_rec, cv_train)
  X_cv_te <- bake_X(cv_rec, cv_test)
  y_cv_tr <- as.numeric(as.character(cv_train$iwise12_imp))
  w_cv_tr <- cv_train[[weight_var]]

  # Elastic Net
  m_enet <- glmnet(X_cv_tr, y_cv_tr, family = "binomial", weights = w_cv_tr,
                   alpha = best_alpha, lambda = best_lambda)
  train_df$enet_oof_prob[idx_test] <- predict(m_enet, X_cv_te, type = "response")[, 1]

  # XGBoost
  m_xgb <- xgb.train(params = xgb_params,
                     data = xgb.DMatrix(X_cv_tr, label = y_cv_tr, weight = w_cv_tr / mean(w_cv_tr)),
                     nrounds = best_nrounds, verbose = 0)
  train_df$xgb_oof_prob[idx_test] <- predict(m_xgb, xgb.DMatrix(X_cv_te))

  # Random Forest
  m_rf <- fit_rf(cv_train, w_cv_tr)
  train_df$rf_oof_prob[idx_test] <- predict(m_rf, cv_test)$predictions[, "1"]

  # Mixed effects (fixed effects only for unseen countries)
  m_me <- tryCatch(fit_me(cv_train), error = function(e) {
    message("    Mixed-effects failed in fold ", f, ": ", e$message); NULL })
  if (!is.null(m_me)) {
    train_df$me_oof_prob[idx_test] <- predict(m_me, newdata = cv_test, type = "response",
                                              re.form = NA, allow.new.levels = TRUE)
  }
}

# CHECK: every training respondent has an OOF prediction (should be 0 NAs;
# ME may be NA if a fold failed to converge)
message("OOF NAs - enet: ", sum(is.na(train_df$enet_oof_prob)), " | xgb: ", sum(is.na(train_df$xgb_oof_prob)),
        " | rf: ", sum(is.na(train_df$rf_oof_prob)), " | me: ", sum(is.na(train_df$me_oof_prob)))

# ── Individual-level metrics by fold, then mean and SD across folds (Table S6) ──
oof_cols <- c("Elastic Net" = "enet_oof_prob", "XGBoost" = "xgb_oof_prob",
              "Random Forest" = "rf_oof_prob", "Mixed-Effects LR" = "me_oof_prob")
lco_results <- map_dfr(1:k_folds, function(f) {
  d <- train_df %>% filter(fold == f)
  y <- as.numeric(as.character(d$iwise12_imp))
  map_dfr(names(oof_cols), function(m) {
    p <- d[[oof_cols[[m]]]]
    if (all(is.na(p))) return(NULL)
    evaluate_model(y, p, m, paste0("Fold_", f))
  })
})
save_table(lco_results, "B1_lco_cv_by_fold")

lco_summary <- lco_results %>%
  group_by(model) %>%
  summarise(across(c(AUC_ROC, AUC_PR, Accuracy, Sensitivity, Specificity, F1, Brier),
                   list(mean = mean, sd = sd)), .groups = "drop")
save_table(lco_summary, "B1_lco_cv_fold_summary")

# Extended LCO CV table (Table S6): extra metrics by fold, mean (SD) across folds, plus a
# "pooled" row computed once on all out-of-fold predictions (70 countries).
lco_extra <- map_dfr(1:k_folds, function(f) {
  d <- train_df %>% filter(fold == f)
  y <- as.numeric(as.character(d$iwise12_imp))
  map_dfr(names(oof_cols), function(m) {
    p <- d[[oof_cols[[m]]]]
    if (all(is.na(p))) return(NULL)
    ok <- !is.na(p)
    extra_metrics(y[ok], p[ok], d[[weight_var]][ok], n_distinct(d$country_code[ok])) %>%
      mutate(model = m, test_set = paste0("Fold_", f))
  })
})
lco_by_fold_ext <- lco_results %>% left_join(lco_extra, by = c("model", "test_set"))
save_table(lco_by_fold_ext, "B1_lco_cv_by_fold_extended")

metric_cols <- c("observed_prev", "mean_predicted", "AUC_ROC", "AUC_PR", "Brier", "Brier_scaled",
                 "Threshold", "Accuracy", "Sensitivity", "Specificity", "F1")
mean_sd <- function(x) sprintf("%.3f (%.3f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
lco_fold_fmt <- lco_by_fold_ext %>%
  group_by(model) %>%
  summarise(summary = "Mean (SD) across folds", across(all_of(metric_cols), mean_sd), .groups = "drop")

y_all <- as.numeric(as.character(train_df$iwise12_imp))
lco_pooled <- map_dfr(names(oof_cols), function(m) {
  p <- train_df[[oof_cols[[m]]]]; ok <- !is.na(p)
  if (!any(ok)) return(NULL)
  evaluate_model(y_all[ok], p[ok], m, "Pooled") %>%
    left_join(extra_metrics(y_all[ok], p[ok], train_df[[weight_var]][ok],
                            n_distinct(train_df$country_code[ok])) %>% mutate(model = m),
              by = "model")
})
save_table(lco_pooled, "B1_lco_cv_pooled")
lco_pooled_fmt <- lco_pooled %>%
  transmute(model, summary = "Pooled (all out-of-fold predictions)",
            across(all_of(metric_cols), fmt3))

tableS4_fmt <- bind_rows(lco_fold_fmt, lco_pooled_fmt) %>%
  mutate(model = factor(model, levels = names(oof_cols))) %>%
  arrange(model, summary != "Mean (SD) across folds") %>%
  rename(Model = model, Summary = summary, `Observed prevalence` = observed_prev,
         `Mean predicted` = mean_predicted, `AUC-ROC` = AUC_ROC, `AUC-PR` = AUC_PR,
         `Scaled Brier` = Brier_scaled)
save_table(tableS4_fmt, "B1_tableS06_lco_cv_metrics")


# ══════════════════════════════════════════════════════════════════════════════
# 9. COUNTRY-LEVEL PREDICTED vs OBSERVED (out-of-sample, 76 countries) — Tables S7/S8 input
# ══════════════════════════════════════════════════════════════════════════════
# Individual probabilities are averaged within country using the survey weights
# and compared with the weighted observed prevalence. All 76 observed countries
# are shown, as for the country-level model (B4): 70 predicted in 5-fold LCO CV
# + the 6 hold-out countries predicted by the final model.

country_pred_train <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final, m49_region) %>%
  summarise(n         = n(),
            observed  = weighted.mean(outcome_num,   .data[[weight_var]]),
            pred_enet = weighted.mean(enet_oof_prob, .data[[weight_var]]),
            pred_xgb  = weighted.mean(xgb_oof_prob,  .data[[weight_var]]),
            pred_rf   = weighted.mean(rf_oof_prob,   .data[[weight_var]]),
            pred_me   = weighted.mean(me_oof_prob,   .data[[weight_var]], na.rm = TRUE),
            .groups   = "drop") %>%
  mutate(country_code = as.character(country_code), oos_source = "LCO CV (70)")

country_pred_all <- bind_rows(
  country_pred_train,
  holdout_country_pred %>%
    transmute(country_code, region_final, m49_region, n, observed,
              pred_enet = prev_enet_individual, pred_xgb = prev_xgb_individual,
              pred_rf = prev_rf_individual, pred_me = prev_me_individual,
              oos_source = "Hold-out (6)")
) %>%
  mutate(across(starts_with("pred_"), ~ .x - observed, .names = "diff_{.col}"))
save_table(country_pred_all, "B1_country_pred_vs_obs_oof")
save_data(country_pred_all, "country_pred_individual_oos")   # used by B7 (Figure 5)

# R² (squared correlation) and MAE for LCO CV (70), hold-out (6) and all (76),
# as for the country-level model in B4.
# Highest R² / lowest MAE = best fit to between-country variation.
# (R² on 6 hold-out countries is unstable — read it alongside MAE.)
validation_sets <- list("LCO CV (70)" = country_pred_all$oos_source == "LCO CV (70)",
                        "Hold-out (6)" = country_pred_all$oos_source == "Hold-out (6)",
                        "All (76)"     = rep(TRUE, nrow(country_pred_all)))
country_fit <- expand_grid(set = names(validation_sets),
                           tibble(model = names(oof_cols), col = c("pred_enet", "pred_xgb", "pred_rf", "pred_me"))) %>%
  rowwise() %>%
  mutate(n_countries = sum(validation_sets[[set]] & !is.na(country_pred_all[[col]])),
         R2  = cor(country_pred_all$observed[validation_sets[[set]]],
                   country_pred_all[[col]][validation_sets[[set]]], use = "complete.obs")^2,
         MAE = mean(abs(country_pred_all[[col]] - country_pred_all$observed)[validation_sets[[set]]], na.rm = TRUE)) %>%
  ungroup() %>% select(-col)
save_table(country_fit, "B1_lco_country_R2_MAE")

# Extended country-level metrics (same as the country models in B4: R2, MAE,
# RMSE, mean bias = predicted - observed, calibration slope/intercept from
# regressing observed on predicted). Used by B4 for Tables S7/S8.
country_fit_ext <- expand_grid(set = names(validation_sets),
                               tibble(model = names(oof_cols), col = c("pred_enet", "pred_xgb", "pred_rf", "pred_me"))) %>%
  rowwise() %>%
  mutate(metrics = list({
    d <- country_pred_all[validation_sets[[set]], ]
    obs <- d$observed; pred <- d[[col]]; ok <- !is.na(pred)
    cal <- coef(lm(obs[ok] ~ pred[ok]))
    tibble(n_countries = sum(ok), n_respondents = sum(d$n[ok]),
           R2 = cor(obs[ok], pred[ok])^2, MAE = mean(abs(pred[ok] - obs[ok])),
           RMSE = sqrt(mean((pred[ok] - obs[ok])^2)), mean_bias = mean(pred[ok] - obs[ok]),
           calibration_intercept = unname(cal[1]), calibration_slope = unname(cal[2]))
  })) %>%
  ungroup() %>% select(-col) %>% unnest(metrics) %>%
  mutate(model_level = "Individual level")
save_table(country_fit_ext, "B1_country_level_fit_extended")
save_data(country_fit_ext, "individual_country_fit")   # used by B4 (Tables S7/S8)

# ── Figure S10: all four models, coloured by M49 region ──
country_pred_long <- country_pred_all %>%
  pivot_longer(cols = c(pred_enet, pred_xgb, pred_rf, pred_me), names_to = "model", values_to = "predicted") %>%
  mutate(model = factor(recode(model, pred_enet = "Elastic Net", pred_xgb = "XGBoost",
                               pred_rf = "Random Forest", pred_me = "Mixed-Effects LR"),
                        levels = names(oof_cols)),
         m49_region = factor(m49_region, levels = m49_region_levels))

ax_max <- max(country_pred_long$observed, na.rm = TRUE) * 1.1
p_all_models <- ggplot(country_pred_long, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = m49_region, size = n, shape = oos_source), alpha = 0.7) +
  geom_text(aes(label = country_code), size = 2, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~ model, ncol = 2) +
  scale_colour_manual(values = m49_col) +
  scale_shape_manual(values = c("LCO CV (70)" = 16, "Hold-out (6)" = 17)) +
  coord_equal(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  labs(title = "Predicted vs Observed Prevalence — All Models (out-of-sample, 76 countries)",
       x = "Predicted Prevalence", y = "Observed Prevalence",
       colour = "M49 region", size = "Sample size", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")
save_fig(p_all_models, "B1_figS10_pred_vs_obs_all_models.png", width = 12, height = 10, dpi = 300)

# ── Figure S9: ENET and XGB only, coloured by MODEL region ──
p_enet_xgb <- country_pred_long %>%
  filter(model %in% c("Elastic Net", "XGBoost")) %>%
  ggplot(aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = region_final, size = n, shape = oos_source), alpha = 0.7) +
  geom_text(aes(label = country_code), size = 2.5, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~ model) +
  scale_shape_manual(values = c("LCO CV (70)" = 16, "Hold-out (6)" = 17)) +
  coord_equal(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  labs(title = "Predicted vs Observed Water Insecurity Prevalence",
       subtitle = "Out-of-sample predictions: 70 countries in LCO CV + 6 hold-out countries (weighted)",
       x = "Observed Prevalence", y = "Predicted Prevalence",
       colour = "Model region", size = "Sample size", shape = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")
save_fig(p_enet_xgb, "B1_figS09_pred_vs_obs_enet_xgb.png", width = 12, height = 7, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 10. ADDITIONAL DIAGNOSTICS (supplementary checks)
# ══════════════════════════════════════════════════════════════════════════════

# Within-country AUC (XGB OOF): which countries are hardest to rank individuals in
country_auc <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final) %>%
  filter(n() >= 30, n_distinct(iwise12_imp) == 2) %>%
  summarise(n = n(), prev = mean(outcome_num),
            auc = as.numeric(auc(roc(outcome_num, xgb_oof_prob, quiet = TRUE))),
            .groups = "drop") %>%
  arrange(auc)
save_table(country_auc, "B1_country_level_auc_xgb")

# Predicted-probability distributions and calibration in the hold-out countries.
# Calibration: points close to the dashed 1:1 line = well calibrated.
prob_dist_plot <- function(y, p, title) {
  tibble(actual = factor(y, levels = c(0, 1), labels = c("Not water insecure", "Water insecure")),
         predicted_prob = p) %>%
    ggplot(aes(x = predicted_prob, fill = actual)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
    labs(title = title, x = "Predicted Probability", y = "Density", fill = "Actual Outcome") +
    theme_minimal()
}
calibration_plot <- function(y, p, title) {
  tibble(actual = y, predicted = p) %>%
    mutate(bin = cut(predicted, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(mean_predicted = mean(predicted), mean_actual = mean(actual), n = n(), .groups = "drop") %>%
    ggplot(aes(x = mean_predicted, y = mean_actual)) +
    geom_point(aes(size = n), colour = "#E41A1C") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_size_continuous(range = c(2, 8)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "Mean Predicted Probability", y = "Observed Proportion", size = "n obs") +
    theme_minimal()
}
p_dist <- (prob_dist_plot(y_test_focal, enet_prob_focal, "India Holdout (Elastic Net)") +
             prob_dist_plot(y_test_other, enet_prob_other, "Other Countries Holdout (Elastic Net)")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(p_dist, "B1_prob_distribution_holdouts.png", width = 14, height = 5, dpi = 300)

p_cal <- (calibration_plot(y_test_focal, enet_prob_focal, "Calibration - India (Elastic Net)") +
            calibration_plot(y_test_other, enet_prob_other, "Calibration - Other Countries (Elastic Net)")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(p_cal, "B1_calibration_holdouts.png", width = 12, height = 6, dpi = 300)

# Performance by survey year and by model region (ENET OOF).
# Look for: similar AUC across years/regions; a much lower AUC flags where the
# model generalises less well.
perf_by <- function(group_var) {
  train_df %>%
    mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= 50, n_distinct(iwise12_imp) == 2) %>%
    summarise(n_obs = n(), n_countries = n_distinct(country_code),
              prevalence = round(mean(outcome_num), 4),
              auc = round(as.numeric(auc(roc(outcome_num, enet_oof_prob, quiet = TRUE))), 4),
              .groups = "drop")
}
year_performance   <- perf_by("YEAR_WAVE")
region_performance <- perf_by("region_final")
save_table(year_performance,   "B1_sensitivity_year_performance")
save_table(region_performance, "B1_sensitivity_region_performance")

region_calibration <- country_pred_all %>%
  group_by(region_final) %>%
  filter(n() >= 3) %>%
  summarise(n_countries = n(),
            R2 = round(cor(observed, pred_enet)^2, 4),
            MAE = round(mean(abs(pred_enet - observed)), 4),
            mean_bias = round(mean(pred_enet - observed), 4),   # >0 = over-prediction
            .groups = "drop")
save_table(region_calibration, "B1_sensitivity_region_calibration")


# ══════════════════════════════════════════════════════════════════════════════
# 11. FEATURE IMPORTANCE: ENET COEFFICIENTS (Figure 6) AND SHAP (Figure S15)
# ══════════════════════════════════════════════════════════════════════════════

feature_type <- c(
  "Age" = "Sociodemographic", "Household size" = "Sociodemographic", "Female" = "Sociodemographic",
  "Secondary education" = "Sociodemographic", "College education" = "Sociodemographic",
  "Married/Partnered" = "Sociodemographic", "Separated/Widowed" = "Sociodemographic",
  "Underemployed" = "Sociodemographic", "Unemployed" = "Sociodemographic",
  "Out of workforce" = "Sociodemographic", "Towns/semi-dense" = "Sociodemographic",
  "Cities" = "Sociodemographic", "Urbanicity missing" = "Sociodemographic",
  "Urbanicity self-reported" = "Sociodemographic", "Population density (log)" = "Sociodemographic",
  "Income quintile 2" = "Economic", "Income quintile 3" = "Economic", "Income quintile 4" = "Economic",
  "Income quintile 5 (richest)" = "Economic", "Income inadequate" = "Economic",
  "Getting by on income" = "Economic", "GNI per capita (log)" = "Economic",
  "Water dissatisfied" = "Water & WASH", "Basic drinking water (%)" = "Water & WASH",
  "Water stress (log)" = "Water & WASH", "WASH mortality" = "Water & WASH",
  "Mean annual rainfall" = "Climate & environment",
  "ND-GAIN index (high score = low vulnerability)" = "Climate & environment",
  "Life expectancy" = "Health", "Fertility rate" = "Health", "Under-5 mortality" = "Health",
  "Corruption index (low score = high corruption)" = "Governance",
  "Income missing" = "Missingness", "Urbanicity missing flag" = "Missingness",
  "Eastern Africa" = "Region", "Western Africa" = "Region", "Middle Africa" = "Region",
  "Southern Africa" = "Region", "Latin America" = "Region", "MENA" = "Region",
  "Southern Asia" = "Region", "N. America/Europe/Australia" = "Region"
)
type_levels <- c("Water & WASH", "Economic", "Sociodemographic", "Climate & environment",
                 "Health", "Governance", "Region", "Missingness")
get_type <- function(x) ifelse(x %in% names(feature_type), feature_type[x], "Other")

# ── Figure 6: top-20 ENET coefficients, grouped with brackets ──
p_coef_data <- enet_coefs %>%
  head(20) %>%
  mutate(type = factor(get_type(feature_label), levels = type_levels),
         # Figure 6 only: "Water & WASH" shown as "Water" (grouped SHAP figure keeps "Water & WASH")
         type = fct_recode(type, "Water" = "Water & WASH"),
         direction = ifelse(coefficient > 0, "Increases predicted risk", "Decreases predicted risk")) %>%
  arrange(type, desc(abs(coefficient))) %>%
  mutate(feature_label = fct_inorder(feature_label), y_pos = row_number())

group_boundaries <- p_coef_data %>%
  group_by(type) %>%
  summarise(y_min = min(y_pos), y_max = max(y_pos), y_mid = mean(y_pos), .groups = "drop")

max_coef  <- max(abs(p_coef_data$coefficient))
bracket_x <- max_coef * 1.08

p_fig6 <- ggplot(p_coef_data, aes(x = coefficient, y = feature_label, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_segment(data = group_boundaries, aes(x = bracket_x, xend = bracket_x, y = y_min - 0.3, yend = y_max + 0.3),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_segment(data = group_boundaries, aes(x = bracket_x, xend = bracket_x - max_coef * 0.02, y = y_max + 0.3, yend = y_max + 0.3),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_segment(data = group_boundaries, aes(x = bracket_x, xend = bracket_x - max_coef * 0.02, y = y_min - 0.3, yend = y_min - 0.3),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_text(data = group_boundaries, aes(x = bracket_x + max_coef * 0.03, y = y_mid, label = type),
            inherit.aes = FALSE, hjust = 0, size = 2.5, fontface = "bold", colour = "grey30") +
  scale_fill_manual(values = c("Increases predicted risk" = "darkred", "Decreases predicted risk" = "navy")) +
  coord_cartesian(clip = "off") +
  labs(title = "Elastic Net coefficients", x = "Coefficient (log-odds scale)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.margin = margin(5, 120, 5, 5))
save_fig(p_fig6, "B1_fig6_enet_coefficients_bracketed.png", width = 10, height = 8, dpi = 300)

# ── SHAP values for XGBoost (random sample of 5,000 training respondents) ──
set.seed(seed)
X_sample <- X_train[sample(nrow(X_train), min(5000, nrow(X_train))), ]

shap <- shapviz(xgb_model, X_pred = X_sample)
colnames(shap$S) <- colnames(shap$X) <- relabel(colnames(shap$S))

save_fig(sv_importance(shap, kind = "beeswarm", max_display = 15) +
           labs(title = "SHAP Feature Importance - XGBoost") + theme_minimal(),
         "B1_figS15_shap_beeswarm_xgb.png", width = 10, height = 8, dpi = 300)
save_fig(sv_importance(shap, kind = "bar", max_display = 15) +
           labs(title = "Mean |SHAP| - XGBoost") + theme_minimal(),
         "B1_shap_bar_xgb.png", width = 8, height = 6, dpi = 300)

# Grouped SHAP beeswarm (top 20 by mean |SHAP|, grouped by feature type)
shap_df <- as.data.frame(shap$S) %>% mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = "feature", values_to = "shap_value") %>%
  left_join(as.data.frame(shap$X) %>% mutate(row_id = row_number()) %>%
              pivot_longer(-row_id, names_to = "feature", values_to = "feat_value"),
            by = c("row_id", "feature")) %>%
  mutate(type = factor(get_type(feature), levels = type_levels))

feature_order <- shap_df %>%
  group_by(feature, type) %>%
  summarise(mean_shap = mean(abs(shap_value)), .groups = "drop") %>%
  slice_max(mean_shap, n = 20) %>%
  arrange(type, desc(mean_shap)) %>%
  mutate(y_pos = row_number())

shap_plot_df <- shap_df %>%
  inner_join(feature_order %>% select(feature, y_pos), by = "feature") %>%
  mutate(feature = fct_reorder(feature, y_pos)) %>%
  group_by(feature) %>%
  mutate(feat_scaled = (feat_value - min(feat_value, na.rm = TRUE)) /
           (max(feat_value, na.rm = TRUE) - min(feat_value, na.rm = TRUE) + 1e-10)) %>%
  ungroup()

group_bounds <- feature_order %>%
  group_by(type) %>%
  summarise(y_min = min(y_pos), y_max = max(y_pos), y_mid = mean(y_pos), .groups = "drop")

set.seed(seed)
shap_plot_sample <- shap_plot_df %>% group_by(feature) %>% slice_sample(n = 500)   # groups with < 500 rows keep all rows %>% ungroup()
max_shap  <- max(abs(shap_plot_sample$shap_value), na.rm = TRUE)
bracket_x <- max(abs(shap_plot_df$shap_value), na.rm = TRUE) * 1.15

p_shap_grouped <- ggplot(shap_plot_sample, aes(x = shap_value, y = feature, colour = feat_scaled)) +
  geom_jitter(height = 0.25, width = 0, size = 0.5, alpha = 0.6) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
  scale_colour_gradient(low = "#0072B2", high = "#D55E00", name = "Feature value",
                        labels = c("Low", "High"), breaks = c(0, 1)) +
  geom_segment(data = group_bounds, aes(x = bracket_x, xend = bracket_x, y = y_min - 0.3, yend = y_max + 0.3),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_segment(data = group_bounds, aes(x = bracket_x, xend = bracket_x - max_shap * 0.02, y = y_max + 0.3, yend = y_max + 0.3),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_segment(data = group_bounds, aes(x = bracket_x, xend = bracket_x - max_shap * 0.02, y = y_min - 0.3, yend = y_min - 0.3),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_text(data = group_bounds, aes(x = bracket_x + max_shap * 0.03, y = y_mid, label = type),
            inherit.aes = FALSE, hjust = 0, size = 2.5, fontface = "bold", colour = "grey30") +
  coord_cartesian(clip = "off") +
  labs(title = "XGBoost SHAP values", x = "SHAP value (impact on prediction)", y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.margin = margin(5, 120, 5, 5))
save_fig(p_shap_grouped, "B1_figS15_shap_beeswarm_grouped.png", width = 10, height = 8, dpi = 300)

# ENET contributions in SHAP form (coefficient x standardised value) for comparison
sv_enet <- shapviz(object = sweep(X_sample, 2, coef(enet_final)[-1, 1], "*"),
                   X = as.data.frame(X_sample), baseline = coef(enet_final)[1, 1])
colnames(sv_enet$S) <- colnames(sv_enet$X) <- relabel(colnames(sv_enet$S))
save_fig(sv_importance(sv_enet, kind = "beeswarm", max_display = 15) +
           labs(title = "SHAP Feature Importance - Elastic Net") + theme_minimal(),
         "B1_shap_beeswarm_enet.png", width = 10, height = 8, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 12. SENSITIVITY: EXCLUDING WATER DISSATISFACTION (Figures S16, S17)
# ══════════════════════════════════════════════════════════════════════════════
# Water dissatisfaction is subjective and related to (but not the same as) water
# insecurity. Refit ENET and XGB without it and compare.
# Look for: how much AUC / country R² drop, and which predictors gain weight.

all_predictors_nowd <- setdiff(all_predictors, "waterdissatisfied_imp")

rec_nowd     <- make_recipe(train_df, all_predictors_nowd)
X_train_nowd <- bake_X(rec_nowd, train_df)
X_focal_nowd <- bake_X(rec_nowd, test_focal_df)
X_other_nowd <- bake_X(rec_nowd, test_other_df)

set.seed(seed)
cv_nowd   <- cv.glmnet(X_train_nowd, y_train, family = "binomial", weights = w_train,
                       alpha = best_alpha, nfolds = 5, type.measure = "auc")
enet_nowd <- glmnet(X_train_nowd, y_train, family = "binomial", weights = w_train,
                    alpha = best_alpha, lambda = cv_nowd$lambda.min)

dtrain_nowd <- xgb.DMatrix(X_train_nowd, label = y_train, weight = w_train_scaled)
set.seed(seed)
xgb_cv_nowd <- xgb.cv(params = xgb_params, data = dtrain_nowd, nrounds = 1000, nfold = 5,
                      early_stopping_rounds = 30, verbose = 0)
best_nrounds_nowd <- xgb_cv_nowd$evaluation_log$iter[which.max(xgb_cv_nowd$evaluation_log$test_auc_mean)]
xgb_nowd <- xgb.train(params = xgb_params, data = dtrain_nowd, nrounds = best_nrounds_nowd, verbose = 0)

sensitivity_wd <- bind_rows(
  evaluate_model(y_test_focal, enet_prob_focal, "Elastic Net (full)",  "India"),
  evaluate_model(y_test_focal, predict(enet_nowd, X_focal_nowd, type = "response")[, 1], "Elastic Net (no WD)", "India"),
  evaluate_model(y_test_other, enet_prob_other, "Elastic Net (full)",  "Other"),
  evaluate_model(y_test_other, predict(enet_nowd, X_other_nowd, type = "response")[, 1], "Elastic Net (no WD)", "Other"),
  evaluate_model(y_test_focal, xgb_prob_focal, "XGBoost (full)", "India"),
  evaluate_model(y_test_focal, predict(xgb_nowd, xgb.DMatrix(X_focal_nowd)), "XGBoost (no WD)", "India"),
  evaluate_model(y_test_other, xgb_prob_other, "XGBoost (full)", "Other"),
  evaluate_model(y_test_other, predict(xgb_nowd, xgb.DMatrix(X_other_nowd)), "XGBoost (no WD)", "Other")
)
save_table(sensitivity_wd, "B1_sensitivity_noWD_holdout")

# OOF predictions without WD (same folds as Section 8)
train_df$enet_nowd_oof <- NA_real_
train_df$xgb_nowd_oof  <- NA_real_
for (f in 1:k_folds) {
  idx_test  <- which(train_df$fold == f)
  idx_train <- which(train_df$fold != f)
  cv_rec  <- make_recipe(train_df[idx_train, ], all_predictors_nowd)
  X_cv_tr <- bake_X(cv_rec, train_df[idx_train, ])
  X_cv_te <- bake_X(cv_rec, train_df[idx_test, ])
  y_cv_tr <- as.numeric(as.character(train_df$iwise12_imp[idx_train]))
  w_cv_tr <- train_df[[weight_var]][idx_train]

  m_enet <- glmnet(X_cv_tr, y_cv_tr, family = "binomial", weights = w_cv_tr,
                   alpha = best_alpha, lambda = cv_nowd$lambda.min)
  train_df$enet_nowd_oof[idx_test] <- predict(m_enet, X_cv_te, type = "response")[, 1]

  m_xgb <- xgb.train(params = xgb_params,
                     data = xgb.DMatrix(X_cv_tr, label = y_cv_tr, weight = w_cv_tr / mean(w_cv_tr)),
                     nrounds = best_nrounds_nowd, verbose = 0)
  train_df$xgb_nowd_oof[idx_test] <- predict(m_xgb, xgb.DMatrix(X_cv_te))
}

country_pred_nowd <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final) %>%
  summarise(observed       = weighted.mean(outcome_num,   .data[[weight_var]]),
            pred_enet_full = weighted.mean(enet_oof_prob, .data[[weight_var]]),
            pred_enet_nowd = weighted.mean(enet_nowd_oof, .data[[weight_var]]),
            pred_xgb_full  = weighted.mean(xgb_oof_prob,  .data[[weight_var]]),
            pred_xgb_nowd  = weighted.mean(xgb_nowd_oof,  .data[[weight_var]]),
            .groups = "drop")
save_table(country_pred_nowd, "B1_sensitivity_noWD_country")
save_table(tibble(model = c("pred_enet_full", "pred_enet_nowd", "pred_xgb_full", "pred_xgb_nowd")) %>%
             rowwise() %>%
             mutate(R2  = cor(country_pred_nowd$observed, country_pred_nowd[[model]])^2,
                    MAE = mean(abs(country_pred_nowd[[model]] - country_pred_nowd$observed))) %>%
             ungroup(),
           "B1_sensitivity_noWD_country_R2_MAE")

# Figure S16: ENET coefficients without WD
enet_nowd_coefs <- coef(enet_nowd) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>% rename(coefficient = s0) %>%
  filter(coefficient != 0, feature != "(Intercept)") %>%
  mutate(feature_label = relabel(feature)) %>%
  arrange(desc(abs(coefficient)))

p_coef_nowd <- enet_nowd_coefs %>%
  head(15) %>%
  mutate(feature_label = fct_reorder(feature_label, abs(coefficient)),
         direction = ifelse(coefficient > 0, "Increases predicted risk", "Decreases predicted risk")) %>%
  ggplot(aes(x = coefficient, y = feature_label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Increases predicted risk" = "darkred", "Decreases predicted risk" = "navy")) +
  labs(title = "Elastic Net coefficients (excl. water dissatisfaction)",
       x = "Coefficient (log-odds scale)", y = NULL, fill = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
save_fig(p_coef_nowd, "B1_figS16_enet_coef_noWD.png", width = 9, height = 7, dpi = 300)

# Figure S17: SHAP without WD
set.seed(seed)
X_sample_nowd <- X_train_nowd[sample(nrow(X_train_nowd), min(5000, nrow(X_train_nowd))), ]
shap_nowd <- shapviz(xgb_nowd, X_pred = X_sample_nowd)
colnames(shap_nowd$S) <- colnames(shap_nowd$X) <- relabel(colnames(shap_nowd$S))
p_shap_nowd <- sv_importance(shap_nowd, kind = "beeswarm", max_display = 15) +
  labs(title = "XGBoost SHAP values (excl. water dissatisfaction)") + theme_minimal()
save_fig(p_shap_nowd, "B1_figS17_shap_noWD.png", width = 10, height = 8, dpi = 300)
save_fig(p_coef_nowd | p_shap_nowd, "B1_figS16_S17_combined_noWD.png", width = 16, height = 8, dpi = 300)

# Coefficient shift when WD is removed (largest shifts = predictors that
# partly carry the WD signal)
coef_shift <- enet_coefs_full %>% select(feature, coef_full = coefficient) %>%
  left_join(coef(enet_nowd) %>% as.matrix() %>% as.data.frame() %>%
              rownames_to_column("feature") %>% rename(coef_nowd = s0), by = "feature") %>%
  mutate(coef_nowd = replace_na(coef_nowd, 0), shift = coef_nowd - coef_full,
         feature_label = relabel(feature)) %>%
  filter(coef_full != 0 | coef_nowd != 0) %>%
  arrange(desc(abs(shift)))
save_table(coef_shift, "B1_coefficient_shift_noWD")


# ══════════════════════════════════════════════════════════════════════════════
# 13. SAVE MODELS AND DATA FOR LATER SCRIPTS
# ══════════════════════════════════════════════════════════════════════════════

save_model(rec,        "ind_recipe")
save_model(enet_final, "ind_enet_final")
save_model(enet_nowd,  "ind_enet_nowd")
invisible(xgb.save(xgb_model, file.path(dir_models, "ind_xgb_model.json")))

# Ranges of numeric predictors in the training data (for extrapolation checks in B2)
numeric_ranges <- train_df %>%
  summarise(across(all_of(c("age_imp", "hhsize", ctry_predictors_ind)),
                   list(min = ~ min(., na.rm = TRUE), max = ~ max(., na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("variable", ".value"), names_pattern = "(.*)_(min|max)$")

save_model(list(best_alpha = best_alpha, best_lambda = best_lambda,
                best_nrounds = best_nrounds, xgb_params = xgb_params,
                all_predictors = all_predictors, train_colnames = colnames(X_train),
                train_countries = levels(droplevels(train_df$country_code)),
                numeric_ranges = numeric_ranges, icc = icc),
           "ind_settings")

save_data(train_df, "train_df_oof")

stop_timer()
