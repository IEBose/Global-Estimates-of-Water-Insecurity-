################################################################################
#  B3_Individual_XGB_predictions.R — TIER 2 PREDICTIONS WITH XGBOOST (SENSITIVITY)
#
#  Purpose : Predict the same 85 Tier 2 countries with the individual-level
#            XGBoost model, to show what the estimates would have been with the
#            non-linear model. Used in Table S10 (Prev xgb individual) and for
#            flagging model disagreement (Table S9). NOT used in the main estimates.
#
#  Depends : B1 (recipe + XGB model), B2 (prepared new_data with ENET predictions)
#  Inputs  : data/new_data_pred.rds, models/ind_recipe.rds, models/ind_xgb_model.json
#  Outputs : tables/B3_tier2_country_prevalence_enet_vs_xgb.csv
#            data/tier2_country_xgb.rds
#
#  How to read the agreement statistics:
#    R²  (squared correlation of country prevalences, ENET vs XGB): close to 1
#        = the two models rank/space countries similarly.
#    MAE (mean absolute difference): average gap between the two models'
#        country prevalences (0.04 = 4 percentage points).
#    Countries with the largest |difference| are the least certain (Table S9).
################################################################################

source("00_setup.R")
start_timer("B3")

new_data  <- load_data("new_data_pred")
rec       <- load_model("ind_recipe")
xgb_model <- xgb.load(file.path(dir_models, "ind_xgb_model.json"))

X_new <- bake(rec, new_data) %>% select(-iwise12_imp) %>% as.matrix()
new_data$pred_xgb <- predict(xgb_model, xgb.DMatrix(data = X_new))

tier2_compare <- new_data %>%
  group_by(country_code, m49_region) %>%
  summarise(prev_enet_individual = weighted.mean(pred_prob, .data[[weight_var]]),
            prev_xgb_individual  = weighted.mean(pred_xgb,  .data[[weight_var]]),
            .groups = "drop") %>%
  mutate(country_code = as.character(country_code),
         diff_xgb_minus_enet = prev_xgb_individual - prev_enet_individual) %>%
  arrange(desc(abs(diff_xgb_minus_enet)))
save_table(tier2_compare, "B3_tier2_country_prevalence_enet_vs_xgb")

message("ENET vs XGB (Tier 2 countries): R² = ",
        round(cor(tier2_compare$prev_enet_individual, tier2_compare$prev_xgb_individual)^2, 3),
        ", MAE = ", round(mean(abs(tier2_compare$diff_xgb_minus_enet)), 3))

save_data(tier2_compare, "tier2_country_xgb")

stop_timer()
