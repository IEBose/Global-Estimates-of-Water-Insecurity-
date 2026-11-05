################################################################################
#  B6_Country_XGB_predictions.R — COUNTRY-LEVEL XGBOOST PREDICTIONS (SENSITIVITY)
#
#  Purpose : Predict prevalence for all 194 countries with the country-level
#            XGBoost model (out-of-sample prediction for the 76 observed countries).
#            Used in Table S10 (Prev xgb country) and to flag model disagreement
#            for Tier 3 countries (Table S9). NOT used in the main estimates.
#
#  Depends : A4, B4, B5 (tier assignment)
#  Inputs  : data/wb_imputed.rds, data/country_df_oos.rds, data/country_preds_enet.rds,
#            models/country_recipe.rds, country_xgb_model.json, country_settings.rds
#  Outputs : data/country_preds_xgb.rds, tables/B6_country_predictions_xgb_194.csv
#
#  How to read: for Tier 3 countries, R² near 1 and a small MAE between ENET and
#  XGB mean the two models agree; the largest |differences| are the least
#  certain estimates.
################################################################################

source("00_setup.R")
start_timer("B6")

wb                 <- load_data("wb_imputed")
country_df_oos     <- load_data("country_df_oos")
country_preds_enet <- load_data("country_preds_enet")
rec_c              <- load_model("country_recipe")
settings_c         <- load_model("country_settings")
xgb_c              <- xgb.load(file.path(dir_models, "country_xgb_model.json"))

all_pred <- wb %>%
  mutate(region_final = factor(region_final, levels = settings_c$region_levels),
         logit_prev   = 0)
X_all <- bake(rec_c, all_pred) %>% select(-logit_prev) %>% as.matrix()
stopifnot(identical(colnames(X_all), settings_c$train_colnames))

all_pred$pred_xgb_country_full <- inv_logit(predict(xgb_c, xgb.DMatrix(data = X_all)))

country_preds_xgb <- all_pred %>%
  select(iso3, pred_xgb_country_full) %>%
  left_join(country_df_oos %>% select(Country_Code, oos_xgb), by = c("iso3" = "Country_Code")) %>%
  mutate(prev_xgb_country = ifelse(!is.na(oos_xgb), oos_xgb, pred_xgb_country_full)) %>%
  left_join(country_preds_enet %>% select(iso3, tier, m49_region, prev_enet_country), by = "iso3") %>%
  mutate(diff_xgb_minus_enet = prev_xgb_country - prev_enet_country)

save_table(country_preds_xgb, "B6_country_predictions_xgb_194")
save_data(country_preds_xgb, "country_preds_xgb")

tier3 <- country_preds_xgb %>% filter(tier == "Predicted (country-level)")
message("Tier 3 (", nrow(tier3), " countries) ENET vs XGB: R² = ",
        round(cor(tier3$prev_enet_country, tier3$prev_xgb_country)^2, 3),
        ", MAE = ", round(mean(abs(tier3$diff_xgb_minus_enet)), 3))

stop_timer()
