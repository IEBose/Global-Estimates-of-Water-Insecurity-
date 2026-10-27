################################################################################
#  B9_Sensitivity_region_estimates.R — SENSITIVITY TO MODEL-REGION ASSIGNMENT
#
#  Purpose : Check how much Tier 2 estimates depend on the model region
#            assigned to a country. Five countries are re-assigned and re-predicted
#            with the SAME recipe and Elastic Net model (no retraining):
#              Japan, South Korea, Singapore -> "N. America, Europe & Australia"
#              Sudan, Yemen                  -> "Eastern Africa"
#            (Table S11)
#
#  Depends : B1 (recipe + ENET), B2 (prepared new_data with predictions)
#  Inputs  : data/new_data_pred.rds, models/ind_recipe.rds, models/ind_enet_final.rds,
#            models/ind_settings.rds
#  Outputs : tables/B9_tableS9_region_reclassification.csv
#
#  How to read: diff = re-assigned minus original prevalence (proportion scale).
#  Small differences mean the estimate is robust to the region choice.
################################################################################

source("00_setup.R")
start_timer("B9")

new_data   <- load_data("new_data_pred")
rec        <- load_model("ind_recipe")
enet_final <- load_model("ind_enet_final")
settings   <- load_model("ind_settings")

region_changes <- tribble(
  ~country_code, ~new_region,
  "JPN", "N. America, Europe & Australia",
  "KOR", "N. America, Europe & Australia",
  "SGP", "N. America, Europe & Australia",
  "SDN", "Eastern Africa",
  "YEM", "Eastern Africa"
)

sens <- new_data %>%
  mutate(country_code = as.character(country_code)) %>%
  filter(country_code %in% region_changes$country_code) %>%
  left_join(region_changes, by = "country_code") %>%
  mutate(original_region = region_final,
         region_final    = factor(new_region, levels = region_final_levels))

X_sens <- bake(rec, sens) %>% select(-iwise12_imp) %>% as.matrix()
stopifnot(identical(colnames(X_sens), settings$train_colnames))   # same columns as training

sens$pred_prob_sens <- predict(enet_final, X_sens, type = "response")[, 1]

comparison <- sens %>%
  group_by(country_code, original_region, new_region) %>%
  summarise(n = n(),
            orig_prev = weighted.mean(pred_prob,      .data[[weight_var]]),
            sens_prev = weighted.mean(pred_prob_sens, .data[[weight_var]]),
            .groups = "drop") %>%
  mutate(diff = sens_prev - orig_prev, pct_change = diff / orig_prev * 100)

save_table(comparison, "B9_tableS11_region_reclassification")

stop_timer()
