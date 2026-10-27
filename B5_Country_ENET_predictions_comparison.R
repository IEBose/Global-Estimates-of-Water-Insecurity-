################################################################################
#  B5_Country_ENET_predictions_comparison.R — COUNTRY-LEVEL ENET PREDICTIONS
#
#  Purpose : (1) Predict prevalence with the country-level Elastic Net for all
#                194 countries. For the 76 observed countries the OUT-OF-SAMPLE
#                prediction is used (the model never saw that country).
#                The 33 Tier 3 countries (no GWP data) use these predictions in
#                the final estimates.
#            (2) Compare individual-level (Tier 2, B2) and country-level
#                predictions for the 85 countries that have both.
#            (3) Extrapolation check for Tier 3 countries (Table S9 input).
#
#  Depends : A4, B2 (Tier 2 country prevalence), B4 (country models)
#  Inputs  : data/wb_imputed.rds, data/country_df_oos.rds, data/new_data_pred.rds,
#            models/country_recipe.rds, country_enet_model.rds, country_beta_model.rds,
#            country_frac_model.rds, country_settings.rds
#  Outputs : data/country_preds_enet.rds
#            tables/B5_country_predictions_enet_194.csv
#            tables/B5_individual_vs_country_model_85.csv
#            tables/B5_tier3_extrapolation_by_variable.csv, _by_country.csv
#            figures/B5_individual_vs_country_model_85.png
#
#  How to read the individual vs country comparison (85 countries):
#    R² near 1 and small MAE = the country-level model reproduces the
#    individual-level estimates well, supporting its use where individual data
#    are missing. Large differences flag countries for Table S9.
################################################################################

source("00_setup.R")
start_timer("B5")
load_pkgs("betareg")

wb               <- load_data("wb_imputed")
country_df_oos <- load_data("country_df_oos")
new_data         <- load_data("new_data_pred")
rec_c            <- load_model("country_recipe")
enet_c           <- load_model("country_enet_model")
beta_c           <- load_model("country_beta_model")
frac_c           <- load_model("country_frac_model")
settings_c       <- load_model("country_settings")

observed_codes <- unique(as.character(country_df_oos$Country_Code))
tier2_codes    <- unique(as.character(new_data$country_code))


# ══════════════════════════════════════════════════════════════════════════════
# 1. ENET PREDICTIONS FOR ALL 194 COUNTRIES
# ══════════════════════════════════════════════════════════════════════════════

all_pred <- wb %>%
  mutate(region_final = factor(region_final, levels = settings_c$region_levels),
         logit_prev   = 0)                       # placeholder required by the recipe

# CHECK: every country has a model region known to the trained model (expect 0)
message("Countries with no usable model region: ", sum(is.na(all_pred$region_final)))

X_all <- bake(rec_c, all_pred) %>% select(-logit_prev) %>% as.matrix()
stopifnot(identical(colnames(X_all), settings_c$train_colnames))

all_pred$pred_enet_country_full <- inv_logit(predict(enet_c, X_all, type = "response")[, 1])

country_preds_enet <- all_pred %>%
  select(iso3, m49_region, m49_sub_region, pred_enet_country_full) %>%
  left_join(country_df_oos %>% select(Country_Code, oos_enet), by = c("iso3" = "Country_Code")) %>%
  mutate(
    # Observed countries: out-of-sample prediction (B4); others: full-model prediction
    prev_enet_country = ifelse(!is.na(oos_enet), oos_enet, pred_enet_country_full),
    tier = case_when(iso3 %in% observed_codes ~ "Observed",
                     iso3 %in% tier2_codes    ~ "Predicted (individual)",
                     TRUE                     ~ "Predicted (country-level)")
  )

# CHECK: 76 / 85 / 33
print(count(country_preds_enet, tier))
save_table(country_preds_enet, "B5_country_predictions_enet_194")
save_data(country_preds_enet, "country_preds_enet")


# ══════════════════════════════════════════════════════════════════════════════
# 2. INDIVIDUAL-LEVEL vs COUNTRY-LEVEL PREDICTIONS (85 TIER 2 COUNTRIES)
# ══════════════════════════════════════════════════════════════════════════════

tier2_indiv <- new_data %>%
  group_by(country_code) %>%
  summarise(n = n(), prev_individual = weighted.mean(pred_prob, .data[[weight_var]]), .groups = "drop") %>%
  mutate(country_code = as.character(country_code))

country_85 <- all_pred %>%
  filter(iso3 %in% tier2_codes) %>%
  mutate(prev_beta = predict(beta_c, newdata = ., type = "response"),
         prev_frac = predict(frac_c, newdata = ., type = "response")) %>%
  left_join(country_preds_enet %>% select(iso3, prev_enet_country), by = "iso3") %>%
  select(country_code = iso3, region_final, m49_region, prev_beta, prev_frac, prev_enet_country) %>%
  left_join(tier2_indiv, by = "country_code") %>%
  mutate(diff_beta = prev_beta - prev_individual,
         diff_frac = prev_frac - prev_individual,
         diff_enet = prev_enet_country - prev_individual)
save_table(country_85, "B5_individual_vs_country_model_85")

save_table(tibble(country_model = c("prev_beta", "prev_frac", "prev_enet_country")) %>%
             rowwise() %>%
             mutate(R2  = cor(country_85$prev_individual, country_85[[country_model]])^2,
                    MAE = mean(abs(country_85[[country_model]] - country_85$prev_individual))) %>%
             ungroup(),
           "B5_individual_vs_country_agreement")

p_85 <- country_85 %>%
  pivot_longer(c(prev_beta, prev_frac, prev_enet_country), names_to = "country_model", values_to = "prev_country") %>%
  mutate(country_model = recode(country_model, prev_beta = "Beta Regression",
                                prev_frac = "Fractional Logit", prev_enet_country = "Elastic Net (country)"),
         m49_region = factor(m49_region, levels = m49_region_levels)) %>%
  ggplot(aes(x = prev_individual, y = prev_country)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = m49_region), size = 2, alpha = 0.7) +
  facet_wrap(~ country_model) +
  scale_colour_manual(values = m49_col) +
  labs(title = "Individual-level vs Country-level Model Predictions",
       subtitle = "85 countries with individual data",
       x = "Individual-level model prediction", y = "Country-level model prediction",
       colour = "M49 region") +
  theme_minimal() + theme(legend.position = "bottom")
save_fig(p_85, "B5_individual_vs_country_model_85.png", width = 14, height = 6, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 3. EXTRAPOLATION CHECK — TIER 3 COUNTRIES (Table S9 input)
# ══════════════════════════════════════════════════════════════════════════════
# Compares each Tier 3 country's predictor values with the range in the 76
# training countries. More variables outside the range = less certain estimate
# (e.g. Pacific islands with low basic water access).

tier3 <- all_pred %>% filter(!iso3 %in% c(observed_codes, tier2_codes))
train_ranges <- country_df_oos %>%
  summarise(across(all_of(country_vars), list(train_min = min, train_max = max))) %>%
  pivot_longer(everything(), names_to = c("variable", ".value"),
               names_pattern = "(.*)_(train_min|train_max)$")

tier3_long <- tier3 %>%
  select(iso3, all_of(country_vars)) %>%
  pivot_longer(-iso3, names_to = "variable", values_to = "value") %>%
  left_join(train_ranges, by = "variable") %>%
  mutate(below = value < train_min, above = value > train_max, outside = below | above)

save_table(tier3_long %>%
             group_by(variable, train_min, train_max) %>%
             summarise(tier3_min = min(value, na.rm = TRUE), tier3_max = max(value, na.rm = TRUE),
                       n_below = sum(below, na.rm = TRUE), n_above = sum(above, na.rm = TRUE),
                       countries = paste(iso3[outside %in% TRUE], collapse = ", "), .groups = "drop"),
           "B5_tier3_extrapolation_by_variable")

save_table(tier3_long %>%
             group_by(iso3) %>%
             summarise(n_extrapolating = sum(outside, na.rm = TRUE),
                       variables = paste(variable[outside %in% TRUE], collapse = ", "), .groups = "drop") %>%
             arrange(desc(n_extrapolating)),
           "B5_tier3_extrapolation_by_country")

stop_timer()
