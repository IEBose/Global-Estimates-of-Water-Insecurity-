################################################################################
#  B8_Global_Estimates_CIs.R — FINAL COUNTRY, REGIONAL AND GLOBAL ESTIMATES + CIs
#
#  Purpose : Combine the three tiers into one estimate per country, quantify
#            uncertainty, and aggregate to M49 regions, M49 sub-regions, World
#            Bank income groups and the globe (population 15+ weighted).
#              Tier 1  Observed (76)                 survey prevalence
#              Tier 2  Individual-level model (85)   ENET (B2)
#              Tier 3  Country-level model (33)      ENET (B5)
#            All tables and the data behind the figures are produced here;
#            the figures themselves are drawn in C1-C3.
#
#  Uncertainty (95% CI) by tier:
#    Tier 1  Taylor-series linearisation from the survey design (PSU, strata, weights);
#            logit 95% CI with design df, as Stata svy: proportion (prop_ci_logit)
#    Tier 2  sampling (cluster bootstrap of PSUs, 1,000 reps) and model uncertainty
#            (residual bootstrap, below) combined in quadrature ON THE LOGIT SCALE:
#            SE_logit = sqrt(SE_sampling,logit² + SE_model,logit²);
#            CI = inverse logit(logit(estimate) ± 1.96 SE_logit)
#    Tier 3  model uncertainty only (residual bootstrap): 2.5th-97.5th percentiles
#    Residual bootstrap: out-of-sample residuals for all 76 observed countries
#            (70 from 5-fold LCO CV + 6 hold-outs), on the LOGIT scale
#            (logit(observed) - logit(predicted)), are resampled, added to the
#            logit of the predicted prevalence and transformed back (1,000 reps).
#            Draws therefore stay between 0 and 1 without clipping, and
#            intervals narrow towards 0% and 100%. Each tier uses its own model's
#            residuals: Tier 2 = individual-level ENET (B1), Tier 3 = country-level ENET (B4).
#  Regional / sub-regional / income-group / global CIs: 1,000 bootstrap reps;
#    in each, countries are resampled with replacement within each group and
#    each country's prevalence is drawn on the logit scale, from
#    Normal(logit(estimate), SE_logit), and transformed back; the
#    population-weighted mean is taken. CI = 2.5th-97.5th percentiles.
#    The global CI from the regional bootstrap is used for all groupings so
#    that the global figure is identical everywhere.
#
#  Depends : A1, B1, B2, B3, B4, B5, B6
#  Outputs : data/all_tiers.rds (final country table, used by C1-C3)
#            tables/B8_all_country_estimates.csv (B8 intervals; Table S10 with the
#              main UIs is made by B11)
#            tables/B8_tableS09_flagged_countries.csv (Table S9, formatted)
#            tables/B8_uncertainty_flag_inputs.csv, B8_outside_training_range_all.csv,
#            B8_predictors_regional_median_by_country.csv, B8_country_flags.csv
#            data/country_flags.rds (Table S9 reasons per country, used by B11)
#            tables/B8_country_CIs_*.csv, B8_ci_<group>.csv, B8_affected_<group>.csv,
#            B8_tier_contributions_<group>.csv
#
#  What to look out for:
#    * Tier counts 76 / 85 / 33 and no missing population or M49 region.
#    * CI widths by tier: observed < individual-predicted < country-predicted
#      is expected.
#    * Global prevalence is identical in the region and income tables.
################################################################################

source("00_setup.R")
start_timer("B8")
load_pkgs("survey")
options(survey.lonely.psu = "certainty")

n_boot <- 1000

df                   <- load_data("df_clean")
train_df             <- load_data("train_df_oof")
new_data             <- load_data("new_data_pred")
holdout_country_pred <- load_data("holdout_country_pred")
tier2_xgb            <- load_data("tier2_country_xgb")
country_preds_enet   <- load_data("country_preds_enet")
country_preds_xgb    <- load_data("country_preds_xgb")
wb                   <- load_data("wb_clean")

# Income group: FY27 classification (calendar-year 2025 data); missing -> "Unclassified"
wb <- wb %>% mutate(income_group = factor(ifelse(is.na(wb_income_fy27), "Unclassified", wb_income_fy27),
                                          levels = income_levels))


# ══════════════════════════════════════════════════════════════════════════════
# 1. COUNTRY ESTIMATES BY TIER
# ══════════════════════════════════════════════════════════════════════════════

# ── Tier 1: survey-design prevalence and SE (Taylor linearisation) ──
tier1_ci <- observed_country_prevalence(df)          # same function as the B4 outcome
tier1_ci <- tier1_ci %>%
  bind_cols(prop_ci_logit(tier1_ci$prevalence, tier1_ci$se, tier1_ci$df_design)) %>%
  select(-df_design) %>%
  mutate(tier = "Observed", method = "Survey design (logit CI)")

# ── Tier 2: point estimate and sampling uncertainty (cluster bootstrap of PSUs) ──
set.seed(seed)
tier2_boot <- new_data %>%
  group_by(country_code) %>%
  group_modify(function(d, ...) {
    psu_list  <- split(seq_len(nrow(d)), d$psu)
    psu_names <- names(psu_list)
    boot_prev <- replicate(n_boot, {
      idx <- unlist(psu_list[sample(psu_names, length(psu_names), replace = TRUE)], use.names = FALSE)
      weighted.mean(d$pred_prob[idx], d$wgtnorm_using[idx])
    })
    tibble(prevalence = weighted.mean(d$pred_prob, d$wgtnorm_using),
           se_sampling = sd(boot_prev),
           ci_lower_sampling = quantile(boot_prev, 0.025),
           ci_upper_sampling = quantile(boot_prev, 0.975))
  }) %>%
  ungroup() %>%
  mutate(country_code = as.character(country_code))

# ── Model uncertainty: residual bootstrap ──
# Out-of-sample residuals (observed - predicted country prevalence) for all 76
# observed countries: 70 from 5-fold LCO CV + 6 hold-out countries.
#   Tier 2 uses the INDIVIDUAL-level ENET residuals (B1)
#   Tier 3 uses the COUNTRY-level ENET residuals (B4)
obs_prev <- tier1_ci %>% select(country_code, observed = prevalence)

resid_individual <- bind_rows(
  train_df %>%
    group_by(country_code) %>%
    summarise(predicted = weighted.mean(enet_oof_prob, wgtnorm_using), .groups = "drop") %>%
    mutate(country_code = as.character(country_code)),
  holdout_country_pred %>% select(country_code, predicted = prev_enet_individual)
) %>%
  left_join(obs_prev, by = "country_code") %>%
  mutate(residual = observed - predicted, residual_logit = logit_c(observed) - logit_c(predicted))

resid_country <- load_data("country_df_oos") %>%
  transmute(country_code = Country_Code, observed = Mod_high_WI, predicted = oos_enet,
            residual = observed - predicted, residual_logit = logit_c(observed) - logit_c(predicted))

# CHECK: 76 residuals each. Mean near 0 = little systematic bias; SD = typical
# size of the error (expect the country-level SD to be larger)
save_table(bind_rows(resid_individual %>% mutate(model = "Individual ENET"),
                     resid_country    %>% mutate(model = "Country ENET")), "B8_out_of_sample_residuals")
for (r in list(list("Individual", resid_individual), list("Country", resid_country))) {
  message(r[[1]], "-level residuals (", nrow(r[[2]]), " countries): mean = ",
          round(mean(r[[2]]$residual), 4), ", SD = ", round(sd(r[[2]]$residual), 4))
}

# Residuals added on the logit scale (residuals_logit), then transformed back
residual_bootstrap <- function(point_df, residuals_logit) {
  point_df %>%
    crossing(boot_id = 1:n_boot) %>%
    mutate(boot_logit = logit_c(prev) + sample(residuals_logit, n(), replace = TRUE),
           boot_prev  = plogis(boot_logit)) %>%
    group_by(country_code) %>%
    summarise(prevalence     = first(prev),
              ci_lower       = quantile(boot_prev, 0.025),
              ci_upper       = quantile(boot_prev, 0.975),
              se_model       = sd(boot_prev),        # prevalence scale (reported)
              se_model_logit = sd(boot_logit),       # logit scale (used for combining / group draws)
              .groups = "drop")
}

tier3_codes <- country_preds_enet$iso3[country_preds_enet$tier == "Predicted (country-level)"]

# (random-number stream continues from the cluster bootstrap)
tier2_resid <- residual_bootstrap(tier2_boot %>% select(country_code, prev = prevalence),
                                  resid_individual$residual_logit)
tier3_resid <- residual_bootstrap(country_preds_enet %>% filter(iso3 %in% tier3_codes) %>%
                                    select(country_code = iso3, prev = prev_enet_country),
                                  resid_country$residual_logit)

# ── Tier 2 final: sampling + model SE in quadrature, on the logit scale ──
tier2_ci <- tier2_boot %>%
  select(country_code, prevalence, se_sampling) %>%
  left_join(tier2_resid %>% select(country_code, se_model, se_model_logit), by = "country_code") %>%
  mutate(se       = sqrt(se_sampling^2 + se_model^2),                              # prevalence scale (reported)
         se_logit = sqrt(logit_se(prevalence, se_sampling)^2 + se_model_logit^2),
         ci_lower = plogis(logit_c(prevalence) - 1.96 * se_logit),
         ci_upper = plogis(logit_c(prevalence) + 1.96 * se_logit),
         tier = "Predicted (individual)", method = "Combined (sampling + model), logit scale")

tier3_ci <- tier3_resid %>%
  rename(se = se_model, se_logit = se_model_logit) %>%
  mutate(tier = "Predicted (country-level)", method = "Residual bootstrap, logit scale")

# Tier 1: delta-method SE on the logit scale (for the group draws)
tier1_ci <- tier1_ci %>% mutate(se_logit = logit_se(prevalence, se))

country_ci <- bind_rows(tier1_ci, tier2_ci, tier3_ci) %>%
  select(country_code, prevalence, se, se_logit, ci_lower, ci_upper, tier, method,
         any_of(c("se_sampling", "se_model")))
save_table(country_ci, "B8_country_CIs_final")

# CI widths by tier (expected: observed < individual < country-level)
save_table(country_ci %>% mutate(width = ci_upper - ci_lower) %>%
             group_by(tier, method) %>%
             summarise(n = n(), mean_width = mean(width), median_width = median(width),
                       min_width = min(width), max_width = max(width), .groups = "drop"),
           "B8_country_CI_width_by_tier")

# Earlier approach, kept for comparison only (NOT used): Tier 2 sampling-only
# CI and Tier 3 CI = range between country-level ENET and XGB
save_table(bind_rows(
  tier2_boot %>% transmute(country_code, prevalence, ci_lower = ci_lower_sampling,
                           ci_upper = ci_upper_sampling, method = "Cluster bootstrap (sampling only)"),
  country_preds_xgb %>% filter(iso3 %in% tier3_codes) %>%
    transmute(country_code = iso3, prevalence = prev_enet_country,
              ci_lower = pmin(prev_enet_country, prev_xgb_country),
              ci_upper = pmax(prev_enet_country, prev_xgb_country),
              method = "Model disagreement range")
), "B8_country_CIs_previous_method_comparison")


# ══════════════════════════════════════════════════════════════════════════════
# 2. FINAL COUNTRY TABLE (Table S10 input) — all_tiers
# ══════════════════════════════════════════════════════════════════════════════
# Region Model  = region used as a predictor (M49-based rule for all tiers)
# Reported region = M49 region / sub-region (used for all aggregation)

model_region <- bind_rows(
  df       %>% distinct(country_code, region_model = region_final),
  new_data %>% distinct(country_code, region_model = region_final)
) %>% mutate(country_code = as.character(country_code), region_model = as.character(region_model))

ind_preds <- bind_rows(
  train_df %>%
    group_by(country_code) %>%
    summarise(prev_enet_individual = weighted.mean(enet_oof_prob, wgtnorm_using),
              prev_xgb_individual  = weighted.mean(xgb_oof_prob,  wgtnorm_using), .groups = "drop") %>%
    mutate(country_code = as.character(country_code), source = "LCO CV out-of-fold"),
  holdout_country_pred,
  tier2_xgb %>% select(country_code, prev_enet_individual, prev_xgb_individual) %>%
    mutate(source = "Prediction")
)

all_tiers <- wb %>%
  select(country_code = iso3, country_name = un_country_name, m49_region, m49_sub_region,
         income_group, pop_15plus, region_final_m49 = region_final) %>%
  left_join(country_preds_enet %>% select(iso3, tier, prev_enet_country), by = c("country_code" = "iso3")) %>%
  left_join(country_preds_xgb  %>% select(iso3, prev_xgb_country),        by = c("country_code" = "iso3")) %>%
  left_join(model_region, by = "country_code") %>%
  mutate(region_model = coalesce(region_model, as.character(region_final_m49))) %>%
  select(-region_final_m49) %>%
  left_join(tier1_ci %>% select(country_code, prev_observed = prevalence), by = "country_code") %>%
  left_join(ind_preds %>% select(country_code, prev_enet_individual, prev_xgb_individual,
                                 individual_pred_source = source), by = "country_code") %>%
  left_join(country_ci %>% select(country_code, prevalence, se, se_logit, ci_lower, ci_upper, ci_method = method),
            by = "country_code") %>%
  mutate(prev_best = case_when(tier == "Observed"                  ~ prev_observed,
                               tier == "Predicted (individual)"    ~ prev_enet_individual,
                               tier == "Predicted (country-level)" ~ prev_enet_country),
         tier       = factor(tier, levels = tier_levels),
         m49_region = factor(m49_region, levels = m49_region_levels),
         m49_sub_region = factor(m49_sub_region, levels = m49_sub_region_levels))

# CHECKS — expected: 194 rows; 76/85/33; prev_best equals the CI point estimate;
# 0 missing M49 region / population / estimate
message("all_tiers: ", nrow(all_tiers), " countries")
print(count(all_tiers, tier))
message("prev_best differs from CI point estimate: ", sum(abs(all_tiers$prev_best - all_tiers$prevalence) > 1e-9, na.rm = TRUE),
        " | missing M49 region: ", sum(is.na(all_tiers$m49_region)),
        " | missing sub-region: ", sum(is.na(all_tiers$m49_sub_region)),
        " | missing pop 15+: ", sum(is.na(all_tiers$pop_15plus)),
        " | missing estimate: ", sum(is.na(all_tiers$prev_best)))

save_table(all_tiers %>% select(country_code, country_name, region_model, m49_region, m49_sub_region,
                                income_group, tier, prev_observed, prev_enet_individual,
                                prev_xgb_individual, prev_enet_country, prev_xgb_country,
                                prev_best, ci_lower, ci_upper, ci_method, pop_15plus,
                                individual_pred_source),
           "B8_all_country_estimates")
save_data(all_tiers, "all_tiers")


# ══════════════════════════════════════════════════════════════════════════════
# 3. REGIONAL, SUB-REGIONAL, INCOME-GROUP AND GLOBAL ESTIMATES WITH CIs
# ══════════════════════════════════════════════════════════════════════════════

# drop_missing_first = TRUE removes countries with no population / group BEFORE
# resampling (as in the original sub-region code); FALSE drops them after the
# draws (as in the original region and income code). This only changes the
# random draws if some country has a missing population or group.
boot_group_ci <- function(data, group_var, drop_missing_first = FALSE) {
  data <- data %>% mutate(across(all_of(group_var), as.character))  # alphabetical group order
  if (drop_missing_first) data <- data %>% filter(!is.na(pop_15plus), !is.na(.data[[group_var]]))
  set.seed(seed)
  boot <- lapply(1:n_boot, function(b) {
    bd <- data %>%
      group_by(across(all_of(group_var))) %>%
      group_modify(~ slice_sample(.x, n = nrow(.x), replace = TRUE)) %>%
      ungroup() %>%
      mutate(boot_prev = plogis(rnorm(n(), mean = logit_c(prevalence), sd = se_logit))) %>%   # logit scale
      filter(!is.na(pop_15plus))
    bind_rows(
      bd %>% group_by(group = as.character(.data[[group_var]])) %>%
        summarise(prev = weighted.mean(boot_prev, pop_15plus), .groups = "drop"),
      bd %>% summarise(group = "Global", prev = weighted.mean(boot_prev, pop_15plus))
    ) %>% mutate(boot_id = b)
  }) %>% bind_rows()

  point <- bind_rows(
    data %>% filter(!is.na(pop_15plus)) %>% group_by(group = as.character(.data[[group_var]])) %>%
      summarise(prevalence = weighted.mean(prevalence, pop_15plus), total_pop = sum(pop_15plus),
                n_countries = n(), .groups = "drop"),
    data %>% filter(!is.na(pop_15plus)) %>%
      summarise(group = "Global", prevalence = weighted.mean(prevalence, pop_15plus),
                total_pop = sum(pop_15plus), n_countries = n())
  )

  boot %>%
    group_by(group) %>%
    summarise(ci_lower = max(quantile(prev, 0.025), 0), ci_upper = min(quantile(prev, 0.975), 1),
              .groups = "drop") %>%
    right_join(point, by = "group") %>%
    rename(!!group_var := group)
}

# Row order (by tier, then country) matches the original, so the bootstrap
# draws are reproduced
ci_input <- all_tiers %>% select(country_code, prevalence, se, se_logit, tier, pop_15plus,
                                 m49_region, m49_sub_region, income_group) %>%
  arrange(tier, country_code)

ci_region    <- boot_group_ci(ci_input, "m49_region")
ci_subregion <- boot_group_ci(ci_input, "m49_sub_region", drop_missing_first = TRUE)   # Figures S2, S3 (B8 version)
ci_income    <- boot_group_ci(ci_input, "income_group")

# Use the global estimate and CI from the regional bootstrap everywhere
global_row <- ci_region %>% filter(m49_region == "Global")
use_global <- function(tab, group_var) {
  tab %>% mutate(
    ci_lower = ifelse(.data[[group_var]] == "Global", global_row$ci_lower, ci_lower),
    ci_upper = ifelse(.data[[group_var]] == "Global", global_row$ci_upper, ci_upper))
}
ci_subregion <- use_global(ci_subregion, "m49_sub_region")
ci_income    <- use_global(ci_income,    "income_group")

# Number of people 15+ with moderate-to-high water insecurity, with CI
add_affected <- function(tab) {
  tab %>% mutate(affected_millions = prevalence * total_pop / 1e6,
                 affected_lower_m  = ci_lower  * total_pop / 1e6,
                 affected_upper_m  = ci_upper  * total_pop / 1e6,
                 pop_millions      = total_pop / 1e6)
}
ci_region    <- add_affected(ci_region)
ci_subregion <- add_affected(ci_subregion)
ci_income    <- add_affected(ci_income)

save_table(ci_region,    "B8_ci_m49_region")
save_table(ci_subregion, "B8_ci_m49_sub_region")
save_table(ci_income,    "B8_ci_income_group")
save_data(list(region = ci_region, sub_region = ci_subregion, income = ci_income), "group_CIs")

# CHECK: global prevalence identical in region and income tables (should be 0)
message("Global prevalence difference region vs income table: ",
        signif(ci_region$prevalence[ci_region$m49_region == "Global"] -
                 ci_income$prevalence[ci_income$income_group == "Global"], 3))


# ══════════════════════════════════════════════════════════════════════════════
# 4. TIER CONTRIBUTIONS (data for the stacked bars in Figure 2 / S2 / S3)
# ══════════════════════════════════════════════════════════════════════════════
# Contribution of each tier = (sum of tier prevalence x pop) / total group pop;
# the contributions add up to the group prevalence.

tier_contributions <- function(group_var) {
  base <- all_tiers %>% filter(!is.na(pop_15plus)) %>% mutate(affected = prevalence * pop_15plus)
  bind_rows(
    base %>% group_by(group = as.character(.data[[group_var]]), tier),
    base %>% mutate(group = "Global") %>% group_by(group, tier)
  ) %>%
    summarise(n_countries = n(), tier_pop = sum(pop_15plus), tier_affected = sum(affected), .groups = "drop") %>%
    group_by(group) %>%
    mutate(total_pop = sum(tier_pop),
           prev_contribution = tier_affected / total_pop,
           affected_millions = tier_affected / 1e6) %>%
    ungroup() %>%
    rename(!!group_var := group)
}
contrib <- list(region     = tier_contributions("m49_region"),
                sub_region = tier_contributions("m49_sub_region"),
                income     = tier_contributions("income_group"))
save_table(contrib$region,     "B8_tier_contributions_m49_region")
save_table(contrib$sub_region, "B8_tier_contributions_m49_sub_region")
save_table(contrib$income,     "B8_tier_contributions_income_group")
save_data(contrib, "tier_contributions")

# Prevalence by data scenario (observed only / + individual / all tiers)
scenario_prev <- map_dfr(list("Observed only" = "Observed",
                              "Observed + individual" = tier_levels[1:2],
                              "All three tiers" = tier_levels), function(tiers) {
  d <- all_tiers %>% filter(tier %in% tiers, !is.na(pop_15plus))
  bind_rows(d %>% group_by(m49_region = as.character(m49_region)) %>%
              summarise(n_countries = n(), prevalence = weighted.mean(prevalence, pop_15plus), .groups = "drop"),
            d %>% summarise(m49_region = "Global", n_countries = n(),
                            prevalence = weighted.mean(prevalence, pop_15plus)))
}, .id = "scenario")
save_table(scenario_prev, "B8_prevalence_by_data_scenario")


# ══════════════════════════════════════════════════════════════════════════════
# 5. TABLE S9 (countries with higher estimation uncertainty)
# ══════════════════════════════════════════════════════════════════════════════
# Per methods 4.7.5: (a) predictors outside the training range (B2, B5);
# (b) disagreement between ENET and XGB; (c) for Tier 2, disagreement between
# the individual and country-level models; plus (d) pre-2015 GWP data.
# (e) national predictors filled with the model-region median (A4).
# Disagreement: the n_top countries with the largest disagreement are listed.

extrap_t2 <- read_csv(file.path(dir_tables, "B2_tier2_extrapolation_by_country.csv"), show_col_types = FALSE) %>%
  mutate(country_code = as.character(country_code))
extrap_t3 <- read_csv(file.path(dir_tables, "B5_tier3_extrapolation_by_country.csv"), show_col_types = FALSE) %>%
  rename(country_code = iso3)

uncertainty <- all_tiers %>%
  filter(tier != "Observed") %>%
  select(country_code, country_name, tier, m49_region, prev_best, prev_enet_individual,
         prev_xgb_individual, prev_enet_country, prev_xgb_country, ci_lower, ci_upper) %>%
  mutate(diff_enet_xgb = ifelse(tier == "Predicted (individual)",
                                prev_xgb_individual - prev_enet_individual,
                                prev_xgb_country - prev_enet_country),
         diff_country_vs_individual = ifelse(tier == "Predicted (individual)",
                                             prev_enet_country - prev_enet_individual, NA_real_),
         max_abs_disagreement = pmax(abs(diff_enet_xgb), abs(diff_country_vs_individual), na.rm = TRUE)) %>%
  left_join(bind_rows(extrap_t2, extrap_t3) %>%
              rename(n_vars_outside_training_range = n_extrapolating, vars_outside_training_range = variables),
            by = "country_code") %>%
  left_join(new_data %>% distinct(country_code, pre_2015) %>% mutate(country_code = as.character(country_code)),
            by = "country_code") %>%
  arrange(desc(max_abs_disagreement))
save_table(uncertainty, "B8_uncertainty_flag_inputs")

# ── Table S9 (formatted): countries flagged for higher uncertainty, by reason ──
# Rows: (1) pre-2015 GWP data; (2) each national predictor outside the training
# range (individual model: 70 training countries; country model: 76 observed
# countries), listing every out-of-range value unless min_pct_outside > 0;
# (3) predictors filled with the model-region median; (4) the n_top countries
# with the largest model disagreement (as above), grouped by model. Values are
# those the models used (MICE-imputed, A4).
min_pct_outside <- 10   # a value is flagged when it lies outside the training range by at least
                        # 10% of that range (0 = report every out-of-range value; all are
                        # in B8_outside_training_range_all.csv with pct_outside)
n_top           <- 15

var_labels <- c(log_gni_pc_ppp = "GNI per capita (log)", log_water_stress = "Water stress (log)",
                log_pop_density = "Population density (log)", basic_water = "Basic drinking water access",
                life_exp = "Life expectancy", fertility = "Fertility rate", cpi_score = "Corruption score",
                u5_mort = "Under-5 mortality", wash_mort = "WASH mortality", precip_mm = "Precipitation",
                cvi = "Climate vulnerability")

# (1) Pre-2015 survey data
pre <- read_csv(file.path(dir_tables, "A1_new_data_pre2015_countries.csv"), show_col_types = FALSE) %>%
  mutate(country_code = as.character(country_code)) %>%
  arrange(YEAR_WAVE) %>%
  left_join(uncertainty %>% distinct(country_code, country_name), by = "country_code")
rows_pre <- tibble(Reason  = "Gallup World Poll data from before 2015",
                   Country = paste(pre$country_name, collapse = "; "),
                   Values  = paste0(pre$country_code, "-", pre$YEAR_WAVE, collapse = "; "))

# (2) National predictors outside the training range
ind_ranges <- load_model("ind_settings")$numeric_ranges %>%
  filter(variable %in% ctry_name_map) %>%
  transmute(model = "Individual", wb_var = names(ctry_name_map)[match(variable, ctry_name_map)], min, max)
ctry_ranges <- load_data("country_df_oos") %>%
  summarise(across(all_of(country_vars), list(min = ~ min(.x, na.rm = TRUE), max = ~ max(.x, na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("wb_var", ".value"), names_pattern = "(.*)_(min|max)$") %>%
  mutate(model = "Country")

extrap <- uncertainty %>%
  distinct(country_code, country_name, tier) %>%
  mutate(model = ifelse(tier == "Predicted (individual)", "Individual", "Country")) %>%
  left_join(load_data("wb_imputed") %>% select(country_code = iso3, all_of(country_vars)) %>%
              pivot_longer(-country_code, names_to = "wb_var", values_to = "value"), by = "country_code") %>%
  inner_join(bind_rows(ind_ranges, ctry_ranges), by = c("model", "wb_var")) %>%
  mutate(direction   = case_when(value > max ~ "Higher", value < min ~ "Lower"),
         pct_outside = 100 * pmax(value - max, min - value, 0) / (max - min)) %>%   # % of training range
  filter(!is.na(direction))
save_table(extrap %>% arrange(desc(pct_outside)), "B8_outside_training_range_all")

# CHECK: same countries/counts as the B2/B5 extrapolation checks (expect 0 rows)
extrap_mismatch <- extrap %>% count(country_code, name = "n_here") %>%
  full_join(uncertainty %>% select(country_code, n_vars_outside_training_range), by = "country_code") %>%
  filter(coalesce(n_here, 0L) != coalesce(n_vars_outside_training_range, 0))
message("Table S9 extrapolation check - mismatches with B2/B5 (expect 0): ", nrow(extrap_mismatch))

rows_extrap <- extrap %>%
  filter(pct_outside >= min_pct_outside) %>%
  arrange(desc(pct_outside)) %>%
  group_by(model, wb_var, direction) %>%
  summarise(Country = paste0(country_name, " (", country_code, ")", collapse = "; "),
            Values  = paste0("Training data ", ifelse(first(direction) == "Higher", "max ", "min "),
                             round(ifelse(first(direction) == "Higher", first(max), first(min)), 2), "; ",
                             paste0(country_code, " ", round(value, 2), collapse = "; ")),
            .groups = "drop") %>%
  mutate(Reason = paste0(direction, " ", var_labels[wb_var], " than training data (",
                         tolower(model), "-level model)")) %>%
  arrange(desc(model), Reason) %>%
  select(Reason, Country, Values)

# (3) Greatest model disagreement: top n_top overall, grouped by model
rows_div <- uncertainty %>%
  slice_max(max_abs_disagreement, n = n_top, with_ties = FALSE) %>%
  arrange(desc(max_abs_disagreement)) %>%
  group_by(tier) %>%
  summarise(Country = paste0(country_name, " (", country_code, ")", collapse = "; "),
            Values  = paste0(country_code, " ", round(100 * max_abs_disagreement, 1), " pp", collapse = "; "),
            .groups = "drop") %>%
  transmute(Reason = paste0("Top ", n_top, " greatest model disagreement (",
                            ifelse(tier == "Predicted (individual)",
                                   "individual-level model: ENET vs XGBoost, or individual vs country-level",
                                   "country-level model: ENET vs XGBoost"), ")"),
            Country, Values)

# (4) National predictors filled with the model-region median in A4 (not observed):
#     gaps in water stress, precipitation, basic water and under-5 mortality
#     (not imputed by MICE), plus corruption / WASH mortality / climate
#     vulnerability / GNI for countries MICE could not impute because one of
#     those predictors was missing. These values lie inside the training range,
#     so they are not caught by (2).
min_n_regional_median <- 1   # e.g. 2 = only list countries with 2+ such predictors
wb_raw_flags <- wb %>% select(country_code = iso3, all_of(c(ctry_vars_wb, ctry_aux_wb)))
regional_median <- wb_raw_flags %>%
  filter(country_code %in% uncertainty$country_code) %>%
  pivot_longer(-country_code, names_to = "wb_var", values_to = "value") %>%
  group_by(country_code) %>%
  mutate(any_aux_missing = any(is.na(value) & wb_var %in% ctry_aux_wb)) %>%
  ungroup() %>%
  # aux gaps are always median-filled; target gaps only when an aux was missing
  filter(is.na(value), wb_var %in% ctry_aux_wb | any_aux_missing) %>%
  left_join(uncertainty %>% distinct(country_code, country_name), by = "country_code") %>%
  group_by(country_code, country_name) %>%
  summarise(n_vars = n(), vars = paste(var_labels[wb_var], collapse = ", "), .groups = "drop") %>%
  filter(n_vars >= min_n_regional_median) %>%
  arrange(desc(n_vars), country_code)
save_table(regional_median, "B8_predictors_regional_median_by_country")
rows_median <- if (nrow(regional_median) > 0) tibble(
  Reason  = "National predictors imputed with the model-region median (not observed)",
  Country = paste0(regional_median$country_name, " (", regional_median$country_code, ")", collapse = "; "),
  Values  = paste0(regional_median$country_code, ": ", regional_median$vars, collapse = "; ")) else NULL

# (5) Individual-level predictors with NO observed values in a Tier 2 country
#     (A4_vars_fully_missing_by_country.csv): MICE fills the whole variable from
#     respondents in other countries, and the models use one imputation, so
#     these estimates depend more on the imputation (see B2b for the spread
#     across the 5 imputations).
ind_var_labels <- c(age = "age", INCOME_5 = "income quintile", inadincome = "income adequacy",
                    education = "education", waterdissatisfied = "water dissatisfaction",
                    maritalstatus = "marital status", urban_final = "urbanicity",
                    employment = "employment", hhsize = "household size")
f_fully_missing <- file.path(dir_tables, "A4_vars_fully_missing_by_country.csv")
fully_imputed <- if (file.exists(f_fully_missing)) {
  read_csv(f_fully_missing, show_col_types = FALSE) %>%
    mutate(country_code = as.character(country_code)) %>%
    filter(country_code %in% uncertainty$country_code[uncertainty$tier == "Predicted (individual)"]) %>%
    left_join(uncertainty %>% distinct(country_code, country_name), by = "country_code") %>%
    group_by(country_code, country_name) %>%
    summarise(vars = paste(coalesce(ind_var_labels[variable], variable), collapse = ", "), .groups = "drop") %>%
    arrange(country_code)
} else {
  message("A4_vars_fully_missing_by_country.csv not found (run A4): Table S9 row on fully imputed variables skipped")
  tibble(country_code = character(), country_name = character(), vars = character())
}
save_table(fully_imputed, "B8_individual_vars_fully_imputed_by_country")
rows_imputed <- if (nrow(fully_imputed) > 0) tibble(
  Reason  = "Individual-level predictors with no observed values in the country (entirely imputed)",
  Country = paste0(fully_imputed$country_name, " (", fully_imputed$country_code, ")", collapse = "; "),
  Values  = paste0(fully_imputed$country_code, ": ", fully_imputed$vars, collapse = "; ")) else NULL

# Table note (last row), so the table carries its own definition of the flag
rows_note <- tibble(
  Reason  = "Note",
  Country = paste0("A national predictor is flagged as outside the training range when its value lies ",
                   "beyond the minimum or maximum of the training data by at least ", min_pct_outside,
                   "% of the training range (training data: 70 countries for the individual-level model, ",
                   "76 observed countries for the country-level model)."),
  Values  = "")

tableS9_flags <- bind_rows(rows_pre, rows_extrap, rows_median, rows_imputed, rows_div, rows_note)
print(tableS9_flags, n = Inf, width = Inf)
readr::write_excel_csv(tableS9_flags, file.path(dir_tables, "B8_tableS09_flagged_countries.csv"))   # Excel-safe (dashes)

# Per-country summary of the Table S9 reasons (used to flag countries in Table S10, B11)
top_div <- uncertainty %>% slice_max(max_abs_disagreement, n = n_top, with_ties = FALSE)
country_flags <- bind_rows(
  pre %>% transmute(country_code, reason = paste0("Pre-2015 GWP data (", YEAR_WAVE, ")")),
  extrap %>% filter(pct_outside >= min_pct_outside) %>%
    group_by(country_code) %>%
    summarise(reason = paste0("Outside training range: ", paste(var_labels[wb_var], collapse = ", ")), .groups = "drop"),
  regional_median %>% transmute(country_code, reason = paste0("Regional-median predictors: ", vars)),
  fully_imputed %>% transmute(country_code, reason = paste0("Entirely imputed individual predictors: ", vars)),
  top_div %>% transmute(country_code, reason = paste0("Top ", n_top, " model disagreement (",
                                                      round(100 * max_abs_disagreement, 1), " pp)"))
) %>%
  group_by(country_code) %>%
  summarise(table_S9_flags = paste(reason, collapse = "; "), .groups = "drop")
save_table(country_flags, "B8_country_flags")
save_data(country_flags, "country_flags")

stop_timer()
