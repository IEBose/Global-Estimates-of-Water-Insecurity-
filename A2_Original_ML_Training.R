################################################################################
#  A2_Original_ML_Training.R — ORIGINAL ML TRAINING (SINGLE IMPUTATION)
#
#  OPTIONAL — FOR INTEREST ONLY. Not needed to reproduce the paper, not part of
#  the main pipeline run time, and no main-pipeline script reads its outputs.
#
#  Purpose : The original version of the individual-level model training, using
#            SIMPLE (single) imputation instead of MICE: country medians/modes,
#            then regional, then global. Otherwise the same models and
#            validation as B1 (ENET, RF, XGB, mixed effects; India + 5 country
#            hold-out; 5-fold LCO CV; no-water-dissatisfaction sensitivity).
#            Used for comparison with the MICE-based results (see A3).
#
#  Depends : A1 only
#  Inputs  : OUTPUT/data/df_clean.rds
#  Outputs : OUTPUT/tables/A2_*.csv, OUTPUT/figures/A2_*.png,
#            OUTPUT/data/A2_objects.RData (objects needed by A3)
#
#  Run time: ~60-90 min (RF and mixed models in 5-fold LCO CV are the slowest)
#
#  NOTE: data loading, model regions, urbanicity and factor coding are now done
#  once in A1 (identical code), so they have been removed from this script.
#  Because A1 removes the rows with a missing outcome BEFORE imputation, the
#  simple-imputation medians/modes here are based on slightly fewer rows than
#  in the original version (426 respondents).
#
#  HOW TO READ THE METRICS: see the header of B1 (AUC-ROC, AUC-PR, Brier,
#  country-level R² = squared correlation, MAE).
################################################################################

source("00_setup.R")
start_timer("A2")
load_pkgs(c("ranger", "lme4", "shapviz"))
select <- dplyr::select
filter <- dplyr::filter

df <- load_data("df_clean")

#Imputation


get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

#individual-level vars
# Missingness flags (income_missing, urban_missing) are created in A1

cat("\n  Countries entirely missing INCOME_5:\n")
df %>%
  group_by(country_code) %>%
  summarise(pct_missing = mean(is.na(INCOME_5)), .groups = "drop") %>%
  filter(pct_missing == 1) %>%
  print()



# Individual-level variables with missingness
df <- df %>%
  mutate(
    age_imp               = age,
    INCOME_5_imp          = INCOME_5,
    inadincome_imp        = inadincome,
    education_imp         = education,
    maritalstatus_imp     = maritalstatus,
    waterdissatisfied_imp = waterdissatisfied,
    urban_final_imp       = urban_final,
    log_GNI_imp           = log_GNI,
    Corruption_imp        = Corruption,
    Mort_WASH_imp         = wash_mort,
    CVI_imp               = cvi
  )


##############################################
# Country-level imputation (region median)
# Only 1 NA each for log_GNI, Corruption, Mort_WASH, CVI
country_vars_imp_only <- c("log_GNI_imp", "Corruption_imp", "Mort_WASH_imp", "CVI_imp")

df <- df %>%
  group_by(region_final) %>%
  mutate(across(all_of(country_vars_imp_only),
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  ungroup() %>%
  mutate(across(all_of(country_vars_imp_only),
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

#for individual vars use country then 
# age_imp: country -> region 

df <- df %>%
  group_by(country_code) %>%
  mutate(age_imp = ifelse(is.na(age_imp), median(age_imp, na.rm = TRUE), age_imp)) %>%
  ungroup() %>%
  group_by(region_final) %>%
  mutate(age_imp = ifelse(is.na(age_imp), median(age_imp, na.rm = TRUE), age_imp)) %>%
  ungroup() %>%
  mutate(age_imp = ifelse(is.na(age_imp), median(age_imp, na.rm = TRUE), age_imp))

# INCOME_5_imp: country -> region 
df <- df %>%
  mutate(income_num = as.numeric(as.character(INCOME_5_imp)))

df <- df %>%
  group_by(country_code) %>%
  mutate(income_num = ifelse(is.na(income_num),
                             median(income_num, na.rm = TRUE), income_num)) %>%
  ungroup() %>%
  group_by(region_final) %>%
  mutate(income_num = ifelse(is.na(income_num),
                             median(income_num, na.rm = TRUE), income_num)) %>%
  ungroup() %>%
  mutate(income_num = ifelse(is.na(income_num),
                             median(income_num, na.rm = TRUE), income_num)) %>%
  mutate(
    income_num   = pmin(pmax(round(income_num), 1), 5),
    INCOME_5_imp = factor(income_num, levels = 1:5)
  ) %>%
  select(-income_num)


#  Categorical: only vars with missingness
# inadincome, education, waterdissatisfied, maritalstatus, urban_final: country -> region -> global

df <- df %>%
  group_by(country_code) %>%
  mutate(
    inadincome_imp        = ifelse(is.na(inadincome_imp),
                                   get_mode(inadincome_imp), as.character(inadincome_imp)),
    education_imp         = ifelse(is.na(education_imp),
                                   get_mode(education_imp), as.character(education_imp)),
    waterdissatisfied_imp = ifelse(is.na(waterdissatisfied_imp),
                                   get_mode(waterdissatisfied_imp), as.character(waterdissatisfied_imp)),
    maritalstatus_imp     = ifelse(is.na(maritalstatus_imp),
                                   get_mode(maritalstatus_imp), as.character(maritalstatus_imp)),
    urban_final_imp       = ifelse(is.na(urban_final_imp),
                                   get_mode(urban_final_imp), as.character(urban_final_imp))
  ) %>%
  ungroup() %>%
  group_by(region_final) %>%
  mutate(
    inadincome_imp        = ifelse(is.na(inadincome_imp),
                                   get_mode(inadincome_imp), inadincome_imp),
    education_imp         = ifelse(is.na(education_imp),
                                   get_mode(education_imp), education_imp),
    waterdissatisfied_imp = ifelse(is.na(waterdissatisfied_imp),
                                   get_mode(waterdissatisfied_imp), waterdissatisfied_imp),
    maritalstatus_imp     = ifelse(is.na(maritalstatus_imp),
                                   get_mode(maritalstatus_imp), maritalstatus_imp),
    urban_final_imp       = ifelse(is.na(urban_final_imp),
                                   get_mode(urban_final_imp), urban_final_imp)
  ) %>%
  ungroup() %>%
  mutate(
    inadincome_imp        = factor(inadincome_imp),
    education_imp         = factor(education_imp,
                                   levels = c("Elementary", "Secondary", "College")),
    maritalstatus_imp     = factor(maritalstatus_imp,
                                   levels = c("Single", "Married/Partnered",
                                              "Separated/Widowed")),
    waterdissatisfied_imp = factor(waterdissatisfied_imp),
    urban_final_imp       = factor(urban_final_imp,
                                   levels = c("Rural", "Towns/semi-dense", "Cities"))
  )

#remove missing outcome data
n_before <- nrow(df)
df <- df %>% filter(!is.na(iwise12_imp))
cat("  Removed", n_before - nrow(df), "rows with missing outcome.\n") #426


# Verification
individual_vars  <- c("age_imp", "hhsize", "female", "INCOME_5_imp",
                      "inadincome_imp", "education_imp", "maritalstatus_imp",
                      "employment", "waterdissatisfied_imp",
                      "urban_final_imp", "urban_source")
missingness_vars <- c("income_missing", "urban_missing")
country_vars     <- c("log_GNI_imp", "log_water_stress", "log_pop_density",
                      "Drinking_basic", "life_exp", "fertility",
                      "Corruption_imp", "u5_mort", "Mort_WASH_imp",
                      "precip_mm", "CVI_imp")

all_check_vars <- c(outcome_var, weight_var,
                    individual_vars, missingness_vars, country_vars, "region_final")

remaining_na <- df %>%
  select(all_of(all_check_vars)) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0)

if (nrow(remaining_na) == 0) {
  cat("  None - all imputation succeeded.\n")
} else {
  print(remaining_na)
}

# Imputation summary

cat("\n  Imputation summary:\n")

tableS1<- tibble(
  variable = c("age", "INCOME_5", "inadincome", "education",
               "maritalstatus", "waterdissatisfied", "urban_final",
               "log_GNI", "Corruption", "wash_mort", "cvi"),
  n_imputed = c(
    sum(is.na(df$age)), sum(is.na(df$INCOME_5)),
    sum(is.na(df$inadincome)), sum(is.na(df$education)),
    sum(is.na(df$maritalstatus)), sum(is.na(df$waterdissatisfied)),
    sum(is.na(df$urban_final)), sum(is.na(df$log_GNI)),
    sum(is.na(df$Corruption)), sum(is.na(df$wash_mort)),
    sum(is.na(df$cvi))
  ),
  n_total     = nrow(df),
  pct_imputed = round(n_imputed / n_total * 100, 2)
)
print(tableS1)



write_csv(tableS1, file.path(dir_tables, "A2_tableS1_Imputation_summary.csv"))

#summary of predictors

all_predictors <- c(individual_vars, missingness_vars, country_vars, "region_final")

cat("  Individual (", length(individual_vars), "):",
    paste(individual_vars, collapse = ", "), "\n")
cat("  Missingness (", length(missingness_vars), "):",
    paste(missingness_vars, collapse = ", "), "\n")
cat("  Country (", length(country_vars), "):",
    paste(country_vars, collapse = ", "), "\n")
cat("  Region: region_final\n")
cat("  Total:", length(all_predictors), "predictors\n")

#Individual ( 11 ): age_imp, hhsize, female, INCOME_5_imp, inadincome_imp, education_imp, maritalstatus_imp, employment, waterdissatisfied_imp, urban_final_imp, urban_source 
#Missingness ( 2 ): income_missing, urban_missing 
#Country ( 11 ): log_GNI_imp, log_water_stress, log_pop_density, Drinking_basic, life_exp, fertility, Corruption_imp, u5_mort, Mort_WASH_imp, precip_mm, CVI_imp 
#Region: region_final
#Total: 25 predictors



table1_region <- df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(region_final) %>%
  summarise(
    n_countries     = n_distinct(country_code),
    n_obs           = n(),
    prevalence      = round(weighted.mean(outcome_num, wgtnorm_using, na.rm = TRUE), 3),
    mean_age        = round(weighted.mean(age, wgtnorm_using, na.rm = TRUE), 1),
    pct_female      = round(weighted.mean(as.numeric(female == "Female"),
                                          wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_urban       = round(weighted.mean(
      as.numeric(urban_final_imp == "Cities"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_difficult   = round(weighted.mean(
      as.numeric(inadincome == "Difficult"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_waterdissat = round(weighted.mean(
      as.numeric(waterdissatisfied == "Yes"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(prevalence))

# Add global total row
table1_global <- df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  summarise(
    region_final    = "GLOBAL",
    n_countries     = n_distinct(country_code),
    n_obs           = n(),
    prevalence      = round(weighted.mean(outcome_num, wgtnorm_using, na.rm = TRUE), 3),
    mean_age        = round(weighted.mean(age, wgtnorm_using, na.rm = TRUE), 1),
    pct_female      = round(weighted.mean(as.numeric(female == "Female"),
                                          wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_urban       = round(weighted.mean(
      as.numeric(urban_final_imp == "Cities"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_difficult   = round(weighted.mean(
      as.numeric(inadincome == "Difficult"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_waterdissat = round(weighted.mean(
      as.numeric(waterdissatisfied == "Yes"),
      wgtnorm_using, na.rm = TRUE) * 100, 1)
  )

table1 <- bind_rows(table1_region, table1_global)
print(table1, n = Inf, width = Inf)
write_csv(table1, file.path(dir_tables, "A2_table1_sample_characteristics.csv"))

##broader_region in m49

table1C_region <- df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(m49_region) %>%
  summarise(
    n_countries     = n_distinct(country_code),
    n_obs           = n(),
    prevalence      = round(weighted.mean(outcome_num, wgtnorm_using, na.rm = TRUE), 3),
    mean_age        = round(weighted.mean(age, wgtnorm_using, na.rm = TRUE), 1),
    pct_female      = round(weighted.mean(as.numeric(female == "Female"),
                                          wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_urban       = round(weighted.mean(
      as.numeric(urban_final_imp == "Cities"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_difficult   = round(weighted.mean(
      as.numeric(inadincome == "Difficult"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    pct_waterdissat = round(weighted.mean(
      as.numeric(waterdissatisfied == "Yes"),
      wgtnorm_using, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(prevalence))


print(table1C_region, n = Inf, width = Inf)
write_csv(table1C_region, file.path(dir_tables, "A2_table1_sample_characteristics_m49region.csv"))

#################
# For each variable, check which level of imputation was used
check_imputation_level <- function(data, var_name, country_col = "country_code", 
                                   region_col = "region_final") {
  data %>%
    group_by(.data[[country_col]]) %>%
    mutate(
      country_has_data = any(!is.na(.data[[var_name]]))
    ) %>%
    ungroup() %>%
    group_by(.data[[region_col]]) %>%
    mutate(
      region_has_data = any(!is.na(.data[[var_name]]))
    ) %>%
    ungroup() %>%
    filter(is.na(.data[[var_name]])) %>%
    mutate(
      imputation_level = case_when(
        country_has_data  ~ "Country",
        region_has_data   ~ "Region",
        TRUE              ~ "Global"
      )
    ) %>%
    group_by(imputation_level) %>%
    summarise(
      n_obs       = n(),
      n_countries = n_distinct(.data[[country_col]]),
      .groups = "drop"
    ) %>%
    mutate(variable = var_name)
}
########
#to check the names of the imputed countries
check_imputation_level <- function(data, var_name, country_col = "country_code", 
                                   region_col = "region_final") {
  missing_rows <- data %>%
    group_by(.data[[country_col]]) %>%
    mutate(country_has_data = any(!is.na(.data[[var_name]]))) %>%
    ungroup() %>%
    group_by(.data[[region_col]]) %>%
    mutate(region_has_data = any(!is.na(.data[[var_name]]))) %>%
    ungroup() %>%
    filter(is.na(.data[[var_name]])) %>%
    mutate(
      imputation_level = as.character(case_when(
        country_has_data  ~ "Country",
        region_has_data   ~ "Region",
        TRUE              ~ "Global"
      ))
    )
  
  if (nrow(missing_rows) == 0) return(NULL)
  
  missing_rows %>%
    group_by(imputation_level) %>%
    summarise(
      n_obs       = n(),
      n_countries = n_distinct(.data[[country_col]]),
      countries   = paste(unique(.data[[country_col]]), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(variable = var_name)
}

###################
# Check all imputed variables (use ORIGINALS, not _imp)
vars_to_check <- c("age", "INCOME_5", "inadincome", "education",
                   "maritalstatus", "waterdissatisfied", "urban_final")

imputation_levels <- lapply(vars_to_check, function(v) {
  check_imputation_level(df, v)
}) %>% bind_rows()

cat("Individual-level imputation sources:\n")
print(imputation_levels, n = Inf)
#Individual-level imputation sources:
 # imputation_level n_obs n_countries variable         
# Country-level variables
country_vars_check <- c("log_GNI", "Corruption", "wash_mort", "cvi")

country_imp_levels <- lapply(country_vars_check, function(v) {
  df %>%
    distinct(country_code, region_final, .keep_all = TRUE) %>%
    check_imputation_level(v)
}) %>% bind_rows()

cat("\nCountry-level imputation sources:\n")
print(country_imp_levels, n = Inf)
#
#Country-level imputation sources:
#imputation_level n_obs n_countries variable  

#############TRAINING DATA

focal_holdout_country <- "IND"

other_countries   <- setdiff(levels(df$country_code), focal_holdout_country)
n_holdout_other   <- 5
set.seed(241189)  # ensure same holdout countries every run
holdout_other_set <- sample(other_countries, n_holdout_other)

cat("  Focal:", focal_holdout_country, "\n")
cat("  Other:", paste(holdout_other_set, collapse = ", "), "\n")

df <- df %>%
  mutate(split = case_when(
    country_code == focal_holdout_country  ~ "holdout_focal",
    country_code %in% holdout_other_set    ~ "holdout_other",
    TRUE                                   ~ "train"
  ))

train_df      <- df %>% filter(split == "train")
test_focal_df <- df %>% filter(split == "holdout_focal")
test_other_df <- df %>% filter(split == "holdout_other")

cat("  Train:        ", nrow(train_df), "(", n_distinct(train_df$country_code), "countries)\n")
cat("  Test (focal): ", nrow(test_focal_df), "(", n_distinct(test_focal_df$country_code), "countries)\n")
cat("  Test (other): ", nrow(test_other_df), "(", n_distinct(test_other_df$country_code), "countries)\n")


####pre-processing

rec <- recipe(iwise12_imp ~ .,
              data = train_df %>%
                select(all_of(c("iwise12_imp", all_predictors))) %>%
                mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  prep(training = train_df %>%
         mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))

train_baked      <- bake(rec, train_df %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
test_focal_baked <- bake(rec, test_focal_df %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
test_other_baked <- bake(rec, test_other_df %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))

y_train      <- train_baked$iwise12_imp
X_train      <- train_baked %>% select(-iwise12_imp) %>% as.matrix()
w_train      <- train_df[[weight_var]]

y_test_focal <- test_focal_baked$iwise12_imp
X_test_focal <- test_focal_baked %>% select(-iwise12_imp) %>% as.matrix()

y_test_other <- test_other_baked$iwise12_imp
X_test_other <- test_other_baked %>% select(-iwise12_imp) %>% as.matrix()


####EVALUATION METRICS

evaluate_model <- function(y_true, y_prob, model_name, test_name) {
  roc_obj <- roc(y_true, y_prob, quiet = TRUE)
  auc_val <- auc(roc_obj)
  coords  <- coords(roc_obj, "best", ret = c("threshold", "sensitivity",
                                             "specificity"))
  y_pred  <- ifelse(y_prob >= coords$threshold, 1, 0)
  cm <- confusionMatrix(factor(y_pred, levels = c("0","1")),
                        factor(y_true, levels = c("0","1")),
                        positive = "1")
  pr <- pr.curve(scores.class0 = y_prob[y_true == 1],
                 scores.class1 = y_prob[y_true == 0],
                 curve = FALSE)
  tibble(
    model       = model_name,
    test_set    = test_name,
    AUC_ROC     = round(as.numeric(auc_val), 4),
    AUC_PR      = round(pr$auc.integral, 4),
    Accuracy    = round(cm$overall["Accuracy"], 4),
    Sensitivity = round(cm$byClass["Sensitivity"], 4),
    Specificity = round(cm$byClass["Specificity"], 4),
    F1          = round(cm$byClass["F1"], 4),
    Threshold   = round(coords$threshold, 4)
  )
}


# MODEL 1: ELASTIC NET


best_alpha  <- 0.5

cv_fit <- cv.glmnet(X_train, y_train,
                    family   = "binomial",
                    weights  = w_train,
                    alpha    = best_alpha,
                    nfolds   = 5,
                    type.measure = "auc")

best_lambda <- cv_fit$lambda.min
cat("  Alpha:", best_alpha, "| Lambda:", round(best_lambda, 6), "\n")
cat("  CV AUC:", round(max(cv_fit$cvm), 4), "\n")

enet_final <- glmnet(X_train, y_train,
                     family  = "binomial",
                     weights = w_train,
                     alpha   = best_alpha,
                     lambda  = best_lambda)

enet_prob_focal <- predict(enet_final, X_test_focal, type = "response")[,1]
enet_prob_other <- predict(enet_final, X_test_other, type = "response")[,1]

# Coefficients
enet_coefs <- coef(enet_final) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>%
  rename(coefficient = s0) %>%
  filter(coefficient != 0, feature != "(Intercept)") %>%
  arrange(desc(abs(coefficient)))

cat("  Non-zero features:", nrow(enet_coefs), "\n")
cat("  Top 10:\n")
print(head(enet_coefs, 10))

# Coefficient plot
enet_coefs %>%
  head(15) %>%
  mutate(
    feature = fct_reorder(feature, abs(coefficient)),
    direction = ifelse(coefficient > 0, "Increases risk", "Decreases risk")
  ) %>%
  ggplot(aes(x = coefficient, y = feature, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Increases risk" = "#E41A1C",
                               "Decreases risk" = "#377EB8")) +
  labs(title = "Elastic Net - Top 15 Coefficients",
       x = "Coefficient (log-odds scale)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(dir_figures, "A2_enet_coefficients.png"), width = 9, height = 7, dpi = 300)

# Odds ratio plot 
enet_coefs %>%
  head(15) %>%
  mutate(
    odds_ratio = exp(coefficient),
    feature    = fct_reorder(feature, abs(coefficient))
  ) %>%
  ggplot(aes(x = odds_ratio, y = feature)) +
  geom_point(size = 3, colour = "#377EB8") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "Elastic Net - Odds Ratios (Top 15 Predictors)",
       subtitle = "Per 1-SD change in standardised predictors",
       x = "Odds Ratio", y = NULL) +
  theme_minimal()

ggsave(file.path(dir_figures, "A2_enet_odds_ratios.png"), width = 9, height = 7, dpi = 300)


# Lookup table for readable labels
label_lookup <- c(
  # Individual-level
  "waterdissatisfied_imp_Yes"       = "Water dissatisfied",
  "inadincome_imp_Difficult"        = "Income inadequate",
  "age_imp"                         = "Age",
  "hhsize"                          = "Household size",
  "female_Female"                   = "Female",
  "INCOME_5_imp_X2"                 = "Income quintile 2",
  "INCOME_5_imp_X3"                 = "Income quintile 3",
  "INCOME_5_imp_X4"                 = "Income quintile 4",
  "INCOME_5_imp_X5"                 = "Income quintile 5 (richest)",
  "education_imp_Secondary"         = "Secondary education",
  "education_imp_College"           = "College education",
  "maritalstatus_imp_Married.Partnered"  = "Married/Partnered",
  "maritalstatus_imp_Separated.Widowed"  = "Separated/Widowed",
  "employment_Underemployed"        = "Underemployed",
  "employment_Unemployed"           = "Unemployed",
  "employment_Out.of.workforce"     = "Out of workforce",
  "urban_final_imp_Towns.semi.dense" = "Towns/semi-dense",
  "urban_final_imp_Cities"          = "Cities",
  "urban_source_missing"            = "Urbanicity missing",
  "urban_source_selfrep"            = "Urbanicity self-reported",
  "income_missing_X1"               = "Income missing",
  "urban_missing_X1"                = "Urbanicity missing flag",
  "inadincome_imp_Getting.by"       = "Getting by on income",
  
  # Country-level
  "log_GNI_imp"                     = "GNI per capita (log)",
  "log_water_stress"                = "Water stress (log)",
  "log_pop_density"                 = "Population density (log)",
  "Drinking_basic"                  = "Basic drinking water (%)",
  "life_exp"                     = "Life expectancy",
  "fertility"                  = "Fertility rate",
  "Corruption_imp"                  = "Corruption index (low score= high corruption) ",
  "u5_mort"                         = "Under-5 mortality",
  "Mort_WASH_imp"                   = "WASH mortality",
  "precip_mm"                            = "Mean annual rainfall",
  "CVI_imp"                         = "Climate vulnerability Index (high score= low vulnerability)",
  
  # Region dummies
  "region_final_Eastern.Africa"               = "Eastern Africa",
  "region_final_Latin.America"                = "Latin America",
  "region_final_MENA"                         = "MENA",
  "region_final_Middle.Africa"                = "Middle Africa",
  "region_final_N..America..Europe...Australia" = "N. America/Europe/Australia",
  "region_final_Southern.Africa"              = "Southern Africa",
  "region_final_Southern.Asia"                = "Southern Asia",
  "region_final_Western.Africa"               = "Western Africa"
)

# Helper function to relabel
relabel <- function(x) {
  ifelse(x %in% names(label_lookup), label_lookup[x], x)
}

p_coef <- enet_coefs %>%
  head(15) %>%
  mutate(
    feature_label = relabel(feature),
    feature_label = fct_reorder(feature_label, abs(coefficient)),
    direction = ifelse(coefficient > 0, "Increases risk", "Decreases risk")
  ) %>%
  ggplot(aes(x = coefficient, y = feature_label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Increases risk" = "darkred",
                               "Decreases risk" = "darkblue")) +
  labs(title = "Elastic Net coefficients",
       x = "Coefficient (log-odds scale)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_coef)

# Full coefficients with odds ratios for export
enet_coefs_full <- coef(enet_final) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>%
  rename(coefficient = s0) %>%
  filter(feature != "(Intercept)") %>%
  mutate(odds_ratio = round(exp(coefficient), 4),
         coefficient = round(coefficient, 4)) %>%
  arrange(desc(abs(coefficient)))

write_csv(enet_coefs_full, file.path(dir_tables, "A2_enet_coefficients_full.csv"))

# Try a few alpha values with the fast approach
for (a in c(0.25, 0.5, 0.75)) {
  cv_fit <- cv.glmnet(X_train, y_train,
                      family = "binomial",
                      weights = w_train,
                      alpha = a,
                      nfolds = 5,
                      type.measure = "auc")
  cat("Alpha:", a, "| Best CV AUC:", round(max(cv_fit$cvm), 4), "\n")
}
# Sensitivity: alpha- Note this is quite slow
cat("\n  Alpha sensitivity:\n")
alpha_sensitivity <- lapply(c(0, 0.25, 0.5, 0.75, 1), function(a) {
  cv_tmp <- cv.glmnet(X_train, y_train,
                      family = "binomial", weights = w_train,
                      alpha = a, nfolds = 5, type.measure = "auc")
  enet_tmp <- glmnet(X_train, y_train, family = "binomial",
                     weights = w_train, alpha = a, lambda = cv_tmp$lambda.min)
  p_focal <- predict(enet_tmp, X_test_focal, type = "response")[,1]
  p_other <- predict(enet_tmp, X_test_other, type = "response")[,1]
  tibble(
    alpha     = a,
    cv_auc    = round(max(cv_tmp$cvm), 4),
    auc_india = round(as.numeric(auc(roc(y_test_focal, p_focal, quiet = TRUE))), 4),
    auc_other = round(as.numeric(auc(roc(y_test_other, p_other, quiet = TRUE))), 4)
  )
}) %>% bind_rows()
print(alpha_sensitivity)
write_csv(alpha_sensitivity, file.path(dir_tables, "A2_sensitivity_alpha.csv"))


#  MODEL 2: RANDOM FOREST 


rf_train <- train_df %>% select(all_of(c("iwise12_imp", all_predictors)))
rf_test_focal <- test_focal_df %>% select(all_of(c("iwise12_imp", all_predictors)))
rf_test_other <- test_other_df %>% select(all_of(c("iwise12_imp", all_predictors)))

rf_model <- ranger(
  iwise12_imp ~ .,
  data            = rf_train,
  case.weights    = w_train,
  num.trees       = 1000,
  mtry            = floor(sqrt(length(all_predictors))),
  min.node.size   = 10,
  probability     = TRUE,
  importance      = "impurity",
  seed            = 241189,
  num.threads     = parallel::detectCores() - 1
)

cat("  OOB error:", round(rf_model$prediction.error, 4), "\n")

rf_prob_focal <- predict(rf_model, rf_test_focal)$predictions[, "1"]
rf_prob_other <- predict(rf_model, rf_test_other)$predictions[, "1"]

rf_importance <- tibble(
  feature    = names(rf_model$variable.importance),
  importance = rf_model$variable.importance
) %>% arrange(desc(importance))

cat("  Top 10:\n")
print(head(rf_importance, 10))

# Feature importance plot
rf_importance %>%
  head(15) %>%
  mutate(feature = fct_reorder(feature, importance)) %>%
  ggplot(aes(x = importance, y = feature)) +
  geom_col(fill = "#377EB8") +
  labs(title = "Random Forest - Feature Importance (Top 15)",
       x = "Impurity Importance", y = NULL) +
  theme_minimal()

ggsave(file.path(dir_figures, "A2_rf_feature_importance.png"), width = 8, height = 6, dpi = 300)


#MODEL 3: XGBOOST


# Rescale weights to mean = 1 (required for min_child_weight)
w_train_scaled <- w_train / mean(w_train)

dtrain <- xgb.DMatrix(data = X_train, label = y_train, weight = w_train_scaled)
dfocal <- xgb.DMatrix(data = X_test_focal, label = y_test_focal)
dother <- xgb.DMatrix(data = X_test_other, label = y_test_other)

set.seed(241189)

xgb_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  eta              = 0.05,
  max_depth        = 6,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 5,
  gamma            = 1,
  seed    = 241189
)

xgb_cv <- xgb.cv(
  params  = xgb_params,
  data    = dtrain,
  nrounds = 1000,
  nfold   = 5,
  early_stopping_rounds = 30,
  verbose = 0
)

eval_log <- xgb_cv$evaluation_log
best_nrounds <- eval_log$iter[which.max(eval_log$test_auc_mean)]
cat("  Best nrounds:", best_nrounds, "\n")
cat("  Best CV AUC:", round(max(eval_log$test_auc_mean), 4), "\n")

#note without setting seed, best rounds changed so checking how this affected AUC
cat("AUC at round 201:", eval_log$test_auc_mean[201], "\n")

xgb_model <- xgb.train(
  params  = xgb_params,
  data    = dtrain,
  nrounds = best_nrounds,
  verbose = 0
)

xgb_prob_focal <- predict(xgb_model, dfocal)
xgb_prob_other <- predict(xgb_model, dother)

xgb_imp <- xgb.importance(model = xgb_model) %>% as_tibble()
cat("  Top 10:\n")
print(head(xgb_imp %>% select(Feature, Gain, Cover, Frequency), 10))

# Feature importance plot
xgb_imp %>%
  head(15) %>%
  mutate(Feature = fct_reorder(Feature, Gain)) %>%
  ggplot(aes(x = Gain, y = Feature)) +
  geom_col(fill = "#4DAF4A") +
  labs(title = "XGBoost - Feature Importance (Top 15)",
       x = "Gain", y = NULL) +
  theme_minimal()

ggsave(file.path(dir_figures, "A2_xgb_feature_importance.png"), width = 8, height = 6, dpi = 300)


# MODEL 4: MIXED-EFFECTS LOGISTIC REGRESSION 


me_train <- train_df %>%
  mutate(
    outcome_num = as.numeric(as.character(iwise12_imp)),
    wt_scaled   = wgtnorm_using / mean(wgtnorm_using)
  )

# Reduced predictor set for convergence
me_fixed_simple <- paste(c(
  "waterdissatisfied_imp", "inadincome_imp", "age_imp", "hhsize",
  "female", "education_imp", "INCOME_5_imp", "urban_final_imp",
  "Drinking_basic", "u5_mort", "CVI_imp", "log_GNI_imp", "log_pop_density",
  "region_final"
), collapse = " + ")

me_formula <- as.formula(paste("outcome_num ~", me_fixed_simple, "+ (1 | country_code)"))
cat("  Formula:", deparse(me_formula), "\n")

me_model <- suppressWarnings(glmer(
  me_formula,
  data    = me_train,
  family  = binomial(link = "logit"),
  weights = me_train$wt_scaled,
  control = glmerControl(
    optimizer = "nloptwrap",
    optCtrl = list(maxeval = 20000),
    calc.derivs = FALSE
  ),
  nAGQ = 0
))


# ICC
country_var <- as.data.frame(VarCorr(me_model))$vcov[1]
icc <- country_var / (country_var + (pi^2 / 3))
cat("  Country variance:", round(country_var, 4), "\n")
cat("  ICC:", round(icc, 4), "\n")

# Predictions (fixed effects only for held-out countries)
me_prob_focal <- predict(me_model, newdata = test_focal_df,
                         type = "response", re.form = NA,
                         allow.new.levels = TRUE)
me_prob_other <- predict(me_model, newdata = test_other_df,
                         type = "response", re.form = NA,
                         allow.new.levels = TRUE)

# Extract coefficients
me_coefs <- data.frame(
  feature    = rownames(summary(me_model)$coefficients),
  estimate   = summary(me_model)$coefficients[, "Estimate"],
  std_error  = summary(me_model)$coefficients[, "Std. Error"],
  z_value    = summary(me_model)$coefficients[, "z value"],
  p_value    = summary(me_model)$coefficients[, "Pr(>|z|)"],
  odds_ratio = exp(summary(me_model)$coefficients[, "Estimate"])
)
rownames(me_coefs) <- NULL
write_csv(me_coefs, file.path(dir_tables, "A2_mixed_effects_coefficients.csv"))


#  RESULTS 


results <- bind_rows(
  evaluate_model(y_test_focal, enet_prob_focal, "Elastic Net", "Focal (India)"),
  evaluate_model(y_test_other, enet_prob_other, "Elastic Net", "Other countries"),
  evaluate_model(as.numeric(as.character(rf_test_focal$iwise12_imp)),
                 rf_prob_focal, "Random Forest", "Focal (India)"),
  evaluate_model(as.numeric(as.character(rf_test_other$iwise12_imp)),
                 rf_prob_other, "Random Forest", "Other countries"),
  evaluate_model(y_test_focal, xgb_prob_focal, "XGBoost", "Focal (India)"),
  evaluate_model(y_test_other, xgb_prob_other, "XGBoost", "Other countries"),
  evaluate_model(as.numeric(as.character(test_focal_df$iwise12_imp)),
                 me_prob_focal, "Mixed-Effects LR", "Focal (India)"),
  evaluate_model(as.numeric(as.character(test_other_df$iwise12_imp)),
                 me_prob_other, "Mixed-Effects LR", "Other countries")
)

print(results, n = Inf, width = Inf)
write_csv(results, file.path(dir_tables, "A2_model_comparison_results.csv"))

# Brier scores
cat("\n  Brier Scores (India holdout):\n")
cat("  Elastic Net:", round(mean((enet_prob_focal - y_test_focal)^2), 4), "\n")
cat("  Random Forest:", round(mean((rf_prob_focal -
                                      as.numeric(as.character(rf_test_focal$iwise12_imp)))^2), 4), "\n")
cat("  XGBoost:", round(mean((xgb_prob_focal - y_test_focal)^2), 4), "\n")
cat("  Mixed-Effects:", round(mean((me_prob_focal -
                                      as.numeric(as.character(test_focal_df$iwise12_imp)))^2), 4), "\n")


# ROC CURVES


plot_roc_comparison <- function(y_true, prob_list, model_names, title) {
  roc_list <- lapply(seq_along(prob_list), function(i) {
    roc(y_true, prob_list[[i]], quiet = TRUE)
  })
  colours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  plot(roc_list[[1]], col = colours[1], lwd = 2, main = title)
  for (i in 2:length(roc_list)) {
    lines(roc_list[[i]], col = colours[i], lwd = 2)
  }
  legend("bottomright",
         legend = paste0(model_names, " (AUC=",
                         sapply(roc_list, function(r) round(auc(r), 3)), ")"),
         col = colours[1:length(roc_list)], lwd = 2, cex = 0.9)
}


# Plot to screen
par(mfrow = c(1, 2))
plot_roc_comparison(
  y_test_focal,
  list(enet_prob_focal, rf_prob_focal, xgb_prob_focal, me_prob_focal),
  c("Elastic Net", "Random Forest", "XGBoost", "Mixed-Effects LR"),
  "ROC - India Holdout"
)
plot_roc_comparison(
  y_test_other,
  list(enet_prob_other, rf_prob_other, xgb_prob_other, me_prob_other),
  c("Elastic Net", "Random Forest", "XGBoost", "Mixed-Effects LR"),
  "ROC - Other Countries Holdout"
)
par(mfrow = c(1, 1))

#remember to save graph
# ROC curves: India hold-out and other hold-out countries, side by side
png(file.path(dir_figures, "A2_roc_holdouts.png"),
    width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
plot_roc_comparison(
  y_test_focal,
  list(enet_prob_focal, rf_prob_focal, xgb_prob_focal, me_prob_focal),
  c("Elastic Net", "Random Forest", "XGBoost", "Mixed-Effects LR"),
  "ROC - India Holdout"
)
plot_roc_comparison(
  y_test_other,
  list(enet_prob_other, rf_prob_other, xgb_prob_other, me_prob_other),
  c("Elastic Net", "Random Forest", "XGBoost", "Mixed-Effects LR"),
  "ROC - Other Countries Holdout"
)
par(mfrow = c(1, 1))
invisible(dev.off())

# LEAVE-COUNTRY-OUT CROSS-VALIDATION 


train_countries <- unique(as.character(train_df$country_code))
k_folds <- 5
set.seed(241189)  # ensure same folds every run
country_folds <- tibble(country_code = train_countries) %>%
  mutate(fold = sample(rep(1:k_folds, length.out = n())))

train_df <- train_df %>%
  select(-any_of("fold")) %>%
  left_join(country_folds, by = "country_code")

lco_results <- lapply(1:k_folds, function(f) {
  cat("  Fold", f, "of", k_folds, "...\n")
  
  cv_test  <- train_df %>% filter(fold == f)
  cv_train <- train_df %>% filter(fold != f)
  
  cv_rec <- recipe(iwise12_imp ~ .,
                   data = cv_train %>%
                     select(all_of(c("iwise12_imp", all_predictors))) %>%
                     mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = cv_train %>%
           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
  
  cv_train_baked <- bake(cv_rec, cv_train %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
  cv_test_baked  <- bake(cv_rec, cv_test %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
  
  y_cv_train <- cv_train_baked$iwise12_imp
  X_cv_train <- cv_train_baked %>% select(-iwise12_imp) %>% as.matrix()
  w_cv_train <- cv_train[[weight_var]]
  
  y_cv_test  <- cv_test_baked$iwise12_imp
  X_cv_test  <- cv_test_baked %>% select(-iwise12_imp) %>% as.matrix()
  
  # Elastic Net
  cv_enet <- glmnet(X_cv_train, y_cv_train,
                    family = "binomial", weights = w_cv_train,
                    alpha = best_alpha, lambda = best_lambda)
  p_enet  <- predict(cv_enet, X_cv_test, type = "response")[,1]
  
  # XGBoost (rescale weights)
  w_cv_scaled <- w_cv_train / mean(w_cv_train)
  d_cv_train <- xgb.DMatrix(data = X_cv_train, label = y_cv_train,
                            weight = w_cv_scaled)
  d_cv_test  <- xgb.DMatrix(data = X_cv_test, label = y_cv_test)
  cv_xgb     <- xgb.train(params = xgb_params, data = d_cv_train,
                          nrounds = best_nrounds, verbose = 0)
  p_xgb      <- predict(cv_xgb, d_cv_test)
  
  # Random Forest
  rf_train_cv <- cv_train %>% select(all_of(c("iwise12_imp", all_predictors)))
  rf_test_cv  <- cv_test %>% select(all_of(c("iwise12_imp", all_predictors)))
  
  rf_cv <- ranger(
    iwise12_imp ~ .,
    data          = rf_train_cv,
    case.weights  = w_cv_train,
    num.trees     = 1000,
    mtry          = floor(sqrt(length(all_predictors))),
    min.node.size = 10,
    probability   = TRUE,
    seed          = 241189
  )
  p_rf <- predict(rf_cv, rf_test_cv)$predictions[, "1"]
  
  # Mixed-Effects (fixed effects only for prediction)
  me_cv_train <- cv_train %>%
    mutate(outcome_num = as.numeric(as.character(iwise12_imp)),
           wt_scaled   = wgtnorm_using / mean(wgtnorm_using))
  
  me_cv <- tryCatch({
    glmer(
      outcome_num ~ waterdissatisfied_imp + inadincome_imp + age_imp +
        hhsize + female + education_imp + INCOME_5_imp + urban_final_imp +
        Drinking_basic + u5_mort + CVI_imp + log_GNI_imp + log_pop_density +
        region_final + (1 | country_code),
      data    = me_cv_train,
      family  = binomial(link = "logit"),
      weights = me_cv_train$wt_scaled,
      control = glmerControl(optimizer = "nloptwrap",
                             optCtrl = list(maxeval = 20000),
                             calc.derivs = FALSE),
      nAGQ = 0
    )
  }, error = function(e) NULL)
  
  if (!is.null(me_cv)) {
    p_me <- predict(me_cv, newdata = cv_test, type = "response",
                    re.form = NA, allow.new.levels = TRUE)
  } else {
    p_me <- rep(NA, nrow(cv_test))
  }
  bind_rows(
    evaluate_model(y_cv_test, p_enet, "Elastic Net", paste0("Fold_", f)),
    evaluate_model(y_cv_test, p_xgb,  "XGBoost",     paste0("Fold_", f)),
    evaluate_model(as.numeric(as.character(rf_test_cv$iwise12_imp)),
                   p_rf, "Random Forest", paste0("Fold_", f)),
    evaluate_model(as.numeric(as.character(cv_test$iwise12_imp)),
                   p_me, "Mixed-Effects LR", paste0("Fold_", f))
  )
}) %>% bind_rows()

lco_summary <- lco_results %>%
  group_by(model) %>%
  summarise(across(c(AUC_ROC, AUC_PR, Accuracy, Sensitivity, Specificity, F1),
                   list(mean = mean, sd = sd)),
            .groups = "drop")

cat("\n  LCO CV Summary:\n")
print(lco_summary, width = Inf)
write_csv(lco_summary, file.path(dir_tables, "A2_lco_cv_summary_all.csv"))


# OUT-OF-FOLD PREDICTIONS 

train_df$enet_oof_prob <- NA_real_
train_df$xgb_oof_prob  <- NA_real_

for (f in 1:k_folds) {
  cat("  Fold", f, "...\n")
  cv_test_idx  <- which(train_df$fold == f)
  cv_train_idx <- which(train_df$fold != f)
  
  cv_rec_tmp <- recipe(iwise12_imp ~ .,
                       data = train_df[cv_train_idx, ] %>%
                         select(all_of(c("iwise12_imp", all_predictors))) %>%
                         mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = train_df[cv_train_idx, ] %>%
           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
  
  X_tmp_train <- bake(cv_rec_tmp, train_df[cv_train_idx, ] %>%
                        mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    select(-iwise12_imp) %>% as.matrix()
  y_tmp_train <- as.numeric(as.character(train_df$iwise12_imp[cv_train_idx]))
  w_tmp_train <- train_df[[weight_var]][cv_train_idx]
  X_tmp_test  <- bake(cv_rec_tmp, train_df[cv_test_idx, ] %>%
                        mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    select(-iwise12_imp) %>% as.matrix()
  
  # Elastic Net OOF
  cv_enet_tmp <- glmnet(X_tmp_train, y_tmp_train,
                        family = "binomial", weights = w_tmp_train,
                        alpha = best_alpha, lambda = best_lambda)
  train_df$enet_oof_prob[cv_test_idx] <- predict(
    cv_enet_tmp, X_tmp_test, type = "response")[,1]
  
  # XGBoost OOF (rescale weights)
  w_tmp_scaled <- w_tmp_train / mean(w_tmp_train)
  d_tmp <- xgb.DMatrix(data = X_tmp_train, label = y_tmp_train,
                       weight = w_tmp_scaled)
  m_tmp <- xgb.train(params = xgb_params, data = d_tmp,
                     nrounds = best_nrounds, verbose = 0)
  train_df$xgb_oof_prob[cv_test_idx] <- predict(
    m_tmp, xgb.DMatrix(data = X_tmp_test))
}
################
train_df$enet_oof_prob <- NA_real_
train_df$xgb_oof_prob  <- NA_real_
train_df$rf_oof_prob   <- NA_real_
train_df$me_oof_prob   <- NA_real_

for (f in 1:k_folds) {
  cat("  Fold", f, "...\n")
  cv_test_idx  <- which(train_df$fold == f)
  cv_train_idx <- which(train_df$fold != f)
  
  # Recipe (for Elastic Net and XGBoost)
  cv_rec_tmp <- recipe(iwise12_imp ~ .,
                       data = train_df[cv_train_idx, ] %>%
                         select(all_of(c("iwise12_imp", all_predictors))) %>%
                         mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = train_df[cv_train_idx, ] %>%
           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
  
  X_tmp_train <- bake(cv_rec_tmp, train_df[cv_train_idx, ] %>%
                        mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    select(-iwise12_imp) %>% as.matrix()
  y_tmp_train <- as.numeric(as.character(train_df$iwise12_imp[cv_train_idx]))
  w_tmp_train <- train_df[[weight_var]][cv_train_idx]
  X_tmp_test  <- bake(cv_rec_tmp, train_df[cv_test_idx, ] %>%
                        mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    select(-iwise12_imp) %>% as.matrix()
  
  # Elastic Net
  cv_enet_tmp <- glmnet(X_tmp_train, y_tmp_train,
                        family = "binomial", weights = w_tmp_train,
                        alpha = best_alpha, lambda = best_lambda)
  train_df$enet_oof_prob[cv_test_idx] <- predict(
    cv_enet_tmp, X_tmp_test, type = "response")[,1]
  
  # XGBoost
  w_tmp_scaled <- w_tmp_train / mean(w_tmp_train)
  d_tmp <- xgb.DMatrix(data = X_tmp_train, label = y_tmp_train,
                       weight = w_tmp_scaled)
  m_tmp <- xgb.train(params = xgb_params, data = d_tmp,
                     nrounds = best_nrounds, verbose = 0)
  train_df$xgb_oof_prob[cv_test_idx] <- predict(
    m_tmp, xgb.DMatrix(data = X_tmp_test))
  
  # Random Forest
  rf_train_tmp <- train_df[cv_train_idx, ] %>%
    select(all_of(c("iwise12_imp", all_predictors)))
  rf_test_tmp  <- train_df[cv_test_idx, ] %>%
    select(all_of(c("iwise12_imp", all_predictors)))
  
  rf_tmp <- ranger(iwise12_imp ~ ., data = rf_train_tmp,
                   case.weights = w_tmp_train,
                   num.trees = 1000,
                   mtry = floor(sqrt(length(all_predictors))),
                   min.node.size = 10,
                   probability = TRUE, seed = 241189)
  train_df$rf_oof_prob[cv_test_idx] <- predict(rf_tmp, rf_test_tmp)$predictions[, "1"]
  
  # Mixed-Effects
  tryCatch({
    me_train_tmp <- train_df[cv_train_idx, ] %>%
      mutate(outcome_num = as.numeric(as.character(iwise12_imp)),
             wt_scaled = wgtnorm_using / mean(wgtnorm_using))
    
    me_tmp <- glmer(
      outcome_num ~ waterdissatisfied_imp + inadincome_imp + age_imp +
        hhsize + female + education_imp + INCOME_5_imp + urban_final_imp +
        Drinking_basic + u5_mort + CVI_imp + log_GNI_imp + log_pop_density +
        region_final + (1 | country_code),
      data = me_train_tmp,
      family = binomial(link = "logit"),
      weights = me_train_tmp$wt_scaled,
      control = glmerControl(optimizer = "nloptwrap",
                             optCtrl = list(maxeval = 20000),
                             calc.derivs = FALSE),
      nAGQ = 0
    )
    train_df$me_oof_prob[cv_test_idx] <- predict(
      me_tmp, newdata = train_df[cv_test_idx, ],
      type = "response", re.form = NA, allow.new.levels = TRUE)
  }, error = function(e) {
    cat("    Mixed-effects failed fold", f, ":", e$message, "\n")
    train_df$me_oof_prob[cv_test_idx] <<- NA
  })
}

# Verify
cat("NAs - enet:", sum(is.na(train_df$enet_oof_prob)),
    "| xgb:", sum(is.na(train_df$xgb_oof_prob)),
    "| rf:", sum(is.na(train_df$rf_oof_prob)),
    "| me:", sum(is.na(train_df$me_oof_prob)), "\n")

##########################
# COUNTRY-LEVEL PREDICTED VS OBSERVED (WEIGHTED) 



country_pred_train <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final) %>%
  summarise(
    n          = n(),
    observed   = weighted.mean(outcome_num, wgtnorm_using),
    pred_enet  = weighted.mean(enet_oof_prob, wgtnorm_using),
    pred_xgb   = weighted.mean(xgb_oof_prob, wgtnorm_using),
    .groups    = "drop"
  )

# R2 and MAE
cat("  Country-level calibration:\n")
country_pred_train %>%
  summarise(
    enet_R2   = round(cor(observed, pred_enet)^2, 4),
    enet_MAE  = round(mean(abs(pred_enet - observed)), 4),
    xgb_R2    = round(cor(observed, pred_xgb)^2, 4),
    xgb_MAE   = round(mean(abs(pred_xgb - observed)), 4)
  ) %>% print()

#  enet_R2 enet_MAE xgb_R2 xgb_MAE

# Comparison table
country_comparison <- country_pred_train %>%
  mutate(
    diff_enet = round(pred_enet - observed, 4),
    diff_xgb  = round(pred_xgb - observed, 4),
    across(c(observed, pred_enet, pred_xgb), ~ round(., 4))
  ) %>%
  arrange(desc(observed))

write_csv(country_comparison, file.path(dir_tables, "A2_country_predicted_vs_observed_weighted.csv"))

country_pred_train <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final) %>%
  summarise(
    n         = n(),
    observed  = weighted.mean(outcome_num, wgtnorm_using),
    pred_enet = weighted.mean(enet_oof_prob, wgtnorm_using),
    pred_xgb  = weighted.mean(xgb_oof_prob, wgtnorm_using),
    pred_rf   = weighted.mean(rf_oof_prob, wgtnorm_using),
    pred_me   = weighted.mean(me_oof_prob, wgtnorm_using, na.rm = TRUE),
    .groups   = "drop"
  )

cat("Country-level calibration (all 4 models):\n")
tibble(
  model = c("Elastic Net", "XGBoost", "Random Forest", "Mixed-Effects"),
  R2  = c(round(cor(country_pred_train$observed, country_pred_train$pred_enet)^2, 4),
          round(cor(country_pred_train$observed, country_pred_train$pred_xgb)^2, 4),
          round(cor(country_pred_train$observed, country_pred_train$pred_rf)^2, 4),
          round(cor(country_pred_train$observed, country_pred_train$pred_me, 
                    use = "complete.obs")^2, 4)),
  MAE = c(round(mean(abs(country_pred_train$pred_enet - country_pred_train$observed)), 4),
          round(mean(abs(country_pred_train$pred_xgb - country_pred_train$observed)), 4),
          round(mean(abs(country_pred_train$pred_rf - country_pred_train$observed)), 4),
          round(mean(abs(country_pred_train$pred_me - country_pred_train$observed), 
                     na.rm = TRUE), 4))
) %>% print()

write_csv(country_pred_train, file.path(dir_tables, "A2_country_predicted_vs_observed_all_models.csv"))



country_pred_long <- country_pred_train %>%
  left_join(df %>% distinct(country_code, m49_region), by = "country_code") %>%
  pivot_longer(cols = c(pred_enet, pred_xgb, pred_rf, pred_me),
               names_to = "model", values_to = "predicted") %>%
  mutate(model = case_when(
    model == "pred_enet" ~ "Elastic Net",
    model == "pred_xgb"  ~ "XGBoost",
    model == "pred_rf"   ~ "Random Forest",
    model == "pred_me"   ~ "Mixed-Effects LR"
  ),
  model = factor(model, levels = c("Elastic Net", "XGBoost",
                                   "Random Forest", "Mixed-Effects LR")))


p_all_models <- ggplot(country_pred_long, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = m49_region, size = n), alpha = 0.7) +
  geom_text(aes(label = country_code), size = 2, vjust = -0.8,
            check_overlap = TRUE) +
  facet_wrap(~ model, ncol = 2) +
  scale_colour_manual(values = m49_col) +
  coord_equal(xlim = c(0, max(country_pred_long$observed, na.rm = TRUE) * 1.1),
              ylim = c(0, max(country_pred_long$observed, na.rm = TRUE) * 1.1)) +
  labs(title = "Predicted vs Observed Prevalence — All Models (OOF)",
       x = "Predicted Prevalence", y = "Observed Prevalence",
       colour = "Region", size = "Sample size") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")

print(p_all_models)
ggsave(file.path(dir_figures, "A2_predicted_vs_observed_all_models.png"), p_all_models,
       width = 12, height = 10, dpi = 300)

#############################################

# Plot
country_pred_long <- country_pred_train %>%
  pivot_longer(cols = c(pred_enet, pred_xgb),     # ENET and XGB only
               names_to = "model", values_to = "predicted") %>%
  mutate(model = ifelse(model == "pred_enet", "Elastic Net", "XGBoost"))

ggplot(country_pred_long, aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = region_final, size = n), alpha = 0.7) +
  geom_text(aes(label = country_code), size = 2.5, vjust = -0.8,
            check_overlap = TRUE) +
  facet_wrap(~ model) +
  coord_equal(xlim = c(0, max(country_pred_long$observed) * 1.1),
              ylim = c(0, max(country_pred_long$observed) * 1.1)) +
  labs(title = "Predicted vs Observed Water Insecurity Prevalence",
       subtitle = "Out-of-fold predictions, all training countries (weighted)",
       x = "Observed Prevalence", y = "Predicted Prevalence",
       colour = "Region", size = "Sample size") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(file.path(dir_figures, "A2_predicted_vs_observed_weighted.png"), width = 12, height = 7, dpi = 300)


# COUNTRY-LEVEL AUC DIAGNOSTICS 



country_auc <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final) %>%
  filter(n() >= 30, n_distinct(iwise12_imp) == 2) %>%
  summarise(
    n    = n(),
    prev = mean(outcome_num),
    auc  = as.numeric(auc(roc(outcome_num, xgb_oof_prob, quiet = TRUE))),
    .groups = "drop"
  ) %>%
  arrange(auc)

cat("  Hardest to predict:\n")
print(head(country_auc, 10))
cat("\n  Easiest to predict:\n")
print(tail(country_auc, 10))

write_csv(country_auc, file.path(dir_tables, "A2_country_level_auc.csv"))


# PREDICTION PERFORMANCE PLOTS 

# India holdout
p_dist_focal <- tibble(
  actual = factor(y_test_focal,
                  labels = c("Not water insecure", "Water insecure")),
  predicted_prob = enet_prob_focal
) %>%
  ggplot(aes(x = predicted_prob, fill = actual)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  labs(title = "Predicted Probability Distribution - India Holdout (Elastic Net)",
       x = "Predicted Probability", y = "Density", fill = "Actual Outcome") +
  theme_minimal()

# Other countries holdout
p_dist_other <- tibble(
  actual = factor(y_test_other,
                  labels = c("Not water insecure", "Water insecure")),
  predicted_prob = enet_prob_other
) %>%
  ggplot(aes(x = predicted_prob, fill = actual)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  labs(title = "Predicted Probability Distribution - Other Countries Holdout (Elastic Net)",
       x = "Predicted Probability", y = "Density", fill = "Actual Outcome") +
  theme_minimal()

print(p_dist_focal)
print(p_dist_other)

ggsave(file.path(dir_figures, "A2_prob_distribution_india.png"), p_dist_focal, width = 8, height = 5, dpi = 300)
ggsave(file.path(dir_figures, "A2_prob_distribution_other.png"), p_dist_other, width = 8, height = 5, dpi = 300)


# Calibration plot
calibration_plot <- function(y_true, y_prob, title) {
  tibble(actual = y_true, predicted = y_prob) %>%
    mutate(bin = cut(predicted, breaks = seq(0, 1, by = 0.1),
                     include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarise(
      mean_predicted = mean(predicted),
      mean_actual    = mean(actual),
      n = n(), .groups = "drop"
    ) %>%
    ggplot(aes(x = mean_predicted, y = mean_actual)) +
    geom_point(aes(size = n), colour = "#E41A1C") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_size_continuous(range = c(2, 8)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "Mean Predicted Probability",
         y = "Observed Proportion", size = "n obs") +
    theme_minimal()
}
c1 <- calibration_plot(y_test_focal, enet_prob_focal, 
                       "Calibration - India (Elastic Net)")
c2 <- calibration_plot(y_test_other, enet_prob_other, 
                       "Calibration - Other Countries (Elastic Net)")

print(c1)
print(c2)

ggsave(file.path(dir_figures, "A2_calibration_india.png"), c1, width = 6, height = 6, dpi = 300)
ggsave(file.path(dir_figures, "A2_calibration_other.png"), c2, width = 6, height = 6, dpi = 300)

library(patchwork)

# Probability distributions side by side
p_dist_combined <- p_dist_focal + p_dist_other +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_dist_combined)
ggsave(file.path(dir_figures, "A2_prob_distribution_both.png"), p_dist_combined, 
       width = 14, height = 5, dpi = 300)

# Calibration plots side by side
c_combined <- c1 + c2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(c_combined)
ggsave(file.path(dir_figures, "A2_calibration_both.png"), c_combined, 
       width = 12, height = 6, dpi = 300)
#
#####################################
# SENSITIVITY ANALYSES 
####################################

year_performance <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(YEAR_WAVE) %>%
  filter(n() >= 50, n_distinct(iwise12_imp) == 2) %>%
  summarise(
    n_obs       = n(),
    n_countries = n_distinct(country_code),
    prevalence  = mean(outcome_num),
    auc = round(as.numeric(auc(roc(outcome_num, enet_oof_prob, quiet = TRUE))), 4),
    .groups = "drop"
  ) %>%
  arrange(YEAR_WAVE)

print(year_performance, n = Inf)
write_csv(year_performance, file.path(dir_tables, "A2_sensitivity_year_performance.csv"))

ggplot(year_performance, aes(x = YEAR_WAVE, y = auc)) +
  geom_point(aes(size = n_obs), colour = "#377EB8") +
  geom_line(colour = "#377EB8", alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(title = "Elastic Net - AUC by Survey Year",
       x = "Survey Year", y = "AUC-ROC", size = "Observations") +
  theme_minimal()

ggsave(file.path(dir_figures, "A2_sensitivity_auc_by_year.png"), width = 8, height = 5, dpi = 300)


region_performance <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(region_final) %>%
  filter(n() >= 50, n_distinct(iwise12_imp) == 2) %>%
  summarise(
    n_obs       = n(),
    n_countries = n_distinct(country_code),
    prevalence  = round(mean(outcome_num), 4),
    auc = round(as.numeric(auc(roc(outcome_num, enet_oof_prob, quiet = TRUE))), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(auc))

print(region_performance, n = Inf)
write_csv(region_performance, file.path(dir_tables, "A2_sensitivity_region_performance.csv"))

region_performance %>%
  mutate(region_final = fct_reorder(region_final, auc)) %>%
  ggplot(aes(x = auc, y = region_final)) +
  geom_point(aes(size = n_obs, colour = prevalence), alpha = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_gradient(low = "#377EB8", high = "#E41A1C") +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title = "Elastic Net - AUC by Region",
       x = "AUC-ROC", y = NULL, size = "Observations", colour = "Prevalence") +
  theme_minimal()

ggsave(file.path(dir_figures, "A2_sensitivity_auc_by_region.png"), width = 9, height = 6, dpi = 300)

# Country-level calibration by region
region_calibration <- country_pred_train %>%
  group_by(region_final) %>%
  filter(n() >= 3) %>%
  summarise(
    n_countries = n(),
    R2          = round(cor(observed, pred_enet)^2, 4),
    MAE         = round(mean(abs(pred_enet - observed)), 4),
    mean_bias   = round(mean(pred_enet - observed), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(R2))

cat("\n  Country-level calibration by region:\n")
print(region_calibration, n = Inf)
write_csv(region_calibration, file.path(dir_tables, "A2_sensitivity_region_calibration.csv"))


# SHAP VALUES (XGBOOST) 


library(shapviz)

set.seed(241189)
idx <- sample(nrow(X_train), min(5000, nrow(X_train)))
X_sample <- X_train[idx, ]
colnames(X_sample) <- colnames(X_train)

shap <- shapviz(xgb_model, X_pred = X_sample)

# Relabel SHAP columns
old_names<- colnames(shap$S)
new_names <- sapply(old_names, relabel)
colnames(shap$S) <- new_names
colnames(shap$X) <- new_names


p_shap <- sv_importance(shap, kind = "beeswarm", max_display = 15) +
  labs(title = "SHAP Feature Importance - XGBoost") +
  theme_minimal()

print(p_shap)
ggsave(file.path(dir_figures, "A2_shap_beeswarm_xgboost.png"), p_shap, width = 10, height = 8, dpi = 300)

p_shap_bar <- sv_importance(shap, kind = "bar", max_display = 15) +
  labs(title = "Mean |SHAP| - XGBoost") +
  theme_minimal()

print(p_shap_bar)
ggsave(file.path(dir_figures, "A2_shap_bar_xgboost.png"), p_shap_bar, width = 8, height = 6, dpi = 300)



# Manual SHAP for linear model: coefficient * feature value
# (subtract intercept contribution)
enet_coefs_all <- coef(enet_final)[-1, 1]  # drop intercept
shap_matrix <- sweep(X_sample, 2, enet_coefs_all, "*")

# Create shapviz object manually
sv_enet <- shapviz(
  object = shap_matrix,
  X = as.data.frame(X_sample),
  baseline = coef(enet_final)[1, 1]
)
# Relabel
old_names_enet <- colnames(sv_enet$S)
new_names_enet <- sapply(old_names_enet, relabel)
colnames(sv_enet$S) <- new_names_enet
colnames(sv_enet$X) <- new_names_enet

p_shap_enet <- sv_importance(sv_enet, kind = "beeswarm", max_display = 15) +
  labs(title = "SHAP Feature Importance - Elastic Net") +
  theme_minimal()

print(p_shap_enet)

ggsave(file.path(dir_figures, "A2_shap_beeswarm_enet.png"), p_shap_enet, width = 10, height = 8, dpi = 300)

# OUTPUTS

#write_csv(results, file.path(dir_tables, "A2_model_comparison_results.csv"))
write_csv(enet_coefs, file.path(dir_tables, "A2_enet_coefficients.csv"))
write_csv(country_auc, file.path(dir_tables, "A2_country_level_auc.csv"))
xgb.save(xgb_model, file.path(dir_models, "A2_xgb_best_model.xgb"))

####################################
# ══════════════════════════════════════════════════════════════
# SENSITIVITY: Excluding waterdissatisfied
# ══════════════════════════════════════════════════════════════


# Remove waterdissatisfied from predictors
all_predictors_nowd <- all_predictors[!all_predictors %in% "waterdissatisfied_imp"]

cat("  Predictors:", length(all_predictors_nowd), "(removed waterdissatisfied_imp)\n")

# New recipe without waterdissatisfied
rec_nowd <- recipe(iwise12_imp ~ .,
                   data = train_df %>%
                     select(all_of(c("iwise12_imp", all_predictors_nowd))) %>%
                     mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  prep(training = train_df %>%
         mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))

train_baked_nowd <- bake(rec_nowd, train_df %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
focal_baked_nowd <- bake(rec_nowd, test_focal_df %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
other_baked_nowd <- bake(rec_nowd, test_other_df %>%
                           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))

y_train_nowd      <- train_baked_nowd$iwise12_imp
X_train_nowd      <- train_baked_nowd %>% select(-iwise12_imp) %>% as.matrix()
X_focal_nowd      <- focal_baked_nowd %>% select(-iwise12_imp) %>% as.matrix()
X_other_nowd      <- other_baked_nowd %>% select(-iwise12_imp) %>% as.matrix()

# ── Elastic Net (no waterdissatisfied) ──
set.seed(241189)
cv_nowd <- cv.glmnet(X_train_nowd, y_train_nowd,
                     family = "binomial", weights = w_train,
                     alpha = 0.5, nfolds = 5, type.measure = "auc")

enet_nowd <- glmnet(X_train_nowd, y_train_nowd,
                    family = "binomial", weights = w_train,
                    alpha = 0.5, lambda = cv_nowd$lambda.min)

enet_nowd_focal <- predict(enet_nowd, X_focal_nowd, type = "response")[,1]
enet_nowd_other <- predict(enet_nowd, X_other_nowd, type = "response")[,1]

# ── XGBoost (no waterdissatisfied) ──
w_train_scaled_nowd <- w_train / mean(w_train)
dtrain_nowd <- xgb.DMatrix(data = X_train_nowd, label = y_train_nowd,
                           weight = w_train_scaled_nowd)

set.seed(241189)
xgb_cv_nowd <- xgb.cv(params = xgb_params, data = dtrain_nowd,
                      nrounds = 1000, nfold = 5,
                      early_stopping_rounds = 30, verbose = 0)

best_nrounds_nowd <- xgb_cv_nowd$evaluation_log$iter[
  which.max(xgb_cv_nowd$evaluation_log$test_auc_mean)]

xgb_nowd <- xgb.train(params = xgb_params, data = dtrain_nowd,
                      nrounds = best_nrounds_nowd, verbose = 0)

xgb_nowd_focal <- predict(xgb_nowd, xgb.DMatrix(data = X_focal_nowd))
xgb_nowd_other <- predict(xgb_nowd, xgb.DMatrix(data = X_other_nowd))

# ── Compare: with vs without waterdissatisfied ──
sensitivity_wd <- bind_rows(
  evaluate_model(y_test_focal, enet_prob_focal, "Elastic Net (full)", "India"),
  evaluate_model(y_test_focal, enet_nowd_focal, "Elastic Net (no WD)", "India"),
  evaluate_model(y_test_other, enet_prob_other, "Elastic Net (full)", "Other"),
  evaluate_model(y_test_other, enet_nowd_other, "Elastic Net (no WD)", "Other"),
  evaluate_model(y_test_focal, xgb_prob_focal, "XGBoost (full)", "India"),
  evaluate_model(y_test_focal, xgb_nowd_focal, "XGBoost (no WD)", "India"),
  evaluate_model(y_test_other, xgb_prob_other, "XGBoost (full)", "Other"),
  evaluate_model(y_test_other, xgb_nowd_other, "XGBoost (no WD)", "Other")
)

cat("  Holdout results:\n")
sensitivity_wd %>%
  select(model, test_set, AUC_ROC, AUC_PR, F1) %>%
  print(n = Inf)

#   model               test_set AUC_ROC AUC_PR    F1

 

# ── OOF predictions for country-level R2/MAE ──
train_df$enet_nowd_oof <- NA_real_
train_df$xgb_nowd_oof  <- NA_real_

for (f in 1:k_folds) {
  cat("  Fold", f, "...\n")
  cv_test_idx  <- which(train_df$fold == f)
  cv_train_idx <- which(train_df$fold != f)
  
  cv_rec_nowd <- recipe(iwise12_imp ~ .,
                        data = train_df[cv_train_idx, ] %>%
                          select(all_of(c("iwise12_imp", all_predictors_nowd))) %>%
                          mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = train_df[cv_train_idx, ] %>%
           mutate(iwise12_imp = as.numeric(as.character(iwise12_imp))))
  
  X_cv_train <- bake(cv_rec_nowd, train_df[cv_train_idx, ] %>%
                       mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    select(-iwise12_imp) %>% as.matrix()
  y_cv_train <- as.numeric(as.character(train_df$iwise12_imp[cv_train_idx]))
  w_cv_train <- train_df[[weight_var]][cv_train_idx]
  X_cv_test  <- bake(cv_rec_nowd, train_df[cv_test_idx, ] %>%
                       mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))) %>%
    select(-iwise12_imp) %>% as.matrix()
  
  # Elastic Net
  cv_enet_nowd <- glmnet(X_cv_train, y_cv_train,
                         family = "binomial", weights = w_cv_train,
                         alpha = 0.5, lambda = cv_nowd$lambda.min)
  train_df$enet_nowd_oof[cv_test_idx] <- predict(
    cv_enet_nowd, X_cv_test, type = "response")[,1]
  
  # XGBoost
  w_cv_scaled <- w_cv_train / mean(w_cv_train)
  d_cv <- xgb.DMatrix(data = X_cv_train, label = y_cv_train, weight = w_cv_scaled)
  cv_xgb_nowd <- xgb.train(params = xgb_params, data = d_cv,
                           nrounds = best_nrounds_nowd, verbose = 0)
  train_df$xgb_nowd_oof[cv_test_idx] <- predict(
    cv_xgb_nowd, xgb.DMatrix(data = X_cv_test))
}

# Country-level calibration comparison
country_pred_nowd <- train_df %>%
  mutate(outcome_num = as.numeric(as.character(iwise12_imp))) %>%
  group_by(country_code, region_final) %>%
  summarise(
    observed       = weighted.mean(outcome_num, wgtnorm_using),
    pred_enet_full = weighted.mean(enet_oof_prob, wgtnorm_using),
    pred_enet_nowd = weighted.mean(enet_nowd_oof, wgtnorm_using),
    pred_xgb_full  = weighted.mean(xgb_oof_prob, wgtnorm_using),
    pred_xgb_nowd  = weighted.mean(xgb_nowd_oof, wgtnorm_using),
    .groups = "drop"
  )

cat("\n  Country-level calibration:\n")
tibble(
  model = c("Elastic Net (full)", "Elastic Net (no WD)",
            "XGBoost (full)", "XGBoost (no WD)"),
  R2 = c(
    round(cor(country_pred_nowd$observed, country_pred_nowd$pred_enet_full)^2, 4),
    round(cor(country_pred_nowd$observed, country_pred_nowd$pred_enet_nowd)^2, 4),
    round(cor(country_pred_nowd$observed, country_pred_nowd$pred_xgb_full)^2, 4),
    round(cor(country_pred_nowd$observed, country_pred_nowd$pred_xgb_nowd)^2, 4)
  ),
  MAE = c(
    round(mean(abs(country_pred_nowd$pred_enet_full - country_pred_nowd$observed)), 4),
    round(mean(abs(country_pred_nowd$pred_enet_nowd - country_pred_nowd$observed)), 4),
    round(mean(abs(country_pred_nowd$pred_xgb_full - country_pred_nowd$observed)), 4),
    round(mean(abs(country_pred_nowd$pred_xgb_nowd - country_pred_nowd$observed)), 4)
  )
) %>% print()

write_csv(sensitivity_wd, file.path(dir_tables, "A2_sensitivity_no_waterdissatisfied_holdout.csv"))
write_csv(country_pred_nowd, file.path(dir_tables, "A2_sensitivity_no_waterdissatisfied_country.csv"))

#model                  R2    MAE


# ══════════════════════════════════════════════════════════════
# FEATURE IMPORTANCE — WITHOUT WATERDISSATISFIED
# ══════════════════════════════════════════════════════════════

library(patchwork)
library(shapviz)

# Elastic Net coefficients (no WD)
enet_nowd_coefs <- coef(enet_nowd) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>%
  rename(coefficient = s0) %>%
  filter(coefficient != 0, feature != "(Intercept)") %>%
  mutate(feature_label = relabel(feature)) %>%
  arrange(desc(abs(coefficient)))

p_coef_nowd <- enet_nowd_coefs %>%
  head(15) %>%
  mutate(
    feature_label = fct_reorder(feature_label, abs(coefficient)),
    direction = ifelse(coefficient > 0, "Increases risk", "Decreases risk")
  ) %>%
  ggplot(aes(x = coefficient, y = feature_label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Increases risk" = "#E41A1C",
                               "Decreases risk" = "#377EB8")) +
  labs(title = "a) Elastic Net coefficients (excl. water dissatisfaction)",
       x = "Coefficient (log-odds scale)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_coef_nowd

p_coef_nowd <- enet_nowd_coefs %>%
  head(15) %>%
  mutate(
    feature_label = fct_reorder(feature_label, abs(coefficient)),
    direction = ifelse(coefficient > 0, "Increases risk", "Decreases risk")
  ) %>%
  ggplot(aes(x = coefficient, y = feature_label, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Increases risk" = "darkred",
                               "Decreases risk" = "navy")) +
  labs(title = "Elastic Net coefficients (excl. water dissatisfaction)",
       x = "Coefficient (log-odds scale)", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_coef_nowd

# XGBoost SHAP (no WD)
set.seed(241189)
idx_nowd <- sample(nrow(X_train_nowd), min(5000, nrow(X_train_nowd)))
X_sample_nowd <- X_train_nowd[idx_nowd, ]

shap_nowd <- shapviz(xgb_nowd, X_pred = X_sample_nowd)

# Relabel SHAP columns
old_names_nowd <- colnames(shap_nowd$S)
new_names_nowd <- sapply(old_names_nowd, relabel)
colnames(shap_nowd$S) <- new_names_nowd
colnames(shap_nowd$X) <- new_names_nowd

p_shap_nowd <- sv_importance(shap_nowd, kind = "beeswarm", max_display = 15) +
  labs(title = "XGBoost SHAP values (excl. water dissatisfaction)") +
  theme_minimal()

p_shap_nowd

# Side by side
p_importance_nowd <- p_coef_nowd | p_shap_nowd

print(p_importance_nowd)
ggsave(file.path(dir_figures, "A2_feature_importance_no_waterdissatisfied.png"), p_importance_nowd,
       width = 16, height = 8, dpi = 300)

# ══════════════════════════════════════════════════════════════
# COMPARISON: Full vs No WD (coefficient shift)
# ══════════════════════════════════════════════════════════════

# How much do coefficients change?
enet_full_coefs <- coef(enet_final) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>%
  rename(coef_full = s0) %>%
  filter(feature != "(Intercept)")

enet_nowd_coefs_all <- coef(enet_nowd) %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("feature") %>%
  rename(coef_nowd = s0) %>%
  filter(feature != "(Intercept)")

coef_shift <- enet_full_coefs %>%
  left_join(enet_nowd_coefs_all, by = "feature") %>%
  mutate(
    coef_nowd = replace_na(coef_nowd, 0),
    shift     = coef_nowd - coef_full,
    feature_label = relabel(feature)
  ) %>%
  filter(coef_full != 0 | coef_nowd != 0) %>%
  arrange(desc(abs(shift)))

cat("\nBiggest coefficient shifts when waterdissatisfied removed:\n")
coef_shift %>%
  select(feature_label, coef_full, coef_nowd, shift) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  head(10) %>%
  print()

write_csv(coef_shift, file.path(dir_tables, "A2_coefficient_shift_no_waterdissatisfied.csv"))

############
country_pred_train %>%
  left_join(df %>% distinct(country_code, m49_region), by = "country_code") %>%
  group_by(m49_region) %>%
  summarise(
    n_countries = n(),
    mean_prev   = round(mean(observed), 3),
    sd_prev     = round(sd(observed), 3),
    .groups = "drop"
  ) %>%
  arrange(mean_prev)

#####################################################

######################################################
#IMPUTATION SENSITIVITY- see separate script


# ══════════════════════════════════════════════════════════════════════════════
# SAVE OBJECTS NEEDED BY THE OPTIONAL A3 SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

save(df, train_df, test_focal_df, test_other_df, all_predictors, all_predictors_nowd,
     w_train, xgb_params, best_alpha, best_lambda, best_nrounds,
     enet_prob_focal, enet_prob_other, xgb_prob_focal, xgb_prob_other,
     y_test_focal, y_test_other, evaluate_model, k_folds, weight_var,
     file = file.path(dir_data, "A2_objects.RData"))

stop_timer(optional = TRUE)

