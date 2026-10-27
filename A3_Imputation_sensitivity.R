################################################################################
#  A3_Imputation_sensitivity.R — SENSITIVITY OF RESULTS TO THE IMPUTATION METHOD
#
#  OPTIONAL — FOR INTEREST ONLY. Not needed to reproduce the paper, not part of
#  the main pipeline run time, and no main-pipeline script reads its outputs.
#
#  Purpose : Compare simple imputation (A2) with
#              (1) model-based INCOME_5 imputation (ordinal model, GNI donors)
#              (2) MICE (5 imputations) — models refitted on EACH completed
#                  dataset and on the pooled (averaged) predictions
#            Look for: hold-out AUC and country-level R²/MAE that are similar
#            across methods -> results are robust to the imputation choice.
#
#  NOTE: the MICE specification here is NOT identical to the main pipeline (A4):
#  it imputes within the 76 observed countries only, ignores hold-out rows when
#  fitting the imputation model, imputes 3 national variables (A4: 4, incl. GNI),
#  and uses country as a predictor except for variables fully missing in a
#  country. This is the sensitivity analysis reported in Table S12 (hold-out
#  AUC by imputation approach: A3_tableS12_imputation_sensitivity.csv).
#
#  Depends : A1, A2 (single-imputation data and models)
#  Inputs  : OUTPUT/data/A2_objects.RData, new_data_clean.rds, wb_clean.rds
#  Outputs : OUTPUT/tables/A3_*.csv (incl. A3_tableS12_imputation_sensitivity.csv),
#            OUTPUT/figures/A3_*.png
#  Run time: ~30-60 min (individual-level MICE ~20 min + model refits)
################################################################################

source("00_setup.R")
start_timer("A3")
load_pkgs(c("mice", "MASS"))
select <- dplyr::select
filter <- dplyr::filter

load(file.path(dir_data, "A2_objects.RData"))   # df, train/test sets, models, settings from A2
new_data <- load_data("new_data_clean")           # prediction countries, NOT imputed
wb       <- load_data("wb_clean")                 # national data, NOT imputed

# National predictors under their un-imputed individual-data names (from A1)
country_vars_original <- c("log_GNI", "log_water_stress", "log_pop_density",
                           "Drinking_basic", "life_exp", "fertility",
                           "Corruption", "u5_mort", "wash_mort", "precip_mm", "cvi")


#set up datasets
dat_raw <- df %>%
  mutate(
    row_id = row_number(),
    # Individual-level: copy originals into _imp names (with NAs)
    age_imp               = age,
    INCOME_5_imp          = INCOME_5,
    inadincome_imp        = inadincome,
    education_imp         = education,
    waterdissatisfied_imp = waterdissatisfied,
    maritalstatus_imp     = maritalstatus,
    urban_final_imp       = urban_final,
    # Country-level: copy originals into _imp names (with NAs)
    Corruption_imp        = Corruption,
    Mort_WASH_imp         = wash_mort,
    CVI_imp               = cvi
  )

# Also add row_id to df
df <- df %>% mutate(row_id = row_number())

# And dat_simple is the current df (already imputed)
dat_simple <- df

#double check that data is missing in raw but not in simple
cat("dat_raw NAs:\n")
dat_raw %>%
  summarise(across(c(age_imp, INCOME_5_imp, Corruption_imp, Mort_WASH_imp, CVI_imp),
                   ~ sum(is.na(.)))) %>%
  print()

cat("\ndat_simple NAs:\n")
dat_simple %>%
  summarise(across(c(age_imp, INCOME_5_imp, Corruption_imp, Mort_WASH_imp, CVI_imp),
                   ~ sum(is.na(.)))) %>%
  print()

########
train_df$row_id <- df$row_id[df$split == "train"]
test_focal_df$row_id <- df$row_id[df$split == "holdout_focal"]
test_other_df$row_id <- df$row_id[df$split == "holdout_other"]

# Verify
cat("train:", length(train_df$row_id), "=", sum(df$split == "train"), "\n")
cat("focal:", length(test_focal_df$row_id), "=", sum(df$split == "holdout_focal"), "\n")
cat("other:", length(test_other_df$row_id), "=", sum(df$split == "holdout_other"), "\n")

test_ids <- c(test_focal_df$row_id, test_other_df$row_id)

set.seed(241189)
# ══════════════════════════════════════════════════════════════
# VARIABLE DEFINITIONS
# ══════════════════════════════════════════════════════════════

id_var      <- "row_id"
country_var <- "country_code"
region_var  <- "region_final"
outcome_var <- "iwise12_imp"

# Individual-level vars with missingness
ind_vars <- c("age_imp", "INCOME_5_imp", "inadincome_imp", "education_imp",
              "waterdissatisfied_imp", "maritalstatus_imp", "urban_final_imp")

# Country-level vars with missingness
ctry_vars <- c("Corruption_imp", "Mort_WASH_imp", "CVI_imp")

# Complete country-level vars that help impute ctry_vars
ctry_aux <- c("log_GNI_imp", "log_water_stress", "log_pop_density",
              "Drinking_basic", "life_exp", "fertility",
              "u5_mort", "precip_mm")

# Complete individual vars used as auxiliary predictors
aux_vars <- c("female", "hhsize")

# Test IDs
test_ids <- c(test_focal_df$row_id, test_other_df$row_id)

# MICE settings
m     <- 5
maxit <- 10

# Factor levels
edu_levels   <- c("Elementary", "Secondary", "College")
mar_levels   <- c("Single", "Married/Partnered", "Separated/Widowed")
urban_levels <- c("Rural", "Towns/semi-dense", "Cities")

# Set types helper
set_types <- function(d) {
  d %>%
    mutate(
      across(all_of(c(country_var, region_var)), as.factor),
      INCOME_5_imp          = as.numeric(as.character(INCOME_5_imp)),
      inadincome_imp        = factor(inadincome_imp),
      waterdissatisfied_imp = factor(waterdissatisfied_imp),
      education_imp         = factor(education_imp, levels = edu_levels),
      maritalstatus_imp     = factor(maritalstatus_imp, levels = mar_levels),
      urban_final_imp       = factor(urban_final_imp, levels = urban_levels)
    )
}

dat_raw    <- set_types(dat_raw)
dat_simple <- set_types(dat_simple)

# Applied before model fitting
prep_model <- function(d) {
  d %>% mutate(INCOME_5_imp = factor(INCOME_5_imp, levels = 1:5))
}


# ══════════════════════════════════════════════════════════════
# 1. CHECKS
# ══════════════════════════════════════════════════════════════


# 1a. Country-level vars should have at most one distinct value per country
dupes <- dat_raw %>%
  group_by(.data[[country_var]]) %>%
  summarise(across(all_of(ctry_vars), ~ n_distinct(na.omit(.)))) %>%
  filter(if_any(all_of(ctry_vars), ~ . > 1))

if (nrow(dupes) > 0) {
  cat("  WARNING: these country vars have >1 value per country:\n")
  print(dupes)
} else {
  cat("  Country-level vars OK (one value per country).\n")
}

# 1b. Individual vars completely missing within a country
fully_missing <- dat_raw %>%
  group_by(.data[[country_var]]) %>%
  summarise(across(all_of(ind_vars), ~ all(is.na(.)))) %>%
  pivot_longer(-all_of(country_var), names_to = "variable", values_to = "all_na") %>%
  filter(all_na)

cat("  Variables fully missing in some countries:\n")
print(fully_missing)
vars_regional <- unique(fully_missing$variable)

# Missingness overview
cat("\n  Missingness rates:\n")
dat_raw %>%
  summarise(across(all_of(c(ind_vars, ctry_vars)), ~ round(mean(is.na(.)) * 100, 2))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing)) %>%
  print(n = Inf)


# ══════════════════════════════════════════════════════════════
# 2. COUNTRY-LEVEL MICE
# ══════════════════════════════════════════════════════════════


ctry_raw <- dat_raw %>%
  group_by(across(all_of(c(country_var, region_var)))) %>%
  summarise(across(all_of(c(ctry_vars, ctry_aux)), ~ first(na.omit(.))),
            .groups = "drop")

ctry_mice_dat <- dplyr::select(ctry_raw, -all_of(country_var))

meth_c <- make.method(ctry_mice_dat)
meth_c[!names(meth_c) %in% ctry_vars] <- ""

pred_c <- make.predictorMatrix(ctry_mice_dat)

set.seed(241189)
imp_ctry <- mice(ctry_mice_dat, m = m, maxit = maxit,
                 method = meth_c, predictorMatrix = pred_c,
                 printFlag = FALSE)

if (!is.null(imp_ctry$loggedEvents)) {
  cat("  MICE logged events:\n")
  print(imp_ctry$loggedEvents)
} else {
  cat("  Country-level MICE complete. No issues.\n")
}

complete_ctry <- function(k) {
  complete(imp_ctry, k) %>%
    mutate(!!country_var := ctry_raw[[country_var]]) %>%
    select(all_of(c(country_var, ctry_vars)))
}


# ══════════════════════════════════════════════════════════════
# 3. INDIVIDUAL-LEVEL MICE
# ══════════════════════════════════════════════════════════════


ind_dat <- dat_raw %>%
  select(all_of(c(id_var, country_var, region_var, ind_vars, aux_vars)))

meth_i <- make.method(ind_dat)
meth_i[!names(meth_i) %in% ind_vars] <- ""

pred_i <- make.predictorMatrix(ind_dat)
pred_i[, id_var] <- 0
pred_i[, region_var] <- 0
# Variables fully missing in some countries: use region instead of country
if (length(vars_regional) > 0) {
  pred_i[vars_regional, country_var] <- 0
  pred_i[vars_regional, region_var]  <- 1
}

ignore_rows <- ind_dat[[id_var]] %in% test_ids


cat("  Running full MICE (m =", m, ", maxit =", maxit, ")...\n")

set.seed(241189)
imp_ind <- mice(ind_dat, m = m, maxit = maxit,
                method = meth_i, predictorMatrix = pred_i,
                ignore = ignore_rows,
                printFlag = TRUE)

if (!is.null(imp_ind$loggedEvents)) {
  cat("  Logged events:\n")
  print(imp_ind$loggedEvents)
} else {
  cat("  Individual-level MICE complete. No issues.\n")
}


# ══════════════════════════════════════════════════════════════
# 4. COMBINE INTO m COMPLETED DATASETS
# ══════════════════════════════════════════════════════════════


complete_k <- function(k) {
  complete(imp_ind, k) %>%
    left_join(complete_ctry(k), by = country_var)
}

# Verify
test_k <- complete_k(1)
cat("  Completed dataset 1: ", nrow(test_k), "rows,", ncol(test_k), "cols\n")
cat("  NAs remaining:", sum(is.na(test_k %>% select(all_of(c(ind_vars, ctry_vars))))), "\n")


# ══════════════════════════════════════════════════════════════
# 5. COMPARE MICE vs SIMPLE IMPUTATION
# ══════════════════════════════════════════════════════════════


all_imp_vars <- c(ind_vars, ctry_vars)
num_vars <- dat_raw %>% select(all_of(all_imp_vars)) %>% select(where(is.numeric)) %>% names()
cat_vars <- setdiff(all_imp_vars, num_vars)

# Cells that were originally missing
miss_cells <- dat_raw %>%
  select(all_of(c(id_var, all_imp_vars))) %>%
  mutate(across(all_of(all_imp_vars), is.na)) %>%
  pivot_longer(all_of(all_imp_vars), names_to = "variable", values_to = "was_na") %>%
  filter(was_na) %>%
  select(-was_na)

# Build MICE long dataset
mice_long <- map(seq_len(m), function(k) {
  complete_k(k) %>% mutate(.imp = k)
}) %>% list_rbind()

# 5a. Numeric variables
if (length(num_vars) > 0) {
  num_mice <- mice_long %>%
    select(.imp, all_of(c(id_var, country_var, num_vars))) %>%
    pivot_longer(all_of(num_vars), names_to = "variable", values_to = "mice_value") %>%
    semi_join(miss_cells, by = c(id_var, "variable"))
  
  num_simple <- dat_simple %>%
    select(all_of(c(id_var, num_vars))) %>%
    pivot_longer(all_of(num_vars), names_to = "variable", values_to = "simple_value")
  
  num_cmp <- num_mice %>%
    group_by(across(all_of(c(id_var, country_var))), variable) %>%
    summarise(mice_mean = mean(mice_value),
              mice_sd   = sd(mice_value),
              mice_lo   = quantile(mice_value, 0.025),
              mice_hi   = quantile(mice_value, 0.975),
              .groups = "drop") %>%
    left_join(num_simple, by = c(id_var, "variable")) %>%
    mutate(diff = mice_mean - simple_value,
           simple_in_mice_range = between(simple_value, mice_lo, mice_hi))
  
  cat("\n  Numeric variable comparison:\n")
  bind_rows(
    num_cmp %>% filter(variable %in% ind_vars) %>% mutate(level = "individual"),
    num_cmp %>% filter(variable %in% ctry_vars) %>%
      distinct(across(all_of(country_var)), variable, .keep_all = TRUE) %>%
      mutate(level = "country")
  ) %>%
    group_by(level, variable) %>%
    summarise(n_imputed        = n(),
              mean_simple      = round(mean(simple_value), 3),
              mean_mice        = round(mean(mice_mean), 3),
              mean_abs_diff    = round(mean(abs(diff)), 3),
              pct_simple_in_95 = round(100 * mean(simple_in_mice_range), 1),
              .groups = "drop") %>%
    print(n = Inf, width = Inf)
}

# 5b. Categorical variables
if (length(cat_vars) > 0) {
  cat_cmp <- mice_long %>%
    select(.imp, all_of(c(id_var, country_var, cat_vars))) %>%
    mutate(across(all_of(cat_vars), as.character)) %>%
    pivot_longer(all_of(cat_vars), names_to = "variable", values_to = "mice_value") %>%
    semi_join(miss_cells, by = c(id_var, "variable")) %>%
    left_join(dat_simple %>%
                select(all_of(c(id_var, cat_vars))) %>%
                mutate(across(all_of(cat_vars), as.character)) %>%
                pivot_longer(all_of(cat_vars), names_to = "variable",
                             values_to = "simple_value"),
              by = c(id_var, "variable"))
  
  cat("\n  Categorical variable agreement (MICE vs simple):\n")
  cat_cmp %>%
    group_by(variable) %>%
    summarise(n_imputed_cells     = n_distinct(.data[[id_var]]),
              pct_agree_with_mode = round(100 * mean(mice_value == simple_value), 1)) %>%
    print(n = Inf)
}



# Distribution plot
if (length(num_vars) > 0 && any(num_vars %in% ind_vars)) {
  obs_long <- dat_raw %>%
    select(all_of(intersect(ind_vars, num_vars))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    drop_na() %>%
    mutate(source = "Observed")
  
  p_dist_compare <- bind_rows(
    obs_long,
    num_cmp %>% filter(variable %in% ind_vars) %>%
      transmute(variable, value = simple_value, source = "Simple (median)"),
    num_mice %>% filter(variable %in% ind_vars) %>%
      transmute(variable, value = mice_value, source = "MICE")
  ) %>%
    ggplot(aes(value, colour = source)) +
    geom_density() +
    facet_wrap(~ variable, scales = "free") +
    labs(x = NULL, y = "Density", colour = NULL,
         title = "Distribution of imputed values: Simple vs MICE") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p_dist_compare)
  ggsave(file.path(dir_figures, "A3_mice_vs_simple_distributions.png"), p_dist_compare,
         width = 12, height = 8, dpi = 300)
}

p_dist_compare_bar <- # Better comparison for discrete INCOME_5
  bind_rows(
    dat_raw %>%
      filter(!is.na(INCOME_5_imp)) %>%
      transmute(value = as.numeric(as.character(INCOME_5_imp)), source = "Observed"),
    dat_simple %>%
      filter(is.na(dat_raw$INCOME_5_imp)) %>%
      transmute(value = as.numeric(as.character(INCOME_5_imp)), source = "Simple (median)"),
    mice_long %>%
      semi_join(miss_cells %>% filter(variable == "INCOME_5_imp"), by = id_var) %>%
      transmute(value = as.numeric(as.character(INCOME_5_imp)), source = "MICE")
  ) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = factor(value), y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "Income quintile", y = "Proportion", fill = NULL,
       title = "Imputed INCOME_5 distribution: Simple vs MICE") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_dist_compare_bar)
ggsave(file.path(dir_figures, "A3_mice_vs_simple_distributions-income_bar.png"), p_dist_compare_bar,
       width = 12, height = 8, dpi = 300)

# ══════════════════════════════════════════════════════════════
# PART A: MODEL-BASED INCOME_5 IMPUTATION
# ══════════════════════════════════════════════════════════════


gni_obs <- dat_raw %>%
  distinct(country_code, log_GNI_obs = log_GNI_imp)

# Create INCOME_5_mod from raw (un-imputed) values
df$INCOME_5_mod <- dat_raw$INCOME_5_imp
df <- df %>%
  mutate(INCOME_5_mod = factor(as.character(INCOME_5_mod), levels = 1:5))

# Also add to train/test splits
train_df$INCOME_5_mod <- df$INCOME_5_mod[df$row_id %in% train_df$row_id]
test_focal_df$INCOME_5_mod <- df$INCOME_5_mod[df$row_id %in% test_focal_df$row_id]
test_other_df$INCOME_5_mod <- df$INCOME_5_mod[df$row_id %in% test_other_df$row_id]

# Countries with no observed income quintiles
no_inc_ctry <- df %>%
  group_by(country_code) %>%
  summarise(all_na = all(is.na(INCOME_5_mod))) %>%
  filter(all_na) %>%
  pull(country_code) %>%
  as.character()

cat("  Countries with no income data:", length(no_inc_ctry), "\n")
cat("  ", paste(no_inc_ctry, collapse = ", "), "\n")

inc_formula <- ~ inadincome_imp + education_imp + urban_final_imp
k_donors    <- 8

set.seed(241189)

# Countries with income data: ordinal model
train_inc <- df %>%
  filter(!is.na(INCOME_5_mod)) %>%
  mutate(across(where(is.factor), droplevels))

fit_inc <- polr(INCOME_5_mod ~ inadincome_imp + education_imp + urban_final_imp + country_code,
                data = train_inc)

miss <- is.na(df$INCOME_5_mod) & !df$country_code %in% no_inc_ctry

if (sum(miss) > 0) {
  df_miss <- df[miss, ] %>%
    mutate(country_code = factor(country_code, levels = levels(train_inc$country_code)),
           inadincome_imp = factor(inadincome_imp, levels = levels(train_inc$inadincome_imp)),
           education_imp = factor(education_imp, levels = levels(train_inc$education_imp)),
           urban_final_imp = factor(urban_final_imp, levels = levels(train_inc$urban_final_imp)))
  
  p <- t(sapply(1:nrow(df_miss), function(i) {
    predict(fit_inc, newdata = df_miss[i, ], type = "probs")
  }))
  
  df$INCOME_5_mod[miss] <- apply(p, 1, function(pr) sample(colnames(p), 1, prob = pr))
  cat("  Imputed", sum(miss), "rows via ordinal model.\n")
} else {
  cat("  No rows to impute (all handled by no_inc_ctry).\n")
}

# Countries with no income data: GNI-based donors
X <- model.matrix(inc_formula, df)[, -1, drop = FALSE]

for (ctry in no_inc_ctry) {
  target_gni <- gni_obs$log_GNI_obs[gni_obs$country_code == ctry]
  if (is.na(target_gni)) {
    cat("  WARNING:", ctry, "has no observed log_GNI — skipping\n")
    next
  }
  
  donors <- gni_obs %>%
    filter(!country_code %in% no_inc_ctry, !is.na(log_GNI_obs)) %>%
    mutate(gni_diff = abs(log_GNI_obs - target_gni)) %>%
    slice_min(gni_diff, n = k_donors, with_ties = FALSE)
  
  cat("  ", ctry, "donors:", paste(donors$country_code, collapse = ", "), "\n")
  
  fit_donor <- polr(update(inc_formula, INCOME_5_mod ~ . + country_code),
                    data = df %>%
                      filter(country_code %in% donors$country_code) %>%
                      droplevels())
  beta <- coef(fit_donor)[colnames(X)]
  
  rows <- df$country_code == ctry
  eta  <- drop(X[rows, , drop = FALSE] %*% beta) + rlogis(sum(rows))
  df$INCOME_5_mod[rows] <- ntile(eta, 5)
}

#nb: ZWE donors: COG, MMR, CMR, PAK, PSE, SEN, KEN, MRT 

# Propagate to train/test
train_df$INCOME_5_mod <- df$INCOME_5_mod[df$row_id %in% train_df$row_id]
test_focal_df$INCOME_5_mod <- df$INCOME_5_mod[df$row_id %in% test_focal_df$row_id]
test_other_df$INCOME_5_mod <- df$INCOME_5_mod[df$row_id %in% test_other_df$row_id]

# Cross-tab: median vs model-based
cat("\n  Median vs model-based (originally missing cells):\n")
print(table(median = df$INCOME_5_imp[is.na(dat_raw$INCOME_5_imp)],
            model  = df$INCOME_5_mod[is.na(dat_raw$INCOME_5_imp)]))


# ══════════════════════════════════════════════════════════════
# PART B: HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════

run_enet_xgb <- function(train, focal, other, predictors,
                         weights = w_train, params = xgb_params, seed = 241189) {
  
  to_num <- function(d) {
    d %>%
      select(all_of(c("iwise12_imp", predictors))) %>%
      mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))
  }
  
  rec <- recipe(iwise12_imp ~ ., data = to_num(train)) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = to_num(train))
  
  train_baked <- bake(rec, to_num(train))
  y_train_b <- train_baked$iwise12_imp
  X_train_b <- train_baked %>% select(-iwise12_imp) %>% as.matrix()
  X_focal_b <- bake(rec, to_num(focal)) %>% select(-iwise12_imp) %>% as.matrix()
  X_other_b <- bake(rec, to_num(other)) %>% select(-iwise12_imp) %>% as.matrix()
  
  # Elastic Net
  set.seed(seed)
  cv <- cv.glmnet(X_train_b, y_train_b, family = "binomial", weights = weights,
                  alpha = 0.5, nfolds = 5, type.measure = "auc")
  enet <- glmnet(X_train_b, y_train_b, family = "binomial", weights = weights,
                 alpha = 0.5, lambda = cv$lambda.min)
  
  # XGBoost
  w_scaled <- weights / mean(weights)
  dtrain <- xgb.DMatrix(data = X_train_b, label = y_train_b, weight = w_scaled)
  set.seed(seed)
  xgb_cv <- xgb.cv(params = params, data = dtrain, nrounds = 1000, nfold = 5,
                   early_stopping_rounds = 30, verbose = 0)
  best_nr <- xgb_cv$evaluation_log$iter[which.max(xgb_cv$evaluation_log$test_auc_mean)]
  xgb <- xgb.train(params = params, data = dtrain, nrounds = best_nr, verbose = 0)
  
  list(enet_focal = predict(enet, X_focal_b, type = "response")[, 1],
       enet_other = predict(enet, X_other_b, type = "response")[, 1],
       xgb_focal  = predict(xgb, xgb.DMatrix(data = X_focal_b)),
       xgb_other  = predict(xgb, xgb.DMatrix(data = X_other_b)))
}

compare_to_main <- function(res, label) {
  bind_rows(
    evaluate_model(y_test_focal, enet_prob_focal, "Elastic Net (full)", "India"),
    evaluate_model(y_test_focal, res$enet_focal,
                   paste0("Elastic Net (", label, ")"), "India"),
    evaluate_model(y_test_other, enet_prob_other, "Elastic Net (full)", "Other"),
    evaluate_model(y_test_other, res$enet_other,
                   paste0("Elastic Net (", label, ")"), "Other"),
    evaluate_model(y_test_focal, xgb_prob_focal, "XGBoost (full)", "India"),
    evaluate_model(y_test_focal, res$xgb_focal,
                   paste0("XGBoost (", label, ")"), "India"),
    evaluate_model(y_test_other, xgb_prob_other, "XGBoost (full)", "Other"),
    evaluate_model(y_test_other, res$xgb_other,
                   paste0("XGBoost (", label, ")"), "Other")
  )
}


# ══════════════════════════════════════════════════════════════
# 6. ML SENSITIVITY: MODEL-BASED INCOME
# ══════════════════════════════════════════════════════════════


all_predictors_inc <- replace(all_predictors,
                              all_predictors == "INCOME_5_imp", "INCOME_5_mod")

res_inc <- run_enet_xgb(train_df, test_focal_df, test_other_df, all_predictors_inc)

sensitivity_inc <- compare_to_main(res_inc, "model income")

cat("  Holdout results:\n")
sensitivity_inc %>%
  select(model, test_set, AUC_ROC, AUC_PR, F1) %>%
  print(n = Inf)

write_csv(sensitivity_inc, file.path(dir_tables, "A3_sensitivity_model_income_holdout.csv"))


# ══════════════════════════════════════════════════════════════
# 7. ML SENSITIVITY: MICE
# ══════════════════════════════════════════════════════════════


mice_vars <- c(ind_vars, ctry_vars)

# Swap median/mode columns for MICE-imputed columns
swap_in_mice <- function(d, completed) {
  # Match by row_id, keep row order unchanged
  d %>%
    select(-all_of(mice_vars)) %>%
    left_join(completed %>% select(row_id, all_of(mice_vars)), by = "row_id")
}

res_mice <- vector("list", m)

for (k in seq_len(m)) {
  cat("  Imputation", k, "of", m, "\n")
  completed_k <- complete_k(k) %>%
    mutate(INCOME_5_imp = factor(INCOME_5_imp, levels = 1:5))
  
  res_mice[[k]] <- run_enet_xgb(
    swap_in_mice(train_df, completed_k),
    swap_in_mice(test_focal_df, completed_k),
    swap_in_mice(test_other_df, completed_k),
    all_predictors
  )
}

# Each imputation separately
sensitivity_mice_each <- map(seq_len(m), function(k) {
  compare_to_main(res_mice[[k]], paste("MICE", k)) %>%
    filter(str_detect(model, "MICE"))
}) %>% list_rbind()

# Pooled: average predicted probabilities across imputations
res_mice_pooled <- map(set_names(names(res_mice[[1]])), function(nm) {
  rowMeans(sapply(res_mice, `[[`, nm))
})

sensitivity_mice <- bind_rows(
  compare_to_main(res_mice_pooled, "MICE pooled"),
  sensitivity_mice_each
)

cat("  Holdout results:\n")
sensitivity_mice %>%
  select(model, test_set, AUC_ROC, AUC_PR, F1) %>%
  print(n = Inf)

write_csv(sensitivity_mice, file.path(dir_tables, "A3_sensitivity_mice_holdout.csv"))


# ══════════════════════════════════════════════════════════════
# 8. COMBINED SUMMARY TABLE
# ══════════════════════════════════════════════════════════════


all_sensitivity <- bind_rows(
  sensitivity_inc %>% mutate(analysis = "Model-based income"),
  sensitivity_mice %>% mutate(analysis = "MICE")
) %>%
  filter(!str_detect(model, "full")) %>%
  select(analysis, model, test_set, AUC_ROC, AUC_PR, F1)

# Add baseline for comparison
baseline <- bind_rows(
  evaluate_model(y_test_focal, enet_prob_focal, "Elastic Net (baseline)", "India"),
  evaluate_model(y_test_other, enet_prob_other, "Elastic Net (baseline)", "Other"),
  evaluate_model(y_test_focal, xgb_prob_focal, "XGBoost (baseline)", "India"),
  evaluate_model(y_test_other, xgb_prob_other, "XGBoost (baseline)", "Other")
) %>%
  mutate(analysis = "Baseline (median/mode)") %>%
  select(analysis, model, test_set, AUC_ROC, AUC_PR, F1)

combined_summary <- bind_rows(baseline, all_sensitivity)

cat("  Full comparison:\n")
combined_summary %>%
  arrange(test_set, model) %>%
  print(n = Inf, width = Inf)

write_csv(combined_summary, file.path(dir_tables, "A3_sensitivity_imputation_all.csv"))

# ── Table S12: hold-out AUC by imputation approach (formatted) ──
# Rows: median/mode (baseline), model-based income, MICE pooled, MICE range
# across the 5 imputations, and the largest difference from the baseline.
auc <- combined_summary %>%
  mutate(algorithm = str_extract(model, "^(Elastic Net|XGBoost)"),
         approach  = case_when(
           analysis == "Baseline (median/mode)" ~ "Median/mode (baseline)",
           analysis == "Model-based income"     ~ "Model-based income",
           str_detect(model, "MICE pooled")     ~ "MICE, pooled",
           TRUE                                 ~ "MICE, range across 5 imputations"),
         column = paste0(algorithm, ": ", ifelse(test_set == "India", "India", "Other hold-outs")))
max_diff <- auc %>%
  group_by(column) %>%
  summarise(AUC = sprintf("%.3f", max(abs(AUC_ROC - AUC_ROC[approach == "Median/mode (baseline)"]))),
            .groups = "drop") %>%
  mutate(approach = "Largest difference from baseline")
tableS12 <- auc %>%
  group_by(approach, column) %>%
  summarise(AUC = if (n() == 1) sprintf("%.3f", AUC_ROC)
                  else paste0(sprintf("%.3f", min(AUC_ROC)), "–", sprintf("%.3f", max(AUC_ROC))),
            .groups = "drop") %>%
  bind_rows(max_diff) %>%
  pivot_wider(names_from = column, values_from = AUC) %>%
  mutate(approach = factor(approach, levels = c("Median/mode (baseline)", "Model-based income",
                                                "MICE, pooled", "MICE, range across 5 imputations",
                                                "Largest difference from baseline"))) %>%
  arrange(approach) %>%
  rename(`Imputation approach` = approach)
print(tableS12, width = Inf)
readr::write_excel_csv(tableS12, file.path(dir_tables, "A3_tableS12_imputation_sensitivity.csv"))
#####################

# ══════════════════════════════════════════════════════════════
# IMPUTATION COMPARISON PLOTS — ALL VARIABLES
# ══════════════════════════════════════════════════════════════

# Age (continuous — density plot)
p_age <- bind_rows(
  dat_raw %>% filter(!is.na(age_imp)) %>%
    transmute(value = age_imp, source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$age_imp)) %>%
    transmute(value = age_imp, source = "Simple (median)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "age_imp"), by = id_var) %>%
    transmute(value = age_imp, source = "MICE")
) %>%
  ggplot(aes(value, colour = source)) +
  geom_density() +
  labs(x = NULL, y = "Density", colour = NULL, title = "Age") +
  theme_minimal() +
  theme(legend.position = "none")

# INCOME_5 (discrete — bar chart)
p_income <- bind_rows(
  dat_raw %>% filter(!is.na(INCOME_5_imp)) %>%
    transmute(value = as.character(INCOME_5_imp), source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$INCOME_5_imp)) %>%
    transmute(value = as.character(INCOME_5_imp), source = "Simple (median)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "INCOME_5_imp"), by = id_var) %>%
    transmute(value = as.character(INCOME_5_imp), source = "MICE")
) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Income quintile") +
  theme_minimal() +
  theme(legend.position = "none")

# Inadincome (categorical — bar chart)
p_inadincome <- bind_rows(
  dat_raw %>% filter(!is.na(inadincome_imp)) %>%
    transmute(value = as.character(inadincome_imp), source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$inadincome_imp)) %>%
    transmute(value = as.character(inadincome_imp), source = "Simple (mode)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "inadincome_imp"), by = id_var) %>%
    transmute(value = as.character(inadincome_imp), source = "MICE")
) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Income adequacy") +
  theme_minimal() +
  theme(legend.position = "none")

# Education (categorical — bar chart)
p_education <- bind_rows(
  dat_raw %>% filter(!is.na(education_imp)) %>%
    transmute(value = as.character(education_imp), source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$education_imp)) %>%
    transmute(value = as.character(education_imp), source = "Simple (mode)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "education_imp"), by = id_var) %>%
    transmute(value = as.character(education_imp), source = "MICE")
) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Education") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# Waterdissatisfied (categorical — bar chart)
p_waterdiss <- bind_rows(
  dat_raw %>% filter(!is.na(waterdissatisfied_imp)) %>%
    transmute(value = as.character(waterdissatisfied_imp), source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$waterdissatisfied_imp)) %>%
    transmute(value = as.character(waterdissatisfied_imp), source = "Simple (mode)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "waterdissatisfied_imp"), by = id_var) %>%
    transmute(value = as.character(waterdissatisfied_imp), source = "MICE")
) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Water dissatisfied") +
  theme_minimal() +
  theme(legend.position = "none")

# Marital status (categorical — bar chart)
p_marital <- bind_rows(
  dat_raw %>% filter(!is.na(maritalstatus_imp)) %>%
    transmute(value = as.character(maritalstatus_imp), source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$maritalstatus_imp)) %>%
    transmute(value = as.character(maritalstatus_imp), source = "Simple (mode)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "maritalstatus_imp"), by = id_var) %>%
    transmute(value = as.character(maritalstatus_imp), source = "MICE")
) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Marital status") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# Urban (categorical — bar chart)
p_urban <- bind_rows(
  dat_raw %>% filter(!is.na(urban_final_imp)) %>%
    transmute(value = as.character(urban_final_imp), source = "Observed"),
  dat_simple %>% filter(is.na(dat_raw$urban_final_imp)) %>%
    transmute(value = as.character(urban_final_imp), source = "Simple (mode)"),
  mice_long %>%
    semi_join(miss_cells %>% filter(variable == "urban_final_imp"), by = id_var) %>%
    transmute(value = as.character(urban_final_imp), source = "MICE")
) %>%
  count(source, value) %>%
  group_by(source) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = NULL, y = "Proportion", fill = NULL, title = "Urbanicity") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# ══════════════════════════════════════════════════════════════
# COMBINE ALL — shared legend at bottom
# ══════════════════════════════════════════════════════════════

library(patchwork)

p_all_imputation <- (p_age | p_income | p_inadincome) /
  (p_education | p_waterdiss | p_marital) /
  (p_urban | plot_spacer() | plot_spacer()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Need to add legend back — recreate one plot with legend
p_legend_source <- bind_rows(
  tibble(value = "a", pct = 1, source = "Observed"),
  tibble(value = "a", pct = 1, source = "Simple (median/mode)"),
  tibble(value = "a", pct = 1, source = "MICE")
) %>%
  ggplot(aes(x = value, y = pct, fill = source)) +
  geom_col() +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9))

legend_only <- cowplot::get_legend(p_legend_source)

p_final <- (p_age | p_income | p_inadincome) /
  (p_education | p_waterdiss | p_marital) /
  (p_urban | plot_spacer() | plot_spacer()) /
  wrap_elements(legend_only) +
  plot_layout(heights = c(1, 1, 1, 0.1))

print(p_final)
ggsave(file.path(dir_figures, "A3_imputation_comparison_all_vars.png"), p_final,
       width = 16, height = 14, dpi = 300)
# ══════════════════════════════════════════════════════════════

ctry_comparison <- lapply(ctry_vars, function(v) {
  missing_country <- dat_raw %>%
    distinct(country_code, .keep_all = TRUE) %>%
    filter(is.na(.data[[v]])) %>%
    pull(country_code)
  
  if (length(missing_country) == 0) return(NULL)
  
  mice_val <- mice_long %>%
    filter(.data[[country_var]] %in% missing_country) %>%
    distinct(.imp, .data[[country_var]], .keep_all = TRUE) %>%
    summarise(
      mice_mean = round(mean(.data[[v]], na.rm = TRUE), 3),
      mice_lo   = round(quantile(.data[[v]], 0.025, na.rm = TRUE), 3),
      mice_hi   = round(quantile(.data[[v]], 0.975, na.rm = TRUE), 3)
    )
  
  simple_val <- dat_simple %>%
    distinct(country_code, .keep_all = TRUE) %>%
    filter(country_code %in% missing_country) %>%
    pull(v) %>%
    round(3)
  
  tibble(
    variable       = v,
    country        = as.character(missing_country),
    simple         = simple_val,
    mice_mean      = mice_val$mice_mean,
    mice_range     = paste0("(", mice_val$mice_lo, "–", mice_val$mice_hi, ")"),
    simple_in_range = simple_val >= mice_val$mice_lo & simple_val <= mice_val$mice_hi
  )
}) %>% bind_rows()

print(ctry_comparison, width = Inf)


write_csv(ctry_comparison, file.path(dir_tables, "A3_ctry_comparison_imputation.csv"))

# ══════════════════════════════════════════════════════════════
# MISSINGNESS CHECK ACROSS ALL THREE DATASETS
# ══════════════════════════════════════════════════════════════


# Individual-level
cat("\n  Individual-level:\n")
df %>%
  summarise(
    n_obs = n(),
    n_countries = n_distinct(country_code),
    across(c(age, INCOME_5, inadincome, education,
             waterdissatisfied, maritalstatus, urban_final),
           ~ round(mean(is.na(.)) * 100, 2))
  ) %>%
  pivot_longer(-c(n_obs, n_countries), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing)) %>%
  print(n = Inf)

# Country-level
cat("\n  Country-level:\n")
df %>%
  distinct(country_code, .keep_all = TRUE) %>%
  summarise(
    n_countries = n(),
    across(all_of(country_vars_original),
           ~ sum(is.na(.)))
  ) %>%
  pivot_longer(-n_countries, names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / first(n_countries) * 100, 1)) %>%
  arrange(desc(n_missing)) %>%
  print(n = Inf)



# Individual-level
cat("\n  Individual-level:\n")
new_data %>%
  summarise(
    n_obs = n(),
    n_countries = n_distinct(country_code),
    across(c(age, INCOME_5, inadincome, education,
             waterdissatisfied, maritalstatus, urban_final),
           ~ round(mean(is.na(.)) * 100, 2))
  ) %>%
  pivot_longer(-c(n_obs, n_countries), names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing)) %>%
  print(n = Inf)

# Country-level
cat("\n  Country-level:\n")

country_vars_new_data <- c("log_gni_pc_ppp", "log_water_stress", "log_pop_density",
                           "basic_water", "life_exp", "fertility",
                           "cpi_score", "u5_mort", "wash_mort", 
                           "precip_mm", "cvi")

new_data %>%
  distinct(country_code, .keep_all = TRUE) %>%
  summarise(
    n_countries = n(),
    across(all_of(country_vars_new_data),
           ~ sum(is.na(.)))
  ) %>%
  pivot_longer(-n_countries, names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / first(n_countries) * 100, 1)) %>%
  arrange(desc(n_missing)) %>%
  print(n = Inf)



# Individual-level by country (worst offenders)
cat("\n  Individual-level by country (>5% missing):\n")
new_data %>%
  group_by(country_code) %>%
  summarise(
    n = n(),
    across(c(age, INCOME_5, inadincome, education,
             waterdissatisfied, maritalstatus, urban_final),
           ~ round(mean(is.na(.)) * 100, 1)),
    .groups = "drop"
  ) %>%
  mutate(max_missing = pmax(age, INCOME_5, inadincome, education,
                            waterdissatisfied, maritalstatus, urban_final)) %>%
  filter(max_missing > 5) %>%
  arrange(desc(max_missing)) %>%
  print(n = Inf, width = Inf)




# Country-level missingness
cat("\n  Country-level missingness:\n")
wb %>%
  summarise(
    n_countries = n(),
    across(all_of(country_vars_new_data),
           ~ sum(is.na(.)))
  ) %>%
  pivot_longer(-n_countries, names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / first(n_countries) * 100, 1)) %>%
  arrange(desc(n_missing)) %>%
  print(n = Inf)

# Which countries are missing what
cat("\n  Countries with missing data:\n")
wb %>%
  select(iso3, all_of(country_vars_new_data)) %>%
  filter(if_any(all_of(country_vars_new_data), is.na)) %>%
  pivot_longer(-iso3, names_to = "variable", values_to = "value") %>%
  filter(is.na(value)) %>%
  group_by(iso3) %>%
  summarise(
    n_vars_missing = n(),
    variables = paste(variable, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_vars_missing)) %>%
  print(n = Inf, width = Inf)


# ══════════════════════════════════════════════════════════════
# COMBINED SUMMARY TABLE
# ══════════════════════════════════════════════════════════════


summary_table <- bind_rows(
  # df individual
  df %>%
    summarise(across(c(age, INCOME_5, inadincome, education,
                       waterdissatisfied, maritalstatus, urban_final),
                     ~ round(mean(is.na(.)) * 100, 2))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
    mutate(dataset = "Observed (76)", level = "Individual"),
  
  # df country
  df %>%
    distinct(country_code, .keep_all = TRUE) %>%
    summarise(across(all_of(country_vars_original),
                     ~ round(sum(is.na(.)) / n() * 100, 1))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
    mutate(dataset = "Observed (76)", level = "Country"),
  
  # new_data individual
  new_data %>%
    summarise(across(c(age, INCOME_5, inadincome, education,
                       waterdissatisfied, maritalstatus, urban_final),
                     ~ round(mean(is.na(.)) * 100, 2))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
    mutate(dataset = "Predicted (85)", level = "Individual"),
  
  # new_data country
  new_data %>%
    distinct(country_code, .keep_all = TRUE) %>%
    summarise(across(all_of(country_vars_new_data),
                     ~ round(sum(is.na(.)) / n() * 100, 1))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
    mutate(dataset = "Predicted (85)", level = "Country"),
  
  # wb_latest country (all 194)
  wb %>%
    summarise(across(all_of(country_vars_new_data),
                     ~ round(sum(is.na(.)) / n() * 100, 1))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
    mutate(dataset = "All 194", level = "Country")
)

summary_wide <- summary_table %>%
  pivot_wider(names_from = dataset, values_from = pct_missing)

cat("\n  Missingness (%) by dataset:\n")
print(summary_wide, n = Inf, width = Inf)

write_csv(summary_wide, file.path(dir_tables, "A3_missingness_summary_all_datasets.csv"))

wb %>%
  filter(is.na(log_water_stress) | is.na(cpi_score) | 
           is.na(precip_mm) | is.na(wash_mort)) %>%
  select(iso3, log_water_stress, cpi_score, precip_mm, wash_mort) %>%
  print(n = Inf, width = Inf)

#mostly pacific island countries

missing_ctry_codes <- c("AND", "ATG", "BHS", "DMA", "FSM", "KIR", "KNA", 
                        "LIE", "MCO", "MHL", "MNE", "NRU", "PLW", "PSE", 
                        "SLB", "SMR", "SRB", "SYC", "TON", "TUV", "VUT", "WSM")

wb %>%
  mutate(is_missing = iso3 %in% missing_ctry_codes) %>%
  summarise(
    total_pop   = sum(pop_15plus, na.rm = TRUE),
    missing_pop = sum(pop_15plus[is_missing], na.rm = TRUE),
    pct         = round(missing_pop / total_pop * 100, 2)
  )

stop_timer(optional = TRUE)
