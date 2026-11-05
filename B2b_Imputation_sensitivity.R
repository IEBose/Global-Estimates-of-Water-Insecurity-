################################################################################
#  B2b_Imputation_sensitivity.R — SENSITIVITY: HOW MUCH DO TIER 2 ESTIMATES
#                                  DEPEND ON THE SINGLE IMPUTATION?
#
#  Purpose : A4 imputes the individual-level variables 5 times (MICE), but the
#            models use imputation 1 only. This script checks how much that
#            choice matters. It does NOT change any result used in the paper.
#              Part 1  PREDICTION side: the final individual ENET (B1, trained
#                      on imputation 1) predicts the 85 Tier 2 countries from
#                      each of the 5 imputed versions of their data.
#              Part 2  TRAINING side: the final individual ENET is refitted on
#                      each of the 5 imputed versions of the 70 training
#                      countries (same alpha and lambda as B1, no
#                      cross-validation), and predicts the Tier 2 countries from
#                      imputation 1. Coefficients are compared too.
#            Countries with whole individual-level variables missing (e.g. no
#            income data at all) are marked; they are expected to vary most.
#
#  How to read
#    range_pp = largest minus smallest country estimate across the 5
#    imputations, in percentage points. Imputation 1 equals the paper estimate
#    (Part 1) / the B1 model (Part 2); the check message should say 0.
#
#  Depends : A4 (mice_completed_all, df_imputed, new_data_imputed), B1 (model,
#            recipe, settings), B2 (Tier 2 estimates, for the check)
#  Outputs : tables/B2b_tier2_by_prediction_imputation.csv (Part 1)
#            tables/B2b_tier2_by_training_imputation.csv   (Part 2)
#            tables/B2b_enet_coefficients_by_training_imputation.csv (Part 2)
#            tables/B2b_imputation_sensitivity_summary.csv
#  Run time: a few minutes (5 ENET fits with fixed lambda; no cross-validation)
################################################################################

source("00_setup.R")
start_timer("B2b")
load_pkgs(c("glmnet", "recipes"))

m_all      <- load_data("mice_completed_all")$individual   # 5 completed datasets (A4)
df_imp1    <- load_data("df_imputed")                      # imputation 1, as used by B1
nd_imp1    <- load_data("new_data_imputed")                # imputation 1, as used by B2
rec        <- load_model("ind_recipe")
enet_final <- load_model("ind_enet_final")
settings   <- load_model("ind_settings")
tier2_main <- read_csv(file.path(dir_tables, "B2_tier2_country_prevalence_enet.csv"), show_col_types = FALSE) %>%
  mutate(country_code = as.character(country_code))
fully_missing <- read_csv(file.path(dir_tables, "A4_vars_fully_missing_by_country.csv"), show_col_types = FALSE) %>%
  mutate(country_code = as.character(country_code)) %>%
  group_by(country_code) %>%
  summarise(vars_fully_imputed = paste(variable, collapse = ", "), .groups = "drop")

m_imp <- length(m_all)
ind_vars_imp <- c("age", "INCOME_5", "inadincome", "education", "waterdissatisfied",
                  "maritalstatus", "urban_final", "employment", "hhsize")
ind_vars_imp <- intersect(ind_vars_imp, names(m_all[[1]]))
num_outcome  <- function(d) d %>% mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))


# ══════════════════════════════════════════════════════════════════════════════
# 1. BUILD THE IMPUTED VERSIONS (same steps as A4 and B2)
# ══════════════════════════════════════════════════════════════════════════════

# Tier 2 data for imputation k: individual variables from completed dataset k,
# national variables as in the paper (mean of the 5 national imputations), then
# the same renaming as B2
new_data_k <- function(k) {
  comp <- m_all[[k]]
  nd <- nd_imp1
  for (v in ind_vars_imp) nd[[v]] <- comp[[v]][match(nd$row_id, comp$row_id)]
  nd %>%
    mutate(
      INCOME_5          = factor(INCOME_5, levels = 1:5),
      inadincome        = factor(inadincome),
      waterdissatisfied = factor(waterdissatisfied),
      education         = factor(education,     levels = lv_edu),
      maritalstatus     = factor(maritalstatus, levels = lv_mar),
      urban_final       = factor(urban_final,   levels = lv_urban),
      employment        = factor(employment,    levels = lv_emp),
      log_GNI_imp    = log_gni_pc_ppp, Corruption_imp = cpi_score, Mort_WASH_imp = wash_mort,
      CVI_imp        = cvi,            Drinking_basic = basic_water,
      age_imp = age, INCOME_5_imp = INCOME_5, inadincome_imp = inadincome,
      education_imp = education, waterdissatisfied_imp = waterdissatisfied,
      maritalstatus_imp = maritalstatus, urban_final_imp = urban_final,
      iwise12_imp = NA_real_)
}

# Training data for imputation k: the _imp columns (or employment / hhsize,
# which have no _imp version) replaced by completed dataset k, as in A4
train_df_k <- function(k) {
  comp <- m_all[[k]]
  d <- df_imp1
  for (v in ind_vars_imp) {
    target <- if (paste0(v, "_imp") %in% names(d)) paste0(v, "_imp") else v
    d[[target]] <- comp[[v]][match(d$row_id, comp$row_id)]
  }
  d %>%
    mutate(
      INCOME_5_imp          = factor(INCOME_5_imp, levels = 1:5),
      inadincome_imp        = factor(inadincome_imp),
      waterdissatisfied_imp = factor(waterdissatisfied_imp),
      education_imp         = factor(education_imp,     levels = lv_edu),
      maritalstatus_imp     = factor(maritalstatus_imp, levels = lv_mar),
      urban_final_imp       = factor(urban_final_imp,   levels = lv_urban),
      employment            = factor(employment,        levels = lv_emp)) %>%
    filter(split == "train")
}

# Survey-weighted country prevalence of predicted probabilities (as B2)
country_prev <- function(nd, prob) {
  nd %>% mutate(p = prob) %>%
    group_by(country_code = as.character(country_code)) %>%
    summarise(prev = weighted.mean(p, .data[[weight_var]]), .groups = "drop")
}

wide_with_range <- function(long) {
  long %>%
    pivot_wider(names_from = imputation, values_from = prev, names_prefix = "imp_") %>%
    rowwise() %>%
    mutate(mean = mean(c_across(starts_with("imp_"))),
           min  = min(c_across(starts_with("imp_"))),
           max  = max(c_across(starts_with("imp_"))),
           range_pp = 100 * (max - min)) %>%
    ungroup() %>%
    left_join(fully_missing, by = "country_code") %>%
    arrange(desc(range_pp))
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. PART 1 — PREDICTION SIDE (final B1 model, 5 imputed Tier 2 datasets)
# ══════════════════════════════════════════════════════════════════════════════

pred_long <- lapply(seq_len(m_imp), function(k) {
  nd <- new_data_k(k)
  X  <- bake(rec, nd) %>% select(-iwise12_imp) %>% as.matrix()
  stopifnot(identical(colnames(X), settings$train_colnames))
  country_prev(nd, predict(enet_final, X, type = "response")[, 1]) %>% mutate(imputation = k)
}) %>% bind_rows()
by_prediction <- wide_with_range(pred_long)
save_table(by_prediction, "B2b_tier2_by_prediction_imputation")

# CHECK: imputation 1 reproduces the B2 estimates (expect 0)
chk <- by_prediction %>% inner_join(tier2_main %>% select(country_code, prev_enet_individual), by = "country_code")
message("Part 1 check - max |imputation 1 - B2 estimate| (pp): ",
        signif(100 * max(abs(chk$imp_1 - chk$prev_enet_individual)), 3))


# ══════════════════════════════════════════════════════════════════════════════
# 3. PART 2 — TRAINING SIDE (ENET refitted on each imputed training dataset)
# ══════════════════════════════════════════════════════════════════════════════
# Same recipe steps as B1 (dummies, zero-variance removal, standardisation),
# re-estimated on each imputed training set; alpha and lambda fixed at B1's.

make_recipe <- function(train_data, predictors) {
  recipe(iwise12_imp ~ ., data = train_data %>% select(all_of(c("iwise12_imp", predictors))) %>% num_outcome()) %>%
    step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors()) %>%
    prep(training = train_data %>% num_outcome())
}
bake_X <- function(rec, d) bake(rec, d %>% num_outcome()) %>% select(-iwise12_imp) %>% as.matrix()

nd1 <- new_data_k(1)
train_fits <- lapply(seq_len(m_imp), function(k) {
  message("  Training imputation ", k, " of ", m_imp)
  tr    <- train_df_k(k)
  rec_k <- make_recipe(tr, all_predictors)
  fit   <- glmnet(bake_X(rec_k, tr), as.numeric(as.character(tr$iwise12_imp)), family = "binomial",
                  weights = tr[[weight_var]], alpha = settings$best_alpha, lambda = settings$best_lambda)
  X_new <- bake(rec_k, nd1) %>% select(-iwise12_imp) %>% as.matrix()
  list(prev  = country_prev(nd1, predict(fit, X_new, type = "response")[, 1]) %>% mutate(imputation = k),
       coefs = coef(fit) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("feature") %>%
         rename(coefficient = s0) %>% mutate(imputation = k))
})

by_training <- wide_with_range(bind_rows(lapply(train_fits, `[[`, "prev")))
save_table(by_training, "B2b_tier2_by_training_imputation")

coefs_by_training <- bind_rows(lapply(train_fits, `[[`, "coefs")) %>%
  filter(feature != "(Intercept)") %>%
  pivot_wider(names_from = imputation, values_from = coefficient, names_prefix = "imp_") %>%
  rowwise() %>%
  mutate(range = max(c_across(starts_with("imp_"))) - min(c_across(starts_with("imp_")))) %>%
  ungroup() %>%
  mutate(feature_label = relabel(feature)) %>%
  arrange(desc(abs(imp_1)))
save_table(coefs_by_training, "B2b_enet_coefficients_by_training_imputation")

# CHECK: imputation 1 reproduces the B1 model's Tier 2 estimates (expect ~0)
chk2 <- by_training %>% inner_join(tier2_main %>% select(country_code, prev_enet_individual), by = "country_code")
message("Part 2 check - max |imputation 1 refit - B2 estimate| (pp): ",
        signif(100 * max(abs(chk2$imp_1 - chk2$prev_enet_individual)), 3))


# ══════════════════════════════════════════════════════════════════════════════
# 4. SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

summ <- function(tab, part) tibble(
  part = part,
  n_countries = nrow(tab),
  median_range_pp = median(tab$range_pp),
  max_range_pp = max(tab$range_pp),
  country_max_range = tab$country_code[which.max(tab$range_pp)],
  n_range_over_0.5pp = sum(tab$range_pp > 0.5),
  n_range_over_1pp   = sum(tab$range_pp > 1),
  countries_over_1pp = paste(tab$country_code[tab$range_pp > 1], collapse = ", "))
summary_tab <- bind_rows(summ(by_prediction, "Prediction side (Tier 2 data, 5 imputations)"),
                         summ(by_training,   "Training side (model refitted on 5 imputations)"))
save_table(summary_tab, "B2b_imputation_sensitivity_summary")
print(summary_tab, width = Inf)
print(by_prediction %>% select(country_code, imp_1, mean, min, max, range_pp, vars_fully_imputed) %>% head(10))

stop_timer(optional = TRUE)
