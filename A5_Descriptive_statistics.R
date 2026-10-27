################################################################################
#  A5_Descriptive_statistics.R — DESCRIPTIVE STATISTICS FOR THE OBSERVED COUNTRIES
#
#  Purpose : Describe the 76 observed (IWISE) countries:
#              1. Individual-level predictors (and the outcome): survey-weighted
#                 means / percentages with 95% CIs (as Stata `svy: mean x`);
#                 for age and household size also SD, median and IQR (as Stata
#                 `summarize x [aw = projectionwgt_svy_yr], detail`),
#                 overall (and by M49 region if by_region = TRUE),
#                 BEFORE and AFTER MICE imputation.
#              2. National-level predictors: distribution across the 76
#                 countries (one value per country), before and after MICE.
#              3. Number of countries (and respondents) per M49 region and
#                 sub-region and per income group, plus all 194 countries by
#                 data source (tier).
#              4. Population 15+ (2025) covered by each tier.
#
#  Survey design (equivalent to Stata):
#    svyset psu [pweight = projectionwgt_svy_yr], strata(strata) singleunit(certainty)
#    - Strata and PSUs are made unique within country and survey year before
#      pooling, so that strata coded 1, 2, ... in different countries are not
#      merged. If your strata are already unique this changes nothing.
#    - singleunit(certainty) (Gallup guidance): strata with one PSU contribute
#      no variance. (Stata singleunit(centered) = R "adjust" gives slightly
#      larger SEs for the ~14 countries with single-PSU strata.)
#    - Means and Taylor-linearised SEs are calculated directly (svy_mean_fast),
#      which takes seconds; the survey package took ~3 min per variable on the
#      full data. Results are identical: the script checks a few estimates
#      against survey::svymean / svyby (check_against_survey; set to FALSE to
#      skip that check, which takes a few minutes).
#    - If projectionwgt_svy_yr is not in the IWISE .dta, the normalised weight
#      (wgtnorm_using) is used instead and a warning is printed. Within a
#      country both give the same result; pooled (overall / regional) results
#      differ, because normalised weights do not reflect population size.
#
#  Non-imputed vs MICE:
#    Non-imputed = original values from A1 (df_clean); missing values are
#                  excluded variable by variable (available-case).
#    MICE        = values used in the models, from A4 (df_imputed / wb_imputed).
#
#  Depends : A1 (df_clean, new_data_clean, wb_clean), A4 (df_imputed, wb_imputed)
#  Outputs : tables/A5_individual_descriptives_long.csv   (all estimates, overall + by region)
#            tables/A5_tableS01A_individual_descriptives.csv (Table S1A: overall, non-imputed vs MICE)
#            tables/A5_tableS01B_national_descriptives.csv   (Table S1B)
#            tables/A5_tableS03_individual_missingness.csv   (Table S3: % missing per
#              individual variable, observed 76 / predicted 85 / all 161 countries)
#            tables/A5_countries_per_m49_region.csv
#            tables/A5_countries_per_m49_subregion.csv
#            tables/A5_countries_by_tier_m49_region.csv
#            tables/A5_countries_per_income_group.csv
#            tables/A5_population_15plus_by_tier.csv
#
#  What to look out for:
#    * "% missing" in the non-imputed columns shows how much MICE contributed.
#    * Non-imputed and MICE estimates should be close; a large difference means
#      the missing values differ systematically from the observed ones.
#    * 76 observed countries in total; the region counts should add up to 76.
################################################################################

source("00_setup.R")
start_timer("A5")
by_region            <- FALSE  # TRUE to also give estimates by M49 region (svy: mean x, over(region))
check_against_survey <- TRUE   # one cross-check against the survey package (a few minutes)

df_raw <- load_data("df_clean")      # original values (with NAs)
df_mice <- load_data("df_imputed")   # MICE-imputed values (used in the models)
nd      <- load_data("new_data_clean")
wb_raw  <- load_data("wb_clean")
wb_mice <- load_data("wb_imputed")

survey_weight <- "projectionwgt_svy_yr"
if (!survey_weight %in% names(df_raw)) {
  warning(survey_weight, " not found in the IWISE data; using wgtnorm_using instead. ",
          "Pooled (overall / regional) estimates will not be population-weighted.")
  survey_weight <- "wgtnorm_using"
}
message("Survey weight used: ", survey_weight)


# ══════════════════════════════════════════════════════════════════════════════
# 1. INDIVIDUAL-LEVEL PREDICTORS (SURVEY-WEIGHTED)
# ══════════════════════════════════════════════════════════════════════════════
# Each row: display label, column in df_clean (non-imputed), column in
# df_imputed (MICE). Variables that were never imputed use the same column.

ind_desc_vars <- tribble(
  ~label,                           ~raw,                ~mice,
  "Moderate-to-high water insecurity", "iwise12_imp",    "iwise12_imp",
  "Age (years)",                    "age",               "age_imp",
  "Household size",                 "hhsize",            "hhsize",
  "Sex",                            "female",            "female",
  "Income quintile",                "INCOME_5",          "INCOME_5_imp",
  "Feelings about household income", "inadincome",       "inadincome_imp",
  "Education",                      "education",         "education_imp",
  "Marital status",                 "maritalstatus",     "maritalstatus_imp",
  "Employment",                     "employment",        "employment",
  "Dissatisfied with water quality", "waterdissatisfied", "waterdissatisfied_imp",
  "Urbanicity",                     "urban_final",       "urban_final_imp",
  "Source of urbanicity",           "urban_source",      "urban_source",
  "Income missing (flag)",          "income_missing",    "income_missing",
  "Urbanicity missing (flag)",      "urban_missing",     "urban_missing"
)
continuous_vars <- c("age", "age_imp", "hhsize")

# Survey design with strata / PSUs unique within country and survey year
prepare_design <- function(d) {
  d %>%
    filter(!is.na(.data[[survey_weight]])) %>%
    mutate(country_code = as.character(country_code),
           strata_id = paste(country_code, YEAR_WAVE, coalesce(as.character(strata), "1"), sep = "_"),
           psu_id    = paste(strata_id, psu, sep = "_"),
           .w        = .data[[survey_weight]])
}

# Survey-weighted mean and its Taylor-linearised SE, calculated directly.
# Gives the same results as survey::svymean(..., na.rm = TRUE) and svyby() with
# options(survey.lonely.psu = "certainty") (Stata: svy: mean, singleunit(certainty)),
# but in seconds rather than minutes on ~90,000 respondents.
#   domain : TRUE for respondents in the group being estimated (e.g. a region)
#            and with a non-missing value. Respondents outside the domain stay
#            in the design (they count towards the PSUs in each stratum), as in
#            Stata subpop / survey's subset.
#   Strata with a single PSU contribute no variance (certainty units).
svy_mean_fast <- function(y, domain, w, strata, psu) {
  y[!domain] <- 0
  wd  <- ifelse(domain, w, 0)
  tw  <- sum(wd)
  est <- sum(wd * y) / tw
  z   <- wd * (y - est) / tw                       # linearised values
  psu_tot  <- rowsum(z, psu)[, 1]                  # PSU totals
  psu_str  <- strata[match(names(psu_tot), psu)]   # stratum of each PSU
  n_h      <- ave(psu_tot, psu_str, FUN = length)  # PSUs per stratum
  dev      <- psu_tot - ave(psu_tot, psu_str, FUN = mean)
  v        <- sum(ifelse(n_h > 1, n_h / (n_h - 1) * dev^2, 0))
  c(estimate = est, se = sqrt(v))
}

# Weighted percentile, as Stata `summarize x [aw = w], detail`:
# with cumulative weights W(i) of the sorted values and P = W * p, take the
# first i with W(i) > P; if W(i-1) equals P, average x(i-1) and x(i).
stata_wpctile <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  W <- cumsum(w); P <- sum(w) * p
  i <- which(W > P * (1 + 1e-12))[1]
  if (i > 1 && abs(W[i - 1] - P) <= 1e-12 * sum(w)) (x[i - 1] + x[i]) / 2 else x[i]
}

# Weighted SD (as sqrt(coef(survey::svyvar())) and Stata summarize with aweights)
svy_sd_fast <- function(y, domain, w) {
  y <- y[domain]; w <- w[domain]; n <- length(y)
  m <- sum(w * y) / sum(w)
  sqrt(sum(w * (y - m)^2) / sum(w) * n / (n - 1))
}

# Categorical variables are turned into 0/1 indicators (one per level, NA where
# the variable is missing) so that every estimate is a weighted mean.
indicator_columns <- function(d, var) {
  x <- d[[var]]
  if (var %in% continuous_vars) {
    d[[paste0(var, "__value")]] <- as.numeric(x)
    return(list(data = d, cols = paste0(var, "__value"), levels = "Mean"))
  }
  x <- as.character(x)
  lv <- if (is.factor(d[[var]])) levels(d[[var]]) else sort(unique(na.omit(x)))
  binary <- identical(sort(as.character(lv)), c("0", "1"))
  if (binary) lv <- "1"   # 0/1 variables: report % = 1 only
  cols <- paste0(var, "__", seq_along(lv))
  for (i in seq_along(lv)) d[[cols[i]]] <- ifelse(is.na(x), NA_real_, as.numeric(x == lv[i]))
  list(data = d, cols = cols, levels = if (binary) "Yes" else lv)
}

describe_individual <- function(d, version) {
  col_name <- if (version == "Non-imputed") "raw" else "mice"
  vars <- ind_desc_vars[[col_name]]
  missing_cols <- setdiff(vars, names(d))
  if (length(missing_cols) > 0) warning("Not found (", version, "), skipped: ", paste(missing_cols, collapse = ", "))
  keep <- which(vars %in% names(d))

  # Create all 0/1 indicator columns first, then prepare the design variables ONCE
  ics <- list()
  for (k in keep) {
    ic <- indicator_columns(d, vars[k])
    d <- ic$data
    ics[[as.character(k)]] <- ic
  }
  dd <- prepare_design(d)
  n_total <- nrow(dd)
  w <- dd$.w; strata <- dd$strata_id; psu <- dd$psu_id
  region <- as.character(dd$m49_region)
  regions <- sort(unique(na.omit(region)))

  out <- list()
  for (k in keep) {
    var <- vars[k]; ic <- ics[[as.character(k)]]
    message("  ", version, ": ", ind_desc_vars$label[k])
    is_cont   <- var %in% continuous_vars
    n_missing <- sum(is.na(dd[[var]]))
    has_value <- !is.na(dd[[var]])

    # As Stata `svy: mean x`: respondents with x missing are dropped from the
    # design for that variable; regions as `svy: mean x, over(region)`.
    keep_rows <- has_value
    res <- list()
    for (j in seq_along(ic$cols)) {
      y <- dd[[ic$cols[j]]][keep_rows]
      all_rows <- rep(TRUE, length(y))
      o <- svy_mean_fast(y, all_rows, w[keep_rows], strata[keep_rows], psu[keep_rows])
      # Continuous variables: also SD, median and IQR (as Stata summarize, detail)
      dist_stats <- function(in_group) {
        if (!is_cont) return(tibble(sd = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_))
        yy <- y[in_group]; ww <- w[keep_rows][in_group]
        tibble(sd = svy_sd_fast(yy, rep(TRUE, length(yy)), ww),
               median = stata_wpctile(yy, ww, 0.50),
               q1 = stata_wpctile(yy, ww, 0.25), q3 = stata_wpctile(yy, ww, 0.75))
      }
      res[[length(res) + 1]] <- tibble(col = ic$cols[j], group = "Overall",
                                       estimate = o[["estimate"]], se = o[["se"]]) %>%
        bind_cols(dist_stats(all_rows))
      for (r in if (by_region) regions else character(0)) {
        in_r <- region[keep_rows] %in% r
        o <- svy_mean_fast(y, in_r, w[keep_rows], strata[keep_rows], psu[keep_rows])
        res[[length(res) + 1]] <- tibble(col = ic$cols[j], group = r,
                                         estimate = o[["estimate"]], se = o[["se"]]) %>%
          bind_cols(dist_stats(in_r))
      }
    }
    overall <- bind_rows(res); by_reg <- NULL

    out[[length(out) + 1]] <- bind_rows(overall, by_reg) %>%
      mutate(level_order = match(col, ic$cols),
             level = ic$levels[level_order],
             version = version, variable = ind_desc_vars$label[k], column = var,
             type = if (is_cont) "mean" else "percent",
             n_respondents = n_total, n_missing = n_missing,
             pct_missing = 100 * n_missing / n_total) %>%
      select(-col)
  }
  bind_rows(out) %>%
    mutate(ci_lower = estimate - 1.96 * se,
           ci_upper = estimate + 1.96 * se,
           across(c(estimate, se, sd, ci_lower, ci_upper),
                  ~ ifelse(type == "percent", 100 * .x, .x)),
           ci_lower = ifelse(type == "percent", pmax(ci_lower, 0), ci_lower),
           ci_upper = ifelse(type == "percent", pmin(ci_upper, 100), ci_upper))
}

ind_long <- bind_rows(describe_individual(df_raw,  "Non-imputed"),
                      describe_individual(df_mice, "MICE")) %>%
  mutate(group = factor(group, levels = c("Overall", m49_region_levels))) %>%
  arrange(factor(variable, levels = ind_desc_vars$label), level_order, group, desc(version)) %>%
  select(variable, level, type, group, version, estimate, se, ci_lower, ci_upper, sd, median, q1, q3,
         n_respondents, n_missing, pct_missing, column, level_order)
# CHECK: the direct calculation should match the survey package exactly.
# Compares the outcome (overall and by region) and mean age (non-imputed).
if (check_against_survey) {
  load_pkgs("survey")
  options(survey.lonely.psu = "certainty")
  dd  <- prepare_design(df_raw) %>%
    mutate(outcome = as.numeric(as.character(iwise12_imp)), age_num = as.numeric(age))
  des <- svydesign(ids = ~psu_id, strata = ~strata_id, weights = ~.w, data = dd, nest = TRUE)
  des_age <- svydesign(ids = ~psu_id, strata = ~strata_id, weights = ~.w,
                       data = filter(dd, !is.na(age_num)), nest = TRUE)   # Stata drops missing
  ref <- bind_rows(
    tibble(variable = "Moderate-to-high water insecurity", group = "Overall",
           estimate = 100 * coef(svymean(~outcome, des)), se = 100 * SE(svymean(~outcome, des))),
    if (by_region) svyby(~outcome, ~m49_region, des, svymean) %>% as_tibble() %>%
      transmute(variable = "Moderate-to-high water insecurity", group = as.character(m49_region),
                estimate = 100 * outcome, se = 100 * se),
    { m <- svymean(~age_num, des_age)
      tibble(variable = "Age (years)", group = "Overall", estimate = coef(m)[[1]], se = SE(m)[[1]]) })
  cmp <- ind_long %>%
    filter(version == "Non-imputed", variable %in% ref$variable, level %in% c("Yes", "Mean")) %>%
    mutate(group = as.character(group)) %>%
    inner_join(ref, by = c("variable", "group"), suffix = c("", "_survey"))
  max_diff <- max(abs(c(cmp$estimate - cmp$estimate_survey, cmp$se - cmp$se_survey)))
  message("Check against survey package (", nrow(cmp), " estimates): max difference = ", signif(max_diff, 3))
  if (max_diff > 1e-6) warning("Direct calculation differs from the survey package — please report.")
}

save_table(ind_long, "A5_individual_descriptives_long")

# Overall table: non-imputed and MICE side by side, formatted "est (lower–upper)"
fmt_ci <- function(e, l, u, digits) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), e, l, u)
}
ind_overall <- ind_long %>%
  filter(group == "Overall") %>%
  mutate(value = fmt_ci(estimate, ci_lower, ci_upper, 1),
         value = ifelse(type == "mean", paste0(value, "; SD ", sprintf("%.1f", sd)), value),
         med_iqr = ifelse(type == "mean", sprintf("%g (%g-%g)", median, q1, q3), "")) %>%
  select(variable, level, type, version, value, med_iqr, pct_missing) %>%
  pivot_wider(names_from = version, values_from = c(value, med_iqr, pct_missing)) %>%
  select(variable, level, type,
         `Non-imputed: % or mean (95% CI)` = `value_Non-imputed`,
         `Non-imputed: median (IQR)`       = `med_iqr_Non-imputed`,
         `Non-imputed: % missing`          = `pct_missing_Non-imputed`,
         `MICE: % or mean (95% CI)`        = value_MICE,
         `MICE: median (IQR)`              = med_iqr_MICE) %>%
  mutate(`Non-imputed: % missing` = round(`Non-imputed: % missing`, 1))
save_table(ind_overall, "A5_tableS01A_individual_descriptives")   # Table S1A
print(ind_overall, n = Inf, width = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 2. NATIONAL-LEVEL PREDICTORS (ONE VALUE PER COUNTRY, 76 COUNTRIES)
# ══════════════════════════════════════════════════════════════════════════════
# Unweighted summaries across countries, on the scale used in the models.

observed_iso3 <- unique(as.character(df_raw$country_code))

describe_national <- function(wb, version) {
  wb %>%
    filter(iso3 %in% observed_iso3) %>%
    select(iso3, all_of(country_vars)) %>%
    pivot_longer(-iso3, names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(n_countries = sum(!is.na(value)),
              n_missing   = sum(is.na(value)),
              mean   = mean(value, na.rm = TRUE),
              sd     = sd(value, na.rm = TRUE),
              median = median(value, na.rm = TRUE),
              q1     = quantile(value, 0.25, na.rm = TRUE),
              q3     = quantile(value, 0.75, na.rm = TRUE),
              min    = min(value, na.rm = TRUE),
              max    = max(value, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(version = version)
}

national <- bind_rows(describe_national(wb_raw,  "Non-imputed"),
                      describe_national(wb_mice, "MICE")) %>%
  mutate(label = relabel(variable),
         variable = factor(variable, levels = country_vars),
         version = factor(version, levels = c("Non-imputed", "MICE"))) %>%
  arrange(variable, version) %>%
  select(variable, label, version, n_countries, n_missing, mean, sd, median, q1, q3, min, max) %>%
  mutate(across(mean:max, ~ round(.x, 3)))
save_table(national, "A5_tableS01B_national_descriptives")         # Table S1B
print(national, n = Inf, width = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 3. NUMBER OF COUNTRIES PER M49 REGION / SUB-REGION
# ══════════════════════════════════════════════════════════════════════════════

obs_countries <- df_raw %>%
  mutate(country_code = as.character(country_code)) %>%
  group_by(country_code, m49_region, m49_sub_region) %>%
  summarise(n_respondents = n(), .groups = "drop")

per_region <- obs_countries %>%
  group_by(m49_region) %>%
  summarise(n_countries = n(), n_respondents = sum(n_respondents), .groups = "drop") %>%
  mutate(m49_region = factor(m49_region, levels = m49_region_levels)) %>%
  arrange(m49_region) %>%
  bind_rows(tibble(m49_region = "Total", n_countries = sum(.$n_countries),
                   n_respondents = sum(.$n_respondents)))
save_table(per_region, "A5_countries_per_m49_region")
print(per_region)

per_subregion <- obs_countries %>%
  group_by(m49_region, m49_sub_region) %>%
  summarise(n_countries = n(), n_respondents = sum(n_respondents),
            countries = paste(sort(country_code), collapse = ", "), .groups = "drop") %>%
  arrange(factor(m49_region, levels = m49_region_levels), m49_sub_region)
save_table(per_subregion, "A5_countries_per_m49_subregion")

# All 194 countries by data source (tier) and M49 region
tier_by_region <- wb_raw %>%
  transmute(iso3, m49_region,
            tier = case_when(iso3 %in% observed_iso3 ~ "Observed",
                             iso3 %in% unique(as.character(nd$country_code)) ~ "Predicted (individual)",
                             TRUE ~ "Predicted (country-level)")) %>%
  count(m49_region, tier) %>%
  pivot_wider(names_from = tier, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  arrange(factor(m49_region, levels = m49_region_levels))
tier_by_region <- bind_rows(tier_by_region,
                            tier_by_region %>% summarise(m49_region = "Total", across(where(is.numeric), sum)))
save_table(tier_by_region, "A5_countries_by_tier_m49_region")
print(tier_by_region)

# Observed countries per World Bank income group (FY27)
per_income <- wb_raw %>%
  filter(iso3 %in% observed_iso3) %>%
  mutate(wb_income_fy27 = coalesce(as.character(wb_income_fy27), "Unclassified")) %>%
  count(wb_income_fy27, name = "n_countries") %>%
  arrange(desc(n_countries))
per_income <- bind_rows(per_income, tibble(wb_income_fy27 = "Total", n_countries = sum(per_income$n_countries)))
save_table(per_income, "A5_countries_per_income_group")
print(per_income)


# ══════════════════════════════════════════════════════════════════════════════
# 4. POPULATION 15+ COVERED BY EACH TIER
# ══════════════════════════════════════════════════════════════════════════════
# pop_15plus is the 2025 population aged 15+ (the weight used for the estimates).

tier2_iso3 <- unique(as.character(nd$country_code))
tier_of <- function(iso3) case_when(iso3 %in% observed_iso3 ~ "Tier 1 (Observed)",
                                    iso3 %in% tier2_iso3    ~ "Tier 2 (Predicted individual)",
                                    TRUE                    ~ "Tier 3 (Predicted country)")

pop_by_tier <- wb_raw %>%
  mutate(tier = tier_of(iso3)) %>%
  group_by(tier) %>%
  summarise(n_countries = n(), pop_15plus = sum(pop_15plus, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_of_total = 100 * pop_15plus / sum(pop_15plus))
pop_by_tier <- bind_rows(pop_by_tier,
                         tibble(tier = "Total", n_countries = sum(pop_by_tier$n_countries),
                                pop_15plus = sum(pop_by_tier$pop_15plus), pct_of_total = 100)) %>%
  mutate(pop_15plus_millions = round(pop_15plus / 1e6), pct_of_total = round(pct_of_total, 1))
save_table(pop_by_tier, "A5_population_15plus_by_tier")

for (i in seq_len(nrow(pop_by_tier))) {
  cat(sprintf("%-31s %5d countries  %6.0f M  %5.1f%%\n", pop_by_tier$tier[i], pop_by_tier$n_countries[i],
              pop_by_tier$pop_15plus[i] / 1e6, pop_by_tier$pct_of_total[i]))
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. TABLE S3 — MISSING DATA FOR THE INDIVIDUAL-LEVEL PREDICTORS
# ══════════════════════════════════════════════════════════════════════════════
# % of respondents missing each variable BEFORE imputation (unweighted), in the
# observed countries (df_clean), the countries predicted by the individual-level
# model (new_data_clean) and both pooled (respondents, not an average of the two
# columns). The largest value in the "All" column is the figure quoted in the
# methods (no variable missing for more than ~2.3% of respondents).
# The MICE method (A4) is kept in the last column for reference.

var_info <- tribble(
  ~variable,           ~Variable,                  ~`Variable type`,              ~mice_method,
  "age",               "Age",                      "Continuous (years)",          "pmm",
  "hhsize",            "Household size",           "Continuous (count)",          "pmm",
  "INCOME_5",          "Income quintile",          "Ordinal categorical (1–5)",   "pmm (imputed as numeric)",
  "inadincome",        "Inadequate income",        "Binary",                      "polyreg (= logistic for 2 categories)",
  "waterdissatisfied", "Dissatisfied with water",  "Binary",                      "logreg",
  "education",         "Education",                "Categorical (3 levels)",      "polyreg",
  "maritalstatus",     "Marital status",           "Categorical (3 levels)",      "polyreg",
  "employment",        "Employment",               "Categorical (4 levels)",      "polyreg",
  "urban_final",       "Urbanicity",               "Categorical (3 levels)",      "polyreg"
)

# Original (pre-imputation) values under the same names in both datasets
miss_obs  <- df_raw %>%
  select(any_of(c(paste0(var_info$variable, "_imp"), "employment", "hhsize"))) %>%
  rename_with(~ sub("_imp$", "", .x)) %>%
  mutate(across(everything(), as.character))
miss_pred <- nd %>% select(any_of(var_info$variable)) %>% mutate(across(everything(), as.character))

miss_pct <- function(d, label) {
  d %>%
    summarise(across(everything(), ~ round(100 * mean(is.na(.x)), 2))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = label)
}
lab_obs  <- paste0("Observed, 76 countries (n = ", format(nrow(miss_obs), big.mark = ","), "), % missing")
lab_pred <- paste0("Predicted individual-level, 85 countries (n = ", format(nrow(miss_pred), big.mark = ","), "), % missing")
lab_all  <- paste0("All, 161 countries (n = ", format(nrow(miss_obs) + nrow(miss_pred), big.mark = ","), "), % missing")

tableS3 <- var_info %>%
  left_join(miss_pct(miss_obs, lab_obs), by = "variable") %>%
  left_join(miss_pct(miss_pred, lab_pred), by = "variable") %>%
  left_join(miss_pct(bind_rows(miss_obs, miss_pred), lab_all), by = "variable") %>%
  select(Variable, `Variable type`, all_of(c(lab_obs, lab_pred, lab_all)), `MICE method (A4)` = mice_method)

# CHECK: every variable found in both datasets (expect no warning)
missing_vars <- setdiff(var_info$variable, intersect(names(miss_obs), names(miss_pred)))
if (length(missing_vars) > 0) warning("Not found in both datasets: ", paste(missing_vars, collapse = ", "))
readr::write_excel_csv(tableS3, file.path(dir_tables, "A5_tableS03_individual_missingness.csv"))
print(tableS3, width = Inf)
message("Largest % missing (all 161 countries): ", max(tableS3[[lab_all]], na.rm = TRUE), "%")


# CHECK: 76 observed countries, 194 in total
stopifnot(sum(per_region$n_countries[per_region$m49_region != "Total"]) == length(observed_iso3))
message("Observed countries: ", length(observed_iso3),
        " | all countries: ", sum(tier_by_region$Total[tier_by_region$m49_region != "Total"]))

stop_timer()
