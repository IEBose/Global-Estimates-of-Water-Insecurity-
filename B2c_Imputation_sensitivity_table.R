################################################################################
#  B2c_Imputation_sensitivity_table.R — SUPPLEMENT TABLE FOR THE B2b
#                                       IMPUTATION SENSITIVITY ANALYSIS
#
#  Purpose : One table for the supplement from B2b's outputs: for the Tier 2
#            countries with individual-level variables that have NO observed
#            values (imputed entirely by MICE from other countries; flagged in
#            Table S9), the main estimate with its 95% UI and the range of the
#            estimate across the 5 MICE imputations
#              - prediction side (final model, 5 imputed Tier 2 datasets)
#              - training side (model refitted on 5 imputed training datasets)
#            plus one summary row for all other Tier 2 countries.
#            Quick (seconds): reads tables only, no model is re-run.
#            Does NOT change any result used in the paper.
#
#  Depends : B2b (by-imputation tables), A4 (A4_vars_fully_missing_by_country),
#            B11 (Table S10, for the main 95% UIs; run after B11)
#  Outputs : tables/B2c_table_imputation_sensitivity.csv   (supplement table)
#            tables/B2c_imputation_sensitivity_long.csv    (all 85 countries,
#                                                           both sides; check)
#  Not numbered yet (supplement numbering to be decided).
################################################################################

source("00_setup.R")
start_timer("B2c")

f_pred <- file.path(dir_tables, "B2b_tier2_by_prediction_imputation.csv")
f_train <- file.path(dir_tables, "B2b_tier2_by_training_imputation.csv")
f_miss <- file.path(dir_tables, "A4_vars_fully_missing_by_country.csv")
f_s10  <- file.path(dir_tables, "B11_tableS10_all_country_estimates.csv")
for (f in c(f_pred, f_train, f_miss)) if (!file.exists(f)) stop("Not found (run A4 and B2b first): ", f)

by_pred  <- readr::read_csv(f_pred,  show_col_types = FALSE) %>% mutate(country_code = as.character(country_code))
by_train <- readr::read_csv(f_train, show_col_types = FALSE) %>% mutate(country_code = as.character(country_code))

# Variables with no observed values, with readable labels (as in Table S9, B8)
ind_var_labels <- c(age = "age", INCOME_5 = "income quintile", inadincome = "income adequacy",
                    education = "education", waterdissatisfied = "water dissatisfaction",
                    maritalstatus = "marital status", urban_final = "urbanicity",
                    employment = "employment", hhsize = "household size")
fully_missing <- readr::read_csv(f_miss, show_col_types = FALSE) %>%
  mutate(country_code = as.character(country_code),
         label = coalesce(unname(ind_var_labels[variable]), variable)) %>%
  group_by(country_code) %>%
  summarise(vars = paste(label, collapse = ", "), .groups = "drop")

s10 <- if (file.exists(f_s10)) {
  readr::read_csv(f_s10, show_col_types = FALSE) %>%
    mutate(country_code = as.character(country_code)) %>%
    select(country_code, any_of("country_name"), prev_best, ui_lower, ui_upper)
} else {
  message("Main Table S10 not found (run B11): table written without UIs")
  tibble(country_code = character(), prev_best = numeric(), ui_lower = numeric(), ui_upper = numeric())
}

pct <- function(x, d = 1) ifelse(is.na(x), NA_character_, sprintf(paste0("%.", d, "f%%"), 100 * x))

# ══════════════════════════════════════════════════════════════════════════════
# 1. BOTH SIDES TOGETHER, ALL 85 TIER 2 COUNTRIES
# ══════════════════════════════════════════════════════════════════════════════

dat <- by_pred %>%
  select(country_code, pred_min = min, pred_max = max, pred_range_pp = range_pp) %>%
  inner_join(by_train %>% select(country_code, train_min = min, train_max = max, train_range_pp = range_pp),
             by = "country_code") %>%
  left_join(fully_missing, by = "country_code") %>%
  left_join(s10, by = "country_code") %>%
  mutate(flagged   = !is.na(vars),
         within_ui = pmin(pred_min, train_min) >= ui_lower & pmax(pred_max, train_max) <= ui_upper)
if (!"country_name" %in% names(dat)) dat$country_name <- dat$country_code
save_table(dat, "B2c_imputation_sensitivity_long")

# CHECK: 85 Tier 2 countries on both sides
message("Tier 2 countries: ", nrow(dat), " (expect 85) | flagged: ", sum(dat$flagged),
        " (", paste(dat$country_code[dat$flagged], collapse = ", "), ")")

# ══════════════════════════════════════════════════════════════════════════════
# 2. SUPPLEMENT TABLE
# ══════════════════════════════════════════════════════════════════════════════

flagged_rows <- dat %>%
  filter(flagged) %>%
  arrange(desc(pmax(pred_range_pp, train_range_pp))) %>%
  transmute(Country = paste0(country_name, " (", country_code, ")"),
            `Variables with no observed values (imputed by MICE from other countries)` = vars,
            `Main estimate (95% UI)` = ifelse(is.na(prev_best), NA_character_,
                                              paste0(pct(prev_best), " (", pct(ui_lower), "–", pct(ui_upper), ")")),
            `Prediction side: estimate across 5 imputations` = paste0(pct(pred_min), "–", pct(pred_max)),
            `Prediction side: range (pp)` = sprintf("%.2f", pred_range_pp),
            `Training side: estimate across 5 imputations` = paste0(pct(train_min), "–", pct(train_max)),
            `Training side: range (pp)` = sprintf("%.2f", train_range_pp),
            `Within main 95% UI` = ifelse(is.na(within_ui), NA_character_, ifelse(within_ui, "Yes", "No")))

others <- filter(dat, !flagged)
summ_range <- function(r, cc) sprintf("median %.2f; max %.2f (%s)", median(r), max(r), cc[which.max(r)])
other_row <- tibble(
  Country = paste0("All other Tier 2 countries (n = ", nrow(others), ")"),
  `Variables with no observed values (imputed by MICE from other countries)` = "none",
  `Main estimate (95% UI)` = "",
  `Prediction side: estimate across 5 imputations` = "",
  `Prediction side: range (pp)` = summ_range(others$pred_range_pp, others$country_code),
  `Training side: estimate across 5 imputations` = "",
  `Training side: range (pp)` = summ_range(others$train_range_pp, others$country_code),
  `Within main 95% UI` = if (all(is.na(others$within_ui))) "" else
    paste0(sum(others$within_ui, na.rm = TRUE), " of ", nrow(others)))

table_out <- bind_rows(flagged_rows, other_row)
readr::write_excel_csv(table_out, file.path(dir_tables, "B2c_table_imputation_sensitivity.csv"))
print(table_out, n = Inf, width = Inf)

# Numbers for the text
message("All 85 Tier 2 countries - prediction side: median ", sprintf("%.2f", median(dat$pred_range_pp)),
        " pp, max ", sprintf("%.2f", max(dat$pred_range_pp)), " (", dat$country_code[which.max(dat$pred_range_pp)], ")",
        " | training side: median ", sprintf("%.2f", median(dat$train_range_pp)),
        " pp, max ", sprintf("%.2f", max(dat$train_range_pp)), " (", dat$country_code[which.max(dat$train_range_pp)], ")",
        " | > 1 pp on either side: ", sum(dat$pred_range_pp > 1 | dat$train_range_pp > 1),
        " | outside the main 95% UI: ", sum(!dat$within_ui, na.rm = TRUE))

stop_timer(optional = TRUE)
