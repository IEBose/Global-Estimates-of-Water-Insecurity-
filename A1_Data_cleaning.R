################################################################################
#  A1_Data_cleaning.R — LOAD AND CLEAN THE THREE ANALYSIS DATASETS
#
#  Purpose : Load the pre-made input files (no API calls) and create three
#            clean datasets for the rest of the pipeline:
#              df        76 observed countries (individual IWISE + GWP data)
#              new_data  85 prediction countries (GWP data, no IWISE)
#              wb        194 countries, national indicators (NOT imputed)
#
#  Inputs  : World Bank snapshot + population file (names in 00_setup.R)
#            IWISE_working_dataset_*.dta, Gallup_World_Poll_*_noWISE.dta
#  Outputs : OUTPUT/data/df_clean.rds, new_data_clean.rds, wb_clean.rds
#            OUTPUT/tables/A1_*.csv (missingness and check tables)
#            OUTPUT/tables/A1_tableS02_national_indicators_coverage.csv (Table S2,
#              copied from A0's coverage file for the same snapshot)
#
#  What to look out for (messages printed while running):
#    * "76 / 85 / 33 countries" — the three tiers must add to 194.
#    * Missing M49 region or model region — should be none.
#    * tables/A1_model_region_assignment.csv — each survey country's M49
#      sub-region and model region (check that they look sensible).
################################################################################

source("00_setup.R")
start_timer("A1")


# ══════════════════════════════════════════════════════════════════════════════
# 1. NATIONAL DATA (194 COUNTRIES) — FROZEN WORLD BANK SNAPSHOT
# ══════════════════════════════════════════════════════════════════════════════

if (wb_snapshot_stamp == "YYYYMMDD") {
  stop("Set wb_snapshot_stamp in 00_setup.R to the date of the World Bank snapshot used.")
}

wb <- readr::read_csv(file_wb, show_col_types = FALSE)

# ND-GAIN index is used as the climate vulnerability predictor (`cvi`)
wb$cvi <- wb$nd_gain

# Population 15+ for the fixed reporting year (pop_year in 00_setup.R).
# Replaces the 'latest year' pop_15plus written by A0, so that every country
# is weighted by population in the same year.
pop <- readr::read_csv(file_wb_pop, show_col_types = FALSE) %>%
  select(iso3, pop_15plus = all_of(paste0("pop_15plus_", pop_year)))

wb <- wb %>%
  rename(pop_15plus_latest = pop_15plus) %>%
  left_join(pop, by = "iso3")

# Model region from M49 (used for imputation and the country-level model)
wb <- wb %>%
  mutate(region_final = assign_region_final_m49(iso3, m49_region, m49_sub_region,
                                                m49_intermediate_region))

# CHECKS — all should be 0
stopifnot(nrow(wb) == 194, !anyDuplicated(wb$iso3))
message("wb: countries missing M49 region: ", sum(is.na(wb$m49_region)),
        " | sub-region: ", sum(is.na(wb$m49_sub_region)),
        " | model region: ", sum(is.na(wb$region_final)),
        " | pop 15+ (", pop_year, "): ", sum(is.na(wb$pop_15plus)))

# Correlation between candidate national predictors (for checking collinearity;
# |r| > 0.8 suggests two indicators carry largely the same information)
cor_mat <- wb %>%
  select(log_gni_pc_ppp, cvi, fertility, primary_complete, log_pop_density,
         cpi_score, food_insec_modsev, u5_mort, wash_mort, life_exp,
         log_water_stress, precip_mm, basic_water) %>%
  cor(use = "pairwise.complete.obs")
save_table(as_tibble(cor_mat, rownames = "variable"), "A1_national_predictor_correlations")


# ══════════════════════════════════════════════════════════════════════════════
# 2. OBSERVED INDIVIDUAL DATA (76 COUNTRIES)
# ══════════════════════════════════════════════════════════════════════════════

iwise_raw <- read_dta_keep(file_iwise_dta, keep_iwise_required, keep_iwise_optional)

# merge() (rather than left_join) is kept on purpose: it sorts rows by country
# code, and the row order affects the random folds in cv.glmnet / xgb.cv and the
# MICE draws. Changing it would change results slightly.
df <- merge(iwise_raw, wb, by.x = "country_code", by.y = "iso3", all.x = TRUE)

# CHECK: countries that did not match the national data (should be none)
unmatched <- unique(df$country_code[is.na(df$life_exp)])
if (length(unmatched) > 0) warning("No national data for: ", paste(unmatched, collapse = ", "))

# Strip Stata labels so values are plain numbers
df <- df %>% mutate(across(where(haven::is.labelled), ~ as.numeric(haven::zap_labels(.))))

# Drop rows without a survey weight (this removes the 2025 Palestine round;
# the 2022 round is used for Palestine)
n0 <- nrow(df)
df <- df %>% filter(!is.na(wgtnorm_using))
message("df: dropped ", n0 - nrow(df), " rows with missing weight")


# ══════════════════════════════════════════════════════════════════════════════
# 3. PREDICTION INDIVIDUAL DATA (85 COUNTRIES, NO IWISE)
# ══════════════════════════════════════════════════════════════════════════════

nowise_raw <- read_dta_keep(file_nowise_dta, keep_nowise_required, keep_nowise_optional)

new_data <- merge(nowise_raw, wb, by.x = "COUNTRY_ISO3", by.y = "iso3", all.x = TRUE)
new_data$country_code <- new_data$COUNTRY_ISO3

new_data <- new_data %>% mutate(across(where(haven::is.labelled), ~ as.numeric(haven::zap_labels(.))))

# Drop non-UN territories (only the 194 are estimated)
dropped_nonUN <- unique(new_data$country_code[new_data$non_UN %in% 1])
new_data <- new_data %>% filter(non_UN != 1 | is.na(non_UN))
message("new_data: dropped non-UN territories: ", paste(dropped_nonUN, collapse = ", "))

# Flag countries whose latest GWP round is before 2015 (estimates to treat with caution)
new_data <- new_data %>% mutate(pre_2015 = ifelse(YEAR_WAVE < 2015, 1, 0))
pre_2015_countries <- new_data %>% filter(pre_2015 == 1) %>% distinct(country_code, YEAR_WAVE)
message("new_data: countries with pre-2015 data (methods: n = 7): ", nrow(pre_2015_countries))
save_table(pre_2015_countries, "A1_new_data_pre2015_countries")

n0 <- nrow(new_data)
new_data <- new_data %>% filter(!is.na(wgtnorm_using))
message("new_data: dropped ", n0 - nrow(new_data), " rows with missing weight")


# ══════════════════════════════════════════════════════════════════════════════
# 4. MODEL REGION (region_final)
# ══════════════════════════════════════════════════════════════════════════════
# region_final is the 9-category MODEL region used as a predictor. It comes
# with the national data (Section 1), where it is built from UN M49 with
# assign_region_final_m49() in 00_setup.R, so the observed, prediction and
# country-level data all use the same rule. ("MENA" includes Central Asia.)

message("Missing model region — df: ", sum(is.na(df$region_final)),
        " | new_data: ", sum(is.na(new_data$region_final)))

# Record which model region each survey country is in (for checking)
region_assignment <- bind_rows(
  df       %>% distinct(country_code, m49_sub_region, region_final) %>% mutate(data = "observed"),
  new_data %>% distinct(country_code, m49_sub_region, region_final) %>% mutate(data = "prediction")
) %>% arrange(region_final, country_code)
save_table(region_assignment, "A1_model_region_assignment")


# ══════════════════════════════════════════════════════════════════════════════
# 5. URBANICITY
# ══════════════════════════════════════════════════════════════════════════════
# DEGURBA is the primary measure. Where missing (large share in Burkina Faso and
# Rep. Congo), self-reported urbanicity is used: rural -> rural; suburb or small
# town -> towns/semi-dense; large city -> cities. urban_source records which.

harmonise_urban <- function(d) {
  d %>%
    mutate(
      selfrep_harmonised = case_when(urbanscale == 0 ~ 1L,
                                     urbanscale %in% c(1, 2) ~ 2L,
                                     urbanscale == 3 ~ 3L,
                                     TRUE ~ NA_integer_),
      urban_final = case_when(!is.na(degurba)            ~ as.integer(degurba),
                              !is.na(selfrep_harmonised) ~ selfrep_harmonised,
                              TRUE                       ~ NA_integer_),
      urban_final  = factor(urban_final, levels = 1:3, labels = lv_urban),
      urban_source = factor(case_when(!is.na(degurba)            ~ "degurba",
                                      !is.na(selfrep_harmonised) ~ "selfrep",
                                      TRUE                       ~ "missing"),
                            levels = lv_urb_src)
    ) %>%
    select(-selfrep_harmonised)
}
df       <- harmonise_urban(df)
new_data <- harmonise_urban(new_data)


# ══════════════════════════════════════════════════════════════════════════════
# 6. FACTORS, MISSINGNESS FLAGS AND _imp COLUMNS
# ══════════════════════════════════════════════════════════════════════════════

set_factors <- function(d) {
  d %>%
    mutate(
      female            = factor(female,            levels = c(0, 1),    labels = lv_female),
      inadincome        = factor(inadincome,        levels = c(0, 1),    labels = lv_inadinc),
      waterdissatisfied = factor(waterdissatisfied, levels = c(0, 1),    labels = lv_wdis),
      INCOME_5          = factor(INCOME_5,          levels = 1:5),
      education         = factor(education,         levels = c(0, 1, 2), labels = lv_edu),
      maritalstatus     = factor(maritalstatus,     levels = c(0, 1, 2), labels = lv_mar),
      employment        = factor(employment,        levels = 0:3,        labels = lv_emp),
      country_code      = factor(country_code)
    )
}
df       <- set_factors(df) %>% mutate(iwise12_imp = factor(iwise12_imp, levels = c(0, 1)))
new_data <- set_factors(new_data)

# Missingness flags are created HERE, from the original values, before any
# imputation (A4). They are model predictors, so they must be defined the same
# way in the training data (df) and the prediction data (new_data).
add_flags_and_imp_cols <- function(d) {
  d %>%
    mutate(
      income_missing = factor(ifelse(is.na(INCOME_5), 1, 0), levels = c(0, 1)),
      urban_missing  = factor(ifelse(is.na(urban_final), 1, 0), levels = c(0, 1)),
      # _imp columns start as copies of the originals (with NAs); A4 fills them
      age_imp               = age,
      INCOME_5_imp          = INCOME_5,
      inadincome_imp        = inadincome,
      education_imp         = education,
      maritalstatus_imp     = maritalstatus,
      waterdissatisfied_imp = waterdissatisfied,
      urban_final_imp       = urban_final,
      # national predictors under the names used by the individual-level model
      log_GNI_imp    = log_gni_pc_ppp,
      Corruption_imp = cpi_score,
      Mort_WASH_imp  = wash_mort,
      CVI_imp        = cvi,
      Drinking_basic = basic_water,
      # un-imputed copies used by the optional single-imputation scripts (A2/A3)
      log_GNI        = log_gni_pc_ppp,
      Corruption     = cpi_score
    )
}
df       <- add_flags_and_imp_cols(df)
new_data <- add_flags_and_imp_cols(new_data)

# Rows with a missing outcome (respondents missing 4+ IWISE items) are removed.
# Methods: n = 426.
n0 <- nrow(df)
df <- df %>% filter(!is.na(iwise12_imp))
message("df: removed ", n0 - nrow(df), " rows with missing outcome (methods: 426)")


# ══════════════════════════════════════════════════════════════════════════════
# 7. TEST / TRAIN SPLIT (COUNTRY HOLD-OUT)
# ══════════════════════════════════════════════════════════════════════════════
# India (most sampled country) + 5 randomly selected countries are held out as
# the test set; the remaining 70 countries are used for training and LCO CV.
# The seed is set immediately before sampling so the same 5 are drawn each time.

focal_holdout_country <- "IND"
other_countries <- setdiff(levels(df$country_code), focal_holdout_country)
set.seed(seed)
holdout_other_set <- sample(other_countries, 5)

# CHECK: these were the hold-out countries in the analysis. If the list of
# observed countries changes, a different set would be drawn — this warns.
expected_holdout <- c("NAM", "CMR", "SOM", "PSE", "BRA")
if (!setequal(holdout_other_set, expected_holdout)) {
  warning("Hold-out countries differ from the analysis (", paste(expected_holdout, collapse = ", "),
          "): now ", paste(holdout_other_set, collapse = ", "))
}

df <- df %>%
  mutate(split = case_when(country_code == focal_holdout_country ~ "holdout_focal",
                           country_code %in% holdout_other_set   ~ "holdout_other",
                           TRUE                                   ~ "train"))


# ══════════════════════════════════════════════════════════════════════════════
# 8. CHECKS AND MISSINGNESS SUMMARIES
# ══════════════════════════════════════════════════════════════════════════════

obs_codes  <- unique(as.character(df$country_code))
pred_codes <- unique(as.character(new_data$country_code))
tier3_n    <- sum(!wb$iso3 %in% c(obs_codes, pred_codes))

# Expect 76 / 85 / 33 and no overlap between observed and prediction countries
message("Countries — observed: ", length(obs_codes), " | individual predictions: ",
        length(pred_codes), " | country-level only: ", tier3_n,
        " | in both df and new_data: ", length(intersect(obs_codes, pred_codes)))

# Each observed country should have a single survey year (should print 0).
# If not, normalised weights and the Stata projection weights (normalised
# weight x population in the survey year) would give different prevalences.
multi_year_obs <- df %>% distinct(country_code, YEAR_WAVE) %>% count(country_code) %>% filter(n > 1)
message("Observed countries with more than one survey year: ", nrow(multi_year_obs))

# Each prediction country should have a single GWP round (should print 0)
multi_wave <- new_data %>% distinct(country_code, YEAR_WAVE) %>% count(country_code) %>% filter(n > 1)
message("new_data countries with more than one GWP round: ", nrow(multi_wave))

# Missingness by variable (before imputation). Methods: < 2.3% for any variable.
miss_summary <- function(d, label) {
  d %>%
    summarise(across(all_of(c("age", "hhsize", "female", "INCOME_5", "inadincome",
                              "education", "maritalstatus", "employment",
                              "waterdissatisfied", "urban_final")), ~ mean(is.na(.)) * 100)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "pct_missing") %>%
    mutate(dataset = label)
}
save_table(bind_rows(miss_summary(df, "df (observed)"), miss_summary(new_data, "new_data (prediction)")),
           "A1_individual_missingness")

save_table(wb %>%
             summarise(across(all_of(country_vars), ~ sum(is.na(.)))) %>%
             pivot_longer(everything(), names_to = "variable", values_to = "n_countries_missing"),
           "A1_national_missingness")


# ══════════════════════════════════════════════════════════════════════════════
# 9. KEEP ONLY VARIABLES USED LATER, AND SAVE
# ══════════════════════════════════════════════════════════════════════════════

keep_common <- c("country_code", "wgtnorm_using", "projectionwgt_svy_yr", "psu", "strata", "YEAR_WAVE",
                 "region_final", "un_country_name", "m49_region", "m49_sub_region",
                 "m49_intermediate_region", "wb_income_fy27", "pop_15plus",
                 "age", "hhsize", "female", "INCOME_5", "inadincome", "education",
                 "maritalstatus", "employment", "waterdissatisfied", "urban_final",
                 "urban_source", "income_missing", "urban_missing",
                 "age_imp", "INCOME_5_imp", "inadincome_imp", "education_imp",
                 "maritalstatus_imp", "waterdissatisfied_imp", "urban_final_imp",
                 "log_GNI_imp", "Corruption_imp", "Mort_WASH_imp", "CVI_imp",
                 "Drinking_basic", "log_GNI", "Corruption",
                 country_vars)

df       <- df       %>% select(any_of(c(keep_common, "iwise12_imp", "split")))
new_data <- new_data %>% select(any_of(c(keep_common, "pre_2015")))

save_data(df,       "df_clean")
save_data(new_data, "new_data_clean")
save_data(wb,       "wb_clean")      # NOT imputed — C2/C3 use these original values

# Table S2 (national indicators): coverage per indicator from the same World Bank
# snapshot (written by A0): countries with data / missing, most recent year,
# year range and range of raw (un-logged) values. Copied to tables/ so all
# paper tables are in one place; source and definition columns are added by hand.
if (file.exists(file_wb_cov)) {
  file.copy(file_wb_cov, file.path(dir_tables, "A1_tableS02_national_indicators_coverage.csv"), overwrite = TRUE)
} else {
  warning("Coverage file not found (run A0 for this snapshot): ", file_wb_cov)
}

stop_timer()
