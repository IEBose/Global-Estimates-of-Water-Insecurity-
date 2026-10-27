################################################################################
#  00_setup.R — SHARED SETTINGS FOR THE GLOBAL WATER INSECURITY ESTIMATES PIPELINE
#
#  Every script starts with  source("00_setup.R")  so that it can be run on its
#  own (as long as the files it reads from earlier scripts already exist).
#
#  This file defines, in one place:
#    1. Folder locations                    (EDIT THESE ONCE)
#    2. Input file names                    (incl. the frozen World Bank snapshot)
#    3. Packages and the random seed
#    4. Variable keep-lists for the large Stata files
#    5. Predictor names used by the models
#    6. Region helpers (model regions and UN M49 reporting regions)
#    7. Plot palettes, labels and theme
#    8. Small helpers: saving/loading intermediate files, run-time logging
#
#  Nothing in this file runs an analysis.
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# 1. FOLDERS  — EDIT THIS SECTION ONLY
# ══════════════════════════════════════════════════════════════════════════════
# root_dir    : project folder that contains INPUT/ and where OUTPUT/ will be made
# scripts_dir : folder that holds these R scripts (used by run_all.R)

root_dir    <- "path/to/project"            # <- EDIT: folder containing INPUT/ (OUTPUT/ is created here)
scripts_dir <- "path/to/scripts"            # <- EDIT: folder containing these R scripts

input_dir   <- file.path(root_dir, "INPUT")
output_dir  <- file.path(root_dir, "OUTPUT")
# Sensitivity runs (e.g. A4b) set `output_dir_override` in the global environment
# so that re-run scripts write to a separate folder; normal runs are unaffected.
if (exists("output_dir_override", envir = globalenv(), inherits = FALSE)) {
  output_dir <- get("output_dir_override", envir = globalenv())
  message("NOTE: writing to SENSITIVITY output folder: ", output_dir)
}

# Sub-folders for outputs (created automatically)
dir_data    <- file.path(output_dir, "data")     # .rds files passed between scripts
dir_models  <- file.path(output_dir, "models")   # fitted models and recipes
dir_tables  <- file.path(output_dir, "tables")   # .csv tables for the paper / checks
dir_figures <- file.path(output_dir, "figures")  # .png / .pdf figures
dir_logs    <- file.path(output_dir, "logs")     # run times, session info

for (d in c(input_dir, output_dir, dir_data, dir_models, dir_tables,
            dir_figures, dir_logs)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. INPUT FILES
# ══════════════════════════════════════════════════════════════════════════════

# ── World Bank snapshot (written by A0, read by A1) ──
# A0 saves each pull with its download date, e.g. wb_latest_194_20261015.csv,
# and never overwrites an earlier pull. To make the paper reproducible, set the
# date stamp of the pull used for the paper below. Re-running A0 later creates a
# NEW file and does not change the results unless this line is changed.
wb_snapshot_stamp <- "20261014"   # <- EDIT after running A0 (e.g. "20261015")

wb_snapshot_dir <- file.path(input_dir, "wb_snapshots")
file_wb         <- file.path(wb_snapshot_dir, paste0("wb_latest_194_", wb_snapshot_stamp, ".csv"))
file_wb_pop     <- file.path(wb_snapshot_dir, paste0("wb_population_2020_2025_194_", wb_snapshot_stamp, ".csv"))
file_wb_cov     <- file.path(wb_snapshot_dir, paste0("wb_coverage_194_", wb_snapshot_stamp, ".csv"))   # Table S2 (A0)
pop_year        <- 2025           # population year used for all weighting/totals

# ── External spreadsheets used by A0 (non-World Bank sources) ──
file_cpi  <- file.path(input_dir, "external", "CPI2025_ResultsCPI2025.csv")     # Transparency International CPI 2025
file_gain <- file.path(input_dir, "external", "gain.csv")                       # ND-GAIN country index
file_jmp  <- file.path(input_dir, "external", "m49_water_coverage.csv")         # M49 regions, WB income groups, JMP coverage
file_inc  <- file.path(input_dir, "external", "Income_classifications.xlsx")    # WB income classification history

# ── Survey data (Stata) ──
file_iwise_dta  <- file.path(input_dir, "IWISE_working_dataset_20260922.dta")      # 76 observed countries
file_nowise_dta <- file.path(input_dir, "Gallup_World_Poll_Mar26_noWISE.dta")      # 85 countries without IWISE

# ── UN map boundaries (downloaded once by C1 if not present) ──
file_un_bnda_zip <- file.path(input_dir, "maps", "BNDA_simplified.zip")
dir_un_bnda      <- file.path(input_dir, "maps", "BNDA_simplified")
url_un_bnda      <- "https://geoportal.un.org/arcgis/sharing/rest/content/items/e4ee80edac9d4e08b8303522dd4a5fc1/data"


# ══════════════════════════════════════════════════════════════════════════════
# 3. PACKAGES AND SEED
# ══════════════════════════════════════════════════════════════════════════════
# Install any missing packages ONCE before running the pipeline, e.g.
#   install.packages(c("tidyverse", "haven", ...))
# Scripts do not install packages themselves.

core_pkgs <- c("tidyverse", "haven", "readxl", "recipes", "glmnet", "xgboost",
               "pROC", "PRROC", "caret", "patchwork", "scales")

load_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Please install these packages first: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
}
load_pkgs(core_pkgs)

# dplyr versions of select/filter take priority (MASS and others mask them)
select <- dplyr::select
filter <- dplyr::filter

seed <- 241189
set.seed(seed)

# Uncertainty intervals used for the regional / sub-regional / income / global
# figures and tables (C1):
#   "B11" = country (case) bootstrap with model refits (main; needs B11)
#   "B8"  = B8 intervals (countries resampled with replacement)
interval_source <- "B11"
interval_label  <- "95% UI"        # uncertainty interval


# ══════════════════════════════════════════════════════════════════════════════
# 4. VARIABLES TO KEEP FROM THE LARGE STATA FILES
# ══════════════════════════════════════════════════════════════════════════════
# A1 reads only these columns (haven::read_dta(col_select = ...)), which is much
# faster than loading the full files. If a later script needs another variable
# from the survey data, add it here.
# "required" variables: A1 stops with a message if any are not in the file.
# "optional" variables: kept if present, ignored if not.

ind_vars_raw <- c("age", "hhsize", "female", "INCOME_5", "inadincome",
                  "education", "maritalstatus", "employment",
                  "waterdissatisfied", "degurba", "urbanscale")

keep_iwise_required <- c("country_code", "iwise12_imp", "wgtnorm_using",
                         "psu", "strata", "YEAR_WAVE",
                         ind_vars_raw)
keep_iwise_optional <- c("projectionwgt_svy_yr")   # survey projection weight (A5 descriptives)

keep_nowise_required <- c("COUNTRY_ISO3", "wgtnorm_using", "psu", "YEAR_WAVE", "non_UN",
                          ind_vars_raw)
keep_nowise_optional <- c("strata")

# Read only the columns needed from a .dta file, with a clear error if any
# required column is missing (e.g. after a variable has been renamed in Stata).
read_dta_keep <- function(path, required, optional = character(0)) {
  if (!file.exists(path)) stop("File not found: ", path)
  available <- names(haven::read_dta(path, n_max = 0))
  missing   <- setdiff(required, available)
  if (length(missing) > 0) {
    stop("These required variables are not in ", basename(path), ": ",
         paste(missing, collapse = ", "),
         "\n  Check the name in Stata or edit the keep-list in 00_setup.R.")
  }
  keep <- c(required, intersect(optional, available))
  haven::read_dta(path, col_select = dplyr::all_of(keep))
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. MODEL VARIABLES
# ══════════════════════════════════════════════════════════════════════════════

outcome_var <- "iwise12_imp"      # 1 = IWISE score >= 12 (moderate-to-high WI)
weight_var  <- "wgtnorm_using"    # GWP weight normalised within country

# Individual-level model: 25 predictors
#   11 individual + 2 missingness flags + 11 national + 1 region
ind_predictors <- c("age_imp", "hhsize", "female", "INCOME_5_imp",
                    "inadincome_imp", "education_imp", "maritalstatus_imp",
                    "employment", "waterdissatisfied_imp", "urban_final_imp",
                    "urban_source")
missing_flags  <- c("income_missing", "urban_missing")
ctry_predictors_ind <- c("log_GNI_imp", "log_water_stress", "log_pop_density",
                         "Drinking_basic", "life_exp", "fertility",
                         "Corruption_imp", "u5_mort", "Mort_WASH_imp",
                         "precip_mm", "CVI_imp")
all_predictors <- c(ind_predictors, missing_flags, ctry_predictors_ind, "region_final")

# Country-level model: 12 predictors (11 national + region)
# (World Bank file names; MICE-imputed in A4)
country_vars <- c("log_gni_pc_ppp", "log_water_stress", "log_pop_density",
                  "basic_water", "life_exp", "fertility",
                  "cpi_score", "u5_mort", "wash_mort", "precip_mm", "cvi")

# National variables imputed by MICE (A4) and complete auxiliary variables
ctry_vars_wb <- c("cpi_score", "wash_mort", "cvi", "log_gni_pc_ppp")
ctry_aux_wb  <- c("log_water_stress", "log_pop_density", "basic_water",
                  "life_exp", "fertility", "u5_mort", "precip_mm")

# World Bank name -> name used in the individual-level data
ctry_name_map <- c(
  "cpi_score"        = "Corruption_imp",
  "wash_mort"        = "Mort_WASH_imp",
  "cvi"              = "CVI_imp",
  "log_gni_pc_ppp"   = "log_GNI_imp",
  "log_water_stress" = "log_water_stress",
  "log_pop_density"  = "log_pop_density",
  "basic_water"      = "Drinking_basic",
  "life_exp"         = "life_exp",
  "fertility"        = "fertility",
  "u5_mort"          = "u5_mort",
  "precip_mm"        = "precip_mm"
)

# Factor levels (must be identical in training and prediction data)
lv_female  <- c("Male", "Female")
lv_inadinc <- c("Getting by", "Difficult")
lv_wdis    <- c("No", "Yes")
lv_edu     <- c("Elementary", "Secondary", "College")
lv_mar     <- c("Single", "Married/Partnered", "Separated/Widowed")
lv_emp     <- c("Employed", "Underemployed", "Unemployed", "Out of workforce")
lv_urban   <- c("Rural", "Towns/semi-dense", "Cities")
lv_urb_src <- c("degurba", "missing", "selfrep")


# ══════════════════════════════════════════════════════════════════════════════
# 6. REGIONS
# ══════════════════════════════════════════════════════════════════════════════
# Two different region variables are used, for two different purposes:
#
#   region_final  — 9 "model regions", used as a PREDICTOR in the models and in
#                   imputation. Built from UN M49 for ALL data (observed,
#                   prediction and national), so every tier uses the same rule.
#   m49_region / m49_sub_region — UN M49 regions, used ONLY for REPORTING
#                   (all regional tables and figures). These replace the old
#                   `region_report` grouping, which is no longer used anywhere.

region_final_levels <- c("East & SE Asia", "Eastern Africa", "Latin America",
                         "MENA", "Middle Africa", "N. America, Europe & Australia",
                         "Southern Africa", "Southern Asia", "Western Africa")

# Model region from the M49 classification (all 194 countries; all scripts).
# "MENA" here = Northern Africa + Western Asia + Central Asia.
# Pacific islands (except AUS/NZL) -> East & SE Asia.
assign_region_final_m49 <- function(iso3, m49_region, m49_sub_region,
                                    m49_intermediate_region) {
  out <- dplyr::case_when(
    m49_intermediate_region == "Eastern Africa"  ~ "Eastern Africa",
    m49_intermediate_region == "Western Africa"  ~ "Western Africa",
    m49_intermediate_region == "Middle Africa"   ~ "Middle Africa",
    m49_intermediate_region == "Southern Africa" ~ "Southern Africa",
    m49_sub_region == "Eastern Asia"             ~ "East & SE Asia",
    m49_sub_region == "South-eastern Asia"       ~ "East & SE Asia",
    m49_sub_region == "Southern Asia"            ~ "Southern Asia",
    m49_sub_region == "Central Asia"             ~ "MENA",
    m49_sub_region == "Northern Africa"          ~ "MENA",
    m49_sub_region == "Western Asia"             ~ "MENA",
    m49_sub_region == "Northern America"         ~ "N. America, Europe & Australia",
    m49_region == "Europe"                       ~ "N. America, Europe & Australia",
    iso3 %in% c("AUS", "NZL")                    ~ "N. America, Europe & Australia",
    m49_region == "Oceania"                      ~ "East & SE Asia",
    m49_region == "Americas"                     ~ "Latin America",
    TRUE                                         ~ NA_character_
  )
  factor(out, levels = region_final_levels)
}

# M49 reporting levels (order used in all tables and figures)
m49_region_levels <- c("Africa", "Americas", "Asia", "Europe", "Oceania")
m49_sub_region_levels <- c(
  "Northern Africa", "Sub-Saharan Africa",
  "Latin America and the Caribbean", "Northern America",
  "Central Asia", "Eastern Asia", "South-eastern Asia", "Southern Asia", "Western Asia",
  "Eastern Europe", "Northern Europe", "Southern Europe", "Western Europe",
  "Australia and New Zealand", "Melanesia", "Micronesia", "Polynesia"
)

income_levels <- c("Low income", "Lower middle income", "Upper middle income",
                   "High income", "Unclassified")

tier_levels <- c("Observed", "Predicted (individual)", "Predicted (country-level)")


# ══════════════════════════════════════════════════════════════════════════════
# 7. PLOT SETTINGS
# ══════════════════════════════════════════════════════════════════════════════

m49_col <- c(Africa = "#E69F00", Americas = "#0072B2", Asia = "#009E73",
             Europe = "#CC79A7", Oceania = "#56B4E9")

tier_col <- c("Observed"                  = "darkgreen",
              "Predicted (individual)"    = "seagreen3",
              "Predicted (country-level)" = "olivedrab2")

theme_nature <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line         = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      axis.text         = ggplot2::element_text(colour = "black"),
      panel.grid        = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 1,
                                                hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.tag          = ggplot2::element_text(face = "bold", size = base_size + 3),
      legend.key        = ggplot2::element_blank(),
      legend.title      = ggplot2::element_text(size = base_size - 1),
      legend.text       = ggplot2::element_text(size = base_size - 1),
      legend.background = ggplot2::element_blank()
    )
}

# Readable labels for model features (dummy names as created by recipes)
label_lookup <- c(
  "waterdissatisfied_imp_Yes"            = "Water dissatisfied",
  "inadincome_imp_Difficult"             = "Income inadequate",
  "inadincome_imp_Getting.by"            = "Getting by on income",
  "age_imp"                              = "Age",
  "hhsize"                               = "Household size",
  "female_Female"                        = "Female",
  "INCOME_5_imp_X2"                      = "Income quintile 2",
  "INCOME_5_imp_X3"                      = "Income quintile 3",
  "INCOME_5_imp_X4"                      = "Income quintile 4",
  "INCOME_5_imp_X5"                      = "Income quintile 5 (richest)",
  "education_imp_Secondary"              = "Secondary education",
  "education_imp_College"                = "College education",
  "maritalstatus_imp_Married.Partnered"  = "Married/Partnered",
  "maritalstatus_imp_Separated.Widowed"  = "Separated/Widowed",
  "employment_Underemployed"             = "Underemployed",
  "employment_Unemployed"                = "Unemployed",
  "employment_Out.of.workforce"          = "Out of workforce",
  "urban_final_imp_Towns.semi.dense"     = "Towns/semi-dense",
  "urban_final_imp_Cities"               = "Cities",
  "urban_source_missing"                 = "Urbanicity missing",
  "urban_source_selfrep"                 = "Urbanicity self-reported",
  "income_missing_X1"                    = "Income missing",
  "urban_missing_X1"                     = "Urbanicity missing flag",
  # National (individual-level model names)
  "log_GNI_imp"                          = "GNI per capita (log)",
  "log_water_stress"                     = "Water stress (log)",
  "log_pop_density"                      = "Population density (log)",
  "Drinking_basic"                       = "Basic drinking water (%)",
  "life_exp"                             = "Life expectancy",
  "fertility"                            = "Fertility rate",
  "Corruption_imp"                       = "Corruption index (low score = high corruption)",
  "u5_mort"                              = "Under-5 mortality",
  "Mort_WASH_imp"                        = "WASH mortality",
  "precip_mm"                            = "Mean annual rainfall",
  "CVI_imp"                              = "ND-GAIN index (high score = low vulnerability)",
  # National (country-level model names)
  "log_gni_pc_ppp"                       = "GNI per capita (log)",
  "basic_water"                          = "Basic drinking water (%)",
  "cpi_score"                            = "Corruption index (low score = high corruption)",
  "wash_mort"                            = "WASH mortality",
  "cvi"                                  = "ND-GAIN index (high score = low vulnerability)",
  # Model regions (reference level = East & SE Asia)
  "region_final_Eastern.Africa"                  = "Eastern Africa",
  "region_final_Latin.America"                   = "Latin America",
  "region_final_MENA"                            = "MENA",
  "region_final_Middle.Africa"                   = "Middle Africa",
  "region_final_N..America..Europe...Australia"  = "N. America/Europe/Australia",
  "region_final_Southern.Africa"                 = "Southern Africa",
  "region_final_Southern.Asia"                   = "Southern Asia",
  "region_final_Western.Africa"                  = "Western Africa"
)
relabel <- function(x) ifelse(x %in% names(label_lookup), label_lookup[x], x)


# ══════════════════════════════════════════════════════════════════════════════
# 8. HELPERS
# ══════════════════════════════════════════════════════════════════════════════

# ── Save / load intermediate objects (.rds) ──
save_data  <- function(x, name) saveRDS(x, file.path(dir_data,   paste0(name, ".rds")))
save_model <- function(x, name) saveRDS(x, file.path(dir_models, paste0(name, ".rds")))
load_data  <- function(name) {
  f <- file.path(dir_data, paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing input '", name, ".rds'. Run the script that creates it first.")
  readRDS(f)
}
load_model <- function(name) {
  f <- file.path(dir_models, paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing model '", name, ".rds'. Run the script that creates it first.")
  readRDS(f)
}
save_table <- function(x, name) readr::write_csv(x, file.path(dir_tables, paste0(name, ".csv")))
save_fig   <- function(plot, name, ...) {
  ggplot2::ggsave(file.path(dir_figures, name), plot, ...)
}

# ── Run-time logging ──
# Each script calls start_timer("A1") at the top and stop_timer() at the end.
# Times are appended to logs/run_times.csv so that the time for each script
# (and for the whole pipeline) is recorded every time it is run.
start_timer <- function(script) {
  assign(".timer_script", script, envir = .GlobalEnv)
  assign(".timer_start", Sys.time(), envir = .GlobalEnv)
  message("── ", script, " started ", format(Sys.time(), "%H:%M:%S"))
}
stop_timer <- function(optional = FALSE) {
  mins <- as.numeric(difftime(Sys.time(), .timer_start, units = "mins"))
  row  <- tibble::tibble(script = .timer_script, optional = optional,
                         started = format(.timer_start, "%Y-%m-%d %H:%M:%S"),
                         minutes = round(mins, 2))
  f <- file.path(dir_logs, "run_times.csv")
  readr::write_csv(row, f, append = file.exists(f))
  message("── ", .timer_script, " finished in ", round(mins, 1), " minutes")
}

# ── Observed country prevalence with survey-design SE ──
# Survey-weighted prevalence of moderate-to-high water insecurity for each
# observed country, accounting for PSU, strata and weights (Taylor
# linearisation). Used for the country-level model outcome (B4) and Tier 1 (B8).
# The normalised weight is used: within a single country and survey year,
# multiplying every weight by the same constant (e.g. the country's population,
# as in the Stata projection weights) does not change the prevalence or its SE.
observed_country_prevalence <- function(df) {
  options(survey.lonely.psu = "certainty")    # = Stata singleunit(certainty)
  d <- df %>% dplyr::mutate(outcome_num = as.numeric(as.character(iwise12_imp)))
  des <- survey::svydesign(ids = ~psu, strata = ~strata, weights = ~wgtnorm_using, data = d, nest = TRUE)
  survey::svyby(~outcome_num, ~country_code, des, survey::svymean, na.rm = TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(country_code = as.character(country_code), prevalence = outcome_num, se)
}

# ── Inverse logit ──
inv_logit <- function(x) exp(x) / (1 + exp(x))

# ── UN world boundaries (Equal Earth projection), shared by C1 and C2 ──
# Downloads the UN simplified boundaries once into INPUT/maps, then caches the
# projected map in OUTPUT/data/world_map_un.rds.
map_crs <- "+proj=eqearth"
get_un_world_map <- function() {
  cache <- file.path(dir_data, "world_map_un.rds")
  if (file.exists(cache)) return(readRDS(cache))
  if (!dir.exists(dir_un_bnda)) {
    dir.create(dirname(file_un_bnda_zip), recursive = TRUE, showWarnings = FALSE)
    download.file(url_un_bnda, destfile = file_un_bnda_zip, mode = "wb")
    unzip(file_un_bnda_zip, exdir = dir_un_bnda)
  }
  world <- sf::st_read(dir_un_bnda, quiet = TRUE) |>
    sf::st_make_valid() |>
    dplyr::select(iso3_join = iso3cd, name = nam_en, geometry) |>
    dplyr::filter(!is.na(iso3_join)) |>
    sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) |>
    sf::st_transform(map_crs) |>
    sf::st_make_valid()
  saveRDS(world, cache)
  world
}
