################################################################################
#  A0_External_data_merge.R — PULL AND SAVE A DATED WORLD BANK SNAPSHOT
#
#  Purpose : Download national indicators for the 194 countries (193 UN Member
#            States + State of Palestine) from the World Bank API, merge the
#            non-World Bank sources (CPI, ND-GAIN, M49/JMP file, income groups),
#            and save everything as a DATED snapshot.
#
#  Run     : Only when deliberately refreshing the national data. The rest of
#            the pipeline never calls the API: A1 reads the snapshot named in
#            00_setup.R (wb_snapshot_stamp), so re-running A0 later cannot change
#            the paper's results unless that stamp is changed.
#
#  Inputs  : World Bank API (internet), external files listed in 00_setup.R
#  Outputs : INPUT/wb_snapshots/
#              wb_latest_194_<date>.csv              one row per country, latest value
#              wb_coverage_194_<date>.csv            coverage per indicator (Table S2)
#              wb_population_2020_2025_194_<date>.csv population by year (pop 15+)
#              wb_pull_log_<date>.csv                API "last updated" date per indicator
#
#  After running: set wb_snapshot_stamp in 00_setup.R to <date>.
################################################################################
source("00_setup.R")
start_timer("A0")
load_pkgs(c("jsonlite", "readxl"))

stamp <- format(Sys.Date(), "%Y%m%d")
if (!dir.exists(wb_snapshot_dir)) dir.create(wb_snapshot_dir, recursive = TRUE)

# Never overwrite an existing snapshot
out_main <- file.path(wb_snapshot_dir, paste0("wb_latest_194_", stamp, ".csv"))
if (file.exists(out_main)) {
  stop("A snapshot dated ", stamp, " already exists. Delete or rename it first ",
       "if you really want to replace it.")
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. THE 194 COUNTRIES (ISO alpha-3)
# ══════════════════════════════════════════════════════════════════════════════

un194 <- c(
  "AFG", "AGO", "ALB", "AND", "ARE", "ARG", "ARM", "ATG", "AUS", "AUT",
  "AZE", "BDI", "BEL", "BEN", "BFA", "BGD", "BGR", "BHR", "BHS", "BIH",
  "BLR", "BLZ", "BOL", "BRA", "BRB", "BRN", "BTN", "BWA", "CAF", "CAN",
  "CHE", "CHL", "CHN", "CIV", "CMR", "COD", "COG", "COL", "COM", "CPV",
  "CRI", "CUB", "CYP", "CZE", "DEU", "DJI", "DMA", "DNK", "DOM", "DZA",
  "ECU", "EGY", "ERI", "ESP", "EST", "ETH", "FIN", "FJI", "FRA", "FSM",
  "GAB", "GBR", "GEO", "GHA", "GIN", "GMB", "GNB", "GNQ", "GRC", "GRD",
  "GTM", "GUY", "HND", "HRV", "HTI", "HUN", "IDN", "IND", "IRL", "IRN",
  "IRQ", "ISL", "ISR", "ITA", "JAM", "JOR", "JPN", "KAZ", "KEN", "KGZ",
  "KHM", "KIR", "KNA", "KOR", "KWT", "LAO", "LBN", "LBR", "LBY", "LCA",
  "LIE", "LKA", "LSO", "LTU", "LUX", "LVA", "MAR", "MCO", "MDA", "MDG",
  "MDV", "MEX", "MHL", "MKD", "MLI", "MLT", "MMR", "MNE", "MNG", "MOZ",
  "MRT", "MUS", "MWI", "MYS", "NAM", "NER", "NGA", "NIC", "NLD", "NOR",
  "NPL", "NRU", "NZL", "OMN", "PAK", "PAN", "PER", "PHL", "PLW", "PNG",
  "POL", "PRK", "PRT", "PRY", "PSE", "QAT", "ROU", "RUS", "RWA", "SAU",
  "SDN", "SEN", "SGP", "SLB", "SLE", "SLV", "SMR", "SOM", "SRB", "SSD",
  "STP", "SUR", "SVK", "SVN", "SWE", "SWZ", "SYC", "SYR", "TCD", "TGO",
  "THA", "TJK", "TKM", "TLS", "TON", "TTO", "TUN", "TUR", "TUV", "TZA",
  "UGA", "UKR", "URY", "USA", "UZB", "VCT", "VEN", "VNM", "VUT", "WSM",
  "YEM", "ZAF", "ZMB", "ZWE"
)
stopifnot(length(un194) == 194, !anyDuplicated(un194))


# ══════════════════════════════════════════════════════════════════════════════
# 2. WORLD BANK INDICATORS
# ══════════════════════════════════════════════════════════════════════════════
# Names = labels used in the pipeline; values = World Bank indicator codes.
# primary_complete and food_insec_modsev are pulled for Table S2 but are NOT
# used as predictors (excluded because of high missingness).

indicators <- c(
  gni_pc_ppp        = "NY.GNP.PCAP.PP.KD",  # GNI per capita, PPP (constant 2021 int. $)
  life_exp          = "SP.DYN.LE00.IN",     # Life expectancy at birth (years)
  u5_mort           = "SH.DYN.MORT",        # Under-5 mortality (per 1,000 live births)
  water_stress      = "ER.H2O.FWST.ZS",     # Freshwater withdrawal as % of available freshwater
  wash_mort         = "SH.STA.WASH.P5",     # Mortality attributed to unsafe WASH (per 100,000)
  basic_water       = "SH.H2O.BASW.ZS",     # At least basic drinking water services (% pop.)
  precip_mm         = "AG.LND.PRCP.MM",     # Average precipitation (mm per year)
  primary_complete  = "SE.PRM.CMPT.ZS",     # Primary completion rate (not used as predictor)
  fertility         = "SP.DYN.TFRT.IN",     # Fertility rate (births per woman)
  food_insec_modsev = "SN.ITK.MSFI.ZS",     # Moderate/severe food insecurity (not used as predictor)
  pop_density       = "EN.POP.DNST",        # Population density (people per sq. km)
  pop_total         = "SP.POP.TOTL",        # Population, total
  pop_0_14          = "SP.POP.0014.TO"      # Population ages 0-14
)

# Records the API "last updated" date of each indicator (saved as a log)
pull_log <- tibble::tibble(indicator = character(), last_updated = character())

# ── Fetch one indicator ──
# latest = TRUE  -> each country's most recent non-missing year (mrnev=1)
# latest = FALSE -> every year in `years` (e.g. "2020:2025")
wb_fetch <- function(indicator, latest = TRUE, years = NULL, per_page = 20000) {
  url <- paste0(
    "https://api.worldbank.org/v2/country/all/indicator/", indicator,
    "?format=json&per_page=", per_page,
    if (latest) "&mrnev=1",
    if (!latest && !is.null(years)) paste0("&date=", years)
  )
  # The API is sometimes slow: allow up to 5 minutes per request and retry
  # up to 3 times (waiting 10 s between attempts) before giving up.
  old_timeout <- options(timeout = max(300, getOption("timeout")))
  on.exit(options(old_timeout), add = TRUE)
  raw <- NULL
  for (attempt in 1:3) {
    raw <- tryCatch(jsonlite::fromJSON(url, flatten = TRUE),
                    error = function(e) {
                      message("  ", indicator, ": attempt ", attempt, " failed (",
                              conditionMessage(e), ")")
                      NULL
                    })
    if (!is.null(raw)) break
    if (attempt < 3) Sys.sleep(10)
  }
  if (is.null(raw)) stop("World Bank API did not respond for ", indicator,
                         " after 3 attempts. Try again later.")

  if (!is.null(raw[[1]]$message)) {
    stop("API error for ", indicator, ": ", raw[[1]]$message[[1]]$value)
  }
  if (is.null(raw[[2]]) || length(raw[[2]]) == 0) {
    warning("No data returned for ", indicator)
    return(NULL)
  }
  pull_log <<- dplyr::bind_rows(pull_log, tibble::tibble(
    indicator = indicator, last_updated = as.character(raw[[1]]$lastupdated)))

  raw[[2]] |>
    as_tibble() |>
    transmute(iso3    = countryiso3code,
              country = country.value,
              year    = as.integer(date),
              value   = as.numeric(value))
}

# ── Pull all indicators ──
# Values older than min_year are set to missing (cross-section) because data
# before 2015 are not considered current. The coverage table (attribute
# "coverage") records n countries with data, years and missing countries.
wb_pull <- function(indicators, panel = FALSE, years = "2015:2025", min_year = 2015) {

  long <- indicators |>
    map(\(ind) wb_fetch(ind, latest = !panel, years = years)) |>
    compact() |>
    list_rbind(names_to = "variable") |>
    filter(iso3 %in% un194, !is.na(value))

  if (!is.null(min_year)) {
    stale <- long |> filter(year < min_year)
    if (nrow(stale) > 0) {
      message("Values older than ", min_year, " set to missing: ",
              paste0(stale$iso3, " (", stale$variable, " ", stale$year, ")", collapse = ", "))
    }
    long <- long |> filter(year >= min_year)
  }

  if (panel) {
    yr_range <- as.integer(strsplit(years, ":")[[1]])
    grid <- expand_grid(iso3 = un194, year = seq(yr_range[1], yr_range[2]))
    out <- long |>
      select(iso3, year, variable, value) |>
      pivot_wider(names_from = variable, values_from = value) |>
      right_join(grid, by = c("iso3", "year")) |>
      left_join(distinct(long, iso3, country), by = "iso3") |>
      relocate(country, .after = iso3) |>
      arrange(iso3, year)
  } else {
    out <- long |>
      select(iso3, country, variable, value, year) |>
      pivot_wider(names_from = variable, values_from = c(value, year),
                  names_glue = "{.value}_{variable}") |>
      rename_with(\(x) sub("^value_", "", x)) |>
      right_join(tibble(iso3 = un194), by = "iso3") |>
      arrange(iso3)
  }

  attr(out, "coverage") <- long |>
    group_by(variable) |>
    summarise(
      n_with_data      = n_distinct(iso3),
      n_missing        = 194 - n_distinct(iso3),
      most_recent_year = max(year),
      year_range       = paste0(min(year), "-", max(year)),
      value_min        = min(value),
      value_max        = max(value),
      missing          = paste(setdiff(un194, iso3), collapse = ", "),
      .groups = "drop"
    ) |>
    mutate(indicator = indicators[variable],
           missing   = if_else(missing == "", "none", missing)) |>
    relocate(indicator, .after = variable) |>
    arrange(match(variable, names(indicators)))

  out
}

wide     <- wb_pull(indicators, panel = FALSE, min_year = 2015)
coverage <- attr(wide, "coverage")


# ══════════════════════════════════════════════════════════════════════════════
# 3. NON-WORLD BANK SOURCES
# ══════════════════════════════════════════════════════════════════════════════

# ── Corruption Perceptions Index 2025 (3 title rows; 0-100, higher = cleaner) ──
cpi <- read.csv(file_cpi, skip = 3, check.names = FALSE)
cpi <- cpi[, names(cpi) != ""] |>
  as_tibble() |>
  transmute(iso3 = ISO3,
            cpi_score = as.numeric(`CPI 2025 score`),
            year_cpi_score = 2025L) |>
  filter(iso3 %in% un194, !is.na(cpi_score))

# ── ND-GAIN country index: keep most recent year from 2015 onwards ──
# NOTE: ND-GAIN combines vulnerability AND readiness; higher = LESS vulnerable.
# In the pipeline it is called `cvi` ("Climate Vulnerability Index").
gain <- read.csv(file_gain, check.names = FALSE) |>
  as_tibble() |>
  pivot_longer(-c(ISO3, Name), names_to = "year", values_to = "nd_gain") |>
  mutate(year = as.integer(year)) |>
  filter(ISO3 %in% un194, !is.na(nd_gain), year >= 2015) |>
  group_by(ISO3) |>
  slice_max(year, n = 1) |>
  ungroup() |>
  transmute(iso3 = ISO3, nd_gain, year_nd_gain = year)

# ── M49 regions, WB income groups and JMP drinking-water coverage ──
# This file is the source of the M49 reporting regions used throughout the paper.
# Read as character so that M49 codes keep leading zeros.
jmp <- read.csv(file_jmp, colClasses = "character") |>
  as_tibble() |>
  filter(iso3 %in% un194) |>
  transmute(
    iso3,
    un_country_name          = country_name,
    m49_region               = region_name,
    m49_sub_region           = sub_region_name,
    m49_intermediate_region  = na_if(intermediate_region_name, ""),
    wb_income_fy26           = na_if(wb_income_fy26, ""),
    wb_income_fy27           = na_if(wb_income_fy27, ""),
    water_safely_managed      = as.numeric(cov_safely_managed),
    year_water_safely_managed = as.integer(year_safely_managed),
    water_available           = as.numeric(cov_available),
    year_water_available      = as.integer(year_available),
    water_on_premises         = as.numeric(cov_on_premises),
    year_water_on_premises    = as.integer(year_on_premises),
    water_free_contam         = as.numeric(cov_free_contam),
    year_water_free_contam    = as.integer(year_free_contam)
  )

# Same 2015 cut-off for the JMP coverage values
for (v in c("water_safely_managed", "water_available",
            "water_on_premises", "water_free_contam")) {
  yr  <- paste0("year_", v)
  old <- !is.na(jmp[[yr]]) & jmp[[yr]] < 2015
  jmp[[v]][old]  <- NA
  jmp[[yr]][old] <- NA
}

# CHECK: every one of the 194 countries must have an M49 region and sub-region.
# Should print 0. If not, the country is missing from m49_water_coverage.csv.
m49_missing <- setdiff(un194, jmp$iso3[!is.na(jmp$m49_region) & !is.na(jmp$m49_sub_region)])
if (length(m49_missing) > 0) {
  warning("Countries without an M49 region/sub-region: ", paste(m49_missing, collapse = ", "))
}

wide <- wide |>
  left_join(cpi,  by = "iso3") |>
  left_join(gain, by = "iso3") |>
  left_join(jmp,  by = "iso3") |>
  relocate(un_country_name, m49_region, m49_sub_region, m49_intermediate_region,
           wb_income_fy26, wb_income_fy27, .after = country)

# ── World Bank income classification history ──
# Two header rows: fiscal year (FY22...) and the calendar year the
# classification is based on (2020...). income_classifications_2025 = FY27.
inc_raw <- readxl::read_excel(file_inc, sheet = "Country Analytical History",
                              col_names = FALSE, .name_repair = "minimal")
names(inc_raw) <- paste0("v", seq_along(inc_raw))

hdr_row  <- which(inc_raw$v2 == "Data for calendar year :")
cal_year <- suppressWarnings(as.integer(unlist(inc_raw[hdr_row, ])))
want     <- 2020:2025
year_col <- match(want, cal_year)
stopifnot(!anyNA(year_col))          # stops if the file layout changes

inc <- inc_raw[(hdr_row + 1):nrow(inc_raw), c(1, year_col)]
names(inc) <- c("iso3", paste0("income_classifications_", want))

inc_labels <- c(L = "Low income", LM = "Lower middle income",
                UM = "Upper middle income", H = "High income")

inc <- inc |>
  filter(iso3 %in% un194) |>
  mutate(across(starts_with("income_classifications_"),
                \(x) unname(inc_labels[trimws(x)])))   # ".." and blanks -> NA

wide <- wide |> left_join(inc, by = "iso3")

# CHECK: the 2025 classification should equal wb_income_fy27 (from the M49 file).
# Any rows printed here are disagreements between the two sources to resolve.
inc_check <- wide |>
  filter(income_classifications_2025 != wb_income_fy27 |
           xor(is.na(income_classifications_2025), is.na(wb_income_fy27))) |>
  select(iso3, country, income_classifications_2025, wb_income_fy27)
if (nrow(inc_check) > 0) {
  message("income_classifications_2025 disagrees with wb_income_fy27 for ",
          nrow(inc_check), " countries:")
  print(inc_check, n = Inf)
}

# ── Add the external sources to the coverage table ──
ext_cov <- function(df, var, source) {
  have <- df$iso3[!is.na(df[[var]])]
  yrs  <- df[[paste0("year_", var)]][!is.na(df[[var]])]
  vals <- df[[var]][!is.na(df[[var]])]
  tibble(variable = var, indicator = source,
         n_with_data = length(have), n_missing = 194 - length(have),
         most_recent_year = max(yrs),
         year_range = paste0(min(yrs), "-", max(yrs)),
         value_min = min(vals), value_max = max(vals),
         missing = { m <- setdiff(un194, have); if (length(m)) paste(m, collapse = ", ") else "none" })
}
coverage <- bind_rows(
  coverage,
  ext_cov(wide, "cpi_score",            "TI CPI 2025"),
  ext_cov(wide, "nd_gain",              "ND-GAIN index"),
  ext_cov(wide, "water_safely_managed", "M49 water coverage file"),
  ext_cov(wide, "water_available",      "M49 water coverage file"),
  ext_cov(wide, "water_on_premises",    "M49 water coverage file"),
  ext_cov(wide, "water_free_contam",    "M49 water coverage file")
)


# ══════════════════════════════════════════════════════════════════════════════
# 4. DERIVED VARIABLES
# ══════════════════════════════════════════════════════════════════════════════

# Population 15+ (latest year) = total minus 0-14, only when both are from the
# same year. NOTE: the paper uses pop 15+ for a FIXED year (pop_year in
# 00_setup.R, from the population file in Section 5); A1 uses that one.
wide <- wide |>
  mutate(pop_15plus      = if_else(year_pop_total == year_pop_0_14,
                                   pop_total - pop_0_14, NA_real_),
         year_pop_15plus = if_else(!is.na(pop_15plus), year_pop_total, NA_integer_))

# Log transforms for skewed variables (all minimums > 0, so no offset needed)
wide <- wide |>
  mutate(log_gni_pc_ppp   = log(gni_pc_ppp),
         log_water_stress = log(water_stress),
         log_pop_density  = log(pop_density))

# (Files are written at the end, once every API pull has succeeded, so a
#  failed run never leaves a half-finished snapshot behind.)


# ══════════════════════════════════════════════════════════════════════════════
# 5. POPULATION BY YEAR, 2020-2025 (separate file)
# ══════════════════════════════════════════════════════════════════════════════
# One row per country, one column per variable-year (pop_total_2020 ... pop_15plus_2025).

pop_years <- "2020:2025"

pop_by_year <- list(
  pop_total = wb_fetch("SP.POP.TOTL",    latest = FALSE, years = pop_years),
  pop_0_14  = wb_fetch("SP.POP.0014.TO", latest = FALSE, years = pop_years)
) |>
  list_rbind(names_to = "variable") |>
  filter(iso3 %in% un194) |>
  select(iso3, variable, year, value) |>
  pivot_wider(names_from = variable, values_from = value) |>
  mutate(pop_15plus = pop_total - pop_0_14) |>
  pivot_wider(names_from = year,
              values_from = c(pop_total, pop_0_14, pop_15plus),
              names_glue = "{.value}_{year}") |>
  right_join(select(jmp, iso3, country_name = un_country_name), by = "iso3") |>
  relocate(country_name, .after = iso3) |>
  arrange(iso3)

# CHECK: should be 194 countries. Countries missing pop_15plus_2025 are listed;
# these get no population weight in the regional/global estimates.
pop_col <- paste0("pop_15plus_", pop_year)
message("Population file: ", nrow(pop_by_year), " countries; missing ", pop_col, ": ",
        paste(pop_by_year$iso3[is.na(pop_by_year[[pop_col]])], collapse = ", "))

# ══════════════════════════════════════════════════════════════════════════════
# 6. SAVE THE SNAPSHOT (only reached if every API pull succeeded)
# ══════════════════════════════════════════════════════════════════════════════

write.csv(wide,     out_main, row.names = FALSE, na = "")
write.csv(coverage, file.path(wb_snapshot_dir, paste0("wb_coverage_194_", stamp, ".csv")),
          row.names = FALSE)

write.csv(pop_by_year,
          file.path(wb_snapshot_dir, paste0("wb_population_2020_2025_194_", stamp, ".csv")),
          row.names = FALSE, na = "")

# API log: when each World Bank indicator was last updated (for the methods text)
write.csv(pull_log %>% mutate(downloaded = as.character(Sys.Date())),
          file.path(wb_snapshot_dir, paste0("wb_pull_log_", stamp, ".csv")),
          row.names = FALSE)

message("Snapshot saved with stamp ", stamp,
        ". Set wb_snapshot_stamp <- \"", stamp, "\" in 00_setup.R to use it.")
stop_timer()
