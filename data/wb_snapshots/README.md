# National-level data snapshot used for the paper (14 October 2026)

These files are the national-level inputs used for every result in the paper. They were created by `A0_External_data_merge.R` on **14 October 2026**. They contain only public, country-level data (no survey data).

| File | Contents | Used by |
|---|---|---|
| `wb_latest_194_20261014.csv` | 194 UN member states: World Bank indicators (most recent value since 2015, with its year), Corruption Perceptions Index, ND-GAIN index, WHO/UNICEF JMP drinking-water indicators, UN M49 regions, World Bank income groups, and the log-transformed predictors | A1 (all national predictors) |
| `wb_population_2020_2025_194_20261014.csv` | Total, 0–14 and 15+ population, 2020–2025 | A1 (2025 population 15+ for weights and numbers affected) |
| `wb_coverage_194_20261014.csv` | For each indicator: countries with / without data, most recent year, year range, value range | A1 (copied as Table S2) |
| `wb_pull_log_20261014.csv` | World Bank indicator codes, date last updated by the World Bank, date downloaded | record only |

## How to use

1. Copy the three `.csv` data files into `INPUT/wb_snapshots/` (see the main README, Section 2).
2. Leave `wb_snapshot_stamp <- "20261014"` in `00_setup.R` unchanged.
3. **Do not run A0** to reproduce the paper: it downloads the current World Bank data, which will have changed since October 2026, and writes a new snapshot with a new date.

## Sources

All of these data are publicly available from the sources below. If you reuse them, please cite the original sources.

| Source | Variables |
|---|---|
| World Bank, World Development Indicators (<https://data.worldbank.org>) | GNI per capita (PPP), life expectancy, under-5 mortality, water stress, WASH mortality, basic drinking water, precipitation, primary completion, fertility, food insecurity, population density, population |
| World Bank country and lending groups (<https://datahelpdesk.worldbank.org/knowledgebase/articles/906519>) | Income groups (FY20–FY27) |
| Transparency International, Corruption Perceptions Index 2025 (<https://www.transparency.org/en/cpi/2025>) | `cpi_score` |
| Notre Dame Global Adaptation Initiative (ND-GAIN) Country Index (<https://gain.nd.edu/our-work/country-index/>) | `nd_gain` |
| WHO/UNICEF Joint Monitoring Programme for Water Supply, Sanitation and Hygiene (<https://washdata.org>) | Safely managed drinking water and its components (available when needed, on premises, free from contamination) |
| United Nations Statistics Division, M49 standard (<https://unstats.un.org/unsd/methodology/m49/>) | M49 regions and sub-regions, country names |

Changes made by A0: values before 2015 set to missing; most recent available year kept for each indicator; GNI per capita, water stress and population density log-transformed; sources merged by ISO3 code.
