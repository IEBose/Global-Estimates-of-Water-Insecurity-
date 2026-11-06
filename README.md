# Global estimates of water insecurity - R analysis pipeline

R code for the country, regional and global estimates of **moderate-to-high water insecurity** (IWISE score ≥ 12) among people aged 15 and over in **194 countries**: the 193 UN Member States and the State of Palestine.

> **Data.** The survey microdata (IWISE / Gallup World Poll) are not distributed with this code; see [Data availability](#8-data-availability). Do not commit `INPUT/` or `OUTPUT/` to any repository.

---

## 1. What the pipeline does

Each country's estimate comes from one of three tiers. A country gets the highest tier its data allow.

| Tier | Source | Countries | How the estimate is made | Uncertainty |
|---|---|---|---|---|
| 1 | Observed (IWISE in the Gallup World Poll) | 76 | Survey-weighted prevalence (PSU, strata, `wgtnorm_using`) | Taylor linearisation |
| 2 | Individual-level model | 85 | Elastic net (ENET) fitted on individual Gallup responses + national predictors, then predictions averaged to country level | PSU bootstrap (sampling) combined with a residual bootstrap of the individual ENET's out-of-sample errors (logit scale) |
| 3 | Country-level model | 33 | ENET on national predictors only | Residual bootstrap of the country ENET's out-of-sample errors (logit scale) |

For each tier, the 118 countries without observed data split as **85 + 33 = 118**.

The country estimates are then aggregated to **UN M49 regions and sub-regions**, **World Bank income groups** and the **world**. Aggregation uses 2025 population aged 15+ as weights. Uncertainty is reported as **95% uncertainty intervals (UIs)** from a bootstrap in which the models are refitted in every repetition (B11, 1,000 reps); B8 and B10 intervals are sensitivity analyses.

XGBoost, random forest and a mixed-effects logistic model are fitted for comparison. **ENET is the model used for the final estimates.** XGBoost predictions (B3, B6) are a sensitivity analysis. RF and the mixed model appear only in training and validation.

### Regions: model vs reporting

There are two separate region variables. Keep them apart.

- **`region_final`** is a 9-category **model predictor**. For all data (observed, prediction and national), it is built from UN M49 by `assign_region_final_m49()` in `00_setup.R`. "MENA" in this variable means Northern Africa, Western Asia and Central Asia. Pacific islands other than Australia and New Zealand go to East & SE Asia.
- **M49 region / sub-region** is used for **all reporting** (tables, figures, uncertainty intervals).

B9 tests how sensitive the estimates are to the model-region assignment: JPN, KOR and SGP are moved to N. America/Europe/Australia, and SDN and YEM to Eastern Africa.

---

## 2. Folder layout

```
<root_dir>/                         set in 00_setup.R, section 1
├── INPUT/
│   ├── IWISE_working_dataset_20260922.dta       76 observed countries
│   ├── Gallup_World_Poll_Mar26_noWISE.dta       85 countries without IWISE
│   ├── wb_snapshots/                            frozen World Bank pulls (written by A0)
│   │   ├── wb_latest_194_<YYYYMMDD>.csv
│   │   ├── wb_coverage_194_<YYYYMMDD>.csv       = Table S2 (A1 copies it to OUTPUT/tables)
│   │   ├── wb_population_2020_2025_194_<YYYYMMDD>.csv
│   │   └── wb_pull_log_<YYYYMMDD>.csv
│   ├── external/
│   │   ├── CPI2025_ResultsCPI2025.csv           Transparency International CPI
│   │   ├── gain.csv                             ND-GAIN
│   │   ├── m49_water_coverage.csv               M49 regions, income groups, JMP coverage
│   │   └── Income_classifications.xlsx          WB income classification history
│   └── maps/                                    UN BNDA boundaries (downloaded by C1 on first run);
│                                                Natural Earth disputed areas + lakes (downloaded by C5)
└── OUTPUT/                                      created automatically
    ├── data/      .rds files passed between scripts
    ├── models/    fitted models, recipes, settings
    ├── tables/    .csv tables (paper + checks)
    ├── figures/   .png / .pdf figures
    └── logs/      run_times.csv, sessionInfo.txt, sessionInfo_<script>.txt (one per script)
```

The scripts themselves (this folder) can live anywhere. Set `scripts_dir` in `00_setup.R`.

---

## 3. How to run

1. **Edit `00_setup.R`, section 1:** set `root_dir` and `scripts_dir`. Check that the input file names in section 2 match your files.
2. **Install packages.** To use the same package versions as the paper (recommended), run once in the scripts folder:
   ```r
   install.packages("renv")
   renv::restore()   # installs the versions recorded in renv.lock (R 4.5.1)
   ```
   Otherwise install the packages yourself; `load_pkgs()` stops and lists any that are missing.
3. **World Bank snapshot.** **To reproduce the paper, do not run A0.** Copy the three files in `data/wb_snapshots/` (the 14 October 2026 snapshot used for the paper) into `INPUT/wb_snapshots/`; `00_setup.R` already points to them (`wb_snapshot_stamp <- "20261014"`). Only run A0 when you deliberately want new data:
   - Run `A0_External_data_merge.R`. It prints the date stamp of the files it writes.
   - Set `wb_snapshot_stamp <- "YYYYMMDD"` in `00_setup.R` to that stamp.
   - A0 never overwrites an existing snapshot. The paper uses the **October 2026** snapshot, so leave the stamp alone to replicate it.
4. **Run the main pipeline:** `source("run_all.R")` with the working directory set to the scripts folder.
   - It runs A1 → A4 → B1–B11 (with B2b) → C1–C5. Each script runs in a fresh environment.
   - Set `run_A0 <- TRUE` to include A0, or `run_optional <- TRUE` to add A2/A3.
5. **Or run scripts one at a time.** Every script starts with `source("00_setup.R")` and reads its inputs from `OUTPUT/`. Any script can be re-run on its own once the scripts it depends on (table below) have been run at least once.

Run times are appended to `OUTPUT/logs/run_times.csv`. Package versions are written to `OUTPUT/logs/sessionInfo.txt` (end of `run_all.R`) and, by every script, to `OUTPUT/logs/sessionInfo_<script>.txt`.

---

## 4. Scripts

Approximate run times are for a standard Windows laptop (32 GB). **Main pipeline total, excluding A0: about 2 hours** (November 2026 run on the full data: A4 ≈ 1 h, B1 ≈ 30 min, B11 ≈ 20 min; everything else a few minutes). Timings may vary depending on computer and if other tasks are running. A2 (≈ 40 min) and A3 (≈ 35 min) are optional and not counted in the total.

### A - Data preparation

| Script | Purpose | Depends on | Run time | Notes |
|---|---|---|---|---|
| `00_setup.R` | Paths, packages, seed (241189), predictor lists (Table S4), region rules, factor levels, helper functions | — | seconds | Sourced by every script |
| `A0_External_data_merge.R` | Pulls 13 WB indicators through the API (values before 2015 set to NA); merges CPI, ND-GAIN, JMP/M49 and income history; writes a dated snapshot and the 2020–2025 population file | Internet; `INPUT/external/` files | 1–3 min | Run only when refreshing data. A1 reads the saved csv, so later API changes cannot alter the pipeline |
| `A1_Data_cleaning.R` | Creates the three analysis files: `df_clean` (76 observed countries), `new_data_clean` (85) and `wb_clean` (194 countries) | Frozen WB snapshot + population file; 2 .dta files | < 1 min | Reads only the variables needed (`read_dta_keep`). Population = 2025 pop 15+. Assigns `region_final`. Drops missing weights, non-UN territories and missing outcome. Defines the hold-out split. Checks 76/85/33 counts |
| `A2_Original_ML_Training.R` | Original ML training on single (simple) imputation | A1 | ~40 min | **Optional**, for interest only |
| `A3_Imputation_sensitivity.R` | Compares MICE with simple imputation and model-based income imputation (hold-out AUC) | A1, A2 | ~35 min | **Optional**. Writes Table S12 (`A3_tableS12_imputation_sensitivity.csv`) |
| `A4_MICE_Imputation.R` | MICE for country-level variables (194 countries; mean of 5 imputations) and individual-level variables (df + new_data; imputation 1) | A1 | ~1 h | Needs only A1. Also writes the spread across imputations (`A4_country_imputation_spread.csv`). National predictors not imputed by MICE (water stress, precipitation, basic water, under-5 mortality) and countries MICE cannot impute get the model-region median (see A4b) |
| `A4b_Sensitivity_MICE_all_national_vars.R` | Sensitivity: MICE imputes ALL incomplete national variables (incl. water stress, precipitation, basic water, under-5 mortality, which A4 fills with model-region medians), then re-runs B1–B6, B8 into `OUTPUT_sensitivity_MICE_all_national/` and compares country, regional and global estimates and LCO CV metrics with the main run | Main pipeline (A1, A4, B1–B8) | ~1–2 h | **Optional**, not in run_all. Main OUTPUT/ is not changed. Run in a fresh R session |
| `A4c_MICE_sensitivity_table_figure.R` | Supplement table (global, region, income: main vs A4b) and figure (country estimates, main vs A4b, with main 95% UIs for countries with national data gaps); list of countries changing by > 1 pp | A4b, B11 | < 1 min | **Optional.** Run automatically at the end of A4b, or on its own to redraw. Writes to the A4b sensitivity folder |
| `A5_Descriptive_statistics.R` | Descriptive statistics for the 76 observed countries: survey-weighted individual predictors (as Stata `svy: mean x`; overall, optionally by M49 region) and national predictors, each **before and after MICE**; number of countries per M49 region / sub-region and by tier; % missing per individual variable (Table S3) | A1, A4 | < 1 min (mostly a one-off check against the survey package) | Tables S1A, S1B, S3. Design as in Stata `svyset psu [pweight = projectionwgt_svy_yr], strata(strata) singleunit(certainty)`, with strata/PSUs made unique within country and year |

### B - Analysis

| Script | Purpose | Depends on | Run time | Notes |
|---|---|---|---|---|
| `B1_Individual_level_Training.R` | Trains and tests ENET, RF, XGB and mixed models on the 76 observed countries (individual + national predictors) | A4 | ~30 min | Tables S5, S6 (ICC with/without predictors for the S6 footnote); Figs 6, S8–S10, S15–S17. Saves out-of-sample predictions used later for UIs. RF (`num.threads = 1`) and XGBoost (`nthread = 1`, `tree_method = "hist"`) run single-threaded so results are identical on every computer |
| `B2_Individual_ENET_predictions.R` | Predicts the 85 `new_data` countries with ENET | A4, B1 | 1–2 min | Tier 2 estimates. Extrapolation checks by variable and by country |
| `B2b_Imputation_sensitivity.R` | **Sensitivity (no change to results):** Tier 2 estimates from each of the 5 MICE imputations (prediction side), and with the ENET refitted on each of the 5 imputed training datasets (training side, α and λ fixed) | A4, B1, B2 | a few min | `B2b_tier2_by_prediction_imputation.csv`, `B2b_tier2_by_training_imputation.csv`, `B2b_enet_coefficients_by_training_imputation.csv`, `B2b_imputation_sensitivity_summary.csv`. Countries with whole variables imputed are marked |
| `B3_Individual_XGB_predictions.R` | Same as B2 using XGB | B1, B2 | <1 min | Sensitivity |
| `B4_Country_level_Training.R` | Trains and tests country-level models (national predictors only) on the 76 observed countries | A1 (hold-outs), A4, B1 (country folds) | 1–3 min | Tables S7, S8; Figs S14, S18, S19. Same hold-outs and 5 folds as B1. Final model refitted on 76. XGBoost single-threaded. `B4_enet_country_coefficients.csv` lists all coefficients (incl. 0) with exp(coefficient). `B4_prediction_interval_coverage.csv` checks the logit-scale (used) and prevalence-scale intervals |
| `B5_Country_ENET_predictions_comparison.R` | Country ENET predictions for all 194; compares individual vs country models for the 85 | A4, B2, B4 | <1 min | Only the 33 Tier 3 countries enter the final estimates |
| `B6_Country_XGB_predictions.R` | Country XGB predictions for all 194 | A4, B4, B5 | <1 min | Sensitivity |
| `B7_Validation_Graphs.R` | Validation graphs for individual and country models (ENET and XGB) | B1, B4, B5 | ~1 min | Figs 5, S11–S13; calibration statistics |
| `B8_Global_Estimates_CIs.R` | Combines the 3 tiers; country CIs; region, sub-region, income and global estimates with CIs | A1, B1–B6 | 5–15 min | Table S9 (flagged countries; out-of-range values flagged at ≥10% of the training range, `min_pct_outside`; entirely imputed individual variables; note row); country estimates (input to Table S10); B8 intervals (sensitivity, logit scale); tier contributions |
| `B9_Sensitivity_region_estimates.R` | Effect of model-region assignment on the estimates | B1, B2 | <1 min | Table S11 |
| `B10_Sensitivity_CI_without_replacement.R` | **Sensitivity:** group intervals with countries sampled *without* replacement (every country in every repetition; only the logit-scale draws of each country's prevalence vary), compared with B8 | B8 | <1 min | Lower bound (assumes independent country errors). Run before B11 so it appears in B11's comparison |
| `B11_Bootstrap_model_CIs.R` | **Main group 95% UIs.** Individual and country ENET refitted in every repetition (country bootstrap stratified by model region, with PSU resampling), plus survey sampling (Tier 1), PSU resampling (Tier 2) and a resampled out-of-sample country residual on the logit scale (Tiers 2–3); groups aggregated with countries fixed. Residual bootstrap run alongside as a check. Compares all methods (B8, B10, both B11) | A4, B1, B2, B4, B5, B8, B10 | ~20 min (n_boot = 1000) | Writes `group_UIs_B11.rds` (used by C1), Table S10 (all countries with UIs) and Fig S20. Set `n_boot <- 50` to time it first |

### C - Visualisation

| Script | Purpose | Depends on | Run time | Notes |
|---|---|---|---|---|
| `C1_Global_Estimates_Figures_Maps.R` | Bar charts and forest plot of the estimates; alternative maps | B8, B11 (UIs) | 2–5 min | Figs 2, S1, S2, S3. Maps saved as `C1_alt_*` (paper maps are from C5). Stops if B11's UIs are missing or older than B8's estimates |
| `C2_Bivariate_maps.R` | Bivariate maps: WI × water indicators (alternative versions, `C2_alt_*`) | A1 (non-imputed indicators), B8 | 1–2 min | Paper versions are from C5. Uses **observed (non-imputed)** water statistics only |
| `C3_Key_Indicator_correlations.R` | Correlations and GNI-adjusted partial correlations of WI with water indicators, all countries and observed only | A1 (non-imputed indicators), B8 | <1 min | Figs 4, S4, S5, S6 (panels A, B, C, C1–C3). Non-imputed indicators; loess span = 1 |
| `C4_UN_boundary_maps.R` | Alternative maps on UN boundaries (`C4_alt_*`): disputed/undetermined areas outlined with dashed lines and filled with the estimate set in `disputed_rules` | A1, B8 | 2–5 min | Same data, colours and classes as C1/C2/C5. Prints `C4_disputed_polygons_check.csv` |
| `C5_WHO_style_maps.R` | **Paper maps** (Figs 1, 3, S7) on UN boundaries, following the WHO Standard Operating Procedures for maps | A1, B8 | ~2 min | Abyei, the Serbia–Kosovo line and lakes come from Natural Earth (downloaded once to `INPUT/maps`). Saves PNG, PDF and SVG |

`run_all.R` runs everything in order.

---

## 5. Paper outputs → script → file

Numbering follows the manuscript. Figures are written to `OUTPUT/figures/` and tables to `OUTPUT/tables/`. Paper files carry the figure/table number (e.g. `B1_figS08_…`, `B4_tableS07_…`); other files are checks or alternative versions. Formatted tables are saved with `write_excel_csv`, so dashes and symbols display correctly in Excel.

### Main paper

| Output | Content | Script | File |
|---|---|---|---|
| Figure 1 | A: prevalence map (predicted countries hatched); B: data source (tier). UN boundaries, WHO map SOPs | C5 | `C5_fig1_map_prevalence_source.png` (panels: `C5_panels/C5_fig1a_…`, `C5_fig1b_…`) |
| Figure 2 | Prevalence (A1, B1) and number affected (A2, B2) by M49 region (A) and income group (B), with global and 95% UIs | C1 (data B8; UIs B11) | `C1_fig2A_region_with_CIs.png`, `C1_fig2B_income_with_CIs.png`, `C1_fig2_combined.png` |
| Figure 3 | Bivariate maps: WI prevalence × A) water stress, B) no basic drinking water, C) no safely managed drinking water | C5 | `C5_fig3_bivariate_maps.png` (single maps and keys: `C5_panels/C5_bivariate_*`, `C5_key_*`) |
| Figure 4 | Correlations of WI prevalence with A) water stress, B) no basic drinking water, C) no SMDWS, C1) water not available when needed, C2) not on premises, C3) not free from contamination (all 194 countries); each panel gives overall r and, in italics, (observed only r) | C3 | `C3_fig4_correlations.png` |
| Figure 5 | ENET validation, LCO CV (70 countries): pred vs obs by region; decile calibration; error by region; data source by region (countries, population 15+) | B7 | `B7_fig5_validation_individual_enet.png/.pdf` |
| Figure 6 | Final individual ENET coefficients (bracketed; groups labelled e.g. "Water"; legend "Increases / Decreases predicted risk") | B1 | `B1_fig6_enet_coefficients_bracketed.png` |

### Supplementary

| Output | Content | Script | File |
|---|---|---|---|
| Table S1A | Individual-level descriptive statistics, observed countries (survey-weighted; before/after MICE) | A5 | `A5_tableS01A_individual_descriptives.csv`, `A5_individual_descriptives_long.csv` |
| Table S1B | National-level descriptive statistics, observed countries (before/after MICE) | A5 | `A5_tableS01B_national_descriptives.csv` |
| Figure S1 | Forest plot of prevalence by M49 region and globally (95% UIs) | C1 (UIs B11) | `C1_figS01_forest_m49_region.png` |
| Figure S2 | Prevalence by M49 sub-region (95% UIs) | C1 (UIs B11) | `C1_figS02_prev_subregion_CIs.png` |
| Figure S3 | Number affected by M49 sub-region (95% UIs) | C1 (UIs B11) | `C1_figS03_pop_subregion_CIs.png` (S2 and S3 together: `C1_figS02_S03_subregion_with_CIs.png`) |
| Figure S4 | GNI-adjusted (partial) correlations, 6 indicators (A, B, C, C1–C3), all countries | C3 | `C3_figS04_partial_correlations_gni.png` |
| Figure S5 | Correlations, 6 indicators, observed countries only | C3 | `C3_figS05_correlations_observed_only.png` |
| Figure S6 | Partial correlations, 6 indicators, observed countries only | C3 | `C3_figS06_partial_correlations_gni_observed_only.png` |
| Figure S7 | Bivariate map: WI prevalence × water not available when needed | C5 | `C5_figS07_bivariate_availability.png` |
| Table S2 | National-level indicators: countries with/without data, most recent year, year range, range of raw values (source and definitions added by hand) | A0 (copied by A1) | `A1_tableS02_national_indicators_coverage.csv` (= `INPUT/wb_snapshots/wb_coverage_194_<date>.csv`) |
| Table S3 | Individual-level indicators: variable type and % missing (observed 76, predicted 85, all 161 countries) | A5 | `A5_tableS03_individual_missingness.csv` |
| Table S4 | Model predictors | 00_setup | predictor lists (no file) |
| Figure S8 | ROC curves, hold-out countries | B1 | `B1_figS08_roc_holdouts.png` |
| Figure S9 | Pred vs obs (76: LCO CV + hold-outs), individual ENET and XGBoost, by model region | B1 | `B1_figS09_pred_vs_obs_enet_xgb.png` |
| Figure S10 | Pred vs obs (76), all individual-level models | B1 | `B1_figS10_pred_vs_obs_all_models.png` |
| Figure S11 | LCO CV residuals by model region (individual and country-level ENET; regions in geographic order) | B7 | `B7_figS11_residuals_by_model_region.png` (table: `B7_residuals_by_model_region_lco.csv`) |
| Figure S12 | LCO CV residuals by income group (individual and country-level ENET) | B7 | `B7_figS12_residuals_by_income_group.png` |
| Table S5 | Hold-out (test) metrics, individual-level models | B1 | `B1_tableS05_holdout_metrics.csv` |
| Table S6 | LCO CV metrics, individual-level models: mean (SD) across folds and pooled; footnote ICC | B1 | `B1_tableS06_lco_cv_metrics.csv`; ICC: `B1_mixed_effects_ICC.csv` |
| Table S7 | Model comparison, LCO CV (70): individual and country-level models + model-region mean: R², MAE, RMSE, calibration slope/intercept | B4 (individual rows from B1) | `B4_tableS07_model_comparison_LCO_CV.csv` (full: `B4_model_comparison_LCO_CV_full.csv`) |
| Table S8 | Hold-out (6) metrics: MAE, RMSE, mean bias | B4 | `B4_tableS08_holdout_metrics.csv` (full: `B4_model_comparison_holdout_all76_full.csv`) |
| Figure S13 | Validation, country-level ENET (LCO CV, 70; as Figure 5) | B7 | `B7_figS13_validation_country_enet.png` |
| Figure S14 | Pred vs obs, all country-level models | B4 | `B4_figS14_country_cv_all_models.png` |
| Figure S15 | SHAP, individual XGBoost | B1 | `B1_figS15_shap_beeswarm_xgb.png` |
| Figure S16 | Individual ENET coefficients without water dissatisfaction | B1 | `B1_figS16_enet_coef_noWD.png` |
| Figure S17 | SHAP without water dissatisfaction | B1 | `B1_figS17_shap_noWD.png` |
| Figure S18 | Feature importance, country-level models | B4 | `B4_figS18_country_feature_importance.png` |
| Figure S19 | SHAP, country-level XGBoost | B4 | `B4_figS19_shap_country_xgb.png` |
| Table S9 | Countries flagged for higher uncertainty, by reason: pre-2015 data; every national predictor outside the training range; predictors filled with the model-region median; top 15 model disagreement | B8 (inputs A1, B2, B3, B5, B6) | `B8_tableS09_flagged_countries.csv` (inputs: `B8_uncertainty_flag_inputs.csv`, `B8_outside_training_range_all.csv`, `B8_predictors_regional_median_by_country.csv`) |
| Table S10 | All 194 countries: estimates by model, final estimate with 95% UI (observed: survey logit CI, as Stata `svy: proportion`; predicted: B11 bootstrap), Table S9 flags | B11 (estimates B8) | `B11_tableS10_all_country_estimates.csv` |
| Figure S20 | Uncertainty-interval sensitivity: main UIs vs residual bootstrap, countries fixed (B10) and countries resampled (B8), by region and income group | B11 (B8, B10) | `B11_figS20_uncertainty_sensitivity.png` (sub-regions: `B11_group_CI_comparison_sub_region.png`) |
| Table S11 | Region reclassification sensitivity | B9 | `B9_tableS11_region_reclassification.csv` |
| Table S12 | Imputation sensitivity: hold-out AUC by imputation approach | A3 (optional; needs A2) | `A3_tableS12_imputation_sensitivity.csv` |
| Text S1 | Comparison with collected national data | — | not produced by the pipeline |

**Not numbered (kept for reference):** Figure 5 and S13 with the 6 hold-outs (`B7_*_all76.png`); prediction-interval coverage (`B4_prediction_interval_coverage.csv`); baselines (`B4_baseline_oos_predictions.csv`); B8 and B10 interval tables; alternative maps (`C1_alt_*`, `C2_alt_*`, `C4_alt_*`); single correlation plots (`C3_wi_vs_*.png`, `C3_panels/`); A4b sensitivity (separate output folder).

Panels are also saved on their own so they can be combined differently: `figures/C1_panels/` (Figure 1 a and b), `figures/C2_panels/` (alternative bivariate maps without their keys, and each key), `figures/C3_panels/` (every correlation panel, raw and GNI-adjusted), `figures/C4_panels/` (alternative UN-boundary maps and keys), `figures/C5_panels/` (paper maps and keys). The ggplot objects are saved in `data/C1_map_plots.rds`, `C2_plots.rds`, `C3_plots.rds`, `C4_map_plots.rds` and `C5_map_plots.rds` for recombining with patchwork.

Each script also writes check tables, prefixed with the script name (e.g. `A1_model_region_assignment.csv`, `B7_calibration_statistics.csv`).

---

## 6. Methods summary

**Outcome.** IWISE score ≥ 12 (`iwise12_imp`), built in Stata.
- Missing items were imputed by single regression imputation within country × year.
- Respondents missing 4 or more items were excluded.

**Missing predictors (A4).**
- Country-level variables (CPI, WASH mortality, ND-GAIN, log GNI): MICE across 194 countries, using the mean of 5 imputations.
- Individual-level variables: MICE on the observed and prediction data combined, with region in place of country and the outcome excluded. Imputation 1 is used.
- Water stress, precipitation, basic water and under-5 mortality are MICE predictors but are not imputed by MICE; a country missing one of them cannot be imputed by MICE (12 countries, mainly Pacific islands and European microstates). These gaps, and the remaining gaps in the four imputed variables, get the model-region median. `A4b` (optional) imputes all of them with MICE as a sensitivity analysis.
- Individual-level gaps MICE cannot fill (Cuba, Djibouti, Guyana income) get the regional median or mode.

**Validation.**
- **Hold-out:** 6 countries are never used in training: IND, plus NAM, CMR, SOM, PSE and BRA drawn with seed 241189.
- **5-fold leave-countries-out CV** on the remaining 70 countries.
- The country-level models use the same hold-outs and folds. The final country model is refitted on all 76.
- **Country-level fit and calibration** (R², MAE, slope, intercept): Figures 5 and S13 show the LCO CV predictions (70 countries); supplementary versions (`_all76`) add the 6 hold-out countries. Tables S7 (LCO CV, 70) and S8 (hold-out, 6) report the model comparison; `B7_calibration_statistics.csv` and the B4 `_full` tables also give all 76.

**Uncertainty intervals (main: B11).** 95% uncertainty intervals (UIs) are the 2.5th and 97.5th percentiles of 1,000 bootstrap repetitions. All country draws are made on the logit scale and transformed back, so they stay between 0 and 1 without clipping (intervals narrow towards 0% and 100%). In each repetition:
- **Models:** the 70 training countries (individual model) and the 76 observed countries (country model) are resampled with replacement within model regions, PSUs are resampled within each selected country (individual model), observed prevalences are drawn from their survey distribution (country model), and both Elastic Net models are refitted with α and λ fixed.
- **Tier 1:** prevalence drawn from Normal(logit(estimate), survey SE on the logit scale) (Taylor linearisation; single-PSU strata as certainty units, as in Gallup guidance). Table S10 gives observed countries the survey 95% CI as in Stata `svy: proportion` (logit interval, design degrees of freedom; same estimate and SE as `svy: mean`).
- **Tier 2:** predictions from the refitted individual model; PSUs resampled within each country; plus a resampled out-of-sample country residual on the logit scale (76 countries).
- **Tier 3:** predictions from the refitted country model plus a resampled out-of-sample country residual on the logit scale (76 countries).
- **Groups:** population-weighted (15+) prevalence of the 194 countries in that repetition; countries are not resampled. Because one refitted model predicts all countries in a repetition, model error is shared across countries.
- `interval_source` in `00_setup.R` switches C1 between B11 (main; labelled "95% UI") and B8 intervals (labelled "95% CI").

**Sensitivity analyses.**
- **B8:** country intervals by tier (Taylor; PSU bootstrap combined with a residual bootstrap; residual bootstrap), then countries resampled *with* replacement within each group (wider: treats the 194 countries as a sample).
- **B10:** as B8 but countries fixed (narrower: assumes independent country errors).
- **B11 residual bootstrap:** new outcomes from fitted values plus resampled residuals instead of resampling countries.

### Reading the checks

| Metric | Meaning | Good |
|---|---|---|
| R² | Squared correlation of predicted vs observed country prevalence; how well the model ranks countries | Closer to 1 |
| MAE / RMSE | Average error in prevalence points (RMSE penalises large errors more) | Smaller |
| AUC | Individual-level discrimination (chance that a random insecure respondent scores higher than a secure one) | 0.5 = chance; >0.7 acceptable; >0.8 good |
| Brier | Mean squared error of individual probabilities | Smaller |
| Calibration slope | <1: high predictions tend to be above, and low predictions below, the observed values; >1: the reverse (high predictions below, low predictions above) | ≈ 1 |
| Calibration intercept | Systematic over- or under-prediction | ≈ 0 |
| ICC (mixed model) | Share of outcome variation between countries (latent scale). B1 reports it with and without predictors; the drop is the share of between-country variance the predictors explain (Table S6 footnote) | Lower with predictors = predictors capture country differences |

Hold-out results rest on only 6 countries, so read them alongside the LCO CV results.

---

## 7. Reproducibility

- **Seed:** `seed <- 241189` is set before every random step (`set.seed`).
- **Same results on every computer:** random forest and XGBoost run on a single thread (with several threads their results depend on the number of CPU cores). ENET, MICE and the survey estimates are single-threaded already, but can still differ very slightly between operating systems or package versions.
- **Package versions:** `renv.lock` records the R version (4.5.1) and the exact version of every package (and its dependencies) used for the paper results; `renv::restore()` installs them (see Section 3, step 2). `sessionInfo_paper_run.txt` is the session record from that run. Every script also writes `OUTPUT/logs/sessionInfo_<script>.txt`, so any rerun leaves its own record to compare against.
- **Frozen data:** the October 2026 WB snapshot is read by date stamp. A0 refuses to overwrite.
- **Logs:** `OUTPUT/logs/run_times.csv`, `sessionInfo.txt` and `sessionInfo_<script>.txt`.
- **Independence:** every script reloads its inputs from `OUTPUT/`; nothing relies on objects left in memory.

---

## 8. Data availability

**The survey data are not included in this repository.** Gallup World Poll (GWP) data, including the Individual Water Insecurity Experiences (IWISE) module, can be made available on reasonable request from the author, Indira Bose (indira.bose@lshtm.ac.uk), or from wise_scales@northwestern.edu.

**The national-level data used for the paper are included** in `data/wb_snapshots/` (World Bank indicators and population, CPI, ND-GAIN, JMP water indicators, M49 regions and income groups, as merged by A0 on 14 October 2026). Use these to reproduce the results; see `data/wb_snapshots/README.md` for the sources to cite.

All other data are publicly available and can be downloaded directly:

- **World Bank World Development Indicators** (including population): downloaded automatically by `A0_External_data_merge.R` through the World Bank API.
- **Corruption Perceptions Index 2025** (Transparency International): <https://www.transparency.org/en/cpi/2025>
- **ND-GAIN Country Index** (Notre Dame Global Adaptation Initiative): <https://gain.nd.edu/our-work/country-index/>
- **WHO/UNICEF JMP drinking water indicators, UN M49 regions and World Bank income classifications:** from their public websites; save them under the file names listed in Section 2 (`INPUT/external/`).
- **UN boundaries for the maps** (UN Geospatial, BNDA simplified): downloaded automatically on first use.
- **Natural Earth** 1:10m disputed areas and 1:50m lakes (public domain; <https://www.naturalearthdata.com>): downloaded automatically by `C5_WHO_style_maps.R`.

To run the code, place the data files in `INPUT/` as shown in Section 2 (including the three files from `data/wb_snapshots/` in `INPUT/wb_snapshots/`) and set `root_dir` and `scripts_dir` in `00_setup.R`. The `.gitignore` prevents data and outputs from being committed (apart from the public snapshot in `data/wb_snapshots/`).

---

## 9. How to cite

If you use this code, please cite the paper:

> Bose I, et al. [Paper title to be added]. [Pre-print server / journal], 2026. doi: [to be added]

The citation details will be updated when the pre-print is available. GitHub's "Cite this repository" button uses `CITATION.cff`.

The code is shared for transparency and reproducibility of the published results. For other uses, please contact the author (indira.bose@lshtm.ac.uk).

