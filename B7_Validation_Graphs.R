################################################################################
#  B7_Validation_Graphs.R — VALIDATION FIGURES (FIGURE 5 AND SUPPLEMENTARY)
#
#  Purpose : Summarise how well the models predict country prevalence in
#            countries they did not see:
#              Panel a  Predicted vs observed country prevalence (M49 region colours)
#              Panel b  Reliability: mean observed vs mean predicted by decile of
#                       predicted prevalence (with bootstrap 95% CI)
#              Panel c  Residuals (predicted - observed) by M49 region
#              Panel d  Source of estimates (observed / individual / country-level)
#                       by share of countries and of population 15+ (Figure 5 only)
#            Versions: individual ENET (Figure 5, main), individual XGB,
#            country ENET (Figure S12), country XGB.
#            Main figures (Figure 5, S12) show the 70 countries predicted in
#            5-fold LCO CV, with statistics for those 70. Supplementary versions
#            ("_all76") add the 6 hold-out countries (triangles), with statistics
#            for all 76; B7_calibration_statistics.csv gives every set
#            (LCO CV 70, all 76, and the 70 / 6 parts of the all-76 versions).
#            Both models use the same folds and hold-outs.
#
#  Depends : B1 (out-of-sample predictions), B4 (out-of-sample predictions), B5 (tiers)
#  Inputs  : data/country_pred_individual_oos.rds, data/country_df_oos.rds,
#            data/country_preds_enet.rds, data/wb_clean.rds
#  Outputs : figures/B7_fig5_validation_individual_enet.png/.pdf (+3-panel version; LCO CV 70)
#            figures/B7_fig5_validation_individual_enet_all76.png (supplement, all 76)
#            figures/B7_figS12_validation_country_enet_all76.png (supplement, all 76)
#            figures/B7_validation_individual_xgb.png, B7_figS12_validation_country_enet.png,
#            B7_validation_country_xgb.png
#            B7_panel_c_model_region_*.png (panel c by model region; for interest)
#            B7_figS11_residuals_by_income_group.png (Figure S11) / B7_residuals_by_income_group.csv
#            tables/B7_calibration_statistics.csv
#
#  How to read the calibration statistics (printed on panel a and saved):
#    slope and intercept from regressing OBSERVED on PREDICTED prevalence.
#    Ideal: slope = 1, intercept = 0. Slope < 1: predictions too extreme
#    (high ones too high, low ones too low); slope > 1: too compressed.
#    R² (squared correlation) and MAE/RMSE as in B1. In panel b, deciles close
#    to the dashed 1:1 line indicate good calibration across the range.
################################################################################

source("00_setup.R")
start_timer("B7")
load_pkgs("cowplot")

ind_oos            <- load_data("country_pred_individual_oos")
country_df_oos   <- load_data("country_df_oos")
country_preds_enet <- load_data("country_preds_enet")
wb                 <- load_data("wb_clean")


# ══════════════════════════════════════════════════════════════════════════════
# 1. FUNCTION: PANELS a-c FROM COUNTRY-LEVEL PREDICTED/OBSERVED DATA
# ══════════════════════════════════════════════════════════════════════════════

calibration_stats <- list()   # filled by build_panels()

build_panels <- function(dat, title_text, region_order = NULL) {
  dat <- dat %>% mutate(residual = predicted - observed,
                        m49_region = factor(m49_region, levels = m49_region_levels),
                        oos_source = factor(oos_source, levels = c("LCO CV (70)", "Hold-out (6)")))
  # Default panel-c order: regions by median residual
  if (is.null(region_order)) {
    region_order <- dat %>% group_by(m49_region) %>%
      summarise(med = median(residual), .groups = "drop") %>%
      arrange(med) %>% pull(m49_region) %>% as.character()
  }

  calc_stats <- function(x, set) {
    cal_fit <- lm(observed ~ predicted, data = x)
    tibble(version = title_text, set = set, n_countries = nrow(x),
           R2 = cor(x$observed, x$predicted)^2,
           RMSE = sqrt(mean(x$residual^2)), MAE = mean(abs(x$residual)),
           slope = unname(coef(cal_fit)[2]), intercept = unname(coef(cal_fit)[1]))
  }
  # Statistics on panel a describe exactly the countries plotted
  plotted_set <- if (all(dat$oos_source == "LCO CV (70)")) "LCO CV" else "All out-of-sample"
  set_label   <- paste0(plotted_set, " (n = ", nrow(dat), ")")
  stats <- calc_stats(dat, set_label)
  extra <- list()
  if (plotted_set != "LCO CV") {   # also record the two parts separately
    for (src in c("LCO CV (70)", "Hold-out (6)")) {
      sub <- filter(dat, oos_source == src)
      if (nrow(sub) >= 3) extra[[src]] <- calc_stats(sub, src)
    }
  }
  calibration_stats[[title_text]] <<- bind_rows(stats, bind_rows(extra))

  stat_label <- sprintf(
    "atop(R^2 == %.2f * ',' ~ RMSE == %.3f * ',' ~ MAE == %.3f, slope == %.2f * ',' ~ intercept == %.2f)",
    stats$R2, stats$RMSE, stats$MAE, stats$slope, stats$intercept)

  axis_max <- min(max(c(dat$observed, dat$predicted)) * 1.1, 1)

  # Panel a — scatter with LOESS trend
  p_a <- ggplot(dat, aes(x = predicted, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.9, se = TRUE, colour = "black",
                fill = "grey75", linewidth = 0.5, alpha = 0.35) +
    geom_point(aes(fill = m49_region, shape = oos_source), colour = "white", stroke = 0.3, size = 3.2, alpha = 0.95) +
    geom_text(data = dat %>% filter(abs(residual) > quantile(abs(residual), 0.85)),
              aes(label = country_code), size = 2.5, vjust = -1, check_overlap = TRUE) +
    scale_fill_manual(values = m49_col, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c("LCO CV (70)" = 21, "Hold-out (6)" = 24), name = NULL) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, axis_max), expand = c(0.01, 0)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, axis_max), expand = c(0.01, 0)) +
    annotate("text", x = 0.02, y = axis_max * 0.97, hjust = 0, vjust = 1,
             label = stat_label, parse = TRUE, size = 2.7, lineheight = 0.95) +
    coord_equal() +
    labs(title = title_text, subtitle = set_label,
         x = "Predicted prevalence", y = "Observed prevalence", tag = "a") +
    guides(fill = guide_legend(override.aes = list(size = 2.8, shape = 21), ncol = 1),
           shape = if (plotted_set == "LCO CV") "none" else
             guide_legend(override.aes = list(size = 2.8, fill = "grey40"), ncol = 1)) +
    theme_nature() +
    theme(legend.position = "right", legend.key.height = unit(10, "pt"), legend.margin = margin(l = 0))

  # Panel b — reliability by decile of predicted prevalence, bootstrap 95% CI
  set.seed(seed)
  rel <- dat %>%
    mutate(bin = ntile(predicted, 10)) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(predicted), mean_obs = mean(observed), n = n(),
              boot = list(replicate(1000, mean(sample(observed, n(), replace = TRUE)))),
              .groups = "drop") %>%
    mutate(lo = vapply(boot, quantile, probs = 0.025, numeric(1)),
           hi = vapply(boot, quantile, probs = 0.975, numeric(1)))

  p_b <- ggplot(rel, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4, colour = "grey30") +
    geom_line(colour = "#222222", linewidth = 0.4) +
    geom_point(shape = 21, fill = "#009E73", colour = "white", stroke = 0.4, size = 2.8) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, axis_max), expand = c(0.01, 0)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, axis_max), expand = c(0.01, 0)) +
    coord_equal() +
    labs(title = "Reliability by decile of predicted prevalence",
         x = "Mean predicted prevalence (decile)", y = "Mean observed prevalence (95% CI)", tag = "b") +
    theme_nature()

  # Panel c — residuals by M49 region; axis labels show the number of countries
  region_n <- dat %>% count(m49_region)
  region_labels <- setNames(paste0(region_n$m49_region, " (n=", region_n$n, ")"),
                            as.character(region_n$m49_region))
  p_c <- dat %>%
    mutate(m49_region = factor(m49_region, levels = region_order)) %>%
    ggplot(aes(x = m49_region, y = residual, fill = m49_region)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey20", linewidth = 0.3, alpha = 0.65) +
    geom_jitter(shape = 21, colour = "white", stroke = 0.25, width = 0.15, size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = m49_col, guide = "none") +
    scale_x_discrete(labels = region_labels) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = "Prediction residuals by region", x = NULL, y = "Predicted − observed", tag = "c") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  list(p_a = p_a, p_b = p_b, p_c = p_c)
}

three_panel <- function(panels) {
  # Panel a's long title would run into panel b's, so it becomes one title for the whole figure
  fig_title    <- panels$p_a$labels$title
  fig_subtitle <- panels$p_a$labels$subtitle
  p_a <- panels$p_a + labs(title = NULL, subtitle = NULL)
  p_b <- panels$p_b + labs(title = "Reliability by decile\nof predicted prevalence")
  (p_a | p_b) / (plot_spacer() | panels$p_c | plot_spacer()) +
    plot_layout(heights = c(1.3, 0.7), widths = c(0.08, 0.84, 0.08)) +
    plot_annotation(title = fig_title, subtitle = fig_subtitle,
                    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
                                  plot.subtitle = element_text(size = 9, hjust = 0.5))) &
    theme(plot.margin = margin(6, 10, 6, 6))
}

# Individual model: 76 countries (70 LCO CV + 6 hold-outs), from B1
country_oof <- function(pred_col) {
  ind_oos %>%
    transmute(country_code, m49_region, oos_source, observed, predicted = .data[[pred_col]], n)
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. FIGURE 5: INDIVIDUAL ELASTIC NET, WITH PANEL d
# ══════════════════════════════════════════════════════════════════════════════

# Panel d — share of countries and of population 15+ by data source and M49 region
coverage_df <- country_preds_enet %>%
  select(country_code = iso3, m49_region, tier) %>%
  left_join(wb %>% select(iso3, pop_15plus), by = c("country_code" = "iso3")) %>%
  mutate(tier = factor(tier, levels = rev(tier_levels)))

coverage_long <- bind_rows(
  coverage_df %>% count(m49_region, tier, name = "value") %>% mutate(metric = "Countries"),
  coverage_df %>% group_by(m49_region, tier) %>% summarise(value = sum(pop_15plus, na.rm = TRUE), .groups = "drop") %>%
    mutate(metric = "Population"),
  coverage_df %>% count(tier, name = "value") %>% mutate(m49_region = "Global", metric = "Countries"),
  coverage_df %>% group_by(tier) %>% summarise(value = sum(pop_15plus, na.rm = TRUE), .groups = "drop") %>%
    mutate(m49_region = "Global", metric = "Population")
) %>%
  group_by(metric, m49_region) %>%
  mutate(pct = value / sum(value)) %>%
  ungroup()
save_table(coverage_long, "B7_fig5d_coverage_by_tier")

# Regions ordered by share of population needing prediction
region_coverage_order <- coverage_long %>%
  filter(metric == "Population", tier != "Observed", m49_region != "Global") %>%
  group_by(m49_region) %>% summarise(pct_predicted = sum(pct), .groups = "drop") %>%
  arrange(desc(pct_predicted)) %>% pull(m49_region)
coverage_long <- coverage_long %>% mutate(m49_region = factor(m49_region, levels = c(region_coverage_order, "Global")))

panel_d_plot <- function(metric_name, title, ylab, tag = NULL) {
  ggplot(coverage_long %>% filter(metric == metric_name), aes(x = m49_region, y = pct, fill = tier)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = length(region_coverage_order) + 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = tier_col, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = ylab, tag = tag) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none",
          plot.title = element_text(margin = margin(t = 0, b = 4)))
}
legend_only <- cowplot::get_legend(
  ggplot(coverage_long, aes(x = m49_region, y = pct, fill = tier)) + geom_col() +
    scale_fill_manual(values = tier_col, name = NULL) + theme_void() +
    theme(legend.position = "top", legend.text = element_text(size = 9), legend.key.size = unit(8, "pt")))
p_d <- wrap_elements(legend_only) /
  (panel_d_plot("Countries", "N countries", "Proportion", "d") | panel_d_plot("Population", "Total population", NULL)) +
  plot_layout(heights = c(0.15, 1))

lco_only <- function(d) filter(d, oos_source == "LCO CV (70)")

# Figure 5 (main): LCO CV predictions only (70 countries)
panels_v1 <- build_panels(lco_only(country_oof("pred_enet")),
                          "Leave-country-out cross-validation (Elastic Net, individual)",
                          region_order = region_coverage_order)   # panel c matches panel d order

fig5 <- (panels_v1$p_a | panels_v1$p_b) / (panels_v1$p_c | p_d) +
  plot_layout(heights = c(1, 1)) & theme(plot.margin = margin(8, 10, 8, 8))
save_fig(fig5, "B7_fig5_validation_individual_enet.png", width = 300, height = 260, units = "mm", dpi = 600)
save_fig(fig5, "B7_fig5_validation_individual_enet.pdf", width = 300, height = 260, units = "mm", device = cairo_pdf)

save_fig(three_panel(panels_v1), "B7_fig5_validation_individual_enet_3panel.png",
         width = 190, height = 200, units = "mm", dpi = 600)

# Supplementary version: all 76 out-of-sample countries (70 LCO CV + 6 hold-outs as triangles)
panels_v1_all <- build_panels(country_oof("pred_enet"),
                              "Out-of-sample validation incl. hold-outs (Elastic Net, individual)",
                              region_order = region_coverage_order)
fig5_all <- (panels_v1_all$p_a | panels_v1_all$p_b) / (panels_v1_all$p_c | p_d) +
  plot_layout(heights = c(1, 1)) & theme(plot.margin = margin(8, 10, 8, 8))
save_fig(fig5_all, "B7_fig5_validation_individual_enet_all76.png", width = 300, height = 260, units = "mm", dpi = 600)


# ══════════════════════════════════════════════════════════════════════════════
# 3. SUPPLEMENTARY VERSIONS
# ══════════════════════════════════════════════════════════════════════════════

# Individual XGBoost (LCO CV, 70 countries)
save_fig(three_panel(build_panels(lco_only(country_oof("pred_xgb")),
                                  "Leave-country-out cross-validation (XGBoost, individual)")),
         "B7_validation_individual_xgb.png", width = 190, height = 200, units = "mm", dpi = 600)

# Country-level models (5-fold LCO CV + hold-outs)
country_oos <- function(pred_col) {
  country_df_oos %>%
    transmute(country_code = Country_Code, m49_region, oos_source, observed = Mod_high_WI, predicted = .data[[pred_col]]) %>%
    filter(!is.na(predicted))
}
# Figure S12: country-level ENET, LCO CV (70) as Figure 5; all-76 version alongside
save_fig(three_panel(build_panels(lco_only(country_oos("oos_enet")),
                                  "Leave-country-out cross-validation (Elastic Net, country-level)")),
         "B7_figS12_validation_country_enet.png", width = 190, height = 200, units = "mm", dpi = 600)
save_fig(three_panel(build_panels(country_oos("oos_enet"),
                                  "Out-of-sample validation incl. hold-outs (Elastic Net, country-level)")),
         "B7_figS12_validation_country_enet_all76.png", width = 190, height = 200, units = "mm", dpi = 600)
save_fig(three_panel(build_panels(lco_only(country_oos("oos_xgb")),
                                  "Leave-country-out cross-validation (XGBoost, country-level)")),
         "B7_validation_country_xgb.png", width = 190, height = 200, units = "mm", dpi = 600)

# ══════════════════════════════════════════════════════════════════════════════
# 4. PANEL c BY MODEL REGION (region_final) — for interest only
# ══════════════════════════════════════════════════════════════════════════════
# Same as panel c (predicted − observed, 76 countries) but grouped by the nine
# regions used as a predictor in the models, for the individual and the
# country-level Elastic Net. Regions are ordered by the individual model's
# median residual (same order in both panels); axis labels give the number of
# countries.

model_region_col <- setNames(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                               "#D55E00", "#CC79A7", "#999999", "#8C564B"), region_final_levels)

panel_c_model_region <- function(dat, title_text, region_order = NULL) {
  dat <- dat %>%
    mutate(residual = predicted - observed, region_final = as.character(region_final),
           oos_source = factor(oos_source, levels = c("LCO CV (70)", "Hold-out (6)")))
  if (is.null(region_order)) {
    region_order <- dat %>% group_by(region_final) %>%
      summarise(med = median(residual), .groups = "drop") %>% arrange(med) %>% pull(region_final)
  }
  region_order <- c(region_order, setdiff(unique(dat$region_final), region_order))
  region_n <- dat %>% count(region_final)
  region_labels <- setNames(paste0(region_n$region_final, " (n=", region_n$n, ")"), region_n$region_final)

  dat %>%
    mutate(region_final = factor(region_final, levels = region_order)) %>%
    ggplot(aes(x = region_final, y = residual, fill = region_final)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey20", linewidth = 0.3, alpha = 0.65) +
    geom_jitter(aes(shape = oos_source), colour = "grey15", stroke = 0.25,
                width = 0.15, height = 0, size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = model_region_col, guide = "none") +
    scale_shape_manual(values = c("LCO CV (70)" = 21, "Hold-out (6)" = 24), name = NULL) +
    scale_x_discrete(labels = region_labels) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title_text, x = NULL, y = "Predicted − observed") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
}

model_region_order <- ind_oos %>%
  group_by(region_final = as.character(region_final)) %>%
  summarise(med = median(pred_enet - observed), .groups = "drop") %>%
  arrange(med) %>% pull(region_final)
p_c_ind_model <- panel_c_model_region(
  ind_oos %>% transmute(region_final, oos_source, observed, predicted = pred_enet),
  "Residuals by model region (Elastic Net, individual)", model_region_order)
p_c_ctry_model <- panel_c_model_region(
  country_df_oos %>% transmute(region_final, oos_source, observed = Mod_high_WI, predicted = oos_enet),
  "Residuals by model region (Elastic Net, country-level)", model_region_order)

save_fig(p_c_ind_model,  "B7_panel_c_model_region_individual_enet.png", width = 180, height = 120, units = "mm", dpi = 600)
save_fig(p_c_ctry_model, "B7_panel_c_model_region_country_enet.png",    width = 180, height = 120, units = "mm", dpi = 600)
save_fig((p_c_ind_model / p_c_ctry_model) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
         "B7_panel_c_model_region_both.png", width = 180, height = 230, units = "mm", dpi = 600)

# ══════════════════════════════════════════════════════════════════════════════
# 5. RESIDUALS BY INCOME GROUP (supplement)
# ══════════════════════════════════════════════════════════════════════════════
# Predicted - observed prevalence by World Bank income group (FY27), LCO CV
# (70 countries), for the individual and country-level Elastic Net. Shows
# whether the models work equally well across income groups (many countries
# without survey data are high-income or small island states).

income_of <- wb %>% transmute(country_code = iso3,
                              income_group = coalesce(as.character(wb_income_fy27), "Unclassified"))
income_cols <- setNames(c("#D55E00", "#E69F00", "#56B4E9", "#0072B2", "#999999"), income_levels)

income_resid <- bind_rows(
  lco_only(country_oof("pred_enet")) %>% mutate(model = "Individual-level Elastic Net"),
  lco_only(country_oos("oos_enet"))  %>% mutate(model = "Country-level Elastic Net")
) %>%
  mutate(country_code = as.character(country_code), residual = predicted - observed) %>%
  left_join(income_of, by = "country_code") %>%
  mutate(model = factor(model, levels = c("Individual-level Elastic Net", "Country-level Elastic Net")))

income_summary <- income_resid %>%
  group_by(model, income_group) %>%
  summarise(n_countries = n(), median_residual = median(residual), mean_bias = mean(residual),
            MAE = mean(abs(residual)), .groups = "drop")
save_table(income_summary, "B7_residuals_by_income_group")

income_n <- income_resid %>% count(model, income_group)
p_income <- income_resid %>%
  mutate(income_group = factor(income_group, levels = intersect(income_levels, unique(income_group)))) %>%
  ggplot(aes(x = income_group, y = residual, fill = income_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "grey20", linewidth = 0.3, alpha = 0.65) +
  geom_jitter(shape = 21, colour = "white", stroke = 0.25, width = 0.15, height = 0, size = 1.8, alpha = 0.9) +
  geom_text(data = income_n, aes(x = income_group, y = -Inf, label = paste0("n=", n)),
            inherit.aes = FALSE, vjust = -0.5, size = 2.7, colour = "grey30") +
  facet_wrap(~ model, ncol = 1) +
  scale_fill_manual(values = income_cols, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Prediction residuals by income group (LCO CV, 70 countries)",
       x = NULL, y = "Predicted - observed") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(p_income, "B7_figS11_residuals_by_income_group.png", width = 160, height = 180, units = "mm", dpi = 600)

save_table(bind_rows(calibration_stats), "B7_calibration_statistics")

stop_timer()
