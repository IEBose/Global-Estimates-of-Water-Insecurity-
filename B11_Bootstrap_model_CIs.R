################################################################################
#  B11_Bootstrap_model_CIs.R — MODEL BOOTSTRAP CIs (COMPARISON OF METHODS)
#
#  Purpose : Confidence intervals for country, regional, sub-regional, income-
#            group and global prevalence in which the models are REFITTED in
#            every bootstrap repetition, so that model error is shared between
#            the countries predicted by the same model (as in draws-based global
#            estimates). Two ways of making the bootstrap data are compared:
#
#   "case"      Country (case) bootstrap. Training countries are resampled with
#               replacement WITHIN each model region (so every region stays in
#               the model), and PSUs are resampled within each selected country
#               (individual model). The country model's observed prevalences are
#               drawn from Normal(estimate, survey SE). Models are refitted.
#   "residual"  Residual bootstrap (the outline: fit, compute residuals, add
#               resampled residuals to the fitted values, refit, repeat).
#                 Country model: new logit prevalence = fitted value + resampled
#                   out-of-sample residual (logit scale).
#                 Individual model (0/1 outcome): new outcomes are simulated from
#                   each respondent's fitted probability after shifting each
#                   country by a resampled out-of-sample country residual (logit
#                   scale) — the binary-outcome version of the same idea.
#
#  In every repetition, for both methods:
#    Tier 1 (observed, 76)        prevalence drawn from Normal(estimate, survey SE)
#    Tier 2 (individual, 85)      refitted individual ENET -> respondent
#                                 probabilities; PSUs resampled within the country
#                                 (sampling uncertainty); + one resampled
#                                 out-of-sample country residual
#    Tier 3 (country-level, 33)   refitted country ENET + one resampled
#                                 out-of-sample country residual
#    Groups                       population-weighted (15+, 2025) means of the
#                                 194 countries in that repetition; countries are
#                                 NOT resampled at this step
#  CI = 2.5th and 97.5th percentiles across repetitions, bounded to [0, 1].
#  Point estimates are those from B8 (unchanged); the bootstrap median is also
#  saved to show any shift.
#
#  Fixed settings: ENET alpha and lambda from B1 / B4 (not re-tuned), the same
#  recipes (standardisation) as the final models, seed 241189.
#
#  Notes for the methods / limitations
#    * Out-of-sample residuals include sampling noise in the observed values and
#      the model's own estimation error, so adding them on top of the refit
#      slightly double counts these: CIs are, if anything, conservative.
#    * The residual bootstrap assumes errors of similar size in all countries;
#      the case bootstrap does not.
#
#  Depends : A4, B1, B2, B4, B5, B8 (B10 for Figure S19; run B10 before B11)
#  Outputs : data/group_UIs_B11.rds (MAIN group 95% UIs, used by C1)
#            tables/B11_ui_<group>.csv
#            tables/B11_tableS10_all_country_estimates.csv (Table S10: 194 countries
#              with 95% UIs; observed = survey CI ± 1.96 SE; Table S9 flags)
#            tables/B11_country_CIs_<method>.csv
#            tables/B11_group_CIs_comparison.csv (B8, B10 if run, case, residual)
#            figures/B11_figS19_uncertainty_sensitivity.png (Figure S19: regions and
#              income groups), B11_group_CI_comparison_sub_region.png
#            data/B11_country_draws_<method>.rds (194 x n_boot matrix of draws)
#  Main method for the group uncertainty intervals (C1 uses group_UIs_B11.rds
#  when interval_source = "B11" in 00_setup.R). B8 and its outputs are unchanged
#  and are kept as a sensitivity analysis.
#
#  Run time: the individual model is refitted on ~80,000 respondents in every
#  repetition. Set n_boot = 50 first to time it; progress and an estimated
#  finishing time are printed every 25 repetitions.
################################################################################

source("00_setup.R")
start_timer("B11")
load_pkgs(c("glmnet", "recipes"))

n_boot      <- 1000                       # 50 for a quick timing run
run_methods <- c("case", "residual")      # either or both


# ══════════════════════════════════════════════════════════════════════════════
# 1. INPUTS
# ══════════════════════════════════════════════════════════════════════════════

all_tiers   <- load_data("all_tiers")                    # B8: tiers, point estimates, SEs, pop
train_df    <- load_data("train_df_oof")                 # B1: 70 training countries
new_data    <- load_data("new_data_pred")                # B2: 85 Tier 2 countries
ind_oos     <- load_data("country_pred_individual_oos")  # B1: 76 out-of-sample country predictions
ctry_oos    <- load_data("country_df_oos")               # B4: 76 out-of-sample country predictions
wb          <- load_data("wb_imputed")
rec_i       <- load_model("ind_recipe")
enet_i      <- load_model("ind_enet_final")
settings_i  <- load_model("ind_settings")
rec_c       <- load_model("country_recipe")
enet_c      <- load_model("country_enet_model")
settings_c  <- load_model("country_settings")

tier_codes <- function(t) as.character(all_tiers$country_code[all_tiers$tier == t])
t1 <- tier_codes("Observed"); t2 <- tier_codes("Predicted (individual)"); t3 <- tier_codes("Predicted (country-level)")
message("Tiers: ", length(t1), " / ", length(t2), " / ", length(t3))

clip01 <- function(x) pmin(pmax(x, 0), 1)
logit  <- function(p) { p <- pmin(pmax(p, 0.001), 0.999); log(p / (1 - p)) }

# ── Individual model: design matrices (baked once with the final recipe) ──
num_outcome <- function(d) d %>% mutate(iwise12_imp = as.numeric(as.character(iwise12_imp)))
X_tr  <- bake(rec_i, train_df %>% num_outcome()) %>% select(-iwise12_imp) %>% as.matrix()
y_tr  <- as.numeric(as.character(train_df$iwise12_imp))
w_tr  <- train_df[[weight_var]]
c_tr  <- as.character(train_df$country_code)
p_tr  <- paste(c_tr, train_df$psu, sep = "_")
reg_tr <- as.character(train_df$region_final)
stopifnot(identical(colnames(X_tr), settings_i$train_colnames))

X_new <- bake(rec_i, new_data %>% mutate(iwise12_imp = NA_real_)) %>% select(-iwise12_imp) %>% as.matrix()
stopifnot(identical(colnames(X_new), settings_i$train_colnames))
c_new <- as.character(new_data$country_code)
p_new <- paste(c_new, new_data$psu, sep = "_")
w_new <- new_data[[weight_var]]

eta_tr_fitted <- predict(enet_i, X_tr, type = "link")[, 1]   # for the residual method

# ── Country model: design matrices ──
ctry_train <- ctry_oos %>%
  mutate(region_final = factor(region_final, levels = settings_c$region_levels), logit_prev = 0)
X_c   <- bake(rec_c, ctry_train) %>% select(-logit_prev) %>% as.matrix()
stopifnot(identical(colnames(X_c), settings_c$train_colnames))
obs_c   <- ctry_train$Mod_high_WI
reg_c   <- as.character(ctry_train$region_final)
se_c    <- all_tiers$se[match(ctry_train$Country_Code, all_tiers$country_code)]

all_pred <- wb %>% filter(iso3 %in% t3) %>%
  mutate(region_final = factor(region_final, levels = settings_c$region_levels), logit_prev = 0)
X_t3 <- bake(rec_c, all_pred) %>% select(-logit_prev) %>% as.matrix()
codes_t3 <- all_pred$iso3

eta_c_fitted <- predict(enet_c, X_c, type = "link")[, 1]

# ── Out-of-sample country residuals (76 each) ──
# prevalence scale: added to predictions (step 6); logit scale: residual method
res_i_prev  <- ind_oos$observed - ind_oos$pred_enet
res_i_logit <- logit(ind_oos$observed) - logit(ind_oos$pred_enet)
res_c_prev  <- ctry_oos$Mod_high_WI - ctry_oos$oos_enet
res_c_logit <- logit(ctry_oos$Mod_high_WI) - logit(ctry_oos$oos_enet)

# ── Tier 1 point estimates and survey SEs ──
t1_tab <- all_tiers %>% filter(tier == "Observed")

# ── Tier 2 PSU structure (for resampling PSUs within each country) ──
psu_country_new <- tapply(c_new, p_new, `[`, 1)            # country of each PSU
psu_ids_new     <- names(psu_country_new)
psu_index_new   <- match(p_new, psu_ids_new)
W_psu_new       <- as.numeric(rowsum(w_new, psu_index_new, reorder = TRUE))
psus_by_country <- split(seq_along(psu_ids_new), psu_country_new)[t2]

# ── Training PSU structure (case method) ──
psu_ids_tr   <- unique(p_tr)
psu_index_tr <- match(p_tr, psu_ids_tr)
psu_country_tr <- c_tr[match(psu_ids_tr, p_tr)]
psus_by_country_tr <- split(seq_along(psu_ids_tr), psu_country_tr)
countries_by_region_tr <- split(unique(c_tr), reg_tr[match(unique(c_tr), c_tr)])
countries_by_region_c  <- split(seq_along(reg_c), reg_c)


# ══════════════════════════════════════════════════════════════════════════════
# 2. ONE BOOTSTRAP REPETITION
# ══════════════════════════════════════════════════════════════════════════════

# multiplicity of each element when n items are drawn with replacement from x
draw_counts <- function(x) tabulate(sample.int(length(x), length(x), replace = TRUE), length(x))

refit_individual <- function(method) {
  if (method == "case") {
    # Countries within model regions, then PSUs within each selected country:
    # each row's weight is multiplied by how often its PSU was selected.
    psu_mult <- numeric(length(psu_ids_tr))
    for (reg in names(countries_by_region_tr)) {
      cc <- countries_by_region_tr[[reg]]
      k  <- draw_counts(cc)
      for (j in which(k > 0)) {
        psus <- psus_by_country_tr[[cc[j]]]
        for (rep_j in seq_len(k[j])) psu_mult[psus] <- psu_mult[psus] + draw_counts(psus)
      }
    }
    w_b <- w_tr * psu_mult[psu_index_tr]
    keep <- w_b > 0
    glmnet(X_tr[keep, ], y_tr[keep], family = "binomial", weights = w_b[keep],
           alpha = settings_i$best_alpha, lambda = settings_i$best_lambda)
  } else {
    # Residual (parametric) bootstrap: shift each country by a resampled
    # out-of-sample residual on the logit scale, simulate new 0/1 outcomes
    shift <- setNames(sample(res_i_logit, length(unique(c_tr)), replace = TRUE), unique(c_tr))
    y_b <- rbinom(length(y_tr), 1, plogis(eta_tr_fitted + shift[c_tr]))
    glmnet(X_tr, y_b, family = "binomial", weights = w_tr,
           alpha = settings_i$best_alpha, lambda = settings_i$best_lambda)
  }
}

refit_country <- function(method) {
  if (method == "case") {
    mult <- numeric(length(obs_c))
    for (idx in countries_by_region_c) mult[idx] <- mult[idx] + draw_counts(idx)
    y_b  <- logit(clip01(rnorm(length(obs_c), obs_c, se_c)))      # survey sampling in the outcome
    keep <- mult > 0
    glmnet(X_c[keep, ], y_b[keep], weights = mult[keep],
           alpha = settings_c$best_alpha_c, lambda = settings_c$best_lambda_c)
  } else {
    y_b <- eta_c_fitted + sample(res_c_logit, length(obs_c), replace = TRUE)
    glmnet(X_c, y_b, alpha = settings_c$best_alpha_c, lambda = settings_c$best_lambda_c)
  }
}

one_rep <- function(method) {
  # Tier 1: survey sampling uncertainty
  prev_t1 <- clip01(rnorm(nrow(t1_tab), t1_tab$prevalence, t1_tab$se))

  # Tier 2: refitted individual model + PSU resampling + country residual
  m_i   <- refit_individual(method)
  p_hat <- predict(m_i, X_new, type = "response")[, 1]
  WP_psu <- as.numeric(rowsum(w_new * p_hat, psu_index_new, reorder = TRUE))
  prev_t2 <- vapply(psus_by_country, function(ps) {
    k <- draw_counts(ps)
    sum(k * WP_psu[ps]) / sum(k * W_psu_new[ps])
  }, numeric(1))
  prev_t2 <- clip01(prev_t2 + sample(res_i_prev, length(prev_t2), replace = TRUE))

  # Tier 3: refitted country model + country residual
  m_c <- refit_country(method)
  prev_t3 <- clip01(plogis(predict(m_c, X_t3, type = "link")[, 1]) +
                      sample(res_c_prev, length(codes_t3), replace = TRUE))

  c(setNames(prev_t1, t1_tab$country_code), setNames(prev_t2, names(psus_by_country)),
    setNames(prev_t3, codes_t3))
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. RUN
# ══════════════════════════════════════════════════════════════════════════════

country_order <- as.character(all_tiers$country_code)
draws <- list()
for (method in run_methods) {
  message("── Method: ", method, " (", n_boot, " repetitions)")
  set.seed(seed)
  t0  <- Sys.time()
  mat <- matrix(NA_real_, nrow = length(country_order), ncol = n_boot,
                dimnames = list(country_order, NULL))
  for (b in seq_len(n_boot)) {
    r <- one_rep(method)
    mat[names(r), b] <- r
    if (b %% 25 == 0 || b == n_boot) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf("   %s: %d / %d done (%.1f min; about %.0f min left)",
                      method, b, n_boot, el, el / b * (n_boot - b)))
    }
  }
  stopifnot(!anyNA(mat))
  draws[[method]] <- mat
  save_data(mat, paste0("B11_country_draws_", method))
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. COUNTRY AND GROUP CIs
# ══════════════════════════════════════════════════════════════════════════════

for (method in names(draws)) {
  mat <- draws[[method]]
  country_ci <- all_tiers %>%
    select(country_code, country_name, tier, prevalence, ci_lower_b8 = ci_lower, ci_upper_b8 = ci_upper) %>%
    mutate(boot_median = apply(mat[as.character(country_code), ], 1, median),
           ci_lower    = apply(mat[as.character(country_code), ], 1, quantile, 0.025),
           ci_upper    = apply(mat[as.character(country_code), ], 1, quantile, 0.975))
  save_table(country_ci, paste0("B11_country_CIs_", method))
}

pop <- all_tiers$pop_15plus[match(country_order, all_tiers$country_code)]
group_ci <- function(mat, group_vec, method) {
  ok <- !is.na(pop) & !is.na(group_vec)
  g  <- c(as.character(group_vec[ok]), rep("Global", sum(ok)))
  M  <- rbind(mat[ok, , drop = FALSE], mat[ok, , drop = FALSE])
  P  <- c(pop[ok], pop[ok])
  agg <- rowsum(M * P, g) / as.numeric(rowsum(P, g))
  tibble(group = rownames(agg),
         boot_median = apply(agg, 1, median),
         ci_lower = apply(agg, 1, quantile, 0.025),
         ci_upper = apply(agg, 1, quantile, 0.975),
         method = method)
}
groupings <- list(m49_region     = all_tiers$m49_region[match(country_order, all_tiers$country_code)],
                  m49_sub_region = all_tiers$m49_sub_region[match(country_order, all_tiers$country_code)],
                  income_group   = all_tiers$income_group[match(country_order, all_tiers$country_code)])

boot_groups <- bind_rows(lapply(names(draws), function(m) bind_rows(lapply(names(groupings), function(gn)
  group_ci(draws[[m]], groupings[[gn]], m) %>% mutate(grouping = gn)))))

# B8 (and B10 if it has been run) for comparison; point estimates from B8.
# Labels used in the comparison table and plots (the case bootstrap is the main UI)
lab_main <- "Country bootstrap with model refits (main, B11)"
lab_res  <- "Residual bootstrap (sensitivity, B11)"
lab_b8   <- "Countries resampled with replacement (sensitivity, B8)"
lab_b10  <- "Countries fixed, estimate uncertainty only (sensitivity, B10)"
b8 <- load_data("group_CIs")
b8_long <- bind_rows(
  b8$region     %>% transmute(grouping = "m49_region",     group = as.character(m49_region), prevalence, ci_lower, ci_upper),
  b8$sub_region %>% transmute(grouping = "m49_sub_region", group = as.character(m49_sub_region), prevalence, ci_lower, ci_upper),
  b8$income     %>% transmute(grouping = "income_group",   group = as.character(income_group), prevalence, ci_lower, ci_upper)
) %>% mutate(method = lab_b8)
b10_file <- file.path(dir_tables, "B10_ci_comparison_all.csv")
b10_long <- if (file.exists(b10_file)) {
  readr::read_csv(b10_file, show_col_types = FALSE) %>%
    transmute(grouping, group, prevalence, ci_lower = ci_lower_fixed, ci_upper = ci_upper_fixed,
              method = lab_b10)
} else NULL

method_labels <- c(case = lab_main, residual = lab_res)
comparison <- bind_rows(
  b8_long, b10_long,
  boot_groups %>%
    left_join(b8_long %>% select(grouping, group, prevalence), by = c("grouping", "group")) %>%
    mutate(method = method_labels[method])
) %>%
  mutate(width = ci_upper - ci_lower) %>%
  group_by(grouping, group) %>%
  mutate(width_ratio_vs_B8   = width / width[method == lab_b8][1],
         width_ratio_vs_main = width / width[method == lab_main][1]) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = c(lab_main, lab_res, lab_b10, lab_b8))) %>%
  arrange(grouping, group, method) %>%
  select(grouping, group, method, prevalence, boot_median, ci_lower, ci_upper, width,
         width_ratio_vs_main, width_ratio_vs_B8)
save_table(comparison, "B11_group_CIs_comparison")
print(comparison %>% filter(group == "Global", grouping == "m49_region") %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
      width = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 4b. MAIN UNCERTAINTY INTERVALS (country bootstrap) FOR FIGURES AND TABLE S8
# ══════════════════════════════════════════════════════════════════════════════
# Same structure as B8's group_CIs (point estimates, populations and counts from
# B8; ci_lower / ci_upper replaced by the country-bootstrap 95% UIs), so C1 can
# use either (interval_source in 00_setup.R). The global UI from the regional
# grouping is used for every grouping, as in B8.

if ("case" %in% names(draws)) {
  case_ui <- boot_groups %>% filter(method == "case")
  global_ui <- case_ui %>% filter(grouping == "m49_region", group == "Global")
  to_b8_shape <- function(b8_tab, grouping_name, group_var) {
    ui <- case_ui %>% filter(grouping == grouping_name) %>%
      mutate(ci_lower = ifelse(group == "Global", global_ui$ci_lower, ci_lower),
             ci_upper = ifelse(group == "Global", global_ui$ci_upper, ci_upper)) %>%
      select(group, ci_lower, ci_upper)
    b8_tab %>%
      mutate(across(all_of(group_var), as.character)) %>%
      select(-any_of(c("ci_lower", "ci_upper", "affected_lower_m", "affected_upper_m"))) %>%
      left_join(ui, by = setNames("group", group_var)) %>%
      mutate(affected_lower_m = ci_lower * total_pop / 1e6,
             affected_upper_m = ci_upper * total_pop / 1e6,
             interval = "95% UI, country bootstrap (B11)")
  }
  group_uis <- list(region     = to_b8_shape(b8$region,     "m49_region",     "m49_region"),
                    sub_region = to_b8_shape(b8$sub_region, "m49_sub_region", "m49_sub_region"),
                    income     = to_b8_shape(b8$income,     "income_group",   "income_group"))
  save_data(group_uis, "group_UIs_B11")
  save_table(group_uis$region,     "B11_ui_m49_region")
  save_table(group_uis$sub_region, "B11_ui_m49_sub_region")
  save_table(group_uis$income,     "B11_ui_income_group")

  # Table S10: all 194 countries with 95% UIs.
  #   Observed countries: survey-design 95% CI, estimate ± 1.96 SE (Taylor
  #     linearisation; exact, as in B8 / Stata) — sampling uncertainty only.
  #   Predicted countries: country bootstrap with model refits (case method).
  #   table_S9_flags: the reasons each country is flagged in Table S9 (B8).
  s10 <- readr::read_csv(file.path(dir_tables, "B8_all_country_estimates.csv"), show_col_types = FALSE)
  case_country <- readr::read_csv(file.path(dir_tables, "B11_country_CIs_case.csv"), show_col_types = FALSE)
  flags <- if (file.exists(file.path(dir_data, "country_flags.rds"))) load_data("country_flags") else
    tibble(country_code = character(), table_S9_flags = character())
  s10_ui <- s10 %>%
    left_join(case_country %>% select(country_code, ui_lower_boot = ci_lower, ui_upper_boot = ci_upper),
              by = "country_code") %>%
    mutate(ui_lower  = ifelse(tier == "Observed", ci_lower, ui_lower_boot),
           ui_upper  = ifelse(tier == "Observed", ci_upper, ui_upper_boot),
           ui_method = ifelse(tier == "Observed", "Survey design (estimate ± 1.96 SE)",
                              "Country bootstrap with model refits")) %>%
    left_join(flags, by = "country_code") %>%
    mutate(table_S9_flags = coalesce(table_S9_flags, "")) %>%
    select(-any_of(c("ci_lower", "ci_upper", "ci_method", "individual_pred_source",
                     "ui_lower_boot", "ui_upper_boot"))) %>%
    relocate(ui_lower, ui_upper, ui_method, .after = prev_best)
  save_table(s10_ui, "B11_tableS10_all_country_estimates")
  message("Saved main 95% UIs (group_UIs_B11.rds) and Table S10 with UIs")
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. COMPARISON PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# Main UI first and in black; sensitivity analyses in colour
method_cols <- setNames(c("black", "#009E73", "#D55E00", "#0072B2"), c(lab_main, lab_res, lab_b10, lab_b8))

plot_compare <- function(d, title_text) {
  order <- d %>% filter(group != "Global") %>% group_by(group) %>%
    summarise(p = mean(prevalence, na.rm = TRUE), .groups = "drop") %>% arrange(p) %>% pull(group)
  d %>%
    mutate(group = factor(group, levels = c("Global", order)),
           method = factor(method, levels = names(method_cols))) %>%
    ggplot(aes(x = prevalence, y = group, colour = method)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.5, linewidth = 0.5,
                   position = position_dodge(width = 0.75)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.75)) +
    scale_colour_manual(values = method_cols, name = NULL, drop = TRUE) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title_text, x = "Prevalence of moderate-to-high water insecurity (95% UI)", y = NULL) +
    theme_nature() + theme(legend.position = "bottom") +
    guides(colour = guide_legend(ncol = 1))
}

p_reg <- plot_compare(comparison %>% filter(grouping == "m49_region"), "M49 regions")
p_inc <- plot_compare(comparison %>% filter(grouping == "income_group", group != "Global"), "Income groups")
p_sub <- plot_compare(comparison %>% filter(grouping == "m49_sub_region"), "M49 sub-regions")

save_fig((p_reg / p_inc) + plot_layout(guides = "collect", heights = c(6, 5)) & theme(legend.position = "bottom"),
         "B11_figS19_uncertainty_sensitivity.png", width = 190, height = 240, units = "mm", dpi = 300)
save_fig(p_sub, "B11_group_CI_comparison_sub_region.png", width = 190, height = 260, units = "mm", dpi = 300)

stop_timer()
