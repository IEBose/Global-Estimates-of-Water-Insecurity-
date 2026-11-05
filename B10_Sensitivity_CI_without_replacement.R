################################################################################
#  B10_Sensitivity_CI_without_replacement.R — SENSITIVITY: GROUP CIs WITH
#                                              COUNTRIES FIXED (NO RESAMPLING)
#
#  Purpose : Sensitivity check of the regional, sub-regional, income-group and
#            global 95% CIs from B8.
#              B8           each of 1,000 repetitions resamples countries WITH
#                           replacement within each group, then draws each
#                           country's prevalence on the logit scale,
#                           Normal(logit(estimate), SE_logit).
#              This script  countries are sampled WITHOUT replacement, i.e. every
#                           country appears exactly once in every repetition;
#                           only the (logit-scale) draws of each
#                           country's prevalence vary.
#            The 194 countries are the whole population rather than a sample,
#            so this version reflects only the uncertainty in each country's
#            estimate, not differences between the countries in a group.
#            Same seed, 1,000 repetitions, population 15+ weights and
#            logit-scale draws as B8. B8 and its outputs are not changed.
#
#  How to read
#    width_ratio = width of this CI / width of the B8 CI. Values well below 1
#    show how much of the B8 CI comes from resampling countries. Both versions
#    assume country errors are independent; model errors shared between similar
#    countries would make true intervals wider than this version, so it is best
#    read as a lower bound.
#
#  Depends : B8 (data/all_tiers.rds, data/group_CIs.rds)
#  Outputs : tables/B10_ci_comparison_<group>.csv, tables/B10_ci_comparison_all.csv
#            figures/B10_ci_comparison_m49_region_income.png
#            figures/B10_ci_comparison_m49_sub_region.png
#  Sensitivity analysis; its table is included in B11's comparison if run first.
################################################################################

source("00_setup.R")
start_timer("B10")

n_boot    <- 1000
all_tiers <- load_data("all_tiers")
b8_cis    <- load_data("group_CIs")

# Same input and row order as B8
ci_input <- all_tiers %>%
  select(country_code, prevalence, se, se_logit, tier, pop_15plus, m49_region, m49_sub_region, income_group) %>%
  arrange(tier, country_code)


# ══════════════════════════════════════════════════════════════════════════════
# 1. BOOTSTRAP WITH COUNTRIES FIXED
# ══════════════════════════════════════════════════════════════════════════════

boot_group_ci_fixed <- function(data, group_var) {
  data <- data %>%
    mutate(across(all_of(group_var), as.character)) %>%
    filter(!is.na(pop_15plus), !is.na(.data[[group_var]]))
  set.seed(seed)
  boot <- lapply(1:n_boot, function(b) {
    bd <- data %>%
      mutate(boot_prev = plogis(rnorm(n(), mean = logit_c(prevalence), sd = se_logit)))   # logit scale, as B8
    bind_rows(
      bd %>% group_by(group = .data[[group_var]]) %>%
        summarise(prev = weighted.mean(boot_prev, pop_15plus), .groups = "drop"),
      bd %>% summarise(group = "Global", prev = weighted.mean(boot_prev, pop_15plus))
    )
  }) %>% bind_rows()

  boot %>%
    group_by(group) %>%
    summarise(ci_lower_fixed = max(quantile(prev, 0.025), 0),
              ci_upper_fixed = min(quantile(prev, 0.975), 1), .groups = "drop") %>%
    rename(!!group_var := group)
}

compare_with_b8 <- function(b8_tab, group_var) {
  fixed <- boot_group_ci_fixed(ci_input, group_var)
  b8_tab %>%
    mutate(across(all_of(group_var), as.character)) %>%
    select(all_of(group_var), prevalence, ci_lower_b8 = ci_lower, ci_upper_b8 = ci_upper) %>%
    left_join(fixed, by = group_var) %>%
    mutate(width_b8    = ci_upper_b8 - ci_lower_b8,
           width_fixed = ci_upper_fixed - ci_lower_fixed,
           width_ratio = width_fixed / width_b8)
}

cmp <- list(
  m49_region     = compare_with_b8(b8_cis$region,     "m49_region"),
  m49_sub_region = compare_with_b8(b8_cis$sub_region, "m49_sub_region"),
  income_group   = compare_with_b8(b8_cis$income,     "income_group")
)

for (g in names(cmp)) save_table(cmp[[g]], paste0("B10_ci_comparison_", g))
cmp_all <- bind_rows(lapply(names(cmp), function(g)
  cmp[[g]] %>% rename(group = all_of(g)) %>% mutate(grouping = g, .before = 1)))
save_table(cmp_all, "B10_ci_comparison_all")

# Note: in B8 the global CI comes from the regional bootstrap for every grouping;
# here the global CI also differs slightly between groupings only through the
# random draws, so the regional one is reported for all.
global_fixed <- cmp$m49_region %>% filter(m49_region == "Global")
message("Global prevalence ", round(100 * global_fixed$prevalence, 1), "%: ",
        "B8 CI ", round(100 * global_fixed$ci_lower_b8, 1), "-", round(100 * global_fixed$ci_upper_b8, 1),
        "% | countries fixed ", round(100 * global_fixed$ci_lower_fixed, 1), "-",
        round(100 * global_fixed$ci_upper_fixed, 1), "%")
print(cmp_all %>% mutate(across(where(is.numeric), ~ round(.x, 4))), n = Inf, width = Inf)


# ══════════════════════════════════════════════════════════════════════════════
# 2. COMPARISON PLOTS
# ══════════════════════════════════════════════════════════════════════════════

ci_long <- function(tab, group_var) {
  bind_rows(
    tab %>% transmute(group = .data[[group_var]], prevalence,
                      lower = ci_lower_b8, upper = ci_upper_b8,
                      method = "With replacement (B8)"),
    tab %>% transmute(group = .data[[group_var]], prevalence,
                      lower = ci_lower_fixed, upper = ci_upper_fixed,
                      method = "Without replacement (countries fixed)")
  )
}

plot_compare <- function(d, title_text) {
  order <- d %>% filter(group != "Global") %>% distinct(group, prevalence) %>%
    arrange(prevalence) %>% pull(group)
  d %>%
    mutate(group = factor(group, levels = c("Global", order))) %>%
    ggplot(aes(x = prevalence, y = group, colour = method)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.35, linewidth = 0.6,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = c("With replacement (B8)" = "grey30",
                                   "Without replacement (countries fixed)" = "#D55E00"), name = NULL) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title_text, x = "Prevalence of moderate-to-high water insecurity (95% UI)", y = NULL) +
    theme_nature() +
    theme(legend.position = "bottom")
}

p_region <- plot_compare(ci_long(cmp$m49_region, "m49_region"), "M49 regions")
p_income <- plot_compare(ci_long(cmp$income_group %>% filter(income_group != "Global"), "income_group"),
                         "Income groups")
p_sub    <- plot_compare(ci_long(cmp$m49_sub_region, "m49_sub_region"), "M49 sub-regions")

save_fig((p_region / p_income) + plot_layout(guides = "collect", heights = c(6, 5)) &
           theme(legend.position = "bottom"),
         "B10_ci_comparison_m49_region_income.png", width = 180, height = 200, units = "mm", dpi = 300)
save_fig(p_sub, "B10_ci_comparison_m49_sub_region.png", width = 180, height = 220, units = "mm", dpi = 300)

stop_timer()
