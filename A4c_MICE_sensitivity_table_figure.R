################################################################################
#  A4c_MICE_sensitivity_table_figure.R — SUPPLEMENT TABLE AND FIGURE FOR THE
#                                        A4b SENSITIVITY ANALYSIS
#
#  OPTIONAL — run after A4b (A4b also runs it automatically at the end).
#  Quick (< 1 min): reads A4b's comparison tables and the main Table S10; it
#  does not re-run any model. Main outputs (OUTPUT/) are not changed.
#
#  Outputs (in OUTPUT_sensitivity_MICE_all_national/):
#    tables/A4b_table_group_comparison.csv     global, M49 region and income
#                                              group: main vs sensitivity (%)
#    tables/A4b_table_countries_changed.csv    countries changing by more than
#                                              label_pp, with the main 95% UI
#    figures/A4b_fig_country_comparison.png/.pdf
#              country estimates, main (x) vs sensitivity (y); dashed 1:1 line;
#              countries with national data gaps highlighted; countries changing
#              by more than label_pp labelled; grey bars = main 95% UI of the
#              countries with national data gaps (show_ui)
#  Not numbered yet (supplement numbering to be decided).
################################################################################

source("00_setup.R")
start_timer("A4c")

show_ui  <- TRUE   # grey horizontal bars: main 95% UI of the countries with national data gaps
                   # (the only ones that change; bars for all 194 would hide the points)
label_pp <- 1      # label (and list) countries whose estimate changed by more than this

sens_dir    <- file.path(root_dir, "OUTPUT_sensitivity_MICE_all_national")
sens_tables <- file.path(sens_dir, "tables")
sens_figs   <- file.path(sens_dir, "figures")
dir.create(sens_figs, showWarnings = FALSE, recursive = TRUE)

f_country <- file.path(sens_tables, "A4b_country_comparison.csv")
f_group   <- file.path(sens_tables, "A4b_group_comparison.csv")
f_s10     <- file.path(dir_tables, "B11_tableS10_all_country_estimates.csv")
for (f in c(f_country, f_group)) if (!file.exists(f)) stop("Not found (run A4b first): ", f)

country_cmp <- readr::read_csv(f_country, show_col_types = FALSE)
group_cmp   <- readr::read_csv(f_group,   show_col_types = FALSE)
s10 <- if (file.exists(f_s10)) {
  readr::read_csv(f_s10, show_col_types = FALSE) %>% select(country_code, ui_lower, ui_upper)
} else {
  message("Main Table S10 not found (run B11): figure drawn without UIs")
  show_ui <- FALSE
  tibble(country_code = character(), ui_lower = numeric(), ui_upper = numeric())
}

pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)

# ══════════════════════════════════════════════════════════════════════════════
# 1. TABLE: GLOBAL, M49 REGION, INCOME GROUP
# ══════════════════════════════════════════════════════════════════════════════

group_tab <- group_cmp %>%
  mutate(group_type = case_when(group == "Global"          ~ "Global",
                                grouping == "m49_region"   ~ "M49 region",
                                grouping == "income_group" ~ "Income group")) %>%
  distinct(group_type, group, .keep_all = TRUE) %>%
  mutate(order = case_when(group_type == "Global"       ~ 0,
                           group_type == "M49 region"   ~ 1 + match(group, m49_region_levels) / 100,
                           group_type == "Income group" ~ 2 + match(group, income_levels) / 100)) %>%
  arrange(order) %>%
  transmute(`Group type` = group_type, Group = group,
            `Main analysis` = pct(prev_main),
            `Sensitivity: MICE for all national variables` = pct(prev_sens),
            `Difference (percentage points)` = sprintf("%.2f", diff_pp))
readr::write_excel_csv(group_tab, file.path(sens_tables, "A4b_table_group_comparison.csv"))
print(group_tab, n = Inf, width = Inf)

# ══════════════════════════════════════════════════════════════════════════════
# 2. TABLE: COUNTRIES THAT CHANGED BY MORE THAN label_pp
# ══════════════════════════════════════════════════════════════════════════════

dat <- country_cmp %>%
  left_join(s10, by = "country_code") %>%
  mutate(changed = abs(diff_pp) > label_pp,
         within_ui = prev_sens >= ui_lower & prev_sens <= ui_upper,
         gap = factor(ifelse(had_national_gap, "Missing national data (gaps filled by MICE in the sensitivity analysis)",
                             "Complete national data"),
                      levels = c("Complete national data",
                                 "Missing national data (gaps filled by MICE in the sensitivity analysis)")))

changed_tab <- dat %>%
  filter(changed) %>%
  arrange(desc(abs(diff_pp))) %>%
  transmute(Country = country_name, ISO3 = country_code, Tier = tier,
            `Main analysis (95% UI)` = paste0(pct(prev_main), " (", pct(ui_lower), "–", pct(ui_upper), ")"),
            `Sensitivity` = pct(prev_sens),
            `Difference (percentage points)` = sprintf("%.1f", diff_pp),
            `Within main 95% UI` = ifelse(within_ui, "Yes", "No"))
readr::write_excel_csv(changed_tab, file.path(sens_tables, "A4b_table_countries_changed.csv"))
print(changed_tab, n = Inf, width = Inf)

message("Countries changing by > ", label_pp, " pp: ", sum(dat$changed, na.rm = TRUE),
        " (all with national data gaps: ", all(dat$had_national_gap[dat$changed %in% TRUE]), ")",
        " | max change: ", sprintf("%.1f", max(abs(dat$diff_pp), na.rm = TRUE)), " pp",
        " | outside the main 95% UI: ", sum(!dat$within_ui, na.rm = TRUE))

# ══════════════════════════════════════════════════════════════════════════════
# 3. FIGURE: COUNTRY ESTIMATES, MAIN VS SENSITIVITY
# ══════════════════════════════════════════════════════════════════════════════

gap_cols <- c("Complete national data" = "grey55",
              "Missing national data (gaps filled by MICE in the sensitivity analysis)" = "#D55E00")
lim <- c(0, max(c(dat$prev_main, dat$prev_sens, if (show_ui) dat$ui_upper[dat$had_national_gap]), na.rm = TRUE) * 1.02)

p <- ggplot(dat, aes(x = prev_main, y = prev_sens)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  { if (show_ui) geom_segment(data = dplyr::filter(dat, had_national_gap),
                              aes(x = ui_lower, xend = ui_upper, yend = prev_sens),
                              colour = "#F2C4A6", linewidth = 0.45, na.rm = TRUE) } +
  geom_point(aes(colour = gap), size = 1.8, alpha = 0.9) +
  geom_text(data = dplyr::filter(dat, changed), aes(label = country_code, colour = gap),
            size = 2.6, hjust = -0.25, vjust = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = gap_cols, name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = lim, expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = lim, expand = c(0, 0)) +
  coord_equal() +
  labs(x = "Main analysis (model-region median for national data gaps)",
       y = "Sensitivity analysis (MICE for all national variables)",
       title = "Country prevalence estimates: main vs sensitivity analysis",
       caption = paste0("Dashed line: equal estimates. Labelled: change > ", label_pp, " percentage point.",
                        if (show_ui) "\nHorizontal bars: 95% uncertainty interval of the main estimate (countries with missing national data)." else "")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(sens_figs, "A4b_fig_country_comparison.png"), p,
       width = 170, height = 190, units = "mm", dpi = 600)
ggsave(file.path(sens_figs, "A4b_fig_country_comparison.pdf"), p,
       width = 170, height = 190, units = "mm", device = cairo_pdf)

stop_timer(optional = TRUE)
