################################################################################
#  C3_Key_Indicator_correlations.R — CORRELATIONS WITH KEY WATER INDICATORS
#
#  Purpose : Correlate country water insecurity prevalence with water stress,
#            JMP drinking water indicators and GNI, before and after adjusting
#            for GNI (partial correlation):
#              Figure 4    WI vs A) water stress (log), B) no basic drinking
#                          water, C) no SMDWS, C1) water not available when
#                          needed — all 194 countries
#              Figure S4   GNI-adjusted (partial) correlations, 6 indicators
#                          (A, B, C, C1 + C2 not on premises, C3 not free from
#                          contamination), all 194 countries
#              Figure S5   raw correlations, 6 indicators, observed countries only
#              Figure S6   partial correlations, 6 indicators, observed only
#            Indicators (and log GNI) use ORIGINAL values only (wb_clean, not
#            imputed); countries without the indicator are left out.
#            Prevalence: the final estimate for ALL countries (observed = filled
#            points, predicted = open points).
#
#  Methods wording: "using the original (non-imputed) indicator values, and
#  water insecurity estimates for all countries (observed and predicted)".
#
#  Partial correlation (GNI-adjusted): WI and the indicator are each regressed
#  on log GNI; the correlation of the two sets of residuals is the partial r.
#
#  HOW TO READ
#    r (Pearson)  -1 to 1; sign = direction; |r| ~0.1 weak, ~0.3 moderate, >=0.5 strong.
#    Spearman     rank-based version; robust to outliers / non-linearity.
#    partial r    association remaining after accounting for national income. A
#                 drop from r to partial r means part of the raw association is
#                 explained by GNI.
#    p-values for partial r use n - 3 degrees of freedom.
#
#  Depends : B8 (all_tiers), A1 (wb_clean)
#  Outputs : tables/C3_wi_correlations_all.csv, figures/C3_*.png,
#            figures/C3_panels/C3_panel_*.png (every panel on its own),
#            data/C3_plots.rds (ggplot objects, for recombining panels)
#            Paper figures (section 3): C3_fig4_correlations.png,
#            C3_figS04_partial_correlations_gni.png, C3_figS05_correlations_observed_only.png,
#            C3_figS06_partial_correlations_gni_observed_only.png;
#            also C3_all_correlations.png (6 raw panels, all countries), C3_wi_vs_gni*.png,
#            tables/C3_wi_correlations_observed_only.csv
################################################################################

source("00_setup.R")
start_timer("C3")

# LOESS span: larger = smoother and less sensitive to outliers (0.75 was too
# sensitive to individual countries)
loess_span <- 1

all_tiers <- load_data("all_tiers")
wb        <- load_data("wb_clean")

dat <- all_tiers %>%
  select(country_code, prevalence, tier, m49_region) %>%
  left_join(wb %>% select(iso3, basic_water, log_water_stress, log_gni_pc_ppp, water_safely_managed,
                          water_available, water_free_contam, water_on_premises),
            by = c("country_code" = "iso3")) %>%
  mutate(wi                 = prevalence * 100,
         # % WITHOUT access, so both axes run in the same direction
         gap_basic          = 100 - basic_water,
         gap_safely_managed = 100 - water_safely_managed,
         gap_available      = 100 - water_available,
         gap_on_premises    = 100 - water_on_premises,
         gap_free_contam    = 100 - water_free_contam,
         stress             = log_water_stress,
         gni_adj            = log_gni_pc_ppp,
         source_type        = factor(ifelse(tier == "Observed", "Observed", "Predicted"),
                                     levels = c("Observed", "Predicted")),
         region_name        = factor(m49_region, levels = m49_region_levels)) %>%
  as.data.frame()


# ══════════════════════════════════════════════════════════════════════════════
# 1. CORRELATIONS
# ══════════════════════════════════════════════════════════════════════════════

compute_cor <- function(dat, gap_var, label) {
  d  <- dat[complete.cases(dat$wi, dat[[gap_var]]), ]
  ct <- cor.test(d$wi, d[[gap_var]])
  sp <- cor(d$wi, d[[gap_var]], method = "spearman")

  if (gap_var == "gni_adj") {   # no partial correlation for GNI itself
    p <- d; pc <- NA_real_; pp <- NA_real_
  } else {
    p <- dat[complete.cases(dat$wi, dat[[gap_var]], dat$gni_adj), ]
    p$res_y <- residuals(lm(wi ~ gni_adj, data = p))
    p$res_x <- residuals(lm(reformulate("gni_adj", gap_var), data = p))
    pc <- cor(p$res_x, p$res_y)
    pp <- 2 * pt(-abs(pc * sqrt((nrow(p) - 3) / (1 - pc^2))), nrow(p) - 3)
  }
  list(d = d, p = p, ct = ct, pc = pc,
       row = data.frame(predictor = label, n = nrow(d),
                        r = unname(ct$estimate), ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
                        p = ct$p.value, spearman = sp,
                        n_gni_adj = nrow(p), r_gni_adj = pc, p_gni_adj = pp))
}

run_cors <- function(d) list(
  safely   = compute_cor(d, "gap_safely_managed", "No safely managed service"),
  basic    = compute_cor(d, "gap_basic",          "No basic drinking water"),
  avail    = compute_cor(d, "gap_available",      "Not available when needed"),
  premises = compute_cor(d, "gap_on_premises",    "Not on premises"),
  contam   = compute_cor(d, "gap_free_contam",    "Not free from contamination"),
  stress   = compute_cor(d, "stress",             "Water stress (log)"),
  gni      = compute_cor(d, "gni_adj",            "GNI per capita (log)")
)
cors <- run_cors(dat)
save_table(bind_rows(lapply(cors, `[[`, "row")), "C3_wi_correlations_all")


# ══════════════════════════════════════════════════════════════════════════════
# 2. SCATTER PLOTS
# ══════════════════════════════════════════════════════════════════════════════
# Lines (explained in the legend):
#   solid dark line + grey band  linear fit with 95% CI
#   red dotted line              LOESS trend (span = loess_span)
#   grey dashed line             line of equality (x = y), % panels only
# Points: filled = observed prevalence, open = predicted prevalence.

line_labels <- c(lm    = "Linear fit (grey band: 95% CI)",
                 loess = paste0("LOESS trend (span = ", loess_span, ")"),
                 equal = "Equal percentages (WI = indicator)")

theme_scatter <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"), plot.margin = margin(5, 8, 5, 5))

make_scatter <- function(cr, xvar, xlab, title_text, partial = FALSE, diag_line = FALSE) {
  d     <- if (partial) cr$p else cr$d
  yvar  <- if (partial) "res_y" else "wi"
  xcol  <- if (partial) "res_x" else xvar
  r_val <- if (partial) cr$pc else unname(cr$ct$estimate)
  ylab  <- if (partial) "Water insecurity | GNI (residual)" else "Moderate-to-high water insecurity (%)"
  first_ok   <- which(complete.cases(d[[xcol]], d[[yvar]]))[1]
  legend_pts <- data.frame(x = d[[xcol]][first_ok], y = d[[yvar]][first_ok],
                           region_name = factor(m49_region_levels, levels = m49_region_levels))

  p <- ggplot(d, aes(x = .data[[xcol]], y = .data[[yvar]])) + theme_scatter
  if (diag_line) p <- p + geom_abline(intercept = 0, slope = 1, colour = "grey70", linetype = "dashed", linewidth = 0.5)
  if (partial)   p <- p + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85")

  # Line legend: the same three entries in every panel, so that panels combined
  # into one figure share a single legend
  line_keys <- unname(line_labels[c("lm", "loess", "equal")])
  line_vals <- c("solid", "dotted", "dashed")
  line_cols <- c("grey20", "red", "grey70")

  p +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20", fill = "grey40",
                alpha = 0.20, linewidth = 0.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = loess_span,
                colour = "red", linewidth = 0.8, linetype = "dotted") +
    # Invisible layer that only draws the line legend (identical in every panel)
    geom_path(data = data.frame(x = NA_real_, y = NA_real_, lt = factor(line_keys, levels = line_keys)),
              aes(x = x, y = y, linetype = lt), inherit.aes = FALSE, na.rm = TRUE) +
    scale_linetype_manual(values = setNames(line_vals, line_keys), name = NULL,
                          guide = guide_legend(order = 3, ncol = 1,
                                               override.aes = list(colour = line_cols, linewidth = 0.8))) +
    geom_point(data = ~ subset(., source_type == "Observed"), aes(fill = region_name),
               shape = 21, size = 3.2, stroke = 0.6, colour = "black", show.legend = FALSE) +
    geom_point(data = ~ subset(., source_type == "Predicted"), aes(colour = region_name),
               shape = 21, size = 3.2, stroke = 1.3, fill = "white", show.legend = FALSE) +
    # The region legends are drawn ONLY from these invisible points (one per
    # region, observed and predicted), so the legend is identical in every panel
    # even when a region has no data for an indicator. Otherwise ggplot2 (>= 3.5)
    # draws the legend differently per panel and patchwork repeats it.
    geom_point(data = legend_pts, aes(x = x, y = y, fill = region_name), inherit.aes = FALSE,
               shape = 21, size = 3.2, alpha = 0) +
    geom_point(data = legend_pts, aes(x = x, y = y, colour = region_name), inherit.aes = FALSE,
               shape = 21, size = 3.2, alpha = 0) +
    scale_fill_manual(values = m49_col, name = "Region (Observed)", limits = m49_region_levels, drop = FALSE) +
    scale_colour_manual(values = m49_col, name = "Region (Predicted)", limits = m49_region_levels, drop = FALSE) +
    guides(fill   = guide_legend(order = 1, override.aes = list(colour = "black", stroke = 0.6, alpha = 1)),
           colour = guide_legend(order = 2, override.aes = list(fill = "white", shape = 21, size = 3.2,
                                                                stroke = 1.3, alpha = 1))) +
    labs(x = xlab, y = ylab, title = title_text,
         subtitle = sprintf("%s = %.2f, n = %d", if (partial) "partial r (GNI adjusted)" else "r", r_val, nrow(d)))
}

collect_bottom <- function(p) p + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical")

# Raw correlations
p_stress   <- make_scatter(cors$stress,   "stress",             "Water stress (log)",               "Water stress")
p_basic    <- make_scatter(cors$basic,    "gap_basic",          "No basic drinking water (%)",      "No basic drinking water", diag_line = TRUE)
p_safely   <- make_scatter(cors$safely,   "gap_safely_managed", "No safely managed service (%)",    "No safely managed service", diag_line = TRUE)
p_avail    <- make_scatter(cors$avail,    "gap_available",      "Not available when needed (%)",    "Not available when needed", diag_line = TRUE)
p_premises <- make_scatter(cors$premises, "gap_on_premises",    "Not on premises (%)",              "Not on premises", diag_line = TRUE)
p_contam   <- make_scatter(cors$contam,   "gap_free_contam",    "Not free from contamination (%)",  "Not free from contamination", diag_line = TRUE)
p_gni      <- make_scatter(cors$gni,      "gni_adj",            "GNI per capita (log)",             "GNI per capita")

# GNI-adjusted (partial) correlations
pg <- function(cr, xlab, title) make_scatter(cr, NULL, xlab, title, partial = TRUE)
p_stress_g   <- pg(cors$stress,   "Water stress (log) | GNI",           "Water stress")
p_basic_g    <- pg(cors$basic,    "No basic drinking water | GNI",      "No basic drinking water")
p_safely_g   <- pg(cors$safely,   "No safely managed | GNI",            "No safely managed service")
p_avail_g    <- pg(cors$avail,    "Not available when needed | GNI",    "Not available when needed")
p_premises_g <- pg(cors$premises, "Not on premises | GNI",              "Not on premises")
p_contam_g   <- pg(cors$contam,   "Not free from contamination | GNI",  "Not free from contamination")

save_fig(p_stress,   "C3_wi_vs_water_stress.png",          width = 6.5, height = 6, dpi = 200)
save_fig(p_stress_g, "C3_wi_vs_water_stress_gni_adj.png", width = 6.5, height = 6, dpi = 200)
save_fig(collect_bottom(p_basic | p_safely),     "C3_wi_vs_basic_smdws.png",       width = 12, height = 5.5, dpi = 200)
save_fig(collect_bottom(p_basic_g | p_safely_g), "C3_wi_vs_basic_smdws_gni_adj.png",   width = 12, height = 5.5, dpi = 200)
save_fig(collect_bottom(p_avail | p_avail_g),       "C3_availability_raw_and_gni_adj.png", width = 12, height = 5.5, dpi = 200)
save_fig(collect_bottom(p_premises | p_premises_g), "C3_on_premises_raw_and_gni_adj.png",  width = 12, height = 5.5, dpi = 200)
save_fig(collect_bottom(p_contam | p_contam_g),     "C3_contamination_raw_and_gni_adj.png", width = 12, height = 5.5, dpi = 200)
save_fig(collect_bottom(p_stress_g | p_basic_g | p_safely_g), "C3_gni_adjusted_stress_basic_safely.png",
         width = 16, height = 5.5, dpi = 200)
save_fig(p_gni, "C3_wi_vs_gni.png", width = 6.5, height = 6, dpi = 200)

# ── Every panel saved on its own (to combine differently later) ──
# figures/C3_panels/: one PNG per panel; data/C3_plots.rds: the ggplot objects,
# e.g. readRDS(".../C3_plots.rds") then (p$basic | p$basic_gni) with patchwork.
panels <- list(stress = p_stress, basic = p_basic, safely = p_safely, avail = p_avail,
               premises = p_premises, contam = p_contam, gni = p_gni,
               stress_gni = p_stress_g, basic_gni = p_basic_g, safely_gni = p_safely_g,
               avail_gni = p_avail_g, premises_gni = p_premises_g, contam_gni = p_contam_g)
dir.create(file.path(dir_figures, "C3_panels"), showWarnings = FALSE)
for (nm in names(panels)) {
  save_fig(panels[[nm]], file.path("C3_panels", paste0("C3_panel_", nm, ".png")),
           width = 6.5, height = 6, dpi = 300)
}
save_data(panels, "C3_plots")


# ══════════════════════════════════════════════════════════════════════════════
# 3. PAPER FIGURES: FIGURE 4 AND FIGURES S4-S6 (COMBINED PANELS)
# ══════════════════════════════════════════════════════════════════════════════
# Panels are labelled as in the paper: A water stress, B no basic drinking water,
# C no safely managed (SMDWS), C1 not available when needed, C2 not on premises,
# C3 not free from contamination (C1-C3 are the components of SMDWS).
#   Figure 4  raw correlations, A-C1, all 194 countries
#   Figure S4 GNI-adjusted (partial), A-C3, all countries
#   Figure S5 raw, A-C3, observed countries only
#   Figure S6 GNI-adjusted, A-C3, observed countries only

observed_only <- TRUE   # FALSE to skip Figures S5/S6

panel_tags <- c("A", "B", "C", "C1", "C2", "C3")
tag_panels <- function(pl) Map(function(p, t) p + labs(tag = t) +
                                 theme(plot.tag = element_text(face = "bold", size = 14)),
                               pl, panel_tags[seq_along(pl)])

panel_sets <- function(cr, hide_predicted = FALSE) {
  raw <- list(
    make_scatter(cr$stress,   "stress",             "Water stress (log)",              "Water stress"),
    make_scatter(cr$basic,    "gap_basic",          "No basic drinking water (%)",     "No basic drinking water", diag_line = TRUE),
    make_scatter(cr$safely,   "gap_safely_managed", "No safely managed service (%)",   "No safely managed service", diag_line = TRUE),
    make_scatter(cr$avail,    "gap_available",      "Not available when needed (%)",   "Not available when needed", diag_line = TRUE),
    make_scatter(cr$premises, "gap_on_premises",    "Not on premises (%)",             "Not on premises", diag_line = TRUE),
    make_scatter(cr$contam,   "gap_free_contam",    "Not free from contamination (%)", "Not free from contamination", diag_line = TRUE))
  adj <- list(
    pg(cr$stress,   "Water stress (log) | GNI",          "Water stress"),
    pg(cr$basic,    "No basic drinking water | GNI",     "No basic drinking water"),
    pg(cr$safely,   "No safely managed | GNI",           "No safely managed service"),
    pg(cr$avail,    "Not available when needed | GNI",   "Not available when needed"),
    pg(cr$premises, "Not on premises | GNI",             "Not on premises"),
    pg(cr$contam,   "Not free from contamination | GNI", "Not free from contamination"))
  gni <- make_scatter(cr$gni, "gni_adj", "GNI per capita (log)", "GNI per capita")
  if (hide_predicted) {   # no predicted countries: drop the empty legend
    raw <- lapply(raw, `+`, guides(colour = "none"))
    adj <- lapply(adj, `+`, guides(colour = "none"))
    gni <- gni + guides(colour = "none")
  }
  list(raw = tag_panels(raw), adj = tag_panels(adj), gni = gni)
}
save_grid <- function(pl, file, title, ncol = 3) {
  n_row <- ceiling(length(pl) / ncol)
  save_fig(collect_bottom(wrap_plots(pl, ncol = ncol)) + plot_annotation(title = title),
           file, width = 5.4 * ncol, height = 5 * n_row + 1.5, dpi = 300)
}

# All 194 countries
sets_all <- panel_sets(cors)
save_grid(sets_all$raw[1:4], "C3_fig4_correlations.png",
          "Correlations with water insecurity: all 194 countries (observed and predicted)", ncol = 2)   # Figure 4
save_grid(sets_all$adj, "C3_figS04_partial_correlations_gni.png",
          "GNI-adjusted (partial) correlations: all 194 countries")                                     # Figure S4
save_grid(sets_all$raw, "C3_all_correlations.png",
          "Correlations with water insecurity: all 194 countries (observed and predicted)")             # all 6, not in paper
save_fig(sets_all$gni, "C3_wi_vs_gni_all_countries.png", width = 6.5, height = 6, dpi = 300)

# Observed countries only (Tier 1)
if (observed_only) {
  dat_obs  <- dat %>% filter(source_type == "Observed")
  cors_obs <- run_cors(dat_obs)
  save_table(bind_rows(lapply(cors_obs, `[[`, "row")), "C3_wi_correlations_observed_only")
  sets_obs <- panel_sets(cors_obs, hide_predicted = TRUE)
  obs_lab  <- paste0("observed countries only (n = ", nrow(dat_obs), ")")
  save_grid(sets_obs$raw, "C3_figS05_correlations_observed_only.png",
            paste("Correlations with water insecurity:", obs_lab))                                      # Figure S5
  save_grid(sets_obs$adj, "C3_figS06_partial_correlations_gni_observed_only.png",
            paste("GNI-adjusted (partial) correlations:", obs_lab))                                     # Figure S6
  save_fig(sets_obs$gni + labs(caption = obs_lab), "C3_wi_vs_gni_observed_only.png", width = 6.5, height = 6, dpi = 300)
}

stop_timer()
