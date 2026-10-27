################################################################################
#  C1_Global_Estimates_Figures_Maps.R — MAIN FIGURES FOR THE GLOBAL ESTIMATES
#
#  Purpose : Draw the figures from the final estimates produced in B8
#            (three tiers: observed, individual-level model, country-level model).
#              Figure 2    prevalence and number affected, with 95% UIs,
#                          by M49 region (A1, A2) and by income group (B1, B2)
#              Figure S1   forest plot of prevalence by M49 region
#              Figure S2   prevalence by M49 sub-region (95% UIs)
#              Figure S3   number affected by M49 sub-region (95% UIs)
#              Maps (alternative versions): prevalence and data source. The
#                          paper's Figure 1 is drawn by C4 (UN boundaries with
#                          disputed areas); C1's maps are kept as a backup.
#              (extra)     prevalence by data scenario (observed only / + individual / all)
#
#  Depends : B8
#  Inputs  : data/all_tiers.rds, group_UIs_B11.rds (or group_CIs.rds), tier_contributions.rds,
#            tables/B8_prevalence_by_data_scenario.csv, UN boundaries (downloaded once)
#  Outputs : figures/C1_fig2A_region_with_CIs.png, C1_fig2B_income_with_CIs.png,
#            C1_figS01_forest_m49_region.png, C1_figS02_prev_subregion_CIs.png,
#            C1_figS03_pop_subregion_CIs.png; alternative maps C1_alt_map_*.png
#            (figures/C1_panels/); data/C1_map_plots.rds (ggplot objects)
#
#  What to look out for:
#    * "Countries in estimates but not on the map" should be empty (if not,
#      the ISO3 code differs between the estimates and the UN boundary file).
################################################################################

source("00_setup.R")
start_timer("C1")
load_pkgs("sf")

all_tiers <- load_data("all_tiers")
# Group intervals: country-bootstrap UIs from B11 (main) or B8 (interval_source in 00_setup.R)
cis <- if (interval_source == "B11" && file.exists(file.path(dir_data, "group_UIs_B11.rds"))) {
  load_data("group_UIs_B11")
} else {
  if (interval_source == "B11") warning("group_UIs_B11.rds not found (run B11) - using B8 intervals instead")
  load_data("group_CIs")
}
message("Group intervals used: ", if (!is.null(cis$region$interval)) cis$region$interval[1] else "B8")
contrib   <- load_data("tier_contributions")
scenario  <- read_csv(file.path(dir_tables, "B8_prevalence_by_data_scenario.csv"), show_col_types = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# 1. MAPS (alternative to Figure 1 in C4)
# ══════════════════════════════════════════════════════════════════════════════

prev_colours <- c("<5%" = "#FFFFB2", "5-10%" = "#FED976", "10-20%" = "#FEB24C",
                  "20-30%" = "#FD8D3C", "30-40%" = "#FC4E2A", "40-50%" = "#E31A1C",
                  "50-60%" = "#BD0026", ">60%" = "#800026")

map_data <- all_tiers %>%
  mutate(prev_bucket = cut(prevalence, breaks = c(-Inf, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, Inf),
                           labels = names(prev_colours), right = FALSE),
         tier = as.character(tier))

world <- get_un_world_map()

# CHECK: every estimated country should appear on the map
message("Countries in estimates but not on the map: ",
        paste(setdiff(map_data$country_code, world$iso3_join), collapse = ", "))

world_map <- world %>%
  left_join(map_data %>% select(country_code, prev_bucket, tier, prevalence),
            by = c("iso3_join" = "country_code")) %>%
  mutate(tier = factor(ifelse(is.na(tier), "No data", tier), levels = c(tier_levels, "No data")),
         predicted = tier %in% tier_levels[2:3])

# Diagonal hatching for predicted countries
make_hatch <- function(x, spacing = 1.5, slope = 1) {
  crs <- st_crs(x)
  g   <- st_union(st_make_valid(st_set_crs(st_geometry(x), NA)))
  bb  <- st_bbox(g)
  x0 <- bb[["xmin"]] - spacing; x1 <- bb[["xmax"]] + spacing
  cs <- seq(bb[["ymin"]] - slope * x1, bb[["ymax"]] - slope * x0, by = spacing)
  lines <- lapply(cs, function(cc) st_linestring(cbind(c(x0, x1), slope * c(x0, x1) + cc)))
  out <- suppressWarnings(st_intersection(st_sfc(lines), g))
  st_sf(geometry = st_set_crs(st_collection_extract(out, "LINESTRING"), crs))
}
hatch <- make_hatch(world_map %>% filter(predicted),
                    spacing = as.numeric(diff(st_bbox(world_map)[c("xmin", "xmax")])) / 240)

# Legend key that shows the hatching
observed_key_col <- "#FFFFFF"
draw_key_hatch <- function(data, params, size) {
  col <- if (is.null(data$colour)) "grey15" else data$colour[1]
  box <- grid::rectGrob(gp = grid::gpar(fill = "white", col = "grey40", lwd = 0.6))
  if (is.na(col) || toupper(col) %in% c(toupper(observed_key_col), "WHITE")) return(box)
  cs <- seq(-0.85, 0.85, by = 0.28); xs <- ys <- numeric(0)
  for (cc in cs) {
    if (cc >= 0) { xs <- c(xs, 0, 1 - cc); ys <- c(ys, cc, 1) } else { xs <- c(xs, -cc, 1); ys <- c(ys, 0, 1 + cc) }
  }
  grid::grobTree(box, grid::polylineGrob(x = xs, y = ys, id = rep(seq_along(cs), each = 2),
                                         gp = grid::gpar(col = col, lwd = 0.9)))
}

map_theme <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13, margin = margin(b = 6)),
        legend.position = "bottom", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), legend.key.size = unit(0.55, "cm"),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(8, 8, 8, 8))
map_coord <- coord_sf(crs = map_crs, default_crs = st_crs(4326),
                      xlim = c(-180, 180), ylim = c(-56, 84), expand = FALSE)

p_map_prev <- ggplot() +
  geom_sf(data = world_map, aes(fill = prev_bucket), colour = "grey40", linewidth = 0.12) +
  geom_sf(data = mutate(hatch, key = "Predicted prevalence"), aes(colour = key),
          linewidth = 0.22, key_glyph = draw_key_hatch) +
  geom_sf(data = mutate(hatch[1, ], key = "Observed prevalence"), aes(colour = key),
          alpha = 0, linewidth = 0, key_glyph = draw_key_hatch) +
  scale_fill_manual(values = prev_colours, name = "Prevalence", na.value = "grey90", drop = FALSE,
                    guide = guide_legend(nrow = 2, byrow = FALSE, order = 1)) +
  scale_colour_manual(values = c("Observed prevalence" = observed_key_col, "Predicted prevalence" = "grey15"),
                      limits = c("Observed prevalence", "Predicted prevalence"), drop = FALSE, name = NULL,
                      guide = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  map_coord + labs(title = "a) Estimated water insecurity prevalence") + map_theme +
  theme(legend.box = "horizontal")

p_map_source <- ggplot(world_map) +
  geom_sf(aes(fill = tier), colour = "grey40", linewidth = 0.12) +
  scale_fill_manual(values = c(tier_col, "No data" = "grey90"), name = "Data source") +
  map_coord + labs(title = "b) Data source") + map_theme

save_fig(p_map_prev / p_map_source + plot_layout(heights = c(1, 0.8)),
         "C1_alt_map_prevalence_source.png", width = 11, height = 12, dpi = 300)

# Each map on its own (to combine differently later)
dir.create(file.path(dir_figures, "C1_panels"), showWarnings = FALSE)
save_fig(p_map_prev,   file.path("C1_panels", "C1_alt_map_prevalence.png"), width = 11, height = 6.5, dpi = 300)
save_fig(p_map_source, file.path("C1_panels", "C1_alt_map_data_source.png"), width = 11, height = 5.5, dpi = 300)
save_data(list(map_prevalence = p_map_prev, map_source = p_map_source), "C1_map_plots")


# ══════════════════════════════════════════════════════════════════════════════
# 2. STACKED BARS WITH CIs (Figure 2 A/B)
# ══════════════════════════════════════════════════════════════════════════════
# Bars: contribution of each tier to the population-weighted prevalence (left)
# and to the number of people affected (right). Error bars: 95% UI (B11 or B8,
# see interval_source in 00_setup.R).

plot_stacked <- function(contrib_tab, ci_tab, group_var, group_levels, x_title = NULL) {
  lv <- c(intersect(group_levels, unique(ci_tab[[group_var]])), "Global")
  ct <- contrib_tab %>%
    mutate(tier = factor(tier, levels = rev(tier_levels)),
           grp  = factor(.data[[group_var]], levels = lv))
  ci <- ci_tab %>% mutate(grp = factor(.data[[group_var]], levels = lv))
  vline <- length(lv) - 0.5

  bar_theme <- theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8), legend.position = "none",
          plot.title = element_text(face = "bold"), panel.grid.major.x = element_blank())

  p_prev <- ggplot(ct, aes(x = grp, y = prev_contribution, fill = tier)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    geom_errorbar(data = ci, aes(x = grp, ymin = ci_lower, ymax = ci_upper),
                  inherit.aes = FALSE, width = 0.15, linewidth = 0.2) +
    geom_text(data = ci, aes(x = grp, y = ci_upper,
                             label = paste0(round(prevalence * 100, 1), "%\n(", round(ci_lower * 100, 1),
                                            "-", round(ci_upper * 100, 1), "%)")),
              inherit.aes = FALSE, vjust = -0.3, size = 2.8, fontface = "bold", lineheight = 0.85) +
    geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = tier_col, name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Water Insecurity Prevalence (%)", subtitle = paste("Population-weighted estimate with", interval_label),
         x = x_title, y = "Prevalence") +
    bar_theme

  p_pop <- ggplot(ct, aes(x = grp, y = affected_millions, fill = tier)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    geom_errorbar(data = ci, aes(x = grp, ymin = affected_lower_m, ymax = affected_upper_m),
                  inherit.aes = FALSE, width = 0.15, linewidth = 0.2) +
    geom_text(data = ci, aes(x = grp, y = affected_upper_m,
                             label = paste0(round(affected_millions, 0), "M\n(", round(affected_lower_m, 0),
                                            "-", round(affected_upper_m, 0), "M)")),
              inherit.aes = FALSE, vjust = -0.3, size = 2.8, fontface = "bold", lineheight = 0.85) +
    geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = tier_col, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Population Water Insecure (millions)", subtitle = paste("Adults 15+ with", interval_label),
         x = x_title, y = "Millions") +
    bar_theme

  (p_prev | p_pop) + plot_layout(guides = "collect") &
    theme(legend.position = "top") & guides(fill = guide_legend(reverse = TRUE))
}

# Region order as in the original figures
fig_region_order <- c("Africa", "Asia", "Americas", "Europe", "Oceania")

fig2a <- plot_stacked(contrib$region, cis$region, "m49_region", fig_region_order)
fig2b <- plot_stacked(contrib$income, cis$income, "income_group", income_levels)
save_fig(fig2a, "C1_fig2A_region_with_CIs.png", width = 14, height = 6, dpi = 300)
save_fig(fig2b, "C1_fig2B_income_with_CIs.png", width = 14, height = 6, dpi = 300)
save_fig(fig2a / fig2b, "C1_fig2_combined.png", width = 14, height = 12, dpi = 300)



# ══════════════════════════════════════════════════════════════════════════════
# 3. FIGURES S2 AND S3 — M49 SUB-REGIONS, GROUPED BY PARENT REGION
# ══════════════════════════════════════════════════════════════════════════════
# Sub-regions are grouped (facets) by M49 region, ordered Africa, Asia, Oceania,
# Americas, Europe, and by prevalence (highest first) within each region.

sub_parent_order <- c("Africa", "Asia", "Oceania", "Americas", "Europe")
sub_to_region <- all_tiers %>%
  filter(!is.na(m49_sub_region)) %>%
  distinct(m49_sub_region = as.character(m49_sub_region), parent = as.character(m49_region))

ci_sub <- cis$sub_region %>%
  left_join(sub_to_region, by = "m49_sub_region") %>%
  mutate(parent = ifelse(m49_sub_region == "Global", "Global", parent))

sub_order <- ci_sub %>%
  filter(m49_sub_region != "Global") %>%
  mutate(parent = factor(parent, levels = sub_parent_order)) %>%
  arrange(parent, desc(prevalence)) %>%
  pull(m49_sub_region)
facet_levels <- c(intersect(sub_parent_order, ci_sub$parent), "Global")

ci_sub <- ci_sub %>%
  mutate(m49_sub_region = factor(m49_sub_region, levels = c(sub_order, "Global")),
         region_group   = factor(parent, levels = facet_levels))
ct_sub <- contrib$sub_region %>%
  left_join(sub_to_region, by = "m49_sub_region") %>%
  mutate(parent         = ifelse(m49_sub_region == "Global", "Global", parent),
         m49_sub_region = factor(m49_sub_region, levels = c(sub_order, "Global")),
         region_group   = factor(parent, levels = facet_levels),
         tier           = factor(tier, levels = rev(tier_levels)))

sub_theme <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        strip.text = element_text(size = 7.5, face = "bold", colour = "grey25"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.placement = "outside", panel.spacing.x = unit(2, "pt"),
        legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank())
sub_facets <- facet_grid(cols = vars(region_group), scales = "free_x", space = "free_x", switch = "x")

p_prev_sub <- ggplot(ct_sub, aes(x = m49_sub_region, y = prev_contribution, fill = tier)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_errorbar(data = ci_sub, aes(x = m49_sub_region, ymin = ci_lower, ymax = ci_upper),
                inherit.aes = FALSE, width = 0.15, linewidth = 0.2) +
  geom_text(data = ci_sub, aes(x = m49_sub_region, y = ci_upper,
                               label = paste0(round(prevalence * 100, 1), "%\n(", round(ci_lower * 100, 1),
                                              "-", round(ci_upper * 100, 1), "%)")),
            inherit.aes = FALSE, vjust = -0.3, size = 2.2, fontface = "bold", lineheight = 0.85) +
  sub_facets +
  scale_fill_manual(values = tier_col, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Water Insecurity Prevalence (%)",
       subtitle = paste("Population-weighted estimate by sub-region with", interval_label), x = NULL, y = "Prevalence") +
  sub_theme

p_pop_sub <- ggplot(ct_sub, aes(x = m49_sub_region, y = affected_millions, fill = tier)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_errorbar(data = ci_sub, aes(x = m49_sub_region, ymin = affected_lower_m, ymax = affected_upper_m),
                inherit.aes = FALSE, width = 0.15, linewidth = 0.2) +
  geom_text(data = ci_sub, aes(x = m49_sub_region, y = affected_upper_m,
                               label = paste0(round(affected_millions, 0), "M\n(", round(affected_lower_m, 0),
                                              "-", round(affected_upper_m, 0), "M)")),
            inherit.aes = FALSE, vjust = -0.3, size = 2.2, fontface = "bold", lineheight = 0.85) +
  sub_facets +
  scale_fill_manual(values = tier_col, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Population Water Insecure (millions)",
       subtitle = paste("Adults 15+ by sub-region with", interval_label), x = NULL, y = "Millions") +
  sub_theme

save_fig(p_prev_sub, "C1_figS02_prev_subregion_CIs.png", width = 12, height = 5.5, dpi = 300)
save_fig(p_pop_sub,  "C1_figS03_pop_subregion_CIs.png",  width = 12, height = 5.5, dpi = 300)
save_fig((p_prev_sub | p_pop_sub) + plot_layout(guides = "collect") &
           theme(legend.position = "top") & guides(fill = guide_legend(reverse = TRUE)),
         "C1_figS02_S03_subregion_with_CIs.png", width = 16, height = 6.5, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 4. FIGURE S1 — FOREST PLOT BY M49 REGION
# ══════════════════════════════════════════════════════════════════════════════

forest_data <- cis$region %>%
  mutate(m49_region = factor(m49_region, levels = rev(c(fig_region_order, "Global"))))

p_forest <- ggplot(forest_data, aes(x = prevalence, y = m49_region)) +
  geom_vline(xintercept = forest_data$prevalence[forest_data$m49_region == "Global"],
             linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper), linewidth = 0.8) +
  geom_point(size = 3, colour = "#0072B2") +
  geom_text(aes(label = paste0(round(prevalence * 100, 1), "%")), vjust = -1, size = 3) +
  geom_text(aes(label = paste0("(", round(ci_lower * 100, 1), "-", round(ci_upper * 100, 1), "%)")),
            vjust = 2, size = 2.5, colour = "grey40") +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.05, 0.1))) +
  labs(title = "Water insecurity prevalence by region",
       subtitle = "Population-weighted estimates with 95% uncertainty intervals", x = "Prevalence", y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())
save_fig(p_forest, "C1_figS01_forest_m49_region.png", width = 10, height = 5, dpi = 600)


# ══════════════════════════════════════════════════════════════════════════════
# 5. (EXTRA) PREVALENCE BY DATA SCENARIO
# ══════════════════════════════════════════════════════════════════════════════
# Shows how regional prevalence changes as predicted countries are added.

p_scenario <- scenario %>%
  mutate(scenario = factor(scenario, levels = unique(scenario)),
         m49_region = factor(m49_region, levels = c(fig_region_order, "Global"))) %>%
  ggplot(aes(x = m49_region, y = prevalence, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0("n=", n_countries)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 2.5) +
  geom_vline(xintercept = length(m49_region_levels) + 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = c("navyblue", "steelblue2", "turquoise3"), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Water insecurity prevalence by region",
       subtitle = "Population-weighted estimates across three data tiers", x = NULL, y = "Prevalence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top",
        plot.title = element_text(face = "bold"), panel.grid.major.x = element_blank())
save_fig(p_scenario, "C1_prevalence_by_data_scenario.png", width = 11, height = 6, dpi = 600)

stop_timer()
