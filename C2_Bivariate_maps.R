################################################################################
#  C2_Bivariate_maps.R — BIVARIATE MAPS: WATER INSECURITY x WATER INDICATORS
#
#  Purpose : Show, on one map, the final water insecurity estimate for every
#            country (all three tiers) jointly with one water indicator, using a
#            3 x 3 colour grid (low/medium/high on each variable):
#              x water stress, x lack of basic drinking water, x lack of
#              safely managed drinking water (SMDWS), x water not available
#              when needed.
#            ALTERNATIVE VERSIONS: the paper's Figure 3 (A-C) and Figure S7 are
#            drawn by C4 on UN boundaries with disputed areas; C2 is kept as a
#            backup (same data, colours and classes).
#            Water indicators use ORIGINAL values only (wb_clean, not MICE-imputed);
#            countries without the indicator are shown as "No data".
#
#  Class thresholds (fixed, not quantiles):
#    Water insecurity prevalence:            <10% | 10-<30% | >=30%
#    % without basic / SMDWS / availability: <10% | 10-<30% | >=30%
#    Water stress (% withdrawal, raw scale): <10% | 10-<40% | >=40%
#
#  Depends : B8 (all_tiers), A1 (wb_clean)
#  Outputs : figures/C2_*.png; figures/C2_panels/ (each map without its key, and
#            each 3x3 key on its own); data/C2_plots.rds (ggplot objects)
################################################################################

source("00_setup.R")
start_timer("C2")
load_pkgs("sf")

all_tiers <- load_data("all_tiers")
wb        <- load_data("wb_clean")     # original (non-imputed) indicator values

class3 <- function(x, cut1, cut2) case_when(is.na(x) ~ NA_character_, x < cut1 ~ "1", x < cut2 ~ "2", TRUE ~ "3")

map_data <- all_tiers %>%
  select(country_code, prevalence) %>%
  left_join(wb %>% select(iso3, basic_water, water_stress, water_available, water_safely_managed),
            by = c("country_code" = "iso3")) %>%
  mutate(prev_class    = class3(prevalence * 100, 10, 30),
         water_class   = class3(100 - basic_water, 10, 30),
         stress_class  = class3(water_stress, 10, 40),
         avail_class   = class3(100 - water_available, 10, 30),
         managed_class = class3(100 - water_safely_managed, 10, 30),
         across(c(water_class, stress_class, avail_class, managed_class),
                ~ ifelse(is.na(prev_class) | is.na(.x), "No data", paste0(prev_class, "-", .x)),
                .names = "bivar_{.col}"))

# CHECK: number of countries with "No data" for each indicator
message("Countries without indicator data — basic: ", sum(map_data$bivar_water_class == "No data"),
        " | stress: ", sum(map_data$bivar_stress_class == "No data"),
        " | availability: ", sum(map_data$bivar_avail_class == "No data"),
        " | SMDWS: ", sum(map_data$bivar_managed_class == "No data"))

bivar_colours <- c("1-1" = "#e8e8e8", "2-1" = "#e4acac", "3-1" = "#c85a5a",
                   "1-2" = "#cbc9e2", "2-2" = "#b58aa5", "3-2" = "#9e3a6e",
                   "1-3" = "#9e9ac8", "2-3" = "#7a3d6e", "3-3" = "#5b1a3a",
                   "No data" = "grey70")

world_map <- get_un_world_map()

bivar_map <- function(class_col, title, legend_label) {
  d <- world_map %>% left_join(map_data %>% select(country_code, all_of(class_col)),
                               by = c("iso3_join" = "country_code"))
  p_map <- ggplot(d) +
    geom_sf(aes(fill = .data[[class_col]]), colour = "grey40", linewidth = 0.12) +
    scale_fill_manual(values = bivar_colours, na.value = "grey70", guide = "none") +
    coord_sf(crs = map_crs, default_crs = st_crs(4326), xlim = c(-180, 180), ylim = c(-56, 84), expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(8, 8, 8, 8))

  p_legend <- expand.grid(x = 1:3, y = 1:3) %>%
    mutate(fill = paste0(x, "-", y)) %>%
    ggplot(aes(x = x, y = y, fill = fill)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = bivar_colours, guide = "none") +
    scale_x_continuous(breaks = c(1, 3), labels = c("Low", "High")) +
    scale_y_continuous(breaks = c(1, 3), labels = c("Low", "High")) +
    labs(x = "WI Prevalence →", y = paste0(legend_label, " →")) +
    theme_void() +
    theme(axis.title.x = element_text(size = 7, margin = margin(t = 2)),
          axis.title.y = element_text(size = 7, angle = 90, margin = margin(r = 2)),
          axis.text = element_text(size = 6), plot.margin = margin(2, 2, 2, 2)) +
    coord_equal()

  list(combined = p_map + inset_element(p_legend, left = 0.02, bottom = 0.02, right = 0.18, top = 0.25),
       map = p_map, legend = p_legend)
}

p_stress  <- bivar_map("bivar_stress_class",  "Water insecurity prevalence × Water stress", "Water Stress")
p_basic   <- bivar_map("bivar_water_class",   "Water insecurity prevalence × Lack of basic drinking water access",
                       "No Basic Drinking (BDWS)")
p_managed <- bivar_map("bivar_managed_class", "Water insecurity prevalence × Lack of access to safely managed drinking water",
                       "No SMDWS")
p_avail   <- bivar_map("bivar_avail_class",   "Water insecurity prevalence × Water not available when needed",
                       "Water Not Available")

save_fig(p_stress$combined,  "C2_alt_bivariate_water_stress.png", width = 12, height = 7, dpi = 600)
save_fig(p_basic$combined,   "C2_alt_bivariate_basic_water.png", width = 12, height = 7, dpi = 600)
save_fig(p_managed$combined, "C2_alt_bivariate_smdws.png",         width = 12, height = 7, dpi = 600)
save_fig(p_avail$combined,   "C2_alt_bivariate_availability.png", width = 12, height = 7, dpi = 600)
save_fig(p_basic$combined / p_managed$combined, "C2_alt_bivariate_basic_smdws_combined.png",
         width = 12, height = 12, dpi = 600)

# ── Each map and its 3x3 colour key saved on their own (to combine differently) ──
# figures/C2_panels/: map without the inset key, and the key alone;
# data/C2_plots.rds: the ggplot objects (list of map / legend / combined).
maps <- list(stress = p_stress, basic = p_basic, smdws = p_managed, availability = p_avail)
dir.create(file.path(dir_figures, "C2_panels"), showWarnings = FALSE)
for (nm in names(maps)) {
  save_fig(maps[[nm]]$map,    file.path("C2_panels", paste0("C2_map_", nm, ".png")),
           width = 12, height = 6.5, dpi = 600)
  save_fig(maps[[nm]]$legend, file.path("C2_panels", paste0("C2_key_", nm, ".png")),
           width = 2.2, height = 2.2, dpi = 600)
}
save_data(maps, "C2_plots")

stop_timer()
