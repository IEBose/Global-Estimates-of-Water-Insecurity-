################################################################################
#  C4_UN_boundary_maps.R — ALL MAPS ON UN BOUNDARIES, WITH DISPUTED AREAS
#
#  Purpose : Alternative versions of every map in the paper on the UN Geospatial boundaries
#            (BNDA_simplified), with disputed or undetermined areas outlined
#            with dashed lines and filled with the estimate set in
#            `disputed_rules` (section 1):
#              Figure 1    a) prevalence (predicted countries hatched)
#                          b) data source (tier)
#              Figure 3    bivariate maps: prevalence x A) water stress,
#                          B) lack of basic drinking water, C) lack of SMDWS
#              Figure S7   bivariate: prevalence x water not available when needed
#            ALTERNATIVE versions (C4_alt_*): the paper maps are now drawn by
#            C5 to the WHO map SOPs. Same data, colours and class thresholds
#            as C1 (maps), C2 (bivariate maps) and C5.
#
#  Disputed areas (section 1, `disputed_rules`)
#    * Polygons listed in the rules take the estimate (and, for the bivariate
#      maps, the indicator values) of the country in `data_iso3`, and get a
#      dashed outline if `dashed = TRUE`.
#    * Any other UN polygon with undetermined status (99), and Western Sahara,
#      gets no data (grey) and a dashed outline.
#    * Crimea, Golan Heights, Kosovo, Northern Cyprus and Somaliland have no
#      polygon of their own in the UN file, so they take Ukraine, Syria,
#      Serbia, Cyprus and Somalia.
#    * The objectid values were taken from BNDA_simplified; section 1 prints
#      the UN name of every rule polygon so they can be checked if the file
#      is updated.
#
#  Depends : B8 (all_tiers), A1 (wb_clean)
#  Inputs  : UN BNDA_simplified boundaries (INPUT/maps; downloaded once)
#  Outputs : figures/C4_alt_map_prevalence_source.png, C4_alt_bivariate_maps.png,
#            C4_alt_bivariate_availability.png (alternatives to Figures 1, 3, S7;
#            paper versions from C5); figures/C4_panels/ (each map and key);
#            tables/C4_disputed_polygons_check.csv; data/C4_map_plots.rds
#
#  What to look out for:
#    * "Countries in estimates but not on the map" should be empty.
#    * "Rule polygons not found" should be empty; check the printed UN names
#      match the `area` column of disputed_rules.
################################################################################

source("00_setup.R")
start_timer("C4")
load_pkgs("sf")

all_tiers <- load_data("all_tiers")
wb        <- load_data("wb_clean")     # original (non-imputed) indicator values, as C2


# ══════════════════════════════════════════════════════════════════════════════
# 1. UN BOUNDARIES AND DISPUTED AREAS
# ══════════════════════════════════════════════════════════════════════════════

# Same as get_un_world_map() in 00_setup.R, but keeps the polygon id and UN
# status code (separate cache, so C1/C2 are not affected)
get_un_world_map_v2 <- function() {
  cache <- file.path(dir_data, "world_map_un_v2.rds")
  if (file.exists(cache)) return(readRDS(cache))
  if (!dir.exists(dir_un_bnda)) {
    dir.create(dirname(file_un_bnda_zip), recursive = TRUE, showWarnings = FALSE)
    download.file(url_un_bnda, destfile = file_un_bnda_zip, mode = "wb")
    unzip(file_un_bnda_zip, exdir = dir_un_bnda)
  }
  world <- sf::st_read(dir_un_bnda, quiet = TRUE) |>
    sf::st_make_valid() |>
    dplyr::select(objectid, iso3_un = iso3cd, status = stscod, name = nam_en, geometry) |>
    dplyr::filter(!is.na(iso3_un)) |>
    sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")) |>
    sf::st_transform(map_crs) |>
    sf::st_make_valid()
  saveRDS(world, cache)
  world
}

# objectid  : polygon id in BNDA_simplified
# data_iso3 : whose estimate fills it (NA = no data)
# dashed    : TRUE = dashed outline (disputed / undetermined)
disputed_rules <- tibble::tribble(
  ~area,                    ~objectid, ~data_iso3, ~dashed,
  # Administered in practice by
  "Western Sahara",          67,        "MAR",      TRUE,
  "Jammu and Kashmir",       240,       "IND",      TRUE,   # one UN polygon: Indian- and
                                                            # Pakistani-administered parts
  "Arunachal Pradesh",       249,       "IND",      TRUE,
  "Aksai Chin",              248,       "CHN",      TRUE,
  "Southern Kuril Islands",  254,       "RUS",      TRUE,
  # Follow the UN file's own country code
  "Hala'ib Triangle",        250,       "SDN",      TRUE,
  "Bir Tawil",               251,       "EGY",      TRUE,
  "Ilemi Triangle",          252,       "SSD",      TRUE,
  # UN status 5, administered by China (China's estimate is a mainland sample)
  "Taiwan",                  225,       "CHN",      FALSE,
  "Hong Kong",               94,        "CHN",      FALSE,
  "Macao",                   130,       "CHN",      FALSE
)

apply_disputed_rules <- function(world, rules = disputed_rules) {
  world |>
    dplyr::left_join(dplyr::select(rules, objectid, data_iso3, rule_dashed = dashed), by = "objectid") |>
    dplyr::mutate(
      undetermined = status == 99 | iso3_un == "ESH",
      # rule polygons take the rule's country; other undetermined (99) or
      # Western Sahara polygons never inherit data through their own code
      iso3_join = dplyr::case_when(
        !is.na(rule_dashed) ~ data_iso3,
        undetermined        ~ NA_character_,
        TRUE                ~ iso3_un),
      # dashed: as set in the rules; other undetermined polygons are dashed too
      dashed = dplyr::coalesce(rule_dashed, undetermined)
    ) |>
    dplyr::select(-rule_dashed)
}

world <- get_un_world_map_v2() |> apply_disputed_rules()

# CHECK 1: every rule polygon exists, and its UN name matches the rule (printed)
rule_check <- disputed_rules %>%
  left_join(world %>% st_drop_geometry() %>% select(objectid, un_name = name, iso3_un, status),
            by = "objectid")
# (several disputed polygons have no English name in the UN file, so the check
#  uses the UN code: xjk, xac, xap, xri = Jammu and Kashmir, Aksai Chin,
#  Arunachal Pradesh, Kuril Islands; Hala'ib, Bir Tawil, Ilemi carry SDN, EGY, SSD)
message("Rule polygons not found in the UN file (expect none): ",
        paste(rule_check$area[is.na(rule_check$iso3_un)], collapse = ", "))
print(rule_check, n = Inf, width = Inf)

# CHECK 2: all other undetermined / special polygons (no data, dashed)
other_dashed <- world %>% st_drop_geometry() %>%
  filter(dashed, !objectid %in% disputed_rules$objectid) %>%
  select(objectid, name, iso3_un, status) %>%
  as_tibble()
message("Other undetermined polygons shown as no data with dashed outline: ", nrow(other_dashed))
print(other_dashed, n = Inf)
save_table(bind_rows(rule_check %>% mutate(source = "disputed_rules"),
                     other_dashed %>% rename(un_name = name) %>% mutate(source = "status 99 / ESH, not in rules")),
           "C4_disputed_polygons_check")

# CHECK 3: every estimated country should appear on the map (expect none)
message("Countries in estimates but not on the map: ",
        paste(setdiff(all_tiers$country_code, world$iso3_join), collapse = ", "))


# ══════════════════════════════════════════════════════════════════════════════
# 2. SHARED MAP ELEMENTS
# ══════════════════════════════════════════════════════════════════════════════

map_caption <- paste0(
  "Boundaries: UN Geospatial. Dashed lines mark disputed areas or areas whose final status has not been determined;\n",
  "they are shown with the estimate of the administering country (Hala'ib, Bir Tawil and Ilemi follow the UN country code).\n",
  "Jammu and Kashmir, including the Pakistani-administered part, is shown with India's estimate. Taiwan, Hong Kong and Macao\n",
  "are shown with China's estimate, based on a mainland sample. This does not imply endorsement of any claim or reflect the authors' opinions of any boundaries."
)

map_theme <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13, margin = margin(b = 6)),
        legend.position = "bottom", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), legend.key.size = unit(0.55, "cm"),
        plot.caption = element_text(hjust = 0, size = 7, colour = "grey30"),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(8, 8, 8, 8))
map_coord <- coord_sf(crs = map_crs, default_crs = st_crs(4326),
                      xlim = c(-180, 180), ylim = c(-56, 84), expand = FALSE)

# Dashed outline for disputed areas, drawn last so it sits on top of fills and hatching
disputed_layer <- function(d) geom_sf(data = filter(d, dashed), fill = NA, colour = "grey10",
                                      linetype = "22", linewidth = 0.25, inherit.aes = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# 3. FIGURE 1 — PREVALENCE AND DATA SOURCE (as C1)
# ══════════════════════════════════════════════════════════════════════════════

prev_colours <- c("<5%" = "#FFFFB2", "5-10%" = "#FED976", "10-20%" = "#FEB24C",
                  "20-30%" = "#FD8D3C", "30-40%" = "#FC4E2A", "40-50%" = "#E31A1C",
                  "50-60%" = "#BD0026", ">60%" = "#800026")

map_data <- all_tiers %>%
  mutate(prev_bucket = cut(prevalence, breaks = c(-Inf, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, Inf),
                           labels = names(prev_colours), right = FALSE),
         tier = as.character(tier))

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

p_map_prev <- ggplot() +
  geom_sf(data = world_map, aes(fill = prev_bucket), colour = "grey40", linewidth = 0.12) +
  geom_sf(data = mutate(hatch, key = "Predicted prevalence"), aes(colour = key),
          linewidth = 0.22, key_glyph = draw_key_hatch) +
  geom_sf(data = mutate(hatch[1, ], key = "Observed prevalence"), aes(colour = key),
          alpha = 0, linewidth = 0, key_glyph = draw_key_hatch) +
  disputed_layer(world_map) +
  scale_fill_manual(values = prev_colours, name = "Prevalence", na.value = "grey90", drop = FALSE,
                    guide = guide_legend(nrow = 2, byrow = FALSE, order = 1)) +
  scale_colour_manual(values = c("Observed prevalence" = observed_key_col, "Predicted prevalence" = "grey15"),
                      limits = c("Observed prevalence", "Predicted prevalence"), drop = FALSE, name = NULL,
                      guide = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  map_coord + labs(title = "a) Estimated water insecurity prevalence", caption = map_caption) + map_theme +
  theme(legend.box = "horizontal")

p_map_source <- ggplot(world_map) +
  geom_sf(aes(fill = tier), colour = "grey40", linewidth = 0.12) +
  disputed_layer(world_map) +
  scale_fill_manual(values = c(tier_col, "No data" = "grey90"), name = "Data source") +
  map_coord + labs(title = "b) Data source", caption = map_caption) + map_theme

# Combined figure: caption once, under panel b
save_fig((p_map_prev + labs(caption = NULL)) / p_map_source + plot_layout(heights = c(1, 0.8)),
         "C4_alt_map_prevalence_source.png", width = 11, height = 12.5, dpi = 300)

dir.create(file.path(dir_figures, "C4_panels"), showWarnings = FALSE)
save_fig(p_map_prev,   file.path("C4_panels", "C4_alt_map_prevalence.png"), width = 11, height = 7, dpi = 300)
save_fig(p_map_source, file.path("C4_panels", "C4_alt_map_data_source.png"), width = 11, height = 6, dpi = 300)


# ══════════════════════════════════════════════════════════════════════════════
# 4. BIVARIATE MAPS (as C2)
# ══════════════════════════════════════════════════════════════════════════════
# Thresholds (fixed): prevalence and % without access/availability <10 | 10-<30 | >=30;
# water stress (% withdrawal) <10 | 10-<40 | >=40. Disputed polygons take both
# the prevalence and the indicator of the country in disputed_rules.

class3 <- function(x, cut1, cut2) case_when(is.na(x) ~ NA_character_, x < cut1 ~ "1", x < cut2 ~ "2", TRUE ~ "3")

bivar_data <- all_tiers %>%
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

bivar_colours <- c("1-1" = "#e8e8e8", "2-1" = "#e4acac", "3-1" = "#c85a5a",
                   "1-2" = "#cbc9e2", "2-2" = "#b58aa5", "3-2" = "#9e3a6e",
                   "1-3" = "#9e9ac8", "2-3" = "#7a3d6e", "3-3" = "#5b1a3a",
                   "No data" = "grey70")

bivar_map <- function(class_col, title, legend_label) {
  d <- world %>% left_join(bivar_data %>% select(country_code, all_of(class_col)),
                           by = c("iso3_join" = "country_code"))
  p_map <- ggplot(d) +
    geom_sf(aes(fill = .data[[class_col]]), colour = "grey40", linewidth = 0.12) +
    disputed_layer(d) +
    scale_fill_manual(values = bivar_colours, na.value = "grey70", guide = "none") +
    map_coord +
    labs(title = title, caption = map_caption) +
    map_theme + theme(plot.title = element_text(face = "bold", size = 12))

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

  # key placed bottom-left, above the caption
  list(combined = p_map + inset_element(p_legend, left = 0.02, bottom = 0.10, right = 0.18, top = 0.33),
       map = p_map, legend = p_legend)
}

p_stress  <- bivar_map("bivar_stress_class",  "A) Water insecurity prevalence × Water stress", "Water Stress")
p_basic   <- bivar_map("bivar_water_class",   "B) Water insecurity prevalence × Lack of basic drinking water access",
                       "No Basic Drinking (BDWS)")
p_managed <- bivar_map("bivar_managed_class", "C) Water insecurity prevalence × Lack of access to safely managed drinking water",
                       "No SMDWS")
p_avail   <- bivar_map("bivar_avail_class",   "Water insecurity prevalence × Water not available when needed",
                       "Water Not Available")

# Figure 3: maps A-C stacked, each with its own key; caption once, under map C
with_key <- function(m, caption = FALSE) {
  mm <- if (caption) m$map else m$map + labs(caption = NULL)
  mm + inset_element(m$legend, left = 0.02, bottom = if (caption) 0.10 else 0.02,
                     right = 0.18, top = if (caption) 0.33 else 0.25)
}
save_fig(with_key(p_stress) / with_key(p_basic) / with_key(p_managed, caption = TRUE) +
           plot_layout(heights = c(1, 1, 1.12)),
         "C4_alt_bivariate_maps.png", width = 12, height = 21, dpi = 600)                    # alternative to Figure 3
# Figure S7: water not available when needed
save_fig(p_avail$combined, "C4_alt_bivariate_availability.png", width = 12, height = 7.5, dpi = 600)
# Each Figure 3 map on its own, with key and caption
save_fig(p_stress$combined,  file.path("C4_panels", "C4_alt_bivariate_water_stress.png"), width = 12, height = 7.5, dpi = 600)
save_fig(p_basic$combined,   file.path("C4_panels", "C4_alt_bivariate_basic_water.png"),  width = 12, height = 7.5, dpi = 600)
save_fig(p_managed$combined, file.path("C4_panels", "C4_alt_bivariate_smdws.png"),        width = 12, height = 7.5, dpi = 600)

maps <- list(stress = p_stress, basic = p_basic, smdws = p_managed, availability = p_avail)
for (nm in names(maps)) {
  save_fig(maps[[nm]]$map,    file.path("C4_panels", paste0("C4_map_", nm, ".png")),
           width = 12, height = 7, dpi = 600)
  save_fig(maps[[nm]]$legend, file.path("C4_panels", paste0("C4_key_", nm, ".png")),
           width = 2.2, height = 2.2, dpi = 600)
}
save_data(list(map_prevalence = p_map_prev, map_source = p_map_source, bivariate = maps), "C4_map_plots")

stop_timer()
