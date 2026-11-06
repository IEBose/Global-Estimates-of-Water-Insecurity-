################################################################################
#  C5_WHO_style_maps.R — PAPER MAPS DRAWN TO THE WHO MAP SOPs
#
#  Purpose : All maps in the paper on UN boundaries (BNDA_simplified), drawn
#            following the WHO Standard Operating Procedures for maps
#            (HIS/HSI, November 2014), as requested by a WHO co-author:
#              Figure 1    a) prevalence (predicted estimates hatched)
#                          b) data source
#              Figure 3    bivariate maps: prevalence x A) water stress,
#                          B) lack of basic drinking water, C) lack of SMDWS
#              Figure S7   bivariate: prevalence x water not available when needed
#            Same data, colours and classes as C4 (and C1/C2), which are kept as
#            alternative versions.
#
#  Disputed areas and boundaries (WHO SOP section 2.2.2)
#    * Western Sahara, Jammu and Kashmir, Abyei: grey, "Data not available"
#        - Western Sahara: solid boundary with Morocco
#        - Jammu and Kashmir, Abyei: dotted boundary
#    * Dashed boundaries: Egypt-Sudan (Hala'ib, Bir Tawil), South Sudan-Kenya
#      (Ilemi), Sudan-South Sudan, Israel-West Bank/Gaza, Korean armistice
#      line, Serbia-Kosovo (UNSCR 1244)
#    * Arunachal Pradesh (India), Aksai Chin (China), Southern Kurils (Russia):
#      drawn as part of that country, no line of their own
#    * Hong Kong, Macao, Taiwan: part of China (SOP 2.3.2), no separate boundary
#    * Overseas territories: "Data not available" (they have their own UN code)
#    * Disputed lines have the same colour and width as international
#      boundaries and are drawn on top of every other layer
#    * Large lakes (> ~5,000 km2) in the background colour
#    * Figure 1 legends end with "Data not available"; on the bivariate maps
#      grey ("data not available") is explained in the caption, to keep the
#      small 3x3 key readable. WHO-style disclaimer and projection (Equal
#      Earth) below each map
#  "Data not available" is used for the disputed areas too (rather than the
#  SOP's "Not applicable"), as agreed with co-authors; this is not a WHO
#  publication, so the WHO logo / copyright elements are not used.
#
#  Abyei and the Serbia-Kosovo line are not in BNDA_simplified, so they come
#  from Natural Earth 1:10m disputed areas; lakes from Natural Earth 1:50m.
#  Both are downloaded once to INPUT/maps (paths in 00_setup.R).
#
#  Depends : B8 (all_tiers), A1 (wb_clean)
#  Inputs  : UN BNDA_simplified; Natural Earth disputed areas and lakes
#  Outputs : figures/C5_fig1_map_prevalence_source.png (Figure 1),
#            C5_fig3_bivariate_maps.png (Figure 3),
#            C5_figS07_bivariate_availability.png (Figure S7),
#            figures/C5_panels/ (each map and key on its own),
#            tables/C5_disputed_areas_check.csv, data/C5_map_plots.rds
#
#  What to look out for:
#    * "Countries in estimates but not on the map" should be empty.
#    * "Rule polygons not found" should be empty; check the UN codes printed.
#    * Warnings about Israel and West Bank/Gaza in the same colour class (WHO:
#      never the same colour) or about missing dashed / dotted lines.
#    * Look at the maps zoomed in on the disputed areas before submission.
################################################################################

source("00_setup.R")
start_timer("C5")
load_pkgs("sf")

all_tiers <- load_data("all_tiers")
wb        <- load_data("wb_clean")      # original (non-imputed) indicator values

border_col  <- "grey40"    # international AND disputed boundaries (WHO: same colour)
border_lwd  <- 0.08        # international AND disputed boundaries (WHO: same width)
# Dashed lines are cut into short solid pieces on the map itself (make_dashes),
# so dash and gap lengths are set in km and are the same in PNG and PDF. R's own
# dash patterns cannot be shorter than ~0.25 mm, so a short line such as
# Serbia-Kosovo (~2 mm on the page) showed only 1-3 dashes with large gaps.
dash_len    <- 14000       # dash length (m on the map; ~0.1 mm on an 11-inch world map)
dash_gap    <- 6000        # gap length (m)
# Dotted lines (Jammu and Kashmir, Abyei): R's pattern, dots round
lty_dotted  <- "13"        # 1 unit on, 3 off
bg_col      <- "white"     # map background; lakes use the same colour
lab_no_data <- "Data not available"
col_no_data_prev  <- "#EBEBEB"   # prevalence / data-source maps
col_no_data_bivar <- "grey70"    # bivariate maps (#EBEBEB is too close to class 1-1)
# Switch: TRUE also hatches countries with PREDICTED WI prevalence on the
# bivariate maps (Figure 3, S7), as on Figure 1a; FALSE = no hatching
hatch_bivariate <- TRUE
hatch_lwd <- 0.08         # hatch line width for predicted estimates (Figure 1a and bivariate maps)
# Switch: scale bar off by default. On a world map the scale varies by far more
# than 10% with latitude, so a single scale bar is misleading (WHO SOP: optional)
add_scale_bar <- FALSE
# Vector copies of every map, saved next to each PNG (dashes and thin lines stay
# sharp at any zoom). "pdf" uses cairo_pdf; "svg" needs the svglite package.
# Use character(0) for PNG only.
vector_formats <- c("pdf", "svg")


# ══════════════════════════════════════════════════════════════════════════════
# 1. BOUNDARY LAYERS
# ══════════════════════════════════════════════════════════════════════════════

# UN boundaries with polygon id and status (same cache as C4)
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

# Natural Earth layer: downloaded once into INPUT/maps if not there
get_ne_layer <- function(dir_layer, url) {
  if (!dir.exists(dir_layer) || length(list.files(dir_layer, pattern = "\\.shp$")) == 0) {
    dir.create(dir_layer, recursive = TRUE, showWarnings = FALSE)
    zip <- file.path(dirname(dir_layer), paste0(basename(dir_layer), ".zip"))
    download.file(url, destfile = zip, mode = "wb")
    unzip(zip, exdir = dir_layer)
  }
  sf::st_read(dir_layer, quiet = TRUE) |> sf::st_make_valid() |> sf::st_transform(map_crs)
}

# Rules for special UN polygons (objectids checked against BNDA_simplified)
# data_iso3 : whose estimate fills it (NA = grey, "Data not available")
# role      : "not_available" (grey) | "data"
# outline   : "solid" | "dashed" | "dotted"
who_rules <- tibble::tribble(
  ~area,                    ~objectid, ~data_iso3, ~role,           ~outline,
  "Western Sahara",          67,        NA,         "not_available", "solid",
  "Jammu and Kashmir",       240,       NA,         "not_available", "dotted",
  "Arunachal Pradesh",       249,       "IND",      "data",          "solid",
  "Aksai Chin",              248,       "CHN",      "data",          "solid",
  "Southern Kuril Islands",  254,       "RUS",      "data",          "solid",
  "Hala'ib Triangle",        250,       "SDN",      "data",          "dashed",
  "Bir Tawil",               251,       "EGY",      "data",          "dashed",
  "Ilemi Triangle",          252,       "SSD",      "data",          "dashed",
  "Taiwan",                  225,       "CHN",      "data",          "solid",
  "Hong Kong",               94,        "CHN",      "data",          "solid",
  "Macao",                   130,       "CHN",      "data",          "solid"
)
# Rule areas drawn as part of their country (no line of their own)
merge_ids <- who_rules$objectid[who_rules$role == "data" & who_rules$outline == "solid"]

# Dashed boundaries between two UN units (drawn along their shared edge only)
dashed_pairs <- list(c("ISR", "PSE"), c("PRK", "KOR"), c("SDN", "SSD"))

apply_who_rules <- function(world, rules = who_rules) {
  world |>
    left_join(select(rules, objectid, data_iso3, role, outline), by = "objectid") |>
    mutate(
      role = coalesce(role, if_else(status == 99, "not_available", "data")),
      iso3_join = case_when(
        !is.na(outline)         ~ data_iso3,
        role == "not_available" ~ NA_character_,
        TRUE                    ~ iso3_un),
      outline = coalesce(outline, "solid"),
      # merged into their country: China's SARs and Taiwan, Aksai Chin,
      # Arunachal Pradesh, Southern Kurils (so no boundary is drawn around them)
      unit = case_when(objectid %in% merge_ids | iso3_join %in% "CHN" ~ iso3_join,
                       objectid %in% rules$objectid                   ~ paste0("rule_", objectid),
                       TRUE                                           ~ paste0("un_", objectid))
    ) |>
    group_by(unit) |>
    summarise(across(c(iso3_join, role, outline), first), geometry = st_union(geometry),
              .groups = "drop") |>
    st_make_valid()
}

# Edge of `a` that runs along `b` (1 km tolerance, as simplified neighbours may
# not share every vertex); merged so the dash pattern runs continuously
shared_edge <- function(a, b, tol = 1000) {
  e <- suppressWarnings(st_intersection(st_boundary(st_union(a)), st_buffer(st_union(b), tol)))
  e <- suppressWarnings(st_collection_extract(st_sf(geometry = e), "LINESTRING"))
  if (nrow(e) == 0) return(NULL)
  st_sf(geometry = st_line_merge(st_cast(st_union(st_geometry(e)), "MULTILINESTRING")))
}

world_raw <- get_un_world_map_v2()
world     <- apply_who_rules(world_raw)

# CHECK 1: every rule polygon exists in the UN file
rule_check <- who_rules %>%
  left_join(world_raw %>% st_drop_geometry() %>% select(objectid, un_name = name, iso3_un, status),
            by = "objectid")
message("Rule polygons not found in the UN file (expect none): ",
        if (all(!is.na(rule_check$iso3_un))) "none" else paste(rule_check$area[is.na(rule_check$iso3_un)], collapse = ", "))
print(rule_check, n = Inf, width = Inf)

# Dashed lines: edges of the dashed rule polygons with their neighbours, and
# the shared edges of the country pairs
dashed_polys <- filter(world_raw, objectid %in% who_rules$objectid[who_rules$outline == "dashed"])
others <- filter(world_raw, !objectid %in% dashed_polys$objectid)
dashed_list <- c(
  lapply(split(dashed_polys, dashed_polys$objectid), function(p)
    shared_edge(p, others[lengths(st_is_within_distance(others, p, 1000)) > 0, ])),
  lapply(dashed_pairs, function(pr)
    shared_edge(filter(world_raw, iso3_un == pr[1]), filter(world_raw, iso3_un == pr[2]))))
names(dashed_list) <- c(paste0("objectid_", names(split(dashed_polys, dashed_polys$objectid))),
                        sapply(dashed_pairs, paste, collapse = "-"))
missing_dashed <- names(dashed_list)[vapply(dashed_list, is.null, TRUE)]
if (length(missing_dashed) > 0) warning("No dashed line found for: ", paste(missing_dashed, collapse = ", "))
dashed_lines <- bind_rows(dashed_list) %>% mutate(type = "dashed")

# Dotted lines: full outline of Jammu and Kashmir
dotted_lines <- filter(world, outline == "dotted") %>% st_boundary() %>% st_sf() %>%
  select(geometry) %>% mutate(type = "dotted")

# Natural Earth: Abyei (grey, dotted) and the Serbia-Kosovo line (dashed)
ne_disputed <- get_ne_layer(dir_ne_disputed, url_ne_disputed)
name_col <- intersect(c("BRK_NAME", "brk_name", "NAME", "name"), names(ne_disputed))[1]
if (!all(c("Abyei", "Kosovo") %in% ne_disputed[[name_col]]))
  stop("Abyei / Kosovo not found in the Natural Earth disputed-areas layer (column ", name_col, ")")

abyei <- st_sf(geometry = st_union(ne_disputed[ne_disputed[[name_col]] == "Abyei", ]))
kosovo_poly <- st_union(ne_disputed[ne_disputed[[name_col]] == "Kosovo", ])
# Serbia-Kosovo line: the part of the Natural Earth Kosovo outline that lies
# inside UN Serbia (kept 1 km off Serbia's international border), short pieces
# dropped (Natural Earth and UN borders differ slightly in the south), merged
# into one line so the dashes run continuously, and both ends joined to the
# nearest point on Serbia's UN border so the line has no gaps
serbia      <- st_union(filter(world_raw, iso3_un == "SRB"))
serbia_core <- st_buffer(serbia, -1000)
kl <- suppressWarnings(st_intersection(st_boundary(kosovo_poly), serbia_core)) %>%
  st_sf(geometry = .) %>% st_collection_extract("LINESTRING") %>% suppressWarnings() %>%
  st_cast("LINESTRING") %>% suppressWarnings()
kl <- kl[as.numeric(st_length(kl)) > 20000, ]                      # drop fragments < 20 km
kl <- st_line_merge(st_cast(st_union(st_geometry(kl)), "MULTILINESTRING"))
kl <- st_cast(kl, "LINESTRING")
kl <- kl[which.max(as.numeric(st_length(kl)))]                     # the main Serbia-Kosovo line
xy <- st_coordinates(kl)[, c("X", "Y")]
ends <- st_sfc(st_point(xy[1, ]), st_point(xy[nrow(xy), ]), crs = st_crs(kl))
to_border <- st_nearest_points(ends, st_boundary(serbia))           # short joining segments
kosovo_line <- st_sf(geometry = st_line_merge(st_cast(st_union(c(kl, to_border)), "MULTILINESTRING"))) %>%
  mutate(type = "dashed")
if (nrow(kosovo_line) == 0) warning("Serbia-Kosovo line is empty")

# All disputed lines; no UN line is drawn through Abyei
disputed_lines <- bind_rows(dashed_lines, dotted_lines, kosovo_line) %>%
  st_difference(st_union(abyei)) %>%
  st_collection_extract("LINESTRING") %>%
  suppressWarnings() %>%
  bind_rows(st_boundary(abyei) %>% st_sf() %>% select(geometry) %>% mutate(type = "dotted"))

# Ordinary boundaries: every polygon edge except where a disputed line runs,
# so no solid line sits under a dashed or dotted one
solid_lines <- st_union(st_boundary(st_geometry(world))) %>%
  st_difference(st_union(st_buffer(st_union(st_geometry(disputed_lines)), 2000), st_union(abyei))) %>%
  st_sf(geometry = .)

# Dashed lines as short solid pieces (dash_len on, dash_gap off, along each line)
make_dashes <- function(x, dash, gap) {
  out <- list()
  for (g in st_geometry(st_cast(st_geometry(x), "LINESTRING"))) {
    xy <- st_coordinates(g)[, 1:2, drop = FALSE]
    d  <- c(0, cumsum(sqrt(rowSums(diff(xy)^2))))          # distance along the line
    keep <- !duplicated(d); xy <- xy[keep, , drop = FALSE]; d <- d[keep]
    if (length(d) < 2) next
    for (s in seq(0, max(d), by = dash + gap)) {
      e <- min(s + dash, max(d))
      if (e <= s) next
      p <- rbind(c(approx(d, xy[, 1], s)$y, approx(d, xy[, 2], s)$y),
                 xy[d > s & d < e, , drop = FALSE],
                 c(approx(d, xy[, 1], e)$y, approx(d, xy[, 2], e)$y))
      out[[length(out) + 1]] <- st_linestring(p)
    }
  }
  st_sf(geometry = st_sfc(out, crs = st_crs(x)))
}
dashed_pieces <- make_dashes(filter(disputed_lines, type == "dashed"), dash_len, dash_gap)

# Large lakes (equal-area projection, so st_area is in m2)
lakes <- get_ne_layer(dir_ne_lakes, url_ne_lakes)
lakes <- lakes[as.numeric(st_area(lakes)) > 5e9, ]

save_table(bind_rows(rule_check %>% mutate(source = "UN BNDA_simplified rule"),
                     tibble(area = c("Abyei", "Serbia-Kosovo line"), source = "Natural Earth 1:10m"),
                     tibble(area = names(dashed_list), source = "dashed boundary",
                            found = !vapply(dashed_list, is.null, TRUE))),
           "C5_disputed_areas_check")

# CHECK 2: every estimated country appears on the map
message("Countries in estimates but not on the map: ",
        { m <- setdiff(all_tiers$country_code, world$iso3_join); if (length(m)) paste(m, collapse = ", ") else "none" })


# ══════════════════════════════════════════════════════════════════════════════
# 2. SHARED MAP ELEMENTS
# ══════════════════════════════════════════════════════════════════════════════

map_caption <- paste0(
  "The boundaries and names shown and the designations used on this map do not imply the expression of any opinion whatsoever\n",
  "on the part of the authors concerning the legal status of any country, territory, city or area or of its authorities, or concerning\n",
  "the delimitation of its frontiers or boundaries. Dotted and dashed lines on maps represent approximate border lines for which\n",
  "there may not yet be full agreement.\n",
  "Map boundaries: UN Geospatial; Abyei, Kosovo and lakes: Natural Earth. Projection: Equal Earth."
)

map_theme <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13, margin = margin(b = 6)),
        legend.position = "bottom", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), legend.key.size = unit(0.55, "cm"),
        plot.caption = element_text(hjust = 0, size = 7, colour = "grey30"),
        plot.background = element_rect(fill = bg_col, colour = NA), plot.margin = margin(8, 8, 8, 8))
map_coord <- coord_sf(crs = map_crs, default_crs = st_crs(4326),
                      xlim = c(-180, 180), ylim = c(-56, 84), expand = FALSE)

# Optional scale bar (WHO SOP lists it; only exact at one latitude on a world map)
scale_bar <- if (add_scale_bar && requireNamespace("ggspatial", quietly = TRUE)) {
  ggspatial::annotation_scale(location = "br", width_hint = 0.15, style = "bar",
                              height = unit(0.15, "cm"), text_cex = 0.6)
} else NULL

# Layers on top of the fills: Abyei grey, lakes, ordinary boundaries, then the
# disputed lines (WHO: disputed layers on top of all other layers)
top_layers <- function(abyei_fill) {
  list(
    geom_sf(data = abyei, fill = abyei_fill, colour = NA, inherit.aes = FALSE),
    geom_sf(data = lakes, fill = bg_col, colour = border_col, linewidth = border_lwd, inherit.aes = FALSE),
    geom_sf(data = solid_lines, colour = border_col, linewidth = border_lwd, inherit.aes = FALSE),
    geom_sf(data = dashed_pieces, colour = border_col,
            linewidth = border_lwd, lineend = "butt", inherit.aes = FALSE),
    geom_sf(data = filter(disputed_lines, type == "dotted"), colour = border_col,
            linewidth = border_lwd, linetype = lty_dotted, lineend = "round", inherit.aes = FALSE)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. FIGURE 1 — PREVALENCE AND DATA SOURCE
# ══════════════════════════════════════════════════════════════════════════════

prev_colours <- c("<5%" = "#FFFFB2", "5-10%" = "#FED976", "10-20%" = "#FEB24C",
                  "20-30%" = "#FD8D3C", "30-40%" = "#FC4E2A", "40-50%" = "#E31A1C",
                  "50-60%" = "#BD0026", ">60%" = "#800026")
fill_values <- c(prev_colours, setNames(col_no_data_prev, lab_no_data))

map_data <- all_tiers %>%
  mutate(prev_bucket = cut(prevalence, breaks = c(-Inf, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, Inf),
                           labels = names(prev_colours), right = FALSE),
         tier = as.character(tier))

world_map <- world %>%
  left_join(map_data %>% select(country_code, prev_bucket, tier, prevalence),
            by = c("iso3_join" = "country_code")) %>%
  mutate(
    fill_prev = case_when(role == "not_available" | is.na(prev_bucket) ~ lab_no_data,
                          TRUE ~ as.character(prev_bucket)),
    fill_prev = factor(fill_prev, levels = names(fill_values)),
    fill_tier = case_when(role == "not_available" | is.na(tier) ~ lab_no_data, TRUE ~ tier),
    fill_tier = factor(fill_tier, levels = c(tier_levels, lab_no_data)),
    predicted = tier %in% tier_levels[2:3] & role != "not_available"
  )

# WHO: West Bank and Gaza may never be shown in the same colour as Israel
wm <- st_drop_geometry(world_map)
if (any(wm$fill_prev[wm$iso3_join %in% "ISR"] %in% wm$fill_prev[wm$iso3_join %in% "PSE"]))
  warning("Israel and West Bank and Gaza fall in the same prevalence class; WHO maps never show them in the same colour.")

# Diagonal hatching for predicted estimates (not on the grey areas)
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
hatch_spacing <- as.numeric(diff(st_bbox(world_map)[c("xmin", "xmax")])) / 240
hatch <- make_hatch(world_map %>% filter(predicted), spacing = hatch_spacing)

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
                                         gp = grid::gpar(col = col, lwd = max(hatch_lwd * .pt, 0.5))))
}

p_map_prev <- ggplot() +
  geom_sf(data = world_map, aes(fill = fill_prev), colour = NA) +
  geom_sf(data = mutate(hatch, key = "Predicted prevalence"), aes(colour = key),
          linewidth = hatch_lwd, key_glyph = draw_key_hatch) +
  geom_sf(data = mutate(hatch[1, ], key = "Observed prevalence"), aes(colour = key),
          alpha = 0, linewidth = 0, key_glyph = draw_key_hatch) +
  top_layers(col_no_data_prev) +
  scale_fill_manual(values = fill_values, name = "Prevalence", drop = FALSE,
                    guide = guide_legend(nrow = 3, byrow = FALSE, order = 1)) +
  scale_colour_manual(values = c("Observed prevalence" = observed_key_col, "Predicted prevalence" = "grey15"),
                      limits = c("Observed prevalence", "Predicted prevalence"), drop = FALSE, name = NULL,
                      guide = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  scale_bar +
  map_coord + labs(title = "a) Estimated water insecurity prevalence", caption = map_caption) + map_theme +
  theme(legend.box = "horizontal")

p_map_source <- ggplot() +
  geom_sf(data = world_map, aes(fill = fill_tier), colour = NA) +
  top_layers(col_no_data_prev) +
  scale_fill_manual(values = c(tier_col, setNames(col_no_data_prev, lab_no_data)),
                    name = "Data source", drop = FALSE) +
  scale_bar +
  map_coord + labs(title = "b) Data source", caption = map_caption) + map_theme

# Save a map as PNG plus the vector formats in `vector_formats`
save_map <- function(plot, name, width, height, dpi) {
  save_fig(plot, name, width = width, height = height, dpi = dpi)
  base <- sub("\\.png$", "", name)
  if ("pdf" %in% vector_formats)
    save_fig(plot, paste0(base, ".pdf"), width = width, height = height, device = cairo_pdf)
  if ("svg" %in% vector_formats) {
    if (requireNamespace("svglite", quietly = TRUE)) {
      save_fig(plot, paste0(base, ".svg"), width = width, height = height)
    } else message("svglite not installed: SVG copy of ", name, " skipped (install.packages(\"svglite\"))")
  }
}

# Combined figure: disclaimer once, under panel b
save_map((p_map_prev + labs(caption = NULL)) / p_map_source + plot_layout(heights = c(1, 0.8)),
         "C5_fig1_map_prevalence_source.png", width = 11, height = 12.5, dpi = 600)   # Figure 1
dir.create(file.path(dir_figures, "C5_panels"), showWarnings = FALSE)
save_map(p_map_prev,   file.path("C5_panels", "C5_fig1a_map_prevalence.png"), width = 11, height = 7, dpi = 600)
save_map(p_map_source, file.path("C5_panels", "C5_fig1b_map_data_source.png"), width = 11, height = 6, dpi = 600)


# ══════════════════════════════════════════════════════════════════════════════
# 4. FIGURE 3 AND FIGURE S7 — BIVARIATE MAPS
# ══════════════════════════════════════════════════════════════════════════════
# Same thresholds as C2/C4: prevalence and % without access/availability
# <10 | 10-<30 | >=30; water stress (% withdrawal) <10 | 10-<40 | >=40.
# Disputed grey areas and countries without either value: "Data not available".

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
                ~ ifelse(is.na(prev_class) | is.na(.x), lab_no_data, paste0(prev_class, "-", .x)),
                .names = "bivar_{.col}"))

bivar_colours <- c("1-1" = "#e8e8e8", "2-1" = "#e4acac", "3-1" = "#c85a5a",
                   "1-2" = "#cbc9e2", "2-2" = "#b58aa5", "3-2" = "#9e3a6e",
                   "1-3" = "#9e9ac8", "2-3" = "#7a3d6e", "3-3" = "#5b1a3a",
                   setNames(col_no_data_bivar, lab_no_data))

# Caption note on the bivariate maps when hatching is switched on
bivar_note <- if (hatch_bivariate) "Hatched: predicted water insecurity prevalence." else ""

bivar_map <- function(class_col, title, legend_label) {
  d <- world %>%
    left_join(bivar_data %>% select(country_code, all_of(class_col)), by = c("iso3_join" = "country_code")) %>%
    mutate(fill = ifelse(role == "not_available" | is.na(.data[[class_col]]), lab_no_data, .data[[class_col]]))
  # Optional hatching of predicted WI (only where the map shows data)
  hatch_layer <- if (hatch_bivariate) {
    pred_iso3 <- world_map$iso3_join[world_map$predicted]
    geom_sf(data = make_hatch(filter(d, iso3_join %in% pred_iso3, fill != lab_no_data),
                              spacing = hatch_spacing),
            colour = "grey15", linewidth = hatch_lwd, inherit.aes = FALSE)
  } else NULL
  p_map <- ggplot() +
    geom_sf(data = d, aes(fill = fill), colour = NA) +
    hatch_layer +
    top_layers(col_no_data_bivar) +
    scale_fill_manual(values = bivar_colours, guide = "none") +
    scale_bar +
    map_coord + labs(title = title,
                     caption = paste0("Grey: ", tolower(lab_no_data), ". ", bivar_note, "\n", map_caption)) +
    map_theme + theme(plot.title = element_text(face = "bold", size = 12))

  # 3x3 key. All bivariate maps use `legend_plain` (no "Data not available"
  # tile, which made the small key hard to read); grey is explained in the
  # caption instead. `legend` (with the tile) is kept in case it is needed.
  key <- expand.grid(x = 1:3, y = 1:3) %>% mutate(fill = paste0(x, "-", y))
  p_legend_plain <- ggplot(key, aes(x = x, y = y, fill = fill)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = bivar_colours, guide = "none") +
    scale_x_continuous(breaks = c(1, 3), labels = c("Low", "High")) +
    scale_y_continuous(breaks = c(1, 3), labels = c("Low", "High")) +
    labs(x = "WI Prevalence →", y = paste0(legend_label, " →")) +
    theme_void() +
    theme(axis.title.x = element_text(size = 7, margin = margin(t = 2)),
          axis.title.y = element_text(size = 7, angle = 90, margin = margin(r = 2)),
          axis.text = element_text(size = 6), plot.margin = margin(2, 2, 2, 2)) +
    coord_equal(xlim = c(0.5, 3.5), ylim = c(0.5, 3.5), clip = "off")
  p_legend <- p_legend_plain +
    annotate("tile", x = 1, y = -0.4, fill = col_no_data_bivar, colour = "white", linewidth = 0.5) +
    annotate("text", x = 1.65, y = -0.4, label = lab_no_data, hjust = 0, size = 2.1) +
    coord_equal(xlim = c(0.5, 3.5), ylim = c(-1, 3.5), clip = "off")

  list(combined = p_map + inset_element(p_legend_plain, left = 0.02, bottom = 0.12, right = 0.18, top = 0.36),
       map = p_map, legend = p_legend, legend_plain = p_legend_plain)
}

p_stress  <- bivar_map("bivar_stress_class",  "A) Water insecurity prevalence × Water stress", "Water Stress")
p_basic   <- bivar_map("bivar_water_class",   "B) Water insecurity prevalence × Lack of basic drinking water access",
                       "No Basic Drinking (BDWS)")
p_managed <- bivar_map("bivar_managed_class", "C) Water insecurity prevalence × Lack of access to safely managed drinking water",
                       "No SMDWS")
p_avail   <- bivar_map("bivar_avail_class",   "Water insecurity prevalence × Water not available when needed",
                       "Water Not Available")

# Figure 3: maps A-C stacked, each with its 3x3 key; the caption (grey = data
# not available, hatching, disclaimer, projection) once, under map C
with_key <- function(m, caption = FALSE) {
  mm <- if (caption) m$map else m$map + labs(caption = NULL)
  mm + inset_element(m$legend_plain, left = 0.02, bottom = if (caption) 0.12 else 0.04,
                     right = 0.18, top = if (caption) 0.36 else 0.30)
}
save_map(with_key(p_stress) / with_key(p_basic) / with_key(p_managed, caption = TRUE) +
           plot_layout(heights = c(1, 1, 1.12)),
         "C5_fig3_bivariate_maps.png", width = 12, height = 21, dpi = 400)                    # Figure 3
save_map(p_avail$combined, "C5_figS07_bivariate_availability.png", width = 12, height = 7.5, dpi = 600)   # Figure S7

maps <- list(stress = p_stress, basic = p_basic, smdws = p_managed, availability = p_avail)
for (nm in names(maps)) {
  save_map(maps[[nm]]$combined, file.path("C5_panels", paste0("C5_bivariate_", nm, ".png")),
           width = 12, height = 7.5, dpi = 600)
  save_map(maps[[nm]]$legend_plain, file.path("C5_panels", paste0("C5_key_", nm, ".png")),
           width = 2.2, height = 2.8, dpi = 600)
}
save_data(list(map_prevalence = p_map_prev, map_source = p_map_source, bivariate = maps), "C5_map_plots")

stop_timer()
