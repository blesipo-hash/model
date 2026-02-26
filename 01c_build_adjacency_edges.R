# =============================================================================
# R/01c_build_adjacency_edges.R
# Purpose:
#   Build country adjacency edges (iso3_from, iso3_to) for ICAR/CAR spatial priors.
#   Uses Natural Earth admin-0 boundaries (reproducible global dataset).
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(rnaturalearth)
})

PANEL_IN  <- "data/processed/model_input/sac_panel_merged.csv"
OUT_EDGES <- "data/input/adjacency_edges.csv"
QC_DIR    <- "data/processed/model_input/qc_adj_edges"

dir.create(dirname(OUT_EDGES), recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load iso3 list used by your model (limits edges to only countries you actually model)
if (!file.exists(PANEL_IN)) stop("Missing panel: ", PANEL_IN, call. = FALSE)

panel_iso3 <- read_csv(PANEL_IN, show_col_types = FALSE) %>%
  transmute(iso3 = toupper(str_squish(as.character(iso3)))) %>%
  distinct() %>%
  filter(!is.na(iso3), nchar(iso3) == 3)

if (nrow(panel_iso3) == 0) stop("No iso3 codes found in panel.", call. = FALSE)

# --- Natural Earth admin-0 boundaries (scale = 50 is a good default)
world <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") %>%
  st_make_valid() %>%
  transmute(
    iso3 = toupper(str_squish(iso_a3)),
    name = name_long,
    geom = geometry
  ) %>%
  # Natural Earth uses "-99" for some entries without iso codes; drop them
  filter(!is.na(iso3), iso3 != "-99", nchar(iso3) == 3)

# keep only countries in your model panel
world_in <- world %>% semi_join(panel_iso3, by = "iso3")

# QC: ensure every modeled iso3 has a polygon
missing_shapes <- panel_iso3 %>% anti_join(st_drop_geometry(world_in), by = "iso3")
if (nrow(missing_shapes) > 0) {
  write_csv(missing_shapes, file.path(QC_DIR, "iso3_missing_in_naturalearth.csv"))
  message(
    "NOTE: Dropping iso3 missing in Natural Earth from adjacency build (not fatal). QC: ",
    file.path(QC_DIR, "iso3_missing_in_naturalearth.csv")
  )
}

# Keep only iso3 that have shapes (so adjacency is well-defined)
panel_iso3_ok <- panel_iso3 %>% anti_join(missing_shapes, by = "iso3")

# Re-filter world_in accordingly
world_in <- world %>% semi_join(panel_iso3_ok, by = "iso3")

# --- Define adjacency: queen contiguity via st_touches (shares boundary OR point)
# Note: islands will have degree 0 (expected)
touch_list <- st_touches(world_in)

edges_raw <- tibble(i = seq_len(nrow(world_in))) %>%
  mutate(js = touch_list) %>%
  unnest(js) %>%
  transmute(
    iso3_from = world_in$iso3[i],
    iso3_to   = world_in$iso3[js]
  ) %>%
  filter(iso3_from != iso3_to) %>%
  distinct()

# Make undirected unique edges: store canonical ordering (min,max) to dedupe
edges_undirected <- edges_raw %>%
  transmute(
    iso3_from = pmin(iso3_from, iso3_to),
    iso3_to   = pmax(iso3_from, iso3_to)
  ) %>%
  distinct() %>%
  arrange(iso3_from, iso3_to)

# QC: degree table + isolates
deg <- bind_rows(
  edges_undirected %>% transmute(iso3 = iso3_from),
  edges_undirected %>% transmute(iso3 = iso3_to)
) %>%
  count(iso3, name = "degree") %>%
  right_join(panel_iso3, by = "iso3") %>%
  mutate(degree = replace_na(degree, 0L)) %>%
  arrange(degree, iso3)

write_csv(deg, file.path(QC_DIR, "adj_degree_by_country.csv"))
write_csv(deg %>% filter(degree == 0), file.path(QC_DIR, "adj_isolates.csv"))

# Write edges for Stan-prep script
write_csv(edges_undirected, OUT_EDGES)

message("Wrote edges: ", OUT_EDGES)
message("QC dir: ", QC_DIR)
message("Countries in panel: ", nrow(panel_iso3))
message("Countries with shapes: ", nrow(world_in))
message("Edges (undirected): ", nrow(edges_undirected))
message("Isolates (degree==0): ", sum(deg$degree == 0))
message("DONE.")


