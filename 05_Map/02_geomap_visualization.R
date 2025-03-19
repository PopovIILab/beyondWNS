# Set working directory
main_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(main_dir)

# Install or call libraries
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(dplyr, raster, sf, stars, tmap, terra)

options(scipen = 999)

# Read metadata.tsv
metadata <- read.delim("map_metadata.tsv", header = TRUE, sep = "\t")

# Extract country dynamically
country_usa <- metadata %>%
  filter(grepl("Pseudogymnoascus destructans", AN_OrganismName)) %>%
  pull(Country)

country_antarctica <- metadata %>%
  filter(grepl(
    paste(
      "Thelebolus microsporus",
      "Antarctomyces pellizariae",
      "Antarctomyces psychrotrophicus",
      sep = "|"
    ),
    AN_OrganismName
  )) %>%
  pull(Country)

ne = 'natural_earth_vector/packages/natural_earth_vector.gpkg'
countries = read_sf(ne, 'ne_110m_admin_0_countries')
ocean = read_sf(ne, 'ne_110m_ocean')

lyr = lst(ocean, countries)

lyrp = lapply(lyr, st_transform, crs = "+proj=mill") # Цилиндрическая проекция Миллера

# Load the world climate data grom geodata
temp = geodata::worldclim_global(var = "tavg", res = 10, path = 'data') |>
  st_as_stars() |>
  rename(tavg = 1)

# Calculate average yearly temperature
average_temp = st_apply(temp, c(1, 2), mean, na.rm = TRUE)

# Transform `stars` → `terra`
temp_terra = as(average_temp, "Raster") |> rast()

# Expand coverage to the entire globe
ext(temp_terra) = c(-180, 180, -90, 90)

# Manual assignment of Miller projection
mill_proj = "+proj=mill +lon_0=0 +datum=WGS84 +units=m +no_defs"

# Recalculation into Miller projection
tempp_terra = project(temp_terra, mill_proj)

# Return to stars
tempp = st_as_stars(tempp_terra)

# Filter countries to leave only needed
lyrp$usa <- lyrp$countries %>%
  filter(NAME %in% country_usa)

#lyrp$canada <- lyrp$countries %>%
#filter(NAME %in% c("Canada"))

lyrp$antarctica <- lyrp$countries %>%
  filter(NAME %in% country_antarctica[1])

# Plotting the map
tm_shape(tempp) +
  tm_raster(
    'layer',
    title = 'Average\nYearly\nTemperature, °C',
    colorNA = 'grey',
    textNA = 'NA',
    legend.format = list(text.separator = '–'),
    n = 11,
    midpoint = 0,
    style = 'pretty',
    legend.reverse = TRUE,
    palette = '-RdBu'
  ) +
  tm_shape(lyrp$ocean) +
  tm_fill(fill = 'azure') +
  tm_borders(col = 'steelblue') +
  
  tm_shape(lyrp$usa) +
  tm_dots(
    col = 'brown1',
    size = 1.5,
    border.col = 'black',
    border.lwd = 1.7
  ) +
  tm_text(
    text = "GU_A3",
    size = 1,
    fontface = "bold",
    col = "black",
    #bgcol = "white",
    #bgcol_alpha = 1,
    just = "left",
    xmod = 2,
    shadow = TRUE
  ) +
  
  #tm_shape(lyrp$canada) +
  #tm_dots(
  #col = 'brown1',
  #size = 1.5,
  #border.col = 'black',
  #border.lwd = 1.7
  #) +
  #tm_text(
  #text = "NAME",
  #size = 1,
  #fontface = "bold",
  #col = "black",
  #bgcol = "white",
  #bgcol_alpha = 1,
  #just = "right",
  #xmod = -2.5,
  #shadow = TRUE
  #) +
  
  # Добавляем страны
  tm_shape(lyrp$antarctica) +
  tm_dots(
    col = 'aquamarine2',
    size = 1.5,
    border.col = 'black',
    border.lwd = 1.7
  ) +
  tm_text(
    text = "NAME",
    size = 1,
    fontface = "bold",
    col = "white",
    #bgcol = "white",
    #bgcol_alpha = 1,
    just = "left",
    xmod = 3,
    shadow = TRUE
  ) +
  
  tm_layout(
    legend.position = c('right', 'top'),
    legend.frame = TRUE,
    legend.frame.lwd = 0.2,
    legend.bg.alpha = 0.75,
    legend.bg.color = 'white',
    component.autoscale = FALSE
  ) +
  tm_graticules(
    x = seq(-150, 150, by = 30),
    y = seq(-60, 60, by = 30),
    lwd = 0.2,
    col = "black"
  )

#Create directory to save the image
sub_dir = "imgs"

dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE)

# Save the image
tmap_save(
  filename = "imgs/map.png",
  width = 12,
  height = 8,
  dpi = 600
)

# Special thanks to these guides:
# 1. https://stackoverflow.com/questions/4216753/folder-management-with-r-check-existence-of-directory-and-create-it-if-it-does
# 2. https://tsamsonov.github.io/r-geo-course/11-ThematicMaps.html