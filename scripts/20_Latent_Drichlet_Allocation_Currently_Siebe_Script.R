Latent drichlet allocation structure
```{r data & libraries, include=FALSE}

library(ggplot2) # Plotting
library(dplyr) # Data wrangling
library(fistools) # Google drive function
library(jsonlite) # JSON file writer
library(reticulate) # Python connection
py_require("tqdm") # Python package for progress bar

# See the phyton map structure: https://github.com/FRC-Community-Ecology/PhD_Siebe_Indestege/tree/main/LoRaWAN_performance_testing/Tools/Py1812-main 


# Download the digital terrain model in case it is missing (too large for Github)
fistools::download_gdrive_if_missing(
  "13iuVeMolWuVvhHEMzePaW7Z55ok02Bu1",
  "LoRaWAN_performance_testing/Tools/Py1812-main/data/dtm_height_model.tif")

```

# Create hexagonal grid over national park
We use hexagonal grids to limit corner effects
The LoRa connectivity will be estimated for each single grid cell, hence decreasing cell size increases computation time.

```{r grid creation}

# Make a grid of xxx by xxx meters (hexagonal to limit corner effects)
grid <- st_make_grid(nphk, cellsize = c(30, 30), square = F)

# Visualization
ggplot() +
  geom_sf(data = grid) +
  geom_sf(data = nphk, col = "blue", fill = NA)

# Transform the grid into an sf object
grid_sf <- st_sf(geometry = grid)

# Give each cell of the grid an identicator
grid_sf$id <- 1:nrow(grid_sf)

# Save grid to python data source
st_write(grid_sf, "LoRaWAN_performance_testing/Tools/Py1812-main/data/NPHK_hex_grid.shp", append = F)

```

# Function to transform google maps coordinates directly in the correct format (copy paste the google maps coordinates)

```{r maps function}

google_maps_coords_31370 <- function(lon, lat) {
  return(
    st_coordinates(
      st_transform(
        st_sfc(
          st_point(c(lat, lon)),
          crs = 4326),
        31370)))
}

google_maps_coords_31370(50.979165, 5.633121) # MH
google_maps_coords_31370(50.977806, 5.574365) # Bos As
google_maps_coords_31370(50.997595, 5.667042) # Rode Terril
google_maps_coords_31370(50.947312, 5.654542) # Koeienweide Opgrimbie
google_maps_coords_31370(50.923244, 5.637279) # Vallei van de Ziepbeek

```

# Make a multitude of gateway locations
In order to get insight in the possible locations for gateway placement in the National Park, we will start with simulating the results for a set of systematically distributed locations throughout the study area.

```{r gw locations}

# Make a grid and keep the corner locations as possible gateway locations
nphk_grid <- st_make_grid(nphk, cellsize = c(750,750), square = T, what = "corners") %>%
  st_intersection(nphk) %>%
  st_as_sf() %>%
  mutate(
    id = row_number(),
    x_coord = st_coordinates(.)[,1],
    y_coord = st_coordinates(.)[,2],
    height = 14.0
  )

# Visualization
ggplot() +
  geom_sf(data = nphk_grid) +
  geom_sf(data = nphk, fill = NA)

# Transform the sf object to a list of lists where each sublist represents one gateway with its associated id, x coordinate, y coordinate, and height.
gw_list <- nphk_grid %>%
  st_drop_geometry() %>%
  select(id, x = x_coord, y = y_coord, height) %>%
  purrr::transpose()

```

# Specify all necessary input for ITU R model
Here we specify all the necessary input that the python script needs to fuel the ITU R model. The gateways_list is only necessary when we want to specify specific gateways.
All others are dependent on the hardware of sensors and gateways, and specific conditions.

```{r gateway locations}

# 1. Define a list of Gateways
gateways_list <- list(
  list(id = "GW_MH", x = 238788.9, y = 186192.2, height = 14.0),
  list(id = "GW_rode_terril", x = 241134.8, y = 188283.2, height = 14.0),
  list(id = "GW_vallei", x = 239187.2, y = 179977.2, height = 14.0),
  list(id = "GW_As", x = 234665.6, y = 185972.3, height = 14.0),
  list(id = "GW_kikbeek", x = 240354.4, y = 182661.8, height = 14.0)
)

# 2. Add them to the global settings
simulation_params <- list(
  # Gateway locations and heights
  gateways = gateways_list,
  
  # Sensor specifications
  sensor = list(height = 1, antenna_gain = 2.0, tx_power = 14.0),
  
  # Hardware specifications of gateway
  hardware = list(antenna_gain_gw = 3.0, cable_loss_gw = 2.0, sensitivity_sf7 = -120, sensitivity_sf12 = -130), # Based on Matni et al. 2020
  
  # Other model specifications
  model = list(freq_GHz = 0.868, time_percentage = 50, loc_percentage = 95, distance_interval = 20),
  
  # Output specifications
  output = list(output_name = "scenario_80perc_coverage")
)

# 3. Write the JSON
write_json(simulation_params, "LoRaWAN_performance_testing/Tools/Py1812-main/data/sim_settings.json", auto_unbox = TRUE, pretty = TRUE)

```

# Python script
Now lets run the python script using all of the data prepared here

```{r python ITU-R model, include=FALSE}

# Keep environment from before python
pre_python_env <- ls()

# Define the base folder path clearly
base_path <- "LoRaWAN_performance_testing/Tools/Py1812-main"

# Tell Python to look into the 'src' folder inside that path
py_run_string(paste0("import sys; sys.path.append('", base_path, "/src')"))

# Now source the python file using the full path
source_python(file.path(base_path, "heatmap_builder_R.py"))

# Clean up the environment
rm(list=setdiff(ls(), pre_python_env))
