
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(rgdal)
library(grid)
library(gridExtra)

# Load shapefile data for the desired African country
country_name <- "Rwanda" # Replace with the name of the African country
country_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == country_name)

# Define the coordinates of two points
point1 <- data.frame(lon =29.01502, lat =-2.584748 )
point2 <- data.frame(lon = 29.25639, lat = -1.70278)

# Fetch the elevation data for the country using the 'raster' package
elevation_data <- getData('alt', country = country_name)

# Project the elevation data to match the shapefile
elevation_data_proj <- projectRaster(elevation_data, crs = st_crs(country_sf))

# Convert the raster to a data frame
elevation_df <- as.data.frame(elevation_data_proj, xy = TRUE)
colnames(elevation_df) <- c("lon", "lat", "elevation")

map <- ggplot() +
  geom_sf(data = country_sf) +
  geom_raster(data = elevation_df, aes(x = lon, y = lat, fill = elevation), alpha = 0.5) +
  geom_point(data = point1, aes(x = lon, y = lat), color = "#009E73", size = 3) +
  geom_point(data = point2, aes(x = lon, y = lat), color = "#E7B800", size = 3) +
  scale_fill_gradient(low = "green", high = "brown") +
  theme_bw()












