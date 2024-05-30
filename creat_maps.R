########################################################################
#Programmer: Wenqiao He
#Last update: May 2024
#Purpose: create maps of Africa and the Democratic Republic of the Congo
########################################################################
install.packages("sf")
install.packages("ggplot2")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)


## Plot map of Afirca and highlight the Democratic Republic of the Congo (DR_Congo) and Tanzania
# Load world data
world <- ne_countries(scale = "medium", returnclass = "sf")

# filter the data for Africa
africa <- world[world$continent == "Africa", ]

africa_featurecla <- africa["featurecla"]

# Filter the data for the Democratic Republic of the Congo (DR_Congo) and Tanzania  
tanzania <- subset(world, name == "Tanzania")
DR_Congo <- subset(world, name == "Dem. Rep. Congo")

# Highlight Tanzania in orange and DR_Congo in lightblue
Africa_map <- ggplot(data = africa_featurecla) +
  geom_sf(fill = "lightgrey", color = "white") +
  geom_sf(data = tanzania, fill = "orange", color = "white") +
  geom_sf(data=DR_Congo, fill = "lightblue", color = "white")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Display the map
print(Africa_map)

##Plot map of the Democratic Republic of the Congo (DR_Congo) and highlight Kinshasa, Sud-Kivu, and Bas-Uele

# Get the country data
congo <- ne_countries(scale = "medium", country = "Democratic Republic of the Congo", returnclass = "sf")

# Get administrative boundaries from https://gadm.org/ (we use gadm41_COD_1.shp for the Congo map)
city_boundaries <- st_read("gadm41_COD_1.shp")

# Filter the data for Kinshasa, Sud-Kivu, and Bas-Uele
Kinshasa <- subset(city_boundaries, city_boundaries$NAME_1 == "Kinshasa")
Bas_Uele <- city_boundaries[city_boundaries$NAME_1 == "Bas-Uele", ]
Sud_Kivu <- city_boundaries[city_boundaries$NAME_1 == "Sud-Kivu", ]

# Base map of DRC with city borders
Congo_map <- ggplot(data = congo) +
  geom_sf(fill = "lightgrey", color = "white")+
  geom_sf(data = city_boundaries, fill = NA, color = "white", size = 0.5)

# Highlight Kinshasa, Sud-Kivu, and Bas-Uele
Congo_map <- Congo_map+
  geom_sf(data = Kinshasa, fill = "red", color = "white")+
  geom_sf(data = Bas_Uele, fill = "lightblue", color = "white")+
  geom_sf(data = Sud_Kivu, fill = "lightgreen", color = "white")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Display the map
print(Congo_map)


