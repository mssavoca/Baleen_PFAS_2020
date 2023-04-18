--------------------------------------------------------------------------------
  # Baleen-PFAS map figure
  # Matt Savoca
  # Started on: 4/3/23
--------------------------------------------------------------------------------
  
#install.packages("ggOceanMapsData", repos = c("https://mikkovihtakari.github.io/drat", "https://cloud.r-project.org"))
library(ggOceanMaps)
library(ggOceanMapsData)
library(maps)
library(maptools)
library(usmap)
library(mapdata)
library(ggspatial)

source("Util.R")


states <- map_data('state')


WC_map <- ggplot() + 
  geom_polygon(data= filter(states, region %in% c("california", "oregon")),
               aes(x=long, y=lat, group=group),
               color="black", fill="burlywood", alpha = 0.7) +
  geom_point(data = MapData, aes(x = long, y = lat, color = abbr_binom(SciName)), size = 4) +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  labs(color = "Species",
       x = "Longitude",
       y = "Latitude") +
  xlim(-127.5, -114) +
  theme_minimal(base_size = 16)
WC_map

ggsave("WC_map.pdf")


GOM_map <- ggplot() + 
  geom_polygon(data= filter(states, region %in% c("massachusetts", "maine", "connecticut",
                                                  "rhode island", "new hampshire")),
               aes(x=long, y=lat, group=group),
               color="black", fill="burlywood2", alpha = 0.7) +
  geom_point(data = MapData, aes(x = long, y = lat, color = abbr_binom(SciName)), size = 3.5) +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  labs(color = "Species",
       x = "Longitude",
       y = "Latitude") +
  xlim(-73.75, -66.5) +
  ylim(41, 47.5) +
  theme_minimal(base_size = 16)
GOM_map 

ggsave("GOM_map.pdf")  




interest_region <- c(left=-126,bottom=32,right=-116.5,top=44)
x <- c(interest_region["left"], interest_region["left"], 
       interest_region["right"], interest_region["right"])
y <- c(interest_region["bottom"], interest_region["top"], 
       interest_region["top"], interest_region["bottom"])
df_WC <- data.frame(x, y)

interest_region <- c(left=-74,bottom=40,right=-66,top=46)
x <- c(interest_region["left"], interest_region["left"], 
       interest_region["right"], interest_region["right"])
y <- c(interest_region["bottom"], interest_region["top"], 
       interest_region["top"], interest_region["bottom"])
df_GOM <- data.frame(x, y)


PFAS_baleen_zoomout <- ggOceanMaps::basemap(limits = c(-130, -58, 27, 53), 
                                 lon.interval = 10, lat.interval = 10) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1, "cm"),
                         width = unit(1, "cm")) +
  annotation_scale(location = "br") +
  geom_polygon(aes(x=x, y=y), data=df_GOM, fill = NA, 
               color = "black", size = 0.75) +
  geom_polygon(aes(x=x, y=y), data=df_WC, fill = NA, 
               color = "black", size = 0.75) +
  annotate("text", x = -69.25, y = 43, label = "B", fontface = "bold",
           size = 6, color = "black") +
  annotate("text", x = -123, y = 36, label = "A", fontface = "bold",
           size = 6, color = "black") +
  labs(x = NULL, y = NULL) +
  theme_bw()

PFAS_baleen_zoomout


ggsave("PFAS_baleen_insetmap.pdf")  







# Junk code below here ----

# Inset map for OPC MP RFP
interest_region <- c(left=-122.25,bottom=36.5,right=-121.75,top=37)
x <- c(interest_region["left"], interest_region["left"], 
       interest_region["right"], interest_region["right"])
y <- c(interest_region["bottom"], interest_region["top"], 
       interest_region["top"], interest_region["bottom"])
df_CA <- data.frame(x, y)


MP_OPC_RFP_zoomout <- ggOceanMaps::basemap(limits = c(-124, -119.5, 33.5, 39), 
                                            lon.interval = 2, lat.interval = 2) +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(2.5, "cm"),
                         width = unit(2, "cm"), text_cex = 1.5) +
  annotation_scale(location = "br", line_width = 2,
                   height = unit(.75, "cm"), text_cex = 1.5) +
  geom_polygon(aes(x=x, y=y), data=df_CA, fill = NA, 
               color = "black", size = 1.75) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 22)
MP_OPC_RFP_zoomout

ggsave("MP_OPC_RFP_zoomout.pdf")  

PFAS_baleen_zoomout


ggsave("PFAS_baleen_insetmap.pdf")  


library(sf)
library(ggmap)
library(rgeos)
library(mapproj)
library(ggnewscale)
library(terra)
library(rnaturalearth)
library(prettymapr)
library(ggpubr)



PFASMapData <- readxl::read_xlsx("Baleen-PFAS master sample list.xlsx", sheet = 2) 
PFAS_transform <- usmap_transform(PFASMapData)





states_sf <- st_as_sf(states, coords = c("long", "lat"), crs = 4326)

states_ggocean <- st_transform(states_sf, st_crs(basemap()))

states_sf_lines <- st_cast(states_sf, "LINESTRING")




GOM_map <- basemap(limits = c(-71.75, -67, 39.999, 45), 
        land.col = "burlywood",
        # bathymetry = TRUE, 
        # bathy.style = "contour_grey",
        rotate = TRUE) + 
  geom_point(data = MapData, aes(x = long, y = lat, color = abbr_binom(SciName)), 
             size = 2.75) +
  #geom_polygon(data = map_data, aes(x=long, y=lat, group = group), fill = NA, color = "black") +
  #geom_sf(data = map_data_sf, color = "black") +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  annotation_scale(location = "br") + 
  annotate("text", x = -71.4, y = 42.35,
           label = "Boston",
           size = 4) +
  #annotation_north_arrow(location = "tr", which_north = "true") +
  labs(color = "Species",
       x = "Longitude",
       y = "Latitude") +
  theme(panel.grid = element_blank()) +
  coord_sf(xlim = c(-71.75, -67), ylim = c(39.999, 45))
GOM_map

ggsave("GOM_map.pdf")  



# Create a data frame with the coordinates for the scale bar
scale_df <- data.frame(x = c(-72, -71),
                       y = c(41.25, 41.25),
                       label = c("0", "100 km"))

# Add the scale bar and label to the plot
GOM_map + 
  geom_segment(data = scale_df, aes(x = x, y = y, xend = x + 100000, yend = y),
               arrow = arrow(length = unit(0.02, "npc"), ends = "last", type = "closed"),
               lineend = "round", color = "black", linewidth = 1.5) +
  geom_text(data = scale_df, aes(x = x + 50000, y = y - 0.5, label = label),
            color = "black", size = 4, hjust = 0.5)

params <- scalebarparams(plotunit = "km", extents = c(-74, -66.5, 41, 47.5))

# Plot the map with scale bar
plotscalebar(x = -72, y = 41.5, ht = 0.5, params = params, style = "bar",
             adj = c(0, 0), bar.cols = "black", lwd = 1, linecol = "black")




WC_map <- basemap(limits = c(-125.5, -115.5, 30, 45), 
                  land.col = "burlywood") +
                  # bathymetry = TRUE, 
                  # bathy.style = "contour_grey") + 
  #geom_spatial_polygon(data = states, aes(long, lat),color = "black") +
  geom_point(data = MapData, aes(x = long, y = lat, color = abbr_binom(SciName)), size = 4) +
  scale_colour_manual(values = pal,
                      guide = guide_legend(label.theme = element_text(face = "italic", size = 10))) +
  annotation_scale(location = "bl") + 
  annotate("text", x = c(-122.5,-116.3), y = c(38,32.7),
           label = c("San\nFrancisco","San\nDeigo"),
           size = 3.5) +
  #annotation_north_arrow(location = "tr", which_north = "true") +
  labs(color = "Species",
       x = "Longitude",
       y = "Latitude") +
  theme(panel.grid = element_blank()) 

WC_map

ggsave("WC_map.pdf")











# Create the ggplot object with longitude and latitude as the axes:
ggplot() +
  
  # Add a plain blue background to signify the ocean:
  geom_rect(fill = "blue", xmin = -180, xmax = 180, ymin = -90, ymax = 90, alpha = 0.1) +
  
  # Add the North America map as the background:
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), fill = "gray90", color = "black") +
  
  # Add points for whale strandings and adjust their size and color:
  geom_point(data = MapData, aes(x = Stranding_long, y = Stranding_lat, color = SciName), size = 3) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  
  # Add a title and labels:
  ggtitle("Whale Strandings by Scientific Name") +
  xlab("Longitude") +
  ylab("Latitude")
