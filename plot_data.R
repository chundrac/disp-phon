require(ggplot2)
require(ggrepel)

library(rgeos)
library(rgdal) # must be loaded after rgeos! reversed sequence leads to a full crach
require(geosphere)
library(dplyr)

library(ggrastr)
library(tikzDevice)

library(HDInterval)
library(ggridges)

project_data <-  function(
  df,  # a dataframe with Longitude, Latitude, and  data
  base = "~/Documents/110m_physical", # a base map file path
  base_layer = "ne_110m_land", # name of the layer in the base
  projection = "+proj=eqearth +wktext"
) {
  
  #### Settings:
  proj_setting <- paste(projection, "+lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  #### Base map
  world.sp <- readOGR(dsn = base, layer = base_layer, verbose = F)
  
  # shift central/prime meridian towards west - positive values only
  shift <- 180 + 30
  # create "split line" to split polygons
  WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  split.line = SpatialLines(list(Lines(list(Line(cbind(180 - shift,c(-90,90)))), ID = "line")),
                            proj4string = WGS84)
  
  # intersect line with country polygons
  line.gInt <- gIntersection(split.line, world.sp)
  
  # create a very thin polygon (buffer) out of the intersecting "split line"
  bf <- suppressWarnings(gBuffer(line.gInt, byid = TRUE, width = 0.000001))
  
  # split country polygons using intersecting thin polygon (buffer)
  world.split <- gDifference(world.sp, bf, byid = TRUE)
  
  # transform split country polygons in a data table that ggplot can use
  world.sh.tr.df <- fortify(world.split)
  
  # Shift coordinates
  world.sh.tr.df$long.new <- world.sh.tr.df$long + shift
  world.sh.tr.df$long.new <- ifelse(world.sh.tr.df$long.new  > 180,
                                    world.sh.tr.df$long.new - 360, world.sh.tr.df$long.new)
  world.sh.tr.df[,c('X', 'Y')]  <- project(cbind(world.sh.tr.df$long.new, world.sh.tr.df$lat),
                                           proj = proj_setting)
  
  base_map.df <- subset(world.sh.tr.df, lat > -60 & lat < 85)
  
  #### Graticules
  b.box <- as(raster::extent(-180, 180, -90, 90), "SpatialPolygons")
  
  # assign CRS to box
  proj4string(b.box) <- WGS84
  
  # create graticules/grid lines from box
  grid <- gridlines(b.box,
                    easts  = seq(from = -180, to = 180, by = 20),
                    norths = seq(from = -90,  to = 90,  by = 10))
  
  # transform graticules from SpatialLines to a data frame that ggplot can use
  grid.df <- fortify(SpatialLinesDataFrame(sl = grid, data = data.frame(1:length(grid)),
                                           match.ID = FALSE))
  # assign matrix of projected coordinates as two columns in data table
  grid.df[, c("X","Y")]  <- project(cbind(grid.df$long, grid.df$lat),
                                    proj = gsub("lon_0=0", "lon_0=150", proj_setting))
  
  graticules.df <- subset(grid.df, lat > -60 & lat < 85)
  
  # create labels for graticules
  grid.lbl <- labels(grid, side = 1:4)
  
  # transform labels from SpatialPointsDataFrame to a data table that ggplot can use
  grid.lbl.df <- data.frame(grid.lbl@coords, grid.lbl@data)
  
  # add degree sign and clean up
  grid.lbl.df$labels <- ifelse(grepl("S|W", grid.lbl.df$labels),
                               paste0("-", gsub("\\*degree\\*?([EWSN])?", "", grid.lbl.df$labels), "°"),
                               paste0(gsub("\\*degree\\*?([EWSN])?", "", grid.lbl.df$labels), "°")
  )
  
  # grid.lbl.df$labels <- paste0(grid.lbl.df$labels,"°")
  # grid.lbl.df$labels <- gsub("\\*degree\\*?([EWSN])?", "", grid.lbl.df$labels, perl = T)
  
  # adjust coordinates of labels so that they fit inside the globe
  grid.lbl.df$long <- ifelse(grid.lbl.df$coords.x1 %in% c(-180,180),
                             grid.lbl.df$coords.x1*175/180, grid.lbl.df$coords.x1)
  
  grid.lbl.df$lat <-  ifelse(grid.lbl.df$coords.x2 %in% c(-90,90),
                             grid.lbl.df$coords.x2*60/90, grid.lbl.df$coords.x2)
  grid.lbl.df[, c("X","Y")] <-  project(cbind(grid.lbl.df$long, grid.lbl.df$lat),
                                        proj = gsub("lon_0=0", "lon_0=150", proj_setting))
  grid_label.df <- rbind(subset(grid.lbl.df, pos == 2 & coords.x2 > -70 &
                                  !labels %in% c("-60°", "90°")),
                         subset(grid.lbl.df, pos == 1 & !(abs(coords.x1) == 180 & coords.x2 == -90)))
  
  #### Data
  
  if (any(names(df) %in% c('longitude')) | any(names(df) %in% c('latitude'))) {
    df$Longitude <- df$longitude
    df$Latitude <- df$latitude
  }
  
  if (any(names(df) %in% c('lon')) | any(names(df) %in% c('lat'))) {
    df$Longitude <- df$lon
    df$Latitude <- df$lat
  }
  
  df.sp <- SpatialPointsDataFrame(df[,c("Longitude","Latitude")], df,
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  df.sp.tr <- spTransform(df.sp, CRS(gsub("lon_0=0", "lon_0=150", proj_setting)))
  
  data.df <- data.frame(df.sp.tr)
  
  data.df$X <- data.df$Longitude.1
  data.df$Y <- data.df$Latitude.1
  
  return(list(base_map = base_map.df,
              graticules = graticules.df,
              graticule_labels = grid_label.df,
              data = data.df))
  
}

data.df <- read.csv('character_data.tsv',sep='\t',row.names=1)

data.freq <- data.frame(Glottocode=rownames(data.df),freq=rowSums(data.df))

languages <- read.csv('languages.csv')

merged <- merge(data.freq,languages,by.x='Glottocode',by.y='ID')

merged <- merged[,c('Longitude','Latitude','freq')]

df.g <- project_data(df = merged,
                     base = "~/Documents/110m_physical",
                     base_layer = "ne_110m_land",
                     projection = "+proj=eqearth +wktext"
                     # projection = "+proj=robin"
)


tikzDevice::tikz('data_visualized_final.tex')
ggplot() + theme_void() + 
  
  geom_polygon(data = df.g$base_map,
               aes(x = X, y = Y, group = group),
               fill = 'lightgray',
               color = 'black',
               size = .1
  ) + 
  
  geom_path(data = df.g$graticules,
            aes(x = X, y = Y, group = group),
            linetype = 'dotted',
            colour = 'grey',
            size = .25
  ) + 
  geom_point(data = df.g$data, aes(x = X, y = Y, color = freq, alpha = freq)) + scale_colour_gradient(low='blue',high='red')

dev.off()
