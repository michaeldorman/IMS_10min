###############################################################
# Code to publish interpolated temperature map                #
# based on current observations from the                      #
# Israel Meteorological Service (IMS) website                 #
# By Michael Dorman                                           #
###############################################################

library(shiny)
library(leaflet)
library(magrittr)
library(XML)
library(httr)
library(plyr)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(gstat)
library(automap)

#########################################################################

# UI
ui <- bootstrapPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput("map", width = "100%", height = "100%"),
  absolutePanel(top = 100, left = 10,
    radioButtons(
      "var", 
      label = "Variable", 
      choices = list(
        "Temperature" = "TD",
        "Relative Humidity" = "RH",
        "Rainfall" = "Rain"
        ), 
      selected = "TD"
      ),
    radioButtons(
      "method", 
      label = "Method", 
      choices = list(
        "IDW" = "IDW",
        "OK" = "OK",
        "UK" = "UK"
        ),
      selected = "UK"
      ),
    checkboxInput("markers", "Show markers", FALSE)
  )
)

# SERVER
server <- function(input, output, session) {
  
  # Input data
  url.name = "https://data.gov.il//sites/data.gov.il/files/xml/imslasthour.xml"
  stations = read.csv("metadata10minutesIMS7.csv", stringsAsFactors = FALSE)
  grid = raster("SRTM_Israel_900.tif") %>% aggregate(2)
  names(grid) = "elevation"
  
  # Read XML contents
  url.get = GET(url.name)
  url.content = content(url.get, as = "text")
  dat = 
    xmlParse(url.content) %>% 
    getNodeSet("//Observation") %>% 
    xmlToDataFrame(stringsAsFactors = FALSE)
  
  # Format time
  dat$time_obs = 
    dat$time_obs %>% 
    gsub("T", " ", .) %>% 
    as.POSIXct
  
  # Daylight savings time
  current = Sys.Date()
  dst = seq(as.Date("2015-05-27"), as.Date("2015-10-24"), by = 1)
  if(current %in% dst) { # Daylight saving time in effect
    dat$time_obs = dat$time_obs + 60*60
  }
  
  # Subset most recent observations
  recent = dat$time_obs %>% unique %>% sort(decreasing = TRUE) %>% `[`(1)
  dat = dat[dat$time_obs == recent, ]
  dat$time_obs = dat$time_obs %>% format('%Y-%m-%d %H:%M')
  
  # Convert to point layer
  dat = join(dat, stations, "stn_name")
  dat = dat[!is.na(dat$lon) & !is.na(dat$lat), ]
  coordinates(dat) = ~ lon + lat
  proj4string(dat) = "+proj=longlat +datum=WGS84"
  dat = spTransform(dat, CRS(proj4string(grid)))
  
# Current variable
filtered <- reactive({
    dat$var = dat@data[, input$var] %>% as.numeric
    dat[!is.na(dat$var), ]
  })

# Stations point layer
pnt <- reactive({
  pnt = spTransform(filtered(), CRS("+proj=longlat +datum=WGS84"))
  pnt$label = paste0(pnt$stn_name, ": ", pnt$var)
  pnt
})

# Remove 'Rainfall' button if it is NOT raining
observe({
  if((dat$Rain %>% as.numeric %>% sum(na.rm=TRUE)) < 1) {
  updateRadioButtons(
    session,
    "var", 
    label = "Variable", 
    choices = list(
      "Temperature" = "TD",
      "Relative Humidity" = "RH"
    ), 
    selected = "TD"
  )}
})

# Scale breaks
vals = reactive({
  switch(
      input$var, 
      TD = seq(1, 50, 1),
      RH = seq(0, 100, 5),
      Rain = seq(0, 5, 0.1)
    )
})

# Scale breaks centers
ctrs = reactive({
  vals()[1:(length(vals())-1)] + diff(vals())/2
})

# Scale colours
cols = reactive({
  cols = rainbow(length(ctrs())+5)[1:length(ctrs())]
  if(input$var %in% c("TD", "Rain")) cols = rev(cols)
  # if(input$var == "Rain") cols[1] = NA
  cols
})

# Base map
output$map = renderLeaflet({

  leaflet() %>% 
    addTiles() %>% 
    fitBounds(34.22457, 29.48087, 35.89474, 33.33879)
})

# Add interpolated values
observe({
  
  map = leafletProxy("map", session)
  map %>% 
    clearShapes() %>% 
    clearControls() %>% 
    clearImages()
  
  # Inverse Distance Weighted interpolation
  if(input$method == "IDW") {
    g = gstat(
      formula = var ~ 1, 
      data = filtered())
    z = interpolate(grid, g)
  }
  
  # Ordinary Kriging interpolation
  if(input$method == "OK") {
    v = autofitVariogram(var ~ 1, filtered())
    g = gstat(
      formula = var ~ 1, 
      model = v$var_model, 
      data = filtered())
    z = interpolate(grid, g)
  }
  
  # Universal Kriging interpolation
  if(input$method == "UK") {
    v = autofitVariogram(var ~ elevation, filtered())
    g = gstat(
      formula = var ~ elevation, 
      model = v$var_model, 
      data = filtered())
    z = interpolate(grid, g, xyOnly = FALSE)
  }
  
  # Predicted values range correction
  if(input$var == "RH") {
    z[z<0] = 0
    z[z>100] = 100
  }
  if(input$var == "Rain") {
    z[z<0] = NA
  }
  
  # Reclassify raster
  rcl = matrix(
    c(vals()[1:(length(vals())-1)],
      vals()[2:length(vals())],
      (vals()[1:(length(vals())-1)] + vals()[2:length(vals())]) / 2),
    ncol = 3
  )
  z = reclassify(z, rcl)
  
  # Mask
  z = mask(z, grid)
  
  # Polygons
  pol = rasterToPolygons(z, n = 4, dissolve = TRUE)
  pol = spTransform(pol, CRS("+proj=longlat +datum=WGS84"))
  
  # Contour lines
  line = rasterToContour(z, levels = vals())
  line = spTransform(line, CRS("+proj=longlat +datum=WGS84"))
  line = disaggregate(line)
  
  # Current palette
  current_vals = z %>% values %>% unique %>% sort
  current_cols = cols()[match(current_vals, ctrs())]
  pal = colorNumeric(palette = cols(), domain = vals(), na.color = NA)
  
    map %>% 
      addRasterImage(
        z, 
        colors = pal, 
        opacity = 0.75
      ) %>% 
      addPolygons(data = pol, 
        stroke = FALSE, 
        fillOpacity = 0, 
        popup = ~as.character(layer)
      ) %>% 
      addPolylines(data = line, weight = 1, color = "black") %>% 
      addLegend(
        position = "topright", 
        colors = current_cols %>% substr(1, 7),
        labels = current_vals,
        opacity = 1,
        title = recent %>% format('%d-%b %H:%M')
        )
})

# Show/Hide markers
observe({
  
  map = leafletProxy("map", session)
  map %>% clearMarkers()
  
  if(input$markers)
    map %>% 
      addMarkers(
        data = pnt(), 
        popup = ~label
    )
})
  
}

shinyApp(server = server, ui = ui)