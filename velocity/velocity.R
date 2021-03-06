# Mengru Yuan

library(R.utils)
library(dplyr)
library(ggmap)
library(sf)
library(rgdal)
library(sp)
library(maptools)
library(spdep) 
library(raster)
library(rgeos)
sourceDirectory("R/") 

shpfile = readOGR(dsn = "data", layer = "municipalities_2017_3114", encoding = "UTF-8");

#Load the data
chika = read.csv("data/chika.csv", header = TRUE, stringsAsFactors = FALSE) %>%  mutate(cod_unique_municipio = stringr::str_pad(as.character(cod_unique_municipio), 5, pad = "0"))
zika = read.csv("data/zika.csv", header = TRUE, stringsAsFactors = FALSE) %>%   mutate(cod_unique_municipio = stringr::str_pad(as.character(cod_unique_municipio), 5, pad = "0"))

#------------------------------- FOR ZIKA ---------------------------------------
disease = "zika"; df = zika; order = 5 

#------------------------------- FOR CHIKA --------------------------------------
disease = "chika"; df = chika; order = 5 

#--------------------------APPLY DATASET TO VELOCITY PACKAGE --------------------
ds_new = df %>%dplyr:: select(xcoord,ycoord,week_1st) %>% rename(time=week_1st, X = xcoord, Y= ycoord)

##Use outbreak velocity function to get velocities, predicted data, order and lm model fit
out = outbreak_velocity(ds_new, max.order = order, manual.order = T )
out_lower = outbreak_velocity(ds = ds_new, max.order = order, manual.order = T, value = "lower_ci")
out_upper = outbreak_velocity(ds = ds_new, max.order = order, manual.order = T, value = "upper_ci")

##Make a plot
velocity_vector_map(r=out, shpfile=shpfile,dsid=df$cod_unique_municipio,fid <- "MPIO_CCDGO" )


## generate raster file
getContourPlot(r = out, scale=200, contour.crs = "+proj=tmerc +lat_0=4.596200416666666 +lon_0=-80.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
               save.name= paste0(disease,"_order_", order,"_contour200_", "_",format(Sys.Date(), "%Y%m%d")))

#Generate new dataset for contour plot
df_xygrid = pred_time_on_xygrid(r= out, ds_new = ds_new, bestorder = order, shpfile = shpfile); summary(df_xygrid)

##Make contour plot
plot_contour(r=out, shpfile=shpfile,dsid=df$cod_unique_municipio,fid <- "MPIO_CCDGO", newdf = df_xygrid, bin.width=1, unit.num=1)

##Output data
output <- out; output$velo$id = df$cod_unique_municipio; out$velo$Municipio = df$Municipio; output$velo$week_1st = df$week_1st
output$velo$degree2=ifelse(out$velo$direction.degrees <= 90, (90-out$velo$direction.degrees),(450-out$velo$direction.degrees))
output$velo$X= out$ds$X
output$velo$Y= out$ds$Y
write.csv(output$velo, file = paste0("results/", disease,"_order_", order, "_",format(Sys.Date(), "%Y%m%d"),".csv"))

##Output data for lower ci
output <- out_lower; output$velo$id = df$cod_unique_municipio; out$velo$Municipio = df$Municipio; output$velo$week_1st = df$week_1st
out = out_lower
output$velo$degree2=ifelse(out$velo$direction.degrees <= 90, (90-out$velo$direction.degrees),(450-out$velo$direction.degrees))
output$velo$X= out$ds$X
output$velo$Y= out$ds$Y
write.csv(output$velo, file = paste0("results/", disease,"_order_", order,"_lowerci_", "_",format(Sys.Date(), "%Y%m%d"),".csv"))

##Output data for upper ci
output <- out_upper; output$velo$id = df$cod_unique_municipio; out$velo$Municipio = df$Municipio; output$velo$week_1st = df$week_1st
out = out_upper
output$velo$degree2=ifelse(out$velo$direction.degrees <= 90, (90-out$velo$direction.degrees),(450-out$velo$direction.degrees))
output$velo$X= out$ds$X
output$velo$Y= out$ds$Y
write.csv(output$velo, file = paste0("results/", disease,"_order_", order,"_upperci_", "_",format(Sys.Date(), "%Y%m%d"),".csv"))

