
#' The front-wave velocity vector map
#'
#' Plots the magnitude and direction of the outbreak front-wave velocity
#' @param r The front-wave velocity summary from the \code{\link{outbreak_velocity}}  function
#' @param shpfile A polynomial shapefile object \code{"SpatialPolygonsDataFrame"} from maptools
#' @param fid The name (string) of a unique numeric ID for each field in the shpfile; must correspond to \code{"id"} and must have one element for every value of \code{"id"}
#' @param dsid A vector of the unique numeric ID for each field in the original dataset; must correspond to \code{"fid"}
#' @param newdf A dataset with new XY grid and time variable
#' @param bin.width A parameter used for contour plot
#' @param unit.num A parameter to define the time variable, default=4 as for a month if the time is measured in week

#' @export
plot_contour = function(r, shpfile, fid, dsid, newdf, bin.width=3, unit.num=4){
  fid.name = paste0(fid)
  shp.fort = fortify(shpfile, region=fid)
  r$velo$id = dsid
  shp.map.temp  = merge(shp.fort, r$velo, by='id', all.x=TRUE, order=FALSE)
  shp.map = shp.map.temp[order(shp.map.temp$order), ]
  g = ggplot()+
    # stat_contour(data = newdf, geom = "polygon", #stat_contour is hard in dealing with the close for the polygon, will form weird polygon
    #              aes( x = X, y = Y, z = time, color = ..level.. ), alpha=0.1, binwidth = 4) +
    # scale_color_gradient(low = "red", high = "yellow") +
    geom_raster(data = newdf, aes( x = X, y = Y, fill = time/unit.num), alpha=0.3)+
    scale_fill_gradientn(colours = heat.colors(20),  labels = paste(c("0", "8", "16"), "months", sep=' '),
                         breaks = c(0, 8, 16))+
    geom_contour(data = newdf, aes( x = X, y = Y, z = time), color = "blue", binwidth =bin.width, size = 0.3) +
    geom_polygon(data=shp.map, aes(x=long, y=lat, group = group), col="black", size=0.1, fill = NA ) +
    #geom_segment(data=vectors.mapping, aes(x = x, y = y, xend = x + dx, yend = y + dy), colour='black', arrow = arrow(angle=30,type='open',length = unit(0.3,"cm"))) +
    theme(panel.background = element_rect(fill = 'white', colour = 'grey90')) +
    xlab("m") + ylab("m") +
    theme(legend.position="right", legend.title=element_blank())

  return(g)
}
