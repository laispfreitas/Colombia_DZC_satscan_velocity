

#' The front-wave velocity vector map
#'
#' Plots the magnitude and direction of the outbreak front-wave velocity
#' @param r The front-wave velocity summary from the \code{\link{outbreak_velocity}}  function
#' @param shpfile A polynomial shapefile object \code{"SpatialPolygonsDataFrame"} from maptools
#' @param fid The name (string) of a unique numeric ID for each field in the shpfile; must correspond to \code{"id"} and must have one element for every value of \code{"id"}
#' @param dsid A vector of the unique numeric ID for each field in the original dataset; must correspond to \code{"fid"}
#' @export
velocity_vector_map = function(r, shpfile, fid, dsid) {
  fid.name = paste0(fid)
  shp.fort = fortify(shpfile, region=fid)
  r$velo$id = dsid
  shp.map.temp  = merge(shp.fort, r$velo, by='id', all.x=TRUE, order=FALSE)
  shp.map = shp.map.temp[order(shp.map.temp$order), ]

  vectors.mapping = data.frame(x=r$ds$X, y=r$ds$Y, velocity=r$velo$magnitude.velocity, dx = r$velo$dx, dy = r$velo$dy, dir=r$velo$direction.degrees, id=r$velo$id)

  g = ggplot()  +
    geom_polygon(data=shp.map, aes(x=long, y=lat, group = group, fill=magnitude.velocity), col="white", size=0.1) +
    scale_fill_gradient(low='#F6D2BF', high='#D53534', na.value = "grey70") +
    geom_segment(data=vectors.mapping, aes(x = x, y = y, xend = x + dx, yend = y + dy), colour='black', arrow = arrow(angle=20,type='open',length = unit(0.2,"cm"))) +
    theme(panel.background = element_rect(fill = 'white', colour = 'grey90')) +
    xlab("m") + ylab("m") +
    theme(legend.position="right", legend.title=element_blank())

  return(g)
}

