# clip the xy grid based on the provided shpfile

#' @param dnew The Dataframe providing the predicted date of outbreak and X and Y coordinates and the covariates
#' @param shpfile A polynomial shapefile object \code{"SpatialPolygonsDataFrame"} from maptools
#' @export
#dnew= new.df
clip_xygrid = function(dnew, shpfile) {
  coordinates(dnew)=~X+Y  #to check:spplot(dnew); plot(shpfile); points(dnew)
  proj4string(dnew) = proj4string(shpfile)
  new.df.subset = dnew[shpfile,]
  #to check: plot(shpfile) ; points(new.df.subset)
  df = as.data.frame(new.df.subset@coords)
  df$time = new.df.subset@data[,1]
  return(df)
}
