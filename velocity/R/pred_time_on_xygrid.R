#' Fit polynomial models
#'
#' This function fits the polynomia functions and predict the delay time for a new xy grid
#' used to generate the dataset for the contour map
#' @param ds_new=ds Dataframe providing the date of outbreak and X and Y coordinates
#' @param max.order Integer of highest order polynomial to attempt; defaults to 10
#' @param shpfile A polynomial shapefile object \code{"SpatialPolygonsDataFrame"} from maptools
#' @param r The front-wave velocity summary from the \code{\link{outbreak_velocity}}  function
#' @param bestorder The order of the best performance model
#' @export

pred_time_on_xygrid = function(ds_new, r, bestorder,shpfile, max.order=10) {
  #order = 1:bestorder
  new.df = expand.grid(X = seq((min(r$ds$X)*0.8), (max(r$ds$X)*1.2), length.out = max(r$ds$time, na.rm = T)),
                       Y = seq((min(r$ds$Y)*0.8), (max(r$ds$Y)*1.2), length.out = max(r$ds$time, na.rm = T)))

  new.df$XY = new.df[,'X'] * new.df[,'Y']

  for(i in 2:bestorder) {

    name.x = paste0("X",i)
    name.y = paste0("Y",i)
    new.df[,name.x] = (new.df$X^i)
    new.df[,name.y] = (new.df$Y^i)
  }
  trend.fit = estimate_surfacetrend_models(ds_new, max.order)
  new.df$time = predict(trend.fit[[bestorder]], new.df) # bestorder set to 6
  new.df = new.df[,c("X", "Y", "time")]
  new.df = clip_xygrid(new.df, shpfile)
  return(new.df)
}
