

#' Fit polynomial models
#'
#' This function fits the polynomia functions of the surface trend analysis.
#' @param ds Dataframe providing the date of outbreak and X and Y coordinates
#' @param max.order Integer of highest order polynomial to attempt; defaults to 10
#' @export
estimate_surfacetrend_models = function(ds, max.order=10) {

  order = 1:max.order
  lm.fit = list()
  ntimes = length(ds$time)
  ds.new = ds
  ds.new$XY = ds[,'X'] * ds[,'Y']

  lm.fit[[1]] = lm( time ~ ., data=ds.new)

  for(i in 2:max.order) {

    name.x = paste0("X",i)
    name.y = paste0("Y",i)
    ds.new[,name.x] = (ds.new$X^i)
    ds.new[,name.y] = (ds.new$Y^i)

    lm.fit[[i]] = lm( time ~ ., data=ds.new)
  }
  return(lm.fit)
}








