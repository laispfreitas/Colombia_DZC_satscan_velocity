

#' The front-wave velocity
#'
#' Calculates the front-wave velocity of an outbreak
#' @param ds Dataframe providing the date of outbreak and X and Y coordinates
#' @param max.order Integer of highest order polynomial to attempt; defaults to 10
#' @param measure The measure of model fit, defaults to R-squared: 'r2', can be AIC: 'aic' or BIC: 'bic'
#' @param manual.order If TRUE, the max.order integer is the polynomial model used, regardless of best fit measures; defaults to FALSE
#' @export
outbreak_velocity = function(ds, max.order=10, measure='r2', manual.order=FALSE, value = "point_estimate") {
  if(manual.order == FALSE) {
    trend.fit = estimate_surfacetrend_models(ds, max.order)
    trend.fit.measure = surfacetrend_model_fit(trend.fit)
    trend.bestfit = get_best_model(trend.fit.measure, measure)
    trend.bestfit.pdv = outbreak_rate_functions(m=trend.fit[[trend.bestfit]], measure=measure, x=ds$X, y=ds$Y, order=trend.bestfit, value = value)

    magnitude.velocity = calculate_velocity(trend.bestfit.pdv)
    direction.radians = calculate_vector_direction(trend.bestfit.pdv$d.dx, trend.bestfit.pdv$d.dy )
    direction.degrees = (direction.radians * 180 / pi)

    return(list(
      order=trend.bestfit,
      measure = trend.fit.measure,
      m = trend.fit[[trend.bestfit]],
      ds=ds,
      velo = data.frame(magnitude.velocity=magnitude.velocity,
                        direction.radians=direction.radians,
                        direction.degrees=direction.degrees,
                        dx = trend.bestfit.pdv$d.dx,
                        dy = trend.bestfit.pdv$d.dy )))
  }

  if(manual.order == TRUE) {

    trend.fit = estimate_surfacetrend_models(ds, max.order)
    trend.bestfit.pdv = outbreak_rate_functions(m=trend.fit[[max.order]], measure=measure, x=ds$X, y=ds$Y, order=max.order, value = value)

    magnitude.velocity = calculate_velocity(trend.bestfit.pdv)
    direction.radians = calculate_vector_direction(trend.bestfit.pdv$d.dx, trend.bestfit.pdv$d.dy )
    direction.degrees = (direction.radians * 180 / pi)

    return(list(order=max.order,
                m = trend.fit[[max.order]], ds=ds,
                velo = data.frame(magnitude.velocity=magnitude.velocity,
                                  direction.radians=direction.radians,
                                  direction.degrees=direction.degrees,
                                  dx = trend.bestfit.pdv$d.dx,
                                  dy = trend.bestfit.pdv$d.dy )))
  }

}
