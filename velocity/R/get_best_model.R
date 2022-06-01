
#' Get the best fitting surface trend model
#'
#' Retrieves the order of the best fitting mdoel based on the given measure of model fit
#' @param summary Dataframe returned from function surfacetrend_model_fit
#' @param measure The measure of model fit, defaults to R-squared ('r2'), can be 'aic' or 'bic'
#' @export
get_best_model = function(summary, measure='r2') {
  best.model = summary[which.min(summary[,measure]  ),]
  return(best.model$order)
}
