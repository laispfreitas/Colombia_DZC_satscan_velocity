

#' Assess the model fit
#'
#' This function assess the fit of the surface trend polynomial regression models by AIC, BIC, and R-squared
#' @param lm.fit A list of lm model objects
#' @export
surfacetrend_model_fit = function(lm.fit) {

  get.r2 = function(m) {
    r2 = summary(m)$adj.r.squared
    return(r2)
  }
  summary = data.frame(order=1:length(lm.fit))
  summary$aic = lapply(lm.fit, FUN=AIC)
  summary$bic = lapply(lm.fit, FUN=BIC)
  summary$r2  = lapply(lm.fit, FUN=get.r2)
  return(summary)
}
