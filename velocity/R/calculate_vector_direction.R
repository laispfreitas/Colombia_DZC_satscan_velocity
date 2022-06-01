
#' Front wave velocity direction
#'
#'  Determines the direction of a vector given magnitude in x and y direction
#' @param x The magnitude of a vector in x
#' @param y The magnitude of a vector in y
#' @export
calculate_vector_direction = function(x,y) {
  theta = 0
  for(i in 1:length(x))
    if ( (x[i] > 0) && (y[i] > 0) ) { theta[i] = atan( y[i]/x[i]) } 
  else if ( (x[i] > 0) && (y[i] < 0) ) { theta[i] = 2*pi + atan( y[i]/x[i]) } 
  else if ( (x[i] < 0) && (y[i] < 0) ) { theta[i] = pi + atan( y[i]/x[i]) } 
  else if ( (x[i] < 0) && (y[i] > 0) ) { theta[i] = pi + atan( y[i]/x[i]) } 
  else (theta[i] = NA)
  return(theta)
}

