
#' Calculate the front-wave velocity
#'
#' Calculates the front-wave velocity of an outbreak given the outbreak rate functions
#' @param d The outbreak rate functions returned from outbreak_rate_functions
#' @export
calculate_velocity = function(d) {
  n = length(d[[1]]) ; velo = 0
  for(i in 1:n) { velo[i] = 1 /  sqrt( sum(( d$d.dx[i])^2 , (d$d.dy[i])^2) ) }
  return(velo)
}





