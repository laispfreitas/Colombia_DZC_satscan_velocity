
# velocitymap R package functions
# March 11, 2016

# VELOCITY VECTOR MAPPING ALGORITHM:
# 1. obtain geocoded point data (centroids if aggregates) and one date per point
# 2. response variable is 'time from first case' where first case - first date
# 3. time to event is regressed against X and Y to estimate a continuous surface
# 4. model of best fit is obtained - up to order X (user specified) with some measure of model fit
# 5. take partial derivatives with respect to X, Y for selected best fit model
# 6. find magnitude and direction of disease spread by finding the inner product of the vectors (magnitude = sqrt(x^2 + y^2)) pythagorus,
# 7. velocity is the inverse of the final value

require(ggplot2)
require(reshape2)
require(maptools)
require(ggmap)
require(rgdal)
require(dplyr)
require(grid)
require(raster)
require(rgeos)
require(RColorBrewer)



fitPolyModels = function(ds, max.order) {

  order = 1:max.order
  lm.fit = list()
  ntimes = length(ds$time)
  ds.new = ds
  ds.new$XY = ds[,'X'] * ds[,'Y']

  lm.fit[[1]] = lm( time ~ ., data=ds.new)

  for(i in 2:max.order) {

    # Needs to append a new column each iteration, keeping the old

    name.x = paste0("X",i)
    name.y = paste0("Y",i)
    ds.new[,name.x] = (ds.new$X^i)
    ds.new[,name.y] = (ds.new$Y^i)

    lm.fit[[i]] = lm( time ~ ., data=ds.new)
  }
  return(lm.fit)
}

getFitMeasures = function(lm.fit) {

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

getBestModel = function(summary, measure) {
    best.model = summary[which.min(summary[,measure]  ),]
    return(best.model$order)
}

getPartialDerivatives = function(m, measure, x, y, order) {

# Create the dataset for model order K
# Run the model for order K and save model object
# Write a generalizable way to take partial derivatives for a known order polynomial

  coefs.x = c("X")
  coefs.y = c("Y")
  if (order > 1) {
    coefs.x = c("X", paste0("X",2:order))
    coefs.y = c("Y", paste0("Y",2:order))
  }

  d.dx = m$coef['XY']*(y)
  d.dy = m$coef['XY']*(x)

  for(i in 1:length(coefs.x)) {

    d.dx = d.dx + (i) * m$coef[coefs.x][i] * x^(i - 1)
    d.dy = d.dy + (i) * m$coef[coefs.y][i] * y^(i - 1)
  }

return(list(d.dx=d.dx, d.dy=d.dy))
}

calcVelocity = function(d) {  #where d = object from getPartialDerivatives
  n = length(d[[1]]) ; velo = 0
  for(i in 1:n) { velo[i] = 1 /  sqrt( sum(( d$d.dx[i])^2 , (d$d.dy[i])^2) ) }
  return(velo)
}


getVectorDirection = function(x,y) {
  theta = 0
  for(i in 1:length(x))
    if ( (x[i] > 0) && (y[i] > 0) ) { theta[i] = atan( y[i]/x[i]) } #Q1
  else if ( (x[i] > 0) && (y[i] < 0) ) { theta[i] = 2*pi + atan( y[i]/x[i]) } #Q2
  else if ( (x[i] < 0) && (y[i] < 0) ) { theta[i] = pi + atan( y[i]/x[i]) } #Q3
  else if ( (x[i] < 0) && (y[i] > 0) ) { theta[i] = pi + atan( y[i]/x[i]) } #Q4
  else (theta[i] = NA)
  return(theta)
}


getVelocity = function(ds, max.order, measure, manual.order=FALSE) {
  if(manual.order == FALSE) {
      trend.fit = fitPolyModels(ds, max.order)
      trend.fit.measure = getFitMeasures(trend.fit)
      trend.bestfit = getBestModel(trend.fit.measure, measure)
      trend.bestfit.pdv = getPartialDerivatives(m=trend.fit[[trend.bestfit]], measure=measure, x=ds$X, y=ds$Y, order=trend.bestfit)

      magnitude.velocity = calcVelocity(trend.bestfit.pdv)
      direction.radians = getVectorDirection(trend.bestfit.pdv$d.dx, trend.bestfit.pdv$d.dy )
      direction.degrees = (direction.radians * 180 / pi)

      return(list(
        order=trend.bestfit,
        m = trend.fit[[trend.bestfit]],
        ds=ds,
        velo = data.frame(magnitude.velocity=magnitude.velocity,
                          direction.radians=direction.radians,
                          direction.degrees=direction.degrees,
                          dx = trend.bestfit.pdv$d.dx,
                          dy = trend.bestfit.pdv$d.dy )))
  }

  if(manual.order == TRUE) {

      trend.fit = fitPolyModels(ds, max.order)
      trend.bestfit.pdv = getPartialDerivatives(m=trend.fit[[max.order]], measure=measure, x=ds$X, y=ds$Y, order=max.order)

      magnitude.velocity = calcVelocity(trend.bestfit.pdv)
      direction.radians = getVectorDirection(trend.bestfit.pdv$d.dx, trend.bestfit.pdv$d.dy )
      direction.degrees = (direction.radians * 180 / pi)

      return(list(order=max.order, m = trend.fit[[max.order]], ds=ds,
                  velo = data.frame(magnitude.velocity=magnitude.velocity,
                                    direction.radians=direction.radians,
                                    direction.degrees=direction.degrees,
                                    dx = trend.bestfit.pdv$d.dx,
                                    dy = trend.bestfit.pdv$d.dy )))
  }

}


getContourPlot = function(r, scale, col.cat=4, contour.crs = "+proj=utm +zone=28 +datum=NAD83", save.name='tempcontour') {

  getContourFitData = function(xran,yran) {
    order = r$order
    ds.new =  data.frame(XY = xran*yran, X=xran, Y = yran)

    for(i in 2:order) {
      # Needs to append a new column each iteration, keeping the old
      name.x = paste0("X",i)
      name.y = paste0("Y",i)
      ds.new[,name.x] = (xran^i)
      ds.new[,name.y] = (yran^i)
    }
    return(ds.new)
  }

  by.x = (max(r$ds$X) - min(r$ds$X) ) / (scale - 1)
  by.y = (max(r$ds$Y) - min(r$ds$Y) ) / (scale - 1)

  xran = seq( range(r$ds$X)[1] , range(r$ds$X)[2], by.x) ; xran
  yran = seq( range(r$ds$Y)[1] , range(r$ds$Y)[2], by.y) ; yran

  modelFun = function(X,Y) {
    ds.new = getContourFitData(X,Y)
    temp1 = predict(r$m, newdata=ds.new )
    return(temp1) }

  z = outer(xran, yran, modelFun)

  ds.contour=list()
  ds.contour$xran = xran
  ds.contour$yran = yran
  ds.contour$z = round(z,0)

  raster.contour = raster(ds.contour, crs=contour.crs)
  contour.shp = rasterToPolygons(raster.contour,n=8,dissolve=T)
  writeSpatialShape(contour.shp,paste0(save.name,'.shp'))
}



getVectorMap = function(r, shpfile) {

  shp.fort = fortify(shpfile, region="ID")
  shp.map.temp  = merge(shp.fort, r$velo, by='id', all.x=TRUE, order=FALSE)
  shp.map = shp.map.temp[order(shp.map.temp$order), ]

  vectors.mapping = data.frame(x=r$ds$X, y=r$ds$Y, velocity=r$velo$magnitude.velocity, dx = r$velo$dx, dy = r$velo$dy, dir=r$velo$direction.degrees, id=r$velo$id)

  g = ggplot()  +
  geom_polygon(data=shp.map, aes(x=long, y=lat, group = group, fill=magnitude.velocity), col="white", size=0.1) +
    scale_fill_gradient(low='#F6D2BF', high='#D53534', na.value = "grey70") +
    geom_segment(data=vectors.mapping, aes(x = x, y = y, xend = x + dx, yend = y + dy), colour='black', arrow = arrow(angle=30,type='open',length = unit(0.3,"cm"))) +
    theme(panel.background = element_rect(fill = 'white', colour = 'grey90')) +
    xlab("m") + ylab("m") +
    theme(legend.position="right", legend.title=element_blank())

  return(g)
}






