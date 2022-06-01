#.........................................................................#
#            SATSCAN - DENGUE, ZIKA AND CHIKUNGUNYA IN COLOMBIA           #
#.........................................................................#

# Laís Picinini Freitas
# Last update: 02-Dec-2021

## This code is part of the study: 
# "Spatio-temporal clusters and patterns of spread of dengue, chikungunya, and Zika in Colombia"
# by Laís Picinini Freitas, Mabel Carabali, Mengru Yuan, Gloria I. Jaramillo-Ramirez, 
# Cesar Garcia Balaguera, Berta N. Restrepo, and Kate Zinszer

## Preprint available at https://doi.org/10.1101/2022.03.17.22272536


# Packages ----------------------------------------------------------------
library(rsatscan)
library(sf)
library(tidyverse)
library(colorspace)


# Data --------------------------------------------------------------------

centroids <- read_sf('centroids.shp')
centroids <- centroids %>%
  mutate(lat = unlist(map(centroids$geometry,1)),
         long = unlist(map(centroids$geometry,2)))

munis <- read_sf('municipalities_2017.shp')

load('DZC_satscan_data.RData')


# SaTScan -----------------------------------------------------------------

tmp <- tempdir()
old <- getwd()
setwd(tmp)

# SaTSCan must be installed, it can be downloaded at:
# https://www.satscan.org/download.html 

# Loading satscan and the template with parameters
ss.local <- "~/bin/SaTScan/"  # change path to where SaTScan is installed
ssenv$.ss.params <- readLines('template_satscan.txt') 

# Creating population and centroids files:
write.table(satscan_data[,c(2,1,4)],file='dzc.pop' , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(centroids[,c(3,13,12)],file='dzc.geo' , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")

run_satscan <- function(startdate, enddate, maxpop, maxradius) {
  ss.options(list(CaseFile=ARQ,
                  StartDate=startdate,
                  EndDate=enddate, 
                  PrecisionCaseTimes=3, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic     
                  PopulationFile="dzc.pop",
                  CoordinatesFile="dzc.geo", 
                  CoordinatesType=1, # 0=Cartesian, 1=latitude/longitude
                  AnalysisType=3, # 1=Purely Spatial, 3=Retrospective Space-Time
                  ModelType=0, # 0=Discrete Poisson
                  ScanAreas=1, # 1=High, 2=Low, 3=High and Low
                  TimeAggregationUnits=3, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
                  TimeAggregationLength=7,
                  MaxSpatialSizeInPopulationAtRisk=50,
                  MinimumTemporalClusterSize=14, # 2 weeks
                  MaxTemporalSize=182, # 26 weeks
                  MinimumCasesInHighRateClusters=100,
                  UseDistanceFromCenterOption_Reported='y',
                  MaxSpatialSizeInDistanceFromCenter_Reported=maxradius, # max size radius cluster in km                
                  MaxSpatialSizeInPopulationAtRisk_Reported=maxpop # max population inside cluster
  ))
  
  ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))
  
  modelo <- paste0(ARQ, '_modelo')
  write.ss.prm(tmp, modelo)
  satscan(tmp, modelo, sslocation="~/bin/SaTScan/") # change path to where SaTScan is installed
  
}


## Zika -------------------------------------------------------------------

DISEASE <- 7     
ARQ <- 'zika.cas'

write.table(satscan_data[,c(2,DISEASE,1)],file=ARQ , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")

result_zika <- run_satscan("2015/01/07", "2018/12/26", 20, 150) 

summary(result_zika)
result_zika


# Map with the clusters
rgis <- subset(result_zika$gis) %>% filter(P_VALUE < 0.05)
munis.clus <- munis %>% left_join(rgis, by = c('MPIO_CCDGO' = 'LOC_ID')) 
munis.clus$CL2 <- factor(munis.clus$CLUSTER)

ggplot() + 
  geom_sf(data = munis.clus, aes(fill = CL2), size = 0.1)  + 
  scale_fill_discrete_sequential(palette = 'Viridis', 
                                 na.value = "#EEEEEE",
                                 name = "Cluster") +
  labs(title = "Zika 20% 150 km") +
  theme_void()


## Chikungunya ------------------------------------------------------------

DISEASE <- 6
ARQ <- 'chik.cas'

write.table(satscan_data[,c(2,DISEASE,1)],file=ARQ , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")

result_chik <- run_satscan("2014/01/01", "2018/12/26", 20, 150) 

summary(result_chik)
result_chik


# Map with the clusters
rgis <- subset(result_chik$gis) %>% filter(P_VALUE < 0.05)
munis.clus <- munis %>% left_join(rgis, by = c('MPIO_CCDGO' = 'LOC_ID')) 
munis.clus$CL2 <- factor(munis.clus$CLUSTER)

ggplot() + 
  geom_sf(data = munis.clus, aes(fill = CL2), size = 0.1)  + 
  scale_fill_discrete_sequential(palette = 'Viridis', 
                                 na.value = "#EEEEEE",
                                 name = "Cluster") +
  labs(title = "Chikungunya 20% 150 km") +
  theme_void()


## Dengue -----------------------------------------------------------------

DISEASE <- 5
ARQ <- 'dengue.cas'

write.table(satscan_data[,c(2,DISEASE,1)],file=ARQ , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")

result_dengue <- run_satscan("2014/01/01", "2018/12/26", 20, 150) 

summary(result_dengue)
result_dengue


# Map of the clusters
rgis <- subset(result_dengue$gis) %>% filter(P_VALUE < 0.05)
munis.clus <- munis %>% left_join(rgis, by = c('MPIO_CCDGO' = 'LOC_ID')) 
munis.clus$CL2 <- factor(munis.clus$CLUSTER)

ggplot() + 
  geom_sf(data = munis.clus, aes(fill = CL2), size = 0.1)  + 
  scale_fill_discrete_sequential(palette = 'Viridis', 
                                 na.value = "#EEEEEE",
                                 name = "Cluster") +
  labs(title = "Dengue 20% 150 km") +
  theme_void()


## Multivariate -----------------------------------------------------------

ssenv$.ss.params <- readLines('multiple_datasets.txt')  

ARQ <- 'dengue.cas'
ARQ2 <- 'zika.cas'
ARQ3 <- 'chik.cas'

write.table(satscan_data[,c(2,5,1)],file=ARQ , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(satscan_data[,c(2,7,1)],file=ARQ2 , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(satscan_data[,c(2,6,1)],file=ARQ3 , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")


ss.options(list(CaseFile=ARQ,
                StartDate="2014/01/01",
                EndDate="2018/12/26", 
                PrecisionCaseTimes=3, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic     
                PopulationFile="dzc.pop",
                CoordinatesFile="dzc.geo", 
                CoordinatesType=1, # 0=Cartesian, 1=latitude/longitude
                AnalysisType=3, # 1=Purely Spatial, 3=Retrospective Space-Time
                ModelType=0, # 0=Discrete Poisson
                ScanAreas=1, # 1=High, 2=Low, 3=High and Low
                TimeAggregationUnits=3, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
                TimeAggregationLength=7,
                MaxSpatialSizeInPopulationAtRisk=50,
                MinimumTemporalClusterSize=14, # 2 weeks
                MaxTemporalSize=182, # 26 weeks
                MinimumCasesInHighRateClusters=100,
                UseDistanceFromCenterOption_Reported='y',
                MaxSpatialSizeInDistanceFromCenter_Reported=150, # max size radius cluster in km                
                MaxSpatialSizeInPopulationAtRisk_Reported=20, # max population inside cluster
                MultipleDataSetsPurposeType=0,
                CaseFile2=ARQ2,
                PopulationFile2='dzc.pop',
                CaseFile3=ARQ3,
                PopulationFile3='dzc.pop'
))

ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))
modelo <- paste0(ARQ, '_Modelo')
write.ss.prm(tmp,modelo)
result_dzc <-  satscan(tmp,modelo, sslocation="~/bin/SaTScan/") # change path to where SaTScan is installed

summary(result_dzc)
result_dzc


# Map of the clusters

rgis_dzc <- subset(result_dzc$gis)
munis.clus_dzc <- munis %>% left_join(rgis_dzc, by = c('MPIO_CCDGO' = 'LOC_ID')) 
munis.clus_dzc$CL2 <- factor(munis.clus_dzc$CLUSTER)

ggplot() + 
  geom_sf(data = munis.clus_dzc, aes(fill = CL2), size = 0.1)  + 
  scale_fill_discrete_sequential(palette = 'Viridis', 
                                 na.value = "#EEEEEE",
                                 name = "Cluster") +
  theme_void() 