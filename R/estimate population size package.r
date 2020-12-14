# Estimate population size based on distance sampling and species distribution modelling
# Xinhai Li (xinhai_li_edu@126.com)
#============================================================================================================================
# library(estPopsize)


# Fit distance sampling
#############################################################################################################################
#' Fit data of  distance sampling using the best combination of Key and Adjustment
#'
#' @description This function ranks the 9 combinations of Key and Adjustment in distance sampling
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param survey A data.frame with columns 'species', 'size', 'distance', 'Region.Label', 'Effort'
#'
#' @return
#'
#' @examples
#'
#'  attach(distancesampling)
#'  AICs = distanceSampling(distancesampling)
#'  AICs
#'  ds.kiang <- ds(kiang, key = AICs$Key[2], adjustment = AICs$Adjustment[2], convert.units = 0.001) # truncation=1000
#'  summary(ds.kiang)
#'  SM = summary(ds.kiang)
#'  Average.p = SM[[3]][[14]][1]
#'  plot(ds.kiang, main="Detection rate\n Tibetan wild ass, half-normal detection function")
#'
#' @import Distance
#' @export

distanceSampling = function(survey){
  library(Distance)
  KEY = c("hn", "unif", "hr")
  ADJ = c("cos", "herm", "poly")
  AICs = data.frame(ID=1:9, Key = NA, Adjustment = NA, AIC = NA)
  N = 0
  for (i in 1:length(KEY)){
    for (j in 1:length(ADJ)){
      tryCatch({
      N = N+1
      AICs$Key[N] = KEY[i]
      AICs$Adjustment[N] = ADJ[j]
      ds.fit <- ds(survey, key = KEY[i], adjustment = ADJ[j], convert.units = 0.001)
      AICs$AIC[N] = as.numeric(AIC(ds.fit)[2])
      print(paste('Finished ', round(N/9*100, 2), "%", sep=''))
      }, error = function(e){cat("ERROR: ", conditionMessage(e),"\n")})
    }
  }
  AICs = AICs[order(AICs$AIC), ]
  return(AICs)
}
#############################################################################################################################





# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#' Crop raster layers to fit data in a smaller region
#'
#' @description Crop raster layers to fit data in a smaller region using user defined buffer
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param species A data.frame of species occurrences with columns "Name_E", "Lon", "Lat", "Date", "Count"
#' @param buffer A value of distance (unit: degree) defining the width of buffer zone outside the occurrences
#' @param Envlayers A RasterBrick (a multi-layer raster object)
#'
#' @return
#'
#' @examples
#'  attach(species) # load occurrences data
#'  BioClim <- brick('D:/GIS/CLIMATE/var27.grd')
#'  BioClim = cropLayers(species, buffer=0.5, Envlayers=BioClim)
#'
#' @import raster
#' @export

cropLayers = function(species, buffer, Envlayers){
  lon.r = range(species$Lon) + c(buffer*(-1), buffer)
  lat.r = range(species$Lat) + c(buffer*(-1), buffer)
  newlayers = crop(Envlayers, extent(c(lon.r, lat.r)))
  return(newlayers)
}
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA




# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#' Provide pseudo-absent points and derive environmental variables
#'
#' @description This function provides pseudo-absent points using user defined sample size,
#'  derives values of environmental variables at occurrences and pseudo-absent points.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param species A data.frame of species occurrences with columns "Name_E", "Lon", "Lat", "Date", "Count"
#'
#' @param buffer A value of distance (unit: degree) defining the width of buffer zone outside the occurrences
#'
#' @param absence The number of pseudo-absent points on the longitude side and the latitude side
#'
#' @param Envlayers A RasterBrick (a multi-layer raster object)
#'
#' @return Return a data.frame with values of environmental variables at occurrences and pseudo-absent points.
#'
#' @examples
#'  attach(species) # load occurrences data
#'  BioClim <- brick('D:/GIS/CLIMATE/var27.grd')
#'  Data = getEnvData(species, buffer=0.5, absence=30, Envlayers=BioClim); head(Data)
#'
#' @import raster
#' @export
#'

getEnvData = function(species, buffer, absence, Envlayers){
  species.loc = species[,  c('Lon','Lat')]
  lon = seq(min(species.loc[,1]) - buffer, max(species.loc[,1]) + buffer, length.out = absence)
  lat = seq(min(species.loc[,2]) - buffer, max(species.loc[,2]) + buffer, length.out = absence)
  back = expand.grid(lon, lat)
  bak <- extract(Envlayers, back)
  bak = as.data.frame(bak)
  bak = cbind(Name='absent', Count=0, bak)
  bak = bak[!is.na(bak$elev),]
  bak = cbind(bak, back)
  names(bak)[30:31] = c('Lon','Lat')
  spe = extract(Envlayers, species.loc) #Lon MUST be the first column
  spe = as.data.frame(spe)
  spe = cbind(Name=species$Name_E, Count=species$Count, spe)
  spe = cbind(spe, species.loc)
  ENV = rbind(spe, bak)
  return(ENV)
}
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







# Predict population density
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#' Predict population density across the range of raster layers
#'
#' @description Predict number of individuals at each cell using Random Forest, based on environmental variables
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param ENV A RasterBrick (a multi-layer raster object)
#' @param RF A list object of the model results of Random Forest
#'
#' @return A raster layer with number of individuals at each cell
#'
#' @examples
#'
#'pred = popSize(BioClim, RF)
#'plot(pred)
#'plot(log(1+log(1+pred)), xlab="Longitude", ylab="Latitude", main='',
#'     col=colorRampPalette(c("grey90", "green", "yellow", "red"))(12)) #topo.colors(20)
#'lines(shape)
#'library(calibrate)
#'attach(loc)
#'loc1 = loc[loc$Class==1,]
#'loc2 = loc[loc$Class==2,]
#'loc3 = loc[loc$Class==3,]
#'textxy(loc1$lon, loc1$lat, loc1$Name, cex=0.8);textxy(loc2$lon, loc2$lat, loc2$Name, cex=0.6);textxy(loc3$lon, loc3$lat, loc3$Name, cex=0.4);points(loc1$lon, loc1$lat, pch=16, cex=1, col='brown')
#'points(species$Lon, species$Lat, pch=16, cex=species$Count/20, col=adjustcolor("darkgreen", 0.5))
#' @import randomForest
#' @export
#'

popSize = function(ENV=BioClim, model=RF){
  pop.size <- predict(ENV, model, type="response")
  return(pop.size)
}
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

# shape <- readOGR('d:/GIS/Qinghai/route2017.shp')






#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
#' Derive representative points from survey routes
#'
#' @description This function derives track points, and thinning the points to the density of one point per km.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param shape A GIS shapefile containing survey routes
#'
#' @return A dataframe containing the Lat/Lon of thinned track points.
#'
#' @examples
#'
#'  tracks = trackPoints(shape)
#'
#' @import sp
#' @import rgdal
#' @import maptools
#'
#' @export
#'

trackPoints = function(shape){
  coord = coordinates(shape)# package sp, list of sections
  lonlat = as.data.frame(coord[[1]][[1]])
  names(lonlat) = c('Lon', 'Lat')
  Dist = 0
  for (k in 2:length(coord)){
    LL = as.data.frame(coord[[k]][[1]])
    names(LL) = c('Lon', 'Lat')
    S = LL
    S = cbind(S, dist=NA)
    N = nrow(S)
    for (i in 1:(N-1)){
      S$dist[i+1] = (((S$Lat[i+1]-S$Lat[i])*39946.79/360)^2+
                       ((S$Lon[i+1]-S$Lon[i])*pi*12756.32/360*cos(S$Lat[i]*pi*2/360))^2)^0.5
    }
    dist = sum(S$dist, na.rm=T)
    Dist = dist + Dist # total distance of the survey routes
    lonlat = rbind(lonlat, LL) # total points of the routes
  }
  thin = nrow(lonlat) / Dist
  LonLat = lonlat[seq(1, nrow(lonlat), by = floor(thin)), ]
  return(LonLat)
}
#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP






## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#' Estimate population size at the whole area
#'
#' @description Predict the population size using a species distribution model, and adjust the predicted value
#'  based on real observations by distance sampling
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param pop.pred A raster layer with number of individuals at each cell
#' @param tracks A dataframe containing the Lat/Lon of thinned track points.
#' @param species A data.frame of species occurrences with columns "Name_E", "Lon", "Lat", "Date", "Count"
#' @param detection The average detection during distance sampling
#'
#' @return A vector including the predicted population size by Random Forest, predicted population on the routes,
#' Observed population on the routes, and final adjusted estimation of population size
#'
#' @examples
#'
#'  EST = estPopSize(pred, tracks, species, Average.p); EST
#'
#' @import raster
#' @export
#'

estPopSize = function(pop.pred, tracks, species, detection){
  library(raster)
  pop_ori = cellStats(pop.pred, stat='sum', na.rm=TRUE) # total individuals for the raster
  pre <- extract(pop.pred, tracks) # predicted pop on survey routes
  pop_pre = sum(pre, na.rm=T) # total predicted individuals on the route
  pop_obs = sum(species$Count) # observation
  pop_obs = pop_obs / detection # distance sampling adjustment
  adjust = pop_pre / pop_obs # SDM adjustment
  pop_est = cellStats(pop.pred, stat='sum', na.rm=TRUE) / adjust
  out = c(pop_ori, pop_pre, pop_obs, pop_est)
  names(out) = c("Original prediction", "Predicted population on the routes",
                 "Observed population on the routes", "Final estimation")
  return(out)
}
## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF








## BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
#' Comparing prediction-observation ratio at different zones of the whole area
#'
#' @description Calculating the uncertainty of prediction-observation rate
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param pop.pred A raster layer with number of individuals at each cell
#' @param tracks A dataframe containing the Lat/Lon of thinned track points.
#' @param species A data.frame of species occurrences with columns "Name_E", "Lon", "Lat", "Date", "Count"
#' @param buffer A value of distance (unit: degree) defining the width of buffer zone outside the occurrences
#' @param grid A value defining the number of grids (grid^2) for prediction-observation comparison
#'
#' @return The standard deviation of the observation/prediction ratios at all panels of the whole area
#'
#' @examples
#'
#'  sd.Match = spatialMatch(species, buffer, tracks, pop.pred, 6); sd.Match
#'
#' @import raster
#' @export
#'

spatialMatch = function(species, buffer, tracks, pop.pred, grid=4){
  lon.r = range(species$Lon) + c(buffer*(-1), buffer)
  lat.r = range(species$Lat) + c(buffer*(-1), buffer)
  lon = seq(lon.r[1], lon.r[2], length.out = grid+1)
  lat = seq(lat.r[1], lat.r[2], length.out = grid+1)
  Match = numeric(grid*grid)
  N = 0
  for (i in 1:grid){
    for (j in 1:grid){
      SPE = species[species$Lon > lon[i] & species$Lon < lon[i+1], ]
      SPE = SPE[SPE$Lat > lat[j] & SPE$Lat < lat[j+1], ]
      OBS = sum(SPE$Count, na.rm=T); OBS
      TRK = tracks[tracks$Lon > lon[i] & tracks$Lon < lon[i+1], ]
      TRK = TRK[TRK$Lat > lat[j] & TRK$Lat < lat[j+1], ]
      PRE <- sum(extract(pop.pred, TRK))
      RAT = OBS / PRE
      # PRED = crop(pop.pred, extent(lon[i], lon[i+1], lat[j], lat[j+1])) # whole area
      N = N+1
      Match[N] = RAT
    }
  }
  Match = Match[! Match == Inf]
  Match = sd(Match, na.rm=T)
  return(Match)
}
## BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
#




## RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
#' Calculate the confidence intervals of population size
#'
#' @description Evaluation the confidence intervals of the population estimation taking into account the uncertainty
#'  of distance sampling, species distribution modelling, and prediction-observation matchness
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param pop A value of estimated population size.
#' @param detection 1 / the average detection rate of distance sampling
#' @param SDM 1 / the explained variance by the species distribution model
#' @param Match The standard deviation of the observation/prediction ratios at all panels of the whole area.
#' @param plotrange The plotting range for estimated population size.
#'
#' @return The 95% confidence interval of the estimated population size
#'
#' @examples
#'
#'  CI(10000, 0.2, 0.8, 0.8)
#'  SDM = 1/max(RF$rsq)
#'  detection = 1/Average.p
#'  CI(EST[4], detection, SDM, sd.Match, 100)
#'
#' @export
#'

CI = function(pop, detection, SDM, Match, plotrange=10){
  CI = numeric()
  for (i in 1:10000){
    CI[i] = pop * rnorm(1, 1, Match) * rnorm(1, 1, SDM) * rnorm(1, 1, detection)
  }
  plot(density(CI), xlab='Population size', ylab='Density', main="", xlim=c(pop * (-1), pop*plotrange))
  Q = quantile(CI, probs = c(0.05, 0.95))
  abline(v= Q[1], lty=2); abline(v= Q[2], lty=2)
  return(Q)
}
## RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR



