---
title: "abundanceR: estimating wildlife abundance based on distance sampling and species distribution models"
output: word_document
vignette: >
  %\VignetteIndexEntry{Introduction to mypackage}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


## Introduction

The package abundanceR is designed to estimate wildlife abundance across expansive geographical scales. By integrating Distance Sampling and Species Distribution Models (SDMs), the package leverages field survey data and environmental variables to predict animal abundance in unsurveyed regions. Furthermore, it provides robust functionality for comprehensive uncertainty assessment. It has following key functions:

* Distance Sampling Analysis: Processes survey data to estimate detection functions.
  
* SDM Integration: Incorporates multi-source environmental variables to construct species distribution models, and predict species abundance in every quadrat across the whole study area, including those not surveyed.

* Spatial Correction: Adjusts model prediction results based on empirical data from actual survey routes.

* Uncertainty Quantification: Quantifies and synthesizes three distinct sources of uncertainty: survey-related errors, model performance, and spatial heterogeneity.

* Built-in Data Support: Provides a global terrestrial environmental dataset with a spatial resolution of 1 km².


## Species Applicability

The effectiveness of abundanceR is highly dependent on the ecological characteristics of the target species and their habitat:

* Habitat Requirements: The package is best suited for open environments—such as grasslands and wetlands—where the distance sampling methodology is most viable.
  
* Species Visibility: It is primarily designed for large, conspicuous species that are easily detected during surveys, such as ungulates.

* Distribution Patterns: The model is highly effective for habitat-specific species whose spatial distribution is strictly constrained by environmental variables.

* Limitations: This package is not recommended for elusive or hard-to-detect species (e.g., carnivores) or for habitat generalists whose distribution shows little specificity to environmental gradients.


## Data Preparation

The package abundanceR requires three core types of data: Species Distribution Data, Environmental Variable Data, and Survey Route Data. The specific requirements for each are detailed below.

1. Species Distribution Data
* File Format: CSV or Excel formats are recommended; an R data.frame can also be used directly.
* Mandatory Fields (Field names and case sensitivity must be exact):
  species: Species name (Character).
  size: Group size (Integer), representing the total number of individuals at the observation point.
  distance: Perpendicular distance between the observer and the animal (Numeric, unit: meters).
  Lat: Latitude of the observation point (Numeric, WGS84 coordinate system).
  Lon: Longitude of the observation point (Numeric, WGS84 coordinate system).
* Optional Fields: Includes Side (relative orientation), Elev (elevation in meters), Date, and Time.

2. Environmental Variable Data
* File Format: Multi-band raster data, supporting .grd/.gri, .tif, and other common formats.
* Data Requirements:
  A spatial resolution of 1km² is recommended for consistency with the package's built-in variables.
  All environmental layers must share the same extent, resolution, and coordinate system (WGS84 is recommended).
  Variables can include climate, topography, vegetation, and human footprint indices.
* Built-in Dataset: The package provides a global dataset of 29 environmental variables (var29.grd) at 1km²   resolution.Included Variables: 19 bioclimatic variables, elevation, human footprint index, land cover, wetlands, as well as solar radiation, wind speed, and water vapor pressure for January and July. 
* Data Access Links: Baidu Netdisk (China users): https://pan.baidu.com/s/1noU8A7WcsuYx0MSiQq6CeQ (Code: 1234) 
Google Drive: https://drive.google.com/drive/folders/1bNh4SdikmjrOkgqE5VOVo86SozD2hvmT


3. Survey Route Data
* File Format: A GIS shapefile with polylines of survey routes.

## General Workflow

The first step involves installing and loading the requisite R packages. Ensure that your species distribution CSV file contains the following mandatory fields: species, size, distance, Lat, and Lon. Additionally, the var29.grd environmental dataset must be pre-downloaded.

```{r}
# Install and load core packages 
if (!require(devtools)) install.packages("devtools") 
library(devtools)
install_github("Xinhai-Li/abundanceR", force = TRUE)
library(abundanceR)
```

* Load example survey data, the occurrences of the Tibetan wild ass.
```{r}
data(kiang)
head(kiang)
library(sf)
# download route2017.zip from https://github.com/Xinhai-Li/abundanceR/blob/main/route2017.zip
shape <- st_read('d:/GIS/Qinghai/route2017.shp') # Use your own file location
plot(shape) # survey routes
```

**Fig. 1 The survey routes for the Tibetan wild ass in 2017**

* Keep occurrences within 500 m of the survey routes on both sides to match the 1 km² prediction resolution.
```{r}
kiang = kiang[kiang$distance<=500, ]
mean(kiang$size); sd(kiang$size)
sum(kiang$size[kiang$distance<=500]) # 449 individuals observed
length(kiang$size[kiang$distance<=500]) # 103 occurrences
```

* Process Distance sampling data, fit the detection function.
```{r}
library(Distance)
ds.fit <- ds(kiang, convert_units = 0.001)
set.seed(1)
# Compare model fit of the combinations of all three keys ("hn", "hr", "unif") and adjustment ("cos", "herm", "poly").
AICs = distanceSampling(kiang[kiang$distance<=500, ]) # Take a few seconds
AICs = AICs[!is.na(AICs$AIC),] # in case no AIC value being calculated
# The first row in AICs has the best parameters.
ds.kiang <- ds(kiang, key = AICs$Key[1], adjustment = AICs$Adjustment[1], convert_units = 0.001, truncation=500)
summary(ds.kiang)
SM = summary(ds.kiang)
Average.p = SM$ds$average.p; Average.p # average detection rate across the distance range
survey.uncertainty = 1 - Average.p # survey uncertainty
plot(ds.kiang, main=paste("Key:", AICs$Key[1], "\n", "Adjustment:", AICs$Adjustment[1], sep=" "))
```

**Fig. 2 The detection function of distance sampling for the Tibetan wild ass in 2017**

* Load environmental data.
```{r}
library(raster)
# BioClim <- brick('var29.grd') # var29.grd is the 105 GB file with 29 variables covering the whole world
# You can download the loacl environmental data (22.7 MB) from 
# https://github.com/Xinhai-Li/abundanceR/blob/main/BioClim.zip
BioClim <- brick('D:/GIS/CLIMATE/SJY.grd')
# Crop the global data to fit your study area based on the ranges of the species occurrences
BioClim = cropLayers(kiang, buffer=0.5, Envlayers=BioClim) 
# the unit for buffer is arc degree, and 1 degree is about 100 km for most areas.

BioClim;res(BioClim);nlayers(BioClim);names(BioClim)
layer = BioClim[[27]]
cellStats(layer, stat='mean', na.rm=TRUE);
cellStats(layer, stat='range', na.rm=TRUE)
cellStats(layer, stat='sd', na.rm=TRUE)
par(mfrow=c(3,4)) 
plot(BioClim[[17:27]])
```

**Fig. 3 Representative environmental variables for species distribution modelling in the Three-River-Source National Park**


* Show occurrences of the Tibetan wild ass.
```{r}
par(mfrow=c(1,1))
data(kiang); head(kiang); sum(kiang$size)
plot(BioClim[[20]], col=terrain.colors(20), xlab="Longitude", ylab="Latitude");lines(shape)
points(kiang$Lon, kiang$Lat, pch=16, cex=log(kiang$size)/2, col=adjustcolor("blue", 0.5))
```

**Fig.4 Occurrences of the Tibetan wild ass**


* Derive environmental data.
```{r}
head(kiang)
D = kiang
Data = getEnvData(D, buffer=0.2, absence=30, Envlayers=BioClim) # this buffer should be smaller than that in cropLayers
head(Data)
table(Data$size)
```

* Plot presence and absence data.
```{r}
plot(BioClim[[1]], xlab="Longitude", ylab="Latitude")
points(Data$Lon[Data$Name=="absent"], Data$Lat[Data$Name=="absent"], col='red')
points(Data$Lon[Data$Name=="kiang"], Data$Lat[Data$Name=="kiang"])
```

**Fig.5 The presence data (black circles) and pseudo-absence data (red circles)**

* Conduct species distribution modelling
```{r}
library(randomForest)
no.col = ncol(Data)
#fill null values
Dat.fill <- na.roughfix(Data[,2:(no.col-4)]) # use 27 variables
head(Dat.fill)

set.seed(2)
RF <- randomForest(Dat.fill[,2:ncol(Dat.fill)], Dat.fill[,1], ntree=500, importance=TRUE, na.action=na.roughfix)
RF #
model.uncertainty = 1-max(RF$rsq) # model uncertainty

# Predict animal abundance
pred = popSize(BioClim[[1:27]], RF)
par(mfrow=c(1,1));plot(pred)
points(kiang$Lon, kiang$Lat, pch=16, cex=log(kiang$size)/2, col=adjustcolor("blue", 0.5))
# writeRaster(pred, filename="prediction.tif", format="GTiff", overwrite=TRUE)
```

**Fig.6 Predicted population density of the Tibetan wild ass in the Three-River-Source National Park**


* Display the subdivided study areas.
```{r}
species = kiang
grid=4
buffer = (range(species$Lon)[2]-range(species$Lon)[1])/5
lon.r = range(species$Lon) + c(buffer/1.2*(-1), buffer/1.2)
lat.r = range(species$Lat) + c(buffer*(-1), buffer)
lon = seq(lon.r[1], lon.r[2], length.out = grid+1)
lat = seq(lat.r[1], lat.r[2], length.out = grid+1)

plot(log(1+log(1+pred)), xlab = "Longitude", ylab = "Latitude", main='',
     col = colorRampPalette(c("grey90", "green", "yellow", "red"))(12),
     xlim=lon.r, ylim=lat.r)
for (i in 1:(grid+1)){
  abline(v = lon[i], col="white",lwd=2)
  abline(h = lat[i], col="white",lwd=2)
}
lines(shape, lwd=1.5)
points(species$Lon, species$Lat, pch=16, cex=log(species$size)/2, col=adjustcolor("red", 0.5))
```

**Fig.7 Subdivided study areas for comparing predicted and observed animal densities**


* Zoom in: occurrences and predicted values
```{r}
par(mfrow=c(1,1))
plot(log(1+log(1+pred)), xlab="Longitude", ylab="Latitude", main='',
     col=colorRampPalette(c("grey90", "green", "yellow", "red"))(12),
     xlim=c(93.6,94.2), ylim=c(33.8,34.2)) #topo.colors(20)
lines(shape, lwd=2)
points(Data$Lon[Data$Name=="kiang"], Data$Lat[Data$Name=="kiang"], cex=log(Data$size)/0.7,
       col=adjustcolor("black", 0.5), pch=16)
points(Data$Lon[Data$Name=="kiang"], Data$Lat[Data$Name=="kiang"], cex=log(Data$size)/0.7, col="white", pch=1)
library(calibrate)
textxy(Data$Lon[Data$Name=="kiang"], Data$Lat[Data$Name=="kiang"], Data$size[Data$Name=="kiang"], cex=0.8, col="red")
```

**Fig.8 Local enlarged area showing the occurrences with different groups sizes and the predicted abundance of each pixel**


* Derive points at 1 km interval from the survey routes
```{r}
tracks = trackPoints(shape) # take a few minutes
```

* Estimate animal abundance
```{r}
EST = estPopSize(pred, tracks, kiang, Average.p); EST
```

* Sum animal count and make adjustment
```{r}
# Predicted number of individuals for the study area
(pop_ori = cellStats(pred, stat='sum', na.rm=TRUE)) # 3314850

# Predicted animal count in each quadrat intersected by the survey routes.
pre <- extract(pred, tracks) 

# The predicted number of individuals on the route, 5671
(pop_pre = sum(pre, na.rm=T)) 

# Keep occ within 500m to match quadrat of 1 km
kiang = kiang[kiang$distance<=500, ] 

# Observed number of individuals during the survey, 449
(pop_obs = sum(kiang$size)) 

# Distance sampling adjustment, 601
(pop_obs_adj = pop_obs / Average.p)

# Adjustment rate based real observations
(adjust = pop_pre / pop_obs_adj) # SDM adjustment, 9.4

# Final estimated animal abundance
(pop_est = cellStats(pred, stat='sum', na.rm=TRUE) / adjust) # 351555.2
```


* Calculate confidence intervals
```{r}
adjust.uncertainty = spatialMatch(kiang, tracks, pred, 4); adjust.uncertainty
# argument 4 means the whole area are splitted to 4*4 subregions for comparing predicted-observed ratios

survey.uncertainty = 1 - Average.p
model.uncertainty
adjust.uncertainty

# Confidence intervals
# Include all uncertainties
CI(EST[4], survey.uncertainty, model.uncertainty, adjust.uncertainty) 

# Ignore survey uncertainty
CI(EST[4], 0.001, model.uncertainty, adjust.uncertainty) 

# Ignore model uncertainty
CI(EST[4], survey.uncertainty, 0.001, adjust.uncertainty) 

# Ignore adjustment uncertainty
CI(EST[4], survey.uncertainty, model.uncertainty, 0.001) 

# Ignore both survey and model uncertainty
CI(EST[4], 0.001, 0.001, adjust.uncertainty) 

# Ignore both survey and adjustment uncertainty
CI(EST[4], 0.001, model.uncertainty, 0.001) 

# Ignore both model and adjustment uncertainty
CI(EST[4], survey.uncertainty, 0.001, 0.001)
```
