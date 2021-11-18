
##############################################
### COMPILE WORLDCLIM DATA ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 
# The following takes raw worldClim bioclim data and stacks them into a single CSV
# Input .bil files downloaded from WorldClim, stored in RawData/WorldClim/bio_10m_bil
# Process: stack and save as csv
# output: RawData/WorldClim/bioclim_10min_region.csv

### ###


rm(list=ls())

#set directory to current location of repo

repopath<-"PATH HERE"
setwd(paste0(repopath,"/RawData/WorldClim/bio_10m_bil"))

library(raster)
library(rworldmap)
library(ggplot2)
library(RColorBrewer)

source("../../../code/functions/miscFunctions.R")


########################
###stack the bioclim variables you want
#########################


#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
###BIO3 = Isothermality (BIO2/BIO7) (* 100)
###BIO4 = Temperature Seasonality (standard deviation *100)
###BIO5 = Max Temperature of Warmest Month
###BIO6 = Min Temperature of Coldest Month
####BIO7 = Temperature Annual Range (BIO5-BIO6)
#####BIO8 = Mean Temperature of Wettest Quarter
####BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
####BIO13 = Precipitation of Wettest Month
####BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
####BIO16 = Precipitation of Wettest Quarter
####BIO17 = Precipitation of Driest Quarter
####BIO18 = Precipitation of Warmest Quarter
####BIO19 = Precipitation of Coldest Quarter

#a<-raster("bio1.bil")
#b<-raster("bio2.bil")
#c<-raster("bio3.bil")
#d<-raster("bio4.bil")
e<-raster("bio5.bil")
f<-raster("bio6.bil")
#g<-raster("bio7.bil")
#h<-raster("bio8.bil")
#i<-raster("bio9.bil")
#j<-raster("bio10.bil")
#k<-raster("bio11.bil")
l<-raster("bio12.bil")
#m<-raster("bio13.bil")
#n<-raster("bio14.bil")
o<-raster("bio15.bil")
#p<-raster("bio16.bil")
#q<-raster("bio17.bil")
#r<-raster("bio18.bil")
#s<-raster("bio19.bil")

a1<-coordinates(e)
a2<-as.data.frame(e)
b2<-as.data.frame(f)
l2<-as.data.frame(l)
l3<-as.data.frame(o)

l.w <- stack(e,f,l,o)


#in test run keeping mean temp of coldest, mean temp of warmeest and total annual precipitation

clim_df<-cbind(a1,a2,b2,l2,l3)
head(clim_df)

clim_df2<-na.omit(clim_df)
dim(clim_df2)



##########################################
###for ease of use look up the geographic region and add that too
##########################################


newmap <- getMap(resolution = "coarse")  # different resolutions available

#ben holt's zoogeographic regions map
shape <- shapefile("../../BioGeographicZones/biogeographic_zonesV2.shp")
(shape@data$Name)
plot(shape,col=as.factor(shape@data$Name),axes=T)

#lookup what region is point of world clim data is in, and add a label
yoo<-clim_df2
clookup(yoo)

yoo$region<-rest$Name



########################
###plot to check results
#########################

palette(c("red","blue","pink","green","midnightblue", "yellow", "palegreen2", "orange","brown","grey","maroon"))

plot(yoo$x,yoo$y,col=as.factor(yoo$region),pch=".",cex=1)
yoo_na<-yoo[is.na(yoo$region),]
nrow(yoo_na)/nrow(yoo)

plot(yoo_na$x,yoo_na$y,col="Black",pch=".",cex=2)
unique(yoo$region)


#write out the completed file
numvar<-nlayers(l.w)
write.csv(yoo,paste0("../bioclim_",numvar,"Var_10min_region.csv",row.names=F))




