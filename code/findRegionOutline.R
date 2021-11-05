
#Based on Broenniman et al (2012) PCA method, adapted by Regan to measure niche NND and other dynamics
#Modified by HHakkinen 2017 to process species data. Points are read in, categorised, PCA space is created
#then save the outputs as raster for later summary, as well graphical output

#finds center of niche in climate space, measures locations, distance and native/natur regions and creates df output
#expansion points are compared against the center of the native niche and measured for direction and magnitude

rm(list=ls())

setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/")
setwd("/Users/henry_hakkinen/Documents/Research/RepoCode/nicheExpansion/")

#functions from the original broenniman paper, various versions of input
#source("niche.overlap.functions.R")
#source("occ.prep.functions.R")
source("code/functions/niche_dynamic_functions_Nov11th2016HH.R")
source("code/functions/PCA_functions.R")

library(dplyr)

library(biomod2)
library(ade4)
library(adehabitat)
library(raster)
library(rworldmap)
library(ecospat)

library(circular)
library(CircMLE)
#library(SDMTools)


#library(sp)
#library(gam)
#library(MASS)
#library(mvtnorm)
#library(gbm)
#library(dismo)
#library(rgeos)


#library(rgdal)

#library(cluster)

#library("easyGgplot2")

##########################################################
#Functions for data processing

## Define colours to be used in plots
ntv.col <- rgb(red=145,green=191,blue=219, max=255, alpha=150) ## alpha is for transparances, 0-255
usa.col <- rgb(red=252,green=141,blue=89, max=255, alpha=150)
overlap.col <- rgb(red=255,green=255,blue=191, max=255, alpha=150)
ntv.col.cont <- rgb(red=145,green=191,blue=219, max=255) 
usa.col.cont <- rgb(red=252,green=141,blue=89, max=255)
ntv.col.map <- rgb(red=145,green=191,blue=219, max=255) ## alpha is for transparances, 0-255
usa.col.map <- rgb(red=252,green=141,blue=89, max=255)
overlap.col.map <- c(0, ntv.col, usa.col, rgb(red=255,green=255,blue=191, max=255))
bkg.col <- rgb(red=225,green=225,blue=225, max=255)


###################################
#load files for processing

#load shapefile with zones for range classification
shape <- shapefile("RawData/BiogeographicZones/biogeographic_zonesV2.shp")
#low-res map for plotting purposes
newmap <- getMap(resolution = "coarse")  # different resolutions available
#plot(shape,col=as.factor(shape@data$Name),axes=T)

#load world raster as base for data to be overlaid, and for info on grid cell area (km squared)
world<-raster("RawData/GlobalGrid/Global_10min_grid.grd")


#this is 3 climate variables + NPP
#clim123_ref <- na.exclude(read.csv("RawData/WorldClim/bioclimNPP_10min_region.csv"))

clim123_ref <- na.exclude(read.csv("RawData/WorldClim/bioclim_4Var_10min_region.csv"))


head(clim123_ref)
clim123<-clim123_ref[,c("x","y", "bio5", "bio6", "bio12", "bio15")]

colnames(clim123)<-c("x","y","X1","X2","X3","X4")


#application of the function

sp_list<-unique(clim123_ref$region)
sp_name<-sp_list[1]

yooo<-mapply(niche_cal, sp_list)


#niche_cal(sp_name)

#the actual function to process the data, it's bad form to have it this way round, but saves on scrolling
niche_cal<-function(sp_name){
  
  
  print(sp_name)

  
  native<-clim123[clim123_ref$region==sp_name,]
  natur<-clim123
  
  #if(ncol(native)>2){native<-native[,1:2]}
  #if(ncol(natur)>2){natur<-natur[,1:2]}
  
  native<-native[!is.na(native$x)&!is.na(native$y),]
  natur<-natur[!is.na(natur$x)&!is.na(natur$y),]
  
  #colnames(native)<-c("x","y")
  #colnames(natur)<-c("x","y")
  
  
  #we comment out because we're just pre-desaggregated data
  #occ.sp_native<-ecospat.occ.desaggregation(df=native,colxy=1:2,min.dist=0.16666,plot=F) 
  #occ.sp_natur<-ecospat.occ.desaggregation(df=natur,colxy=1:2,min.dist=0.16666,plot=F) 
  
  occ.sp_native<-native
  occ.sp_natur<-natur
  
  #global climate is treated as equally available, NPP removed
  clim1<-clim123_ref[,c("x","y", "bio5", "bio6", "bio12", "bio15")]
  clim2<-clim123_ref[,c("x","y", "bio5", "bio6", "bio12", "bio15")] 
  
  
  #clim12 is a combination of clim1 (native) and clim2 (naturalised)
  clim12 <- rbind(clim1, clim2)
  

  
  #extract bioclim variables for species occurences
  
  occ.sp1<-native
  occ.sp2<-natur
  
  names(occ.sp1)<-names(clim1)
  names(occ.sp2)<-names(clim1)

      
      ############Below is Regan's adapated method.
      
      #################################################################################################
      ############################## ANALYSIS - selection of parameters ###############################
      #################################################################################################
      
      # selection of the type of analysis.
      # If PROJ =F, the models are calibrated on both ranges.
      # This is a hangover from more complex versions of this script.
      PROJ = F ## Keep this as F (False)
      
      # selection of variables to include in the analyses
      names(clim123)
      Xvar<-c(3:ncol(clim12)) # This means use columns 3 to 5 of the climate files. If you have and want to use more environmental variables, change this to (e.g.) Xvar<-c(3:7)
      nvar<-length(Xvar)
      
      #number of interation for the tests of equivalency and similarity, though this is not used in the current version
      iterations<-100
      
      #resolution of the gridding of the climate space
      R=100
      
      #################################################################################################
      ################### row weigthing and grouping factors for ade4 functions  ######################
      #################################################################################################
      
      
      
      # if PROJ = F
      row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1
      row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
      ##THE following is adjusted from regan's script since we are NOT combining (since environments are the same)
      row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))
      #row.w.occ.PROJT<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
      
      row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
      row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
      
      ##can adjust weighting here depending on climate availability
      row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
      #row.w.env.PROJT<-c(rep(1, nrow(clim1)),rep(0, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
      
      fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))
      
      # global dataset for the analysis and rows for each sub dataset
      data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
      
      row.clim1<-1:nrow(clim1)
      row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
      row.clim12<-1:(nrow(clim1)+nrow(clim2))
      row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
      row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))
      
      
      
      
      #################################################################################################
      #################################### PCA-ENV ####################################################
      #################################################################################################
      
      # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
      
      #fit of the analyse using occurences from both ranges	
      #names(data.env.occ)<-c("Tmax","Tmin","Precip","PrecipSeason")
      #data.env.occ<-data.env.occ[,c(1,3)]
      
      #pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
      
      
      #save global PCA for later use if you like
      save(pca.cal,file=paste0("./IntermediateOutput/findRegionOutline/PCA_contrib_", sp_name))
      #plot the loadings if needed
      ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
      
      
      
      
      # predict the scores on the axes
      scores.clim12<- pca.cal$li[row.clim12,]
      scores.clim1<- pca.cal$li[row.clim1,]
      scores.clim2<- pca.cal$li[row.clim2,]
      scores.sp1<- pca.cal$li[row.sp1,]
      scores.sp2<- pca.cal$li[row.sp2,]
      
      
      #grid.clim.NNU is slightly different to the main plotting functions later (ecospat.grid.clim.dyn), but only in the output. Calculatations are the same
      #return niche matrices both as a matrix and as a raster:
      z2<- grid.clim.NNU2(scores.clim12, scores.clim2, scores.sp2, R) ## Uses the PCA scores of the entire climate space, native climate space, and species distribution to calculate occurrence density on a grid. The entire climate space is only used to set the boundaries of the grid.
      z1<- grid.clim.NNU2(scores.clim12,scores.clim1,scores.sp1,R)
      
      plot(z1$zz)

      writeRaster(z1$zz, paste0("IntermediateOutput/findRegionOutline/4Var/",sp_name,"PCA.tif"),overwrite=T)
      
      
      
      #set threshold for 'occupied' climate
      thresh=0.00
      
      
      #select appropriate occurence density layer
      ntv.distn <- z1["z.uncor"][[1]]
      
      
      #ntv.distn <- z1$z.cor
      
      #remove occurrences that fall under the threshold of occupied climate
      #z.raw is z/Z, so scaled to availability of climate. remove any points from z.uncor that fall under this threshold
      ntv.distn.r.pa<-z1$z.raw
      ntv.distn.r<-ntv.distn*ntv.distn.r.pa
      ntv.distn.dens <- cellStats(ntv.distn.r, stat='sum')
      
      
      ### naturalised distribution
      natur.distn<-z2["z.uncor"][[1]]
      #natur.distn<-z2$z.cor
      
      #remove occurrences that fall under the threshold of occupied climate
      natur.distn.r.pa<-z2$z.raw
      natur.distn.r<-natur.distn*natur.distn.r.pa
      natur.distn.dens <- cellStats(natur.distn.r, stat='sum')
      
      
      #what bits overlap?
      ntv.overlap<-ntv.distn.r.pa+2*natur.distn.r.pa
      
      plotname=paste('IntermediateOutput/findRegionOutline/4Var/',sp_name,'_similarity.pdf', sep='')
      pdf(file=plotname, width=12, height=9)
      par(mfrow=c(2,2))
      #x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
      ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      ecospat.plot.niche(z2,title="PCA-env - naturalised niche",name.axis1="PC1",name.axis2="PC2")
      
      
      #plot.niche.dyn(z1,z2,0.25,int=2, title=paste("Overlap in native and naturalised climatic niches"))  
      plot(z1$zz,main=sp_name)
      
      plot(shape,main=paste("Native and naturalised distributions",sp_name))
      points(native$x,native$y,col="blue",pch=".",cex=4)
      #points(natur$x,natur$y,col="red",pch=".",cex=4)
      
      dev.off()
      
      #natur_PCA<-raster(paste("RawData/RegionPCA/3Var/",sp_name,"PCA.tif",sep=""))
      
      print("finished simulations")

}





