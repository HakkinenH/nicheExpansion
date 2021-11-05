
#Based on Broenniman et al (2012) PCA method, adapted by Regan to measure niche NND and other dynamics
#Modified by HHakkinen 2017 to process species data. Points are read in, categorised, PCA space is created
#then save the outputs as raster for later summary, as well graphical output

#finds center of niche in climate space, measures locations, distance and native/natur regions and creates df output
#cells that species have not filled are isolated and compare to the native niche center

rm(list=ls())

setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/")
#setwd("/Users/henry_hakkinen/Documents/Research/RepoCode/nicheExpansion/")

#functions from the original broenniman paper, various versions of input
#source("niche.overlap.functions.R")
#source("occ.prep.functions.R")
source("code/functions/niche_dynamic_functions_Nov11th2016HH.R")
source("code/functions/PCA_functions.R")

library(dplyr)

library(biomod2)
library(ade4)
#library(adehabitat)
library(raster)
library(rworldmap)

library(ecospat)
library(adehabitatHR)

library(circular)
library(CircMLE)
#library(SDMTools)
library(rgeos)


library(sp)
#library(gam)
#library(MASS)
#library(mvtnorm)
#library(gbm)
#library(dismo)


#library(rgdal)

#library(cluster)

#library("easyGgplot2")

##########################################################
#Functions for data processing



###################################
#load files for processing

#load shapefile with zones for range classification
shape <- shapefile("RawData/BiogeographicZones/biogeographic_zonesV2.shp")
#low-res map for plotting purposes
newmap <- getMap(resolution = "coarse")  # different resolutions available
#plot(shape,col=as.factor(shape@data$Name),axes=T)

#load world raster as base for data to be overlaid, and for info on grid cell area (km squared)
world<-raster("RawData/GlobalGrid/Global_10min_grid.grd")


#global climate space. This is basically the above .grd file but with region assigned to each cell (in csv form)
#it takes a while to calculate this everytime, so I just saved it with region assigned (check it matches shape)
#this is 3 climate variables
#clim123_ref <- na.exclude(read.csv("bioclim_10min_region.csv"))
#clim123<-clim123_ref[,-c(6)]
#colnames(clim123)<-c("x","y","X1","X2","X3")

#this is 3 climate variables + NPP
#clim123_ref <- na.exclude(read.csv("RawData/WorldClim/bioclimNPP_10min_region.csv"))

clim123_ref <- na.exclude(read.csv("RawData/WorldClim/bioclim_4Var_10min_region.csv"))

head(clim123_ref)
clim123<-clim123_ref[,-c(6,7,8)]

colnames(clim123)<-c("x","y","X1","X2","X3")

#previously we've saved a copy of the global PCA (slightly circular but there it is)
#open and extract direction of each component vector

load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecip")
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))
colnames(arr_tab)<-c("x","y")
#table with points which give the direction as compared to the origin (0,0)
#find the angle of direction for each PCA component
arr_tab$angle<-apply(arr_tab,1,rad.ang.center,x2=c(0,0))

plot(0,0,xlim=c(-1,1),ylim=c(-1,1))
points(arr_tab$x,arr_tab$y)
arrows.circular(arr_tab[4,3],col="red")



#list of all species for processing. This is the old total list
#plant_summary<-read.csv("species lists/plant_data_summary0607.csv")
#remove species we ruled out for whatever reason
#plant_list<-plant_summary[plant_summary$Verdict == "Yes",]
#this is the new summary list for informing analogue climate
plant_summary<-read.csv("IntermediateOutput/PCA_find_analogue_byRegion/plant_data_analogue05042019.csv",stringsAsFactors = F)


#########################
#load failures if necessary
########################
getwd()
faillist<-list.files(path="niche_shift_work/failure/",pattern="*.csv")
failsp<-c()

q<-faillist[7]

for(q in faillist){
  a<-unlist(strsplit(q, split="-"))
  
  a_sp<-a[1]
  reg<-a[2:length(a)]
  if(length(reg)>1){reg<-paste(reg,collapse = "-")}
  region<-unlist(strsplit(reg, split=".csv"))[1]
  failsp<-rbind(failsp, c(a_sp,region))
}

failsp<-data.frame(failsp)
#colnames(failsp)<-c("species_name", "region")


##############################

#blank dataframe that will be filled row by row for any species that fail
mergefail.df<-data.frame(sp_name=integer(),native_removed=integer(),naturalised_removed=integer())


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



#when running calculatations we get all of the native/naturalised distribution, this is scaled from 0:1
#there is a choice here. z.uncor is normal density, z.cor is corrected for environment availability
#set this below using cor_sel
cor_sel<-"z.uncor"
thresh_exp_per<-0.1
perc<-"ninetypercentile"


#application of the function
#plant_summary<-plant_summary[1:20,]
plant_summary<-plant_summary[1:nrow(plant_summary),]
head(plant_summary)

head(plant_summary)
yooo<-mapply(niche_cal, plant_summary$species_name,plant_summary$region)


sp_name<-"Agrostis gigantea"
region<-"Sino-Japanese"



#the actual function to process the data, it's bad form to have it this way round, but saves on scrolling
niche_cal<-function(sp_name,region){
  
  
  print(sp_name)
  print(region)
  
  #species lists/plant_dist_data_0802: these are original raw files but we use the deaggregated points for speed
  #species lists/plant_dist_data_0802/desaggregated_data/ is the global disaggregated data
  
  #this is out current dataset
  #deaggregaated, and filtered so only analogue native/naturalised occurrence data is included
  #native (all regions) vs naturalised (the region we are currently looking at only)
  file_natur = paste("IntermediateOutput/PCA_find_analogue_byregion/",region,"/",sp_name,region,"_natur.csv",sep="")
  file_native = paste("IntermediateOutput/PCA_find_analogue_byregion/",region,"/",sp_name,region,"_native.csv",sep="")
  
  
  
  native<-read.csv(file_native)
  natur<-read.csv(file_natur)
  if(ncol(native)>2){native<-native[,1:2]}
  if(ncol(natur)>2){natur<-natur[,1:2]}
  
  native<-native[!is.na(native$x)&!is.na(native$y),]
  natur<-natur[!is.na(natur$x)&!is.na(natur$y),]
  
  colnames(native)<-c("x","y")
  colnames(natur)<-c("x","y")
  
  
  #we comment out because we're just pre-desaggregated data
  #occ.sp_native<-ecospat.occ.desaggregation(df=native,colxy=1:2,min.dist=0.16666,plot=F) 
  #occ.sp_natur<-ecospat.occ.desaggregation(df=natur,colxy=1:2,min.dist=0.16666,plot=F) 
  
  occ.sp_native<-native
  occ.sp_natur<-natur
  
  
  #main 3 climatic variables
  clim_sel<-c(1:5)
  
  
  #global climate is treated as equally available, NPP removed
  clim1<-clim123_ref[,c(1:5)]
  clim2<-clim123_ref[,c(1:5)] 
  
  
  #clim12 is a combination of clim1 (native) and clim2 (naturalised)
  clim12 <- rbind(clim1, clim2)
  
  #filter out occurences with fewer than 5 occurences, PCA cannot be built on fewer.
  if (nrow(occ.sp_native)<=5 | nrow(occ.sp_natur)<=5){
    print ("not enough known occurrences to calculate niche")
  }else{
    
    #extract bioclim variables for species occurences
    
    occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp_native,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=0.16666))
    occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp_natur,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=0.16666))
    
    
    
    if (nrow(occ.sp1)<5 | nrow(occ.sp2)<5){
      print ("Too many data point removals!")
      mergefail.df[nrow(mergefail.df)+1,1]<-c(sp_name)
      print ("not enough known occurrences to calculate niche")
      
    }
    else{
      
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
      
      #pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
      
      #save global PCA for later use if you like
      #save(pca.cal,file="niche_shift_work/PCA_contrib")
      
      
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
      
      
      
      # #test of equivalency taken out due to time delay
      #a<-ecospat.niche.equivalency.test(z1_alt,z2_alt,rep=100)# test of niche equivalency and similarity according to Warren et al. 2008
      #b<-ecospat.niche.similarity.test(z1,z2,rep=100) ## can take some time
      #b2<-ecospat.niche.similarity.test(z2,z1,rep=100) ## can take some time
      
      #we want some expansion info
      dyn.100<-dynamic.index(z1,z2,thresh=0) ## dynamics in all analagous climate space
      exp_sum<-dyn.100[[2]][1]
      
      nnddyn.100<-dynamic.index.nnd(z1,z2,thresh=0) ## dynamics in all analagous climate space
      
      
      print("finished simulations")
      
      
      #take D result (overlap) just in case we want it
      D_result<-round(as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]),3)
      #plot equivalency values, not used currently			
      
      
      #set threshold for 'occupied' climate
      thresh=0.00
      
      
      #select appropriate occurence density layer
      ntv.distn <- z1[cor_sel][[1]]
      
      
      #ntv.distn <- z1$z.cor
      
      #remove occurrences that fall under the threshold of occupied climate
      #z.raw is z/Z, so scaled to availability of climate. remove any points from z.uncor that fall under this threshold
      ntv.distn.r.pa<-z1$z.raw>thresh
      ntv.distn.r<-ntv.distn*ntv.distn.r.pa
      ntv.distn.dens <- cellStats(ntv.distn.r, stat='sum')
      
      
      ### naturalised distribution
      natur.distn<-z2[cor_sel][[1]]
      #natur.distn<-z2$z.cor
      
      #remove occurrences that fall under the threshold of occupied climate
      natur.distn.r.pa<-z2$z.raw>thresh
      natur.distn.r<-natur.distn*natur.distn.r.pa
      natur.distn.dens <- cellStats(natur.distn.r, stat='sum')
      
      
      #what bits overlap?
      ntv.overlap<-ntv.distn.r.pa+2*natur.distn.r.pa
      
      #plot(ntv.overlap)
      stab<-ntv.overlap==3
      exp<-ntv.overlap==2
      res<-ntv.overlap==1
      
      ## Climate in naturalised range
      z2.r<-z2$Z
      
      ## Climate in native range
      z1.r<-z1$Z
      
      ## Overlap between USA and Europe - all climate space
      z1.r100 <- z1.r > 0
      z2.r100 <- z2.r > 0
      overlap100 <- z1.r100 * z2.r100
      
      #in order to calculate what points lie where, what do I need?
      #outline of native occupied niche (100% analogue)
      ntv.distn.r.pa100<-ntv.distn.r.pa*overlap100
      
      #outline of naturalised occuoied niche (100% analogue)
      natur.distn.r.pa100<-natur.distn.r.pa*overlap100
      
      #stab, exp and res in analogue space
      stab100<-stab*overlap100
      exp100<-exp*overlap100
      res100<-res*overlap100
      
      
      ###we have found all areas of overlapping niche space so now we can start to answer questions
      
      #########################1) do species unfill their niche?###################
      
      
      #find all naturalised climate cells that lie in native niche
      
      natur_PCA<-raster(paste("IntermediateOutput/findRegionOutline/3Var/",region,"PCA.tif",sep=""))
      
      #correct extent problem if there is one. Not best practise
      extent(natur_PCA)<-extent(ntv.distn.r.pa100)
      
      #we need information on the center of the native and naturalised niche, what are they like?
      #find center of native niche
      #find center of mass
      ntv_df<-as.data.frame(coordinates(ntv.distn.r))
      ntv_df$layer<-as.data.frame(ntv.distn.r$layer)$layer
      if(nrow(ntv_df[which(ntv_df$layer==0),])>0){ntv_df[which(ntv_df$layer==0),]<-NA}
      ntv_df<-ntv_df[!(is.na(ntv_df$layer)),]
      
      native_center<-COGravity(x=ntv_df$x, y=ntv_df$y, wt=ntv_df$layer)[c(1,3)]
      native_center<-data.frame(x=native_center[1],y=native_center[2])
      
      natur_df<-as.data.frame(coordinates(natur.distn.r))
      natur_df$layer<-as.data.frame(natur.distn.r$layer)$layer
      if(nrow(natur_df[which(natur_df$layer==0),])>0){natur_df[which(natur_df$layer==0),]<-NA}
      natur_df<-natur_df[!(is.na(natur_df$layer)),]
      
      natur_center<-COGravity(x=natur_df$x, y=natur_df$y, wt=natur_df$layer)[c(1,3)]
      natur_center<-data.frame(x=natur_center[1],y=natur_center[2])
      
      
      #we also need info on the potential naturalised niche. Where is the center?
      #take naturalised climate and crop by native niche
      natur_poten <-natur_PCA*ntv.distn.r.pa100


      
      natur_poten_df<-as.data.frame(coordinates(natur_poten))
      natur_poten_df$layer<-as.data.frame(natur_poten$layer)$layer
      if(nrow(natur_poten_df[which(natur_poten_df$layer==0),])>0){natur_poten_df[which(natur_poten_df$layer==0),]<-NA}
      natur_poten_df<-natur_poten_df[!(is.na(natur_poten_df$layer)),]
      
      natur_poten_center<-COGravity(x=natur_poten_df$x, y=natur_poten_df$y, wt=natur_poten_df$layer)[c(1,3)]
      natur_poten_center<-data.frame(x=natur_poten_center[1],y=natur_poten_center[2])
      
      
      
      #an alternative method if more than 1 medoid
      #native_center<-as.data.frame(pam(scores.sp1, 1)$medoids)
      #natur_center<-as.data.frame(pam(scores.sp2, 1)$medoids)
      
      #what is the climate like in native center (center of expansion)
      native_clim_pca<-cbind(clim1,scores.clim1)
      natur_clim_pca<-cbind(clim2,scores.clim2)
      
      
      
      # a<-ggplot(native_clim_pca, aes(Axis1, Axis2)) +
      #   ggtitle("NPP")+
      #   geom_point(aes(colour = NPP))+
      #   scale_colour_gradient(low="blue", high="red")
      # 
      # b<-ggplot(native_clim_pca, aes(Axis1, Axis2)) +
      #   ggtitle("Tmin")+
      #   geom_point(aes(colour = bio6))+
      #   scale_colour_gradient(low="blue", high="red")
      # 
      # c<-ggplot(native_clim_pca, aes(Axis1, Axis2)) +
      #   ggtitle("TMax")+
      #   geom_point(aes(colour = bio5))+
      #   scale_colour_gradient(low="blue", high="red")
      # 
      # d<-ggplot(native_clim_pca, aes(Axis1, Axis2)) +
      #   ggtitle("Precip")+
      #   geom_point(aes(colour = bio12))+
      #   scale_colour_gradient(low="blue", high="red")
      # 
      # 
      # png(file="niche_shift_work/Figures/global_pca_diag.png", width=900, height=700)
      # ggplot2.multiplot(b,c,d, cols=2)
      # dev.off()
      
      # png(file="niche_shift_work/Figures/global_pca_arrows.png", width=900, height=700)
      # par(mfrow=c(2,2))
      # ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      # arrows.circular(npp_dir, shrink=3, length=0.1, x0=native_center$x, y0=native_center$y,col="red")
      # ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      # arrows.circular(tmin_dir, shrink=3, length=0.1, x0=native_center$x, y0=native_center$y,col="red")
      # ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      # arrows.circular(tmax_dir, shrink=3, length=0.1, x0=native_center$x, y0=native_center$y,col="red")
      # ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      # arrows.circular(precip_dir, shrink=3, length=0.1, x0=native_center$x, y0=native_center$y,col="red")
      # dev.off()
      
      
      native_cinfo<-near_point2(native_center,native_clim_pca)
      native_cinfo
      
      circ_model<-NA;circ_rstat<-NA;circ_pvalue<-NA
      exp.mean.dir1<-NA;exp.median.dir1<-NA;exp.ten.dir1<-NA;exp.ninety.dir1<-NA;exp.max.dir1<-NA;
      exp.mean.dir2<-NA;exp.median.dir2<-NA;exp.ten.dir2<-NA;exp.ninety.dir2<-NA;exp.max.dir2<-NA;
      CV<-NA;dir1<-NA;dir2<-NA;circ_altmodel<-NA; conc_dir1<-NA;conc_dir2<-NA;
      tmax_pro<-NA;tmin_pro<-NA;precip_pro<-NA;cold_pro<-NA;dry_pro<-NA;
      
      
      #####################key calculation area #####################
      
      #we need some info on unfilling points
      #slightly different now as we want POTENTIAL cells
      #find all naturalised climate cells that lie in native niche
      
      natur_PCA<-raster(paste("IntermediateOutput/findRegionOutline/3Var/",region,"PCA.tif",sep=""))
      natur_PCA[natur_PCA==0]<-NA
      extent(natur_PCA)<-extent(ntv.distn.r.pa100)
      
      
      
      # #we need some info on stability points
      # #find all stablity cells and their associated density
      stab.distn.r<-stab100*z2$zz
      stab.distn.rKEEP<-stab.distn.r
      #correct for climate availability in the naturalised region
      stab.distn.r<-stab.distn.r/natur_PCA
      
      #rescale 0-1 to make it comparable for other species
      stab.distn.r<-stab.distn.r/cellStats(stab.distn.r,stat="max")
      
      
      #rescale expansion for easier analysis (0-1000)
      stab.distn.r[stab.distn.r==0]<-NA
      
      stab.distn.r2<-stab.distn.r*1000
      stab.distn.r2<-round(stab.distn.r2)
      stab.distn.r2[stab.distn.r2<1]<-NA
      
      
      
      #save information as a dataframe
      stab_df<-as.data.frame(coordinates(stab.distn.r2))
      stab_df$layer_orig<-as.data.frame(stab.distn.r$layer)$layer
      stab_df$layer<-as.data.frame(stab.distn.r2$layer)$layer
      
      stab_df<-stab_df[!is.na(stab_df$layer),]
      
      
      #calculate distance to each stability cell
      stab_df$exp_dist<-unlist(apply(stab_df[,1:2], 1, euc.dist.center, x2=natur_poten_center))
      
      #calculate angle to each cell,
      stab_df$exp_dir<-unlist(apply(stab_df[,1:2], 1, rad.ang.center, x2=natur_poten_center))
      
      ######################################################
      ##sneaky second part with uncorrected numbers
      #######################################################
      stab.distn.r0<-stab.distn.rKEEP
      stab.distn.r0[stab.distn.r0==0]<-NA
      stab.distn.r0<-stab.distn.r0*1000
      stab.distn.r0<-round(stab.distn.r0)
      stab.distn.r0[stab.distn.r0<1]<-NA
      
      stab_raw<-as.data.frame(coordinates(stab.distn.r0))
      stab_raw$layer<-as.data.frame(stab.distn.r0$layer)$layer

      stab_raw<-stab_raw[!is.na(stab_raw$layer),]
      
      
      #we've already constructed the potential naturalised niche above
      natur_poten0<-natur_poten
      natur_poten0[natur_poten0==0]<-NA
      natur_poten0<-natur_poten0*1000
      natur_poten0<-round(natur_poten0)
      natur_poten0[natur_poten0<1]<-NA
      
      natur_pot<-as.data.frame(coordinates(natur_poten0))
      natur_pot$layer<-as.data.frame(natur_poten0$layer)$layer
      natur_pot<-natur_pot[!is.na(natur_pot$layer),]

      #and we have the center

      stab_raw$exp_dir<-unlist(apply(stab_raw[,1:2], 1, rad.ang.center, x2=natur_poten_center))
      natur_pot$exp_dir<-unlist(apply(natur_pot[,1:2], 1, rad.ang.center, x2=natur_poten_center))
      
      
      #specify directions of precip/tmin/max
      precip_dir<-arr_tab["bio12","angle"]
      tmax_dir<-arr_tab["bio5","angle"]
      
      
      stab_precip<-find_slice(stab_raw,precip_dir)
      stab_tmax<-find_slice(stab_raw,tmax_dir)
      stab_dry<-find_slice(stab_raw,precip_dir+pi)
      stab_cold<-find_slice(stab_raw,tmax_dir-pi)
      
      pot_precip<-find_slice(natur_pot,precip_dir)
      pot_tmax<-find_slice(natur_pot,tmax_dir)
      pot_dry<-find_slice(natur_pot,precip_dir+pi)
      pot_cold<-find_slice(natur_pot,tmax_dir-pi)
      

      res_df<-data.frame("sp_name"=sp_name,
                        "region" = region,
                        "precipOccdens" = sum(stab_precip$layer),
                        "precipPotdens" = sum(pot_precip$layer),
                        "precipPro" = sum(stab_precip$layer)/sum(pot_precip$layer),
                        
                        "tmaxOccdens"= sum(stab_tmax$layer),
                        "tmaxPotdens" = sum(pot_tmax$layer),
                        "tmaxPro"= sum(stab_tmax$layer)/sum(pot_tmax$layer),
                        
                        "dryOccdens"= sum(stab_cold$layer),
                        "dryPotdens" = sum(pot_cold$layer),
                        "dryPro"= sum(stab_cold$layer)/sum(pot_cold$layer),
                        
                        "coldOccdens"= sum(stab_dry$layer),
                        "coldPotdens" = sum(pot_dry$layer),
                        "coldPro"= sum(stab_dry$layer)/sum(pot_dry$layer))
      write.csv(res_df, paste0("IntermediateOutput/plant_PCA_nichefilling/findFilling/",sp_name,'_',region,'_fillData.csv'), row.names = F)
      
      
      
      precip_dir_plot<-as.circular(precip_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      tmax_dir_plot<-as.circular(tmax_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      
      dry_dir_plot<-as.circular(precip_dir+pi,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      cold_dir_plot<-as.circular(tmax_dir-pi,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      
    
      plotname=paste('IntermediateOutput/plant_PCA_nichefilling/findFilling/',sp_name,'_',region,'_shiftpie_exp.jpg', sep='')
      png(file=plotname, width=900, height=700)
      
      #let's plot
      par(mfrow=c(2,2))
      
      
      #precip
      plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Precip")
      points(stab_raw$x,stab_raw$y,col="blue")
      
      
      points(natur_poten_center$x,natur_poten_center$y,col="red")
      

      arrows.circular(precip_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(precip_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(precip_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      
      points(stab_precip$x,stab_precip$y,col="red")
      text(2,2,paste("proportion:", round(sum(stab_precip$layer)/sum(pot_precip$layer), digits=2)))
      
      
      
      #tmax
      plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Tmax")
      points(stab_raw$x,stab_raw$y,col="blue")
      
      points(natur_poten_center$x,natur_poten_center$y,col="red")

      arrows.circular(tmax_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(tmax_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(tmax_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      
      points(stab_tmax$x,stab_tmax$y,col="red")
      text(2,2,paste("proportion:", round(sum(stab_tmax$layer)/sum(pot_tmax$layer), digits=2)))
      
      
      #towards DRY
      plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="DRY")
      points(stab_raw$x,stab_raw$y,col="blue")
      
      points(natur_poten_center$x,natur_poten_center$y,col="red")

      arrows.circular(dry_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(dry_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(dry_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      
      points(stab_dry$x,stab_dry$y,col="red")
      text(2,2,paste("proportion:", round(sum(stab_dry$layer)/sum(pot_dry$layer), digits=2)))
      
      
      #towards COLD
      plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="COLD")
      points(stab_raw$x,stab_raw$y,col="blue")
      
      points(natur_poten_center$x,natur_poten_center$y,col="red")

      arrows.circular(cold_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(cold_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      arrows.circular(cold_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      
      points(stab_cold$x,stab_cold$y,col="red")
      text(2,2,paste("proportion:", round(sum(stab_cold$layer)/sum(pot_cold$layer), digits=2)))
      

      dev.off()
      # 
      # #############################################
      # 
      # 
      # #specify directions of precip/tmin/max
      # precip_dir<-arr_tab["bio12","angle"]
      # tmin_dir<-arr_tab["bio6","angle"]
      # tmax_dir<-arr_tab["bio5","angle"]
      # #npp_dir<-arr_tab["NPP","angle"]
      # 
      # arr_tab
      # 
      # 
      # #repeat for stability
      # if(nrow(stab_df)>0){
      #   
      #   ##############2) is there strong expansion in one (or two) directions? Or is it in every direction?##############
      #   #circ_mle comparing the central point of naturalised density with every gridcell of expansion (can I weight by density of gridcell?)
      #   #what model describes it best? Null is every direction at once by not very much
      #   #can we categorise it roughly? 1) no expansion 2) unidirectional expansion 3) bidirectional expansion 4)multidirectional exp
      #   
      #   #calculate distance to each stability cell
      #   stab_df$exp_dist<-unlist(apply(stab_df[,1:2], 1, euc.dist.center, x2=natur_poten_center))
      #   
      #   #calculate angle to each cell, 
      #   stab_df$exp_dir<-unlist(apply(stab_df[,1:2], 1, rad.ang.center, x2=natur_poten_center))
      #   
      #   
      #   
      #   #calculate circular variance http://resources.esri.com/help/9.3/arcgisdesktop/com/gp_toolref/spatial_statistics_tools/how_linear_directional_mean_spatial_statistics_works.htm
      #   CV<-1-(sqrt((sum(sin(stab_df$exp_dir))^2)+(sum(cos(stab_df$exp_dir))^2))/nrow(stab_df))
      #   
      #   
      #   #this is a nice summary table but in order to compensate for varying density of each expansion point
      #   #we are doing a weighting process. Each density point has been rescaled to have a minimum of one
      #   #we then multiply the angle by the weighted number
      #   
      #   exp_dir<-rep(stab_df$exp_dir,stab_df$layer)
      #   #exp_dir<-unfill_df$exp_dir
      #   
      #   
      #   
      #   
      #   circ_res<-circ_mle(exp_dir)
      #   
      #   #what is the best model to describe expansion?
      #   
      #   circ_model<-circ_res$bestmodel
      #   circ_rstat<-circ_res$rt[1]
      #   circ_pvalue<-circ_res$rt[2]
      #   
      #   dAIC<-circ_res$results$deltaAIC
      #   
      #   #select names of other viable models
      #   dAIC2<-rownames(circ_res$results[which(dAIC>0&dAIC<10),])
      #   
      #   if(length(dAIC2)>0){circ_altmodel<-paste(dAIC2,collapse=";")}
      #   
      #   dir1<-circ_res$results[1,2]
      #   dir2<-circ_res$results[1,5]
      #   conc_dir1<-circ_res$results[1,3]
      #   conc_dir2<-circ_res$results[1,6]
      #   
      #   ntv_plot<-ntv.distn.r
      #   ntv_plot[ntv_plot==0]<-NA
      #   
      #   plotname=paste('IntermediateOutput/plant_PCA_nichefilling/bucket_direction/',sp_name,'_',region,'_shift.jpg', sep='')
      #   png(file=plotname, width=900, height=700)
      #   
      #   par(mfrow=c(1,2))
      #   plot(res100,legend=F,main=sp_name,col=c("grey","blue"))
      #   #plot(exp.distn.r2,add=T)
      #   plot(stab.distn.r,add=T)
      #   
      #   points(natur_poten_center,col="red")
      #   arrows(x1=stab_df$x,y1=stab_df$y,
      #          x0=natur_poten_center$x,y0=natur_poten_center$y,
      #          col="red",length=0.05)
      #   
      #   plot_circMLE(exp_dir, circ_res)
      #   
      #   mtext(circ_model)
      #   dev.off()
      #   
      #   #plot breakdown to look at
      #   plotname=paste('IntermediateOutput/plant_PCA_nichefilling/circle_breakdown/',sp_name,'_',region,'.jpg', sep='')
      #   png(file=plotname, width=1000, height=1000, res=75)
      #   
      #   par(mfrow=c(3,4))
      #   
      #   model.names<-row.names(circ_res$results)
      #   for(i in model.names){
      #     plot_circMLE(exp_dir, circ_res, model = i)
      #     mtext(i)
      #   }
      #   
      #   dev.off()
      #   
      #   
      #   
      #   ####################################
      #   #DIRECTION CUT DOWN
      #   ###################################
      #   
      #   ###we have dir1 (and maybe dir2) which are directions of expansion
      #   #we have unfill_df, find all points that are within 45
      #   ###what is the mean expansion distance in this direction?
      #   
      #   #made function for this because we're going to do this several times
      #   #find all points that lie in a 90 degree slice around our primary direction of expansion
      #   exp_df_dir1<-find_slice(stab_df,dir1)
      #   
      #   
      #   ####calculate mean,median,max and 5% 95% of expansion for points within this direction
      #   ####is there any reason to want to know overall mean? 
      #   #this is the length of the line
      #   exp.mean.dir1<-mean(exp_df_dir1$exp_dist)
      #   exp.median.dir1<-median(exp_df_dir1$exp_dist)
      #   exp.ten.dir1<-quantile(exp_df_dir1$exp_dist,0.1)
      #   exp.ninety.dir1<-quantile(exp_df_dir1$exp_dist,0.9)
      #   exp.max.dir1<-max(exp_df_dir1$exp_dist)
      #   if(exp.max.dir1==-Inf){
      #     
      #     print("oh noooooo")
      #     resu<-c(sp_name,"dir1",dir1,exp.max.dir1,nrow(exp_df_dir1))
      #     write.csv(resu,paste("U:/Data/niche_shift_work/failure/",sp_name,"-",region,".csv",sep=""))
      #     
      #   }
      #   
      #   #if direction 2 exists in model, extract the same distance info
      #   if(!is.na(dir2)){
      #     
      #     #find all points that lie in a 90 degree slice around our secondary direction of expansion (if it exists)
      #     exp_df_dir2<-find_slice(stab_df,dir2)
      #     
      #     summary(exp_df_dir2)
      #     
      #     ####calculate mean,median,max and 5% 95% of expansion for points within this direction
      #     ####is there any reason to want to know overall mean? 
      #     #this is the length of the line
      #     exp.mean.dir2<-mean(exp_df_dir2$exp_dist)
      #     exp.median.dir2<-median(exp_df_dir2$exp_dist)
      #     exp.ten.dir2<-quantile(exp_df_dir2$exp_dist,0.1)
      #     exp.ninety.dir2<-quantile(exp_df_dir2$exp_dist,0.9)
      #     exp.max.dir2<-max(exp_df_dir2$exp_dist)
      #     
      #     if(exp.max.dir2==-Inf){
      #       
      #       print("oh noooooo")
      #       resu<-c(sp_name,"dir2",dir2,exp.max.dir2,nrow(exp_df_dir2))
      #       write.csv(resu,paste("U:/Data/niche_shift_work/failure/",sp_name,region,".csv",sep=""))
      #       
      #     }
      #     
      #   }
      #   #find all points that lie in a 90 degree slice around Precip vector (increasing)
      #   
      #   precip_df<-find_slice(stab_df,precip_dir)
      #   #what proportion of exp density is in this slice?
      #   precip_pro<-sum(precip_df$layer_orig)/sum(stab_df$layer_orig)
      #   
      #   #find all points that lie in a 90 degree slice around tmin vector (increasing)
      #   
      #   tmin_df<-find_slice(stab_df,tmin_dir)
      #   #what proportion of exp density is in this slice?
      #   tmin_pro<-sum(tmin_df$layer_orig)/sum(stab_df$layer_orig)
      #   
      #   #find all points that lie in a 90 degree slice around tmax vector (increasing)
      #   
      #   tmax_df<-find_slice(stab_df,tmax_dir)
      #   #what proportion of exp density is in this slice?
      #   tmax_pro<-sum(tmax_df$layer_orig)/sum(stab_df$layer_orig)
      #   
      #   
      #   #find all points that lie in a 90 degree slice around DECREASING precip (inverse)
      #   dry_df<-find_slice(stab_df,precip_dir+pi)
      #   #what proportion of exp density is in this slice?
      #   dry_pro<-sum(dry_df$layer_orig)/sum(stab_df$layer_orig)
      #   
      #   #find all points that lie in a 90 degree slice around DECREASING tmax (inverse)
      #   cold_df<-find_slice(stab_df,tmax_dir-pi)
      #   #what proportion of exp density is in this slice?
      #   cold_pro<-sum(cold_df$layer_orig)/sum(stab_df$layer_orig)
      #   
      #   
      #   
      #   #SECOND ANALYSIS
      #   
      #   #for the precipitation df, find which bits of niche filling lie around the vector of increasing precipitation
      #   precip_df<-find_slice(stab_df,precip_dir)
      #   
      #   #save the summed density of occupied wetter climate
      #   sum(precip_df$layer_orig)
      #   
      #   #find which parts of the naturalised niche fall around this vector as well
      #   
      #   
      #   #find all points that lie in a 90 degree slice around npp vector (deprecated)
      #   
      #   #npp_df<-find_slice(stab_df,npp_dir)
      #   #what proportion of exp density is in this slice?
      #   #npp_pro<-sum(npp_df$layer_orig)/sum(unfill_df$layer_orig)
      #   
      #   #plot that stuff as a summary
      #   
      #   dir1_plot<-as.circular(dir1,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   
      #   precip_dir_plot<-as.circular(precip_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   tmin_dir_plot<-as.circular(tmin_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   tmax_dir_plot<-as.circular(tmax_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   
      #   dry_dir_plot<-as.circular(precip_dir+pi,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   cold_dir_plot<-as.circular(tmax_dir-pi,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   
      #   
      #   #npp_dir_plot<-as.circular(npp_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
      #   
      #   
      #   
      #   
      #   
      #   
      #   plotname=paste('IntermediateOutput/plant_PCA_nichefilling/bucket_direction/',sp_name,'_',region,'_shiftpie_exp.jpg', sep='')
      #   png(file=plotname, width=900, height=700)
      #   
      #   #let's plot
      #   par(mfrow=c(3,2))
      #   
      #   
      #   plot(natur_poten)
      #   points(natur_poten_center,col="red")
      #   
      #   #precip
      #   plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Precip")
      #   points(stab_df$x,stab_df$y,col="blue")
      #   
      #   
      #   points(natur_poten_center$x,natur_poten_center$y,col="red")
      #   
      #   arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   arrows.circular(precip_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(precip_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(precip_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   points(precip_df$x,precip_df$y,col="red")
      #   text(2,2,paste("proportion:", round(precip_pro, digits=2)))
      #   
      #   #tmin
      #   plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Tmin")
      #   points(stab_df$x,stab_df$y,col="blue")
      #   
      #   points(natur_poten_center$x,natur_poten_center$y,col="red")
      #   arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   arrows.circular(tmin_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(tmin_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(tmin_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   points(tmin_df$x,tmin_df$y,col="red")
      #   text(2,2,paste("proportion:", round(tmin_pro, digits=2)))
      #   
      #   #tmax
      #   plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Tmax")
      #   points(stab_df$x,stab_df$y,col="blue")
      #   
      #   points(natur_poten_center$x,natur_poten_center$y,col="red")
      #   arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   arrows.circular(tmax_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(tmax_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(tmax_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   points(tmax_df$x,tmax_df$y,col="red")
      #   text(2,2,paste("proportion:", round(tmax_pro, digits=2)))
      #   
      #   
      #   #towards DRY
      #   plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="DRY")
      #   points(stab_df$x,stab_df$y,col="blue")
      #   
      #   points(natur_poten_center$x,natur_poten_center$y,col="red")
      #   arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   arrows.circular(dry_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(dry_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(dry_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   points(dry_df$x,dry_df$y,col="red")
      #   text(2,2,paste("proportion:", round(dry_pro, digits=2)))
      #   
      #   
      #   #towards COLD
      #   plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="COLD")
      #   points(stab_df$x,stab_df$y,col="blue")
      #   
      #   points(natur_poten_center$x,natur_poten_center$y,col="red")
      #   arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   arrows.circular(cold_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(cold_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   arrows.circular(cold_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   points(cold_df$x,cold_df$y,col="red")
      #   text(2,2,paste("proportion:", round(cold_pro, digits=2)))
      #   
      #   
      #   #NPP - removed
      #   #plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="NPP")
      #   #points(unfill_df$x,unfill_df$y,col="blue")
      #   
      #   #points(natur_poten_center$x,natur_poten_center$y,col="red")
      #   #arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   #arrows.circular(npp_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   #arrows.circular(npp_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   #arrows.circular(npp_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
      #   
      #   #points(npp_df$x,npp_df$y,col="red")
      #   #text(2,2,paste("proportion:", round(npp_pro, digits=2)))
      #   
      #   dev.off()
      #   
      #   writeRaster(stab.distn.r, paste("IntermediateOutput/plant_PCA_nichefilling/shift_PCArasters/",sp_name,'_',region,"_PCAraster.tif",sep=""), options= c("COMPRESSION=LZW","INTERLEAVE=PIXEL"), overwrite=TRUE)
      #   
      # }
      
      ######show plots to summarise these points
      plotname=paste('IntermediateOutput/plant_PCA_nichefilling/shift_PCArasters/',sp_name,'_',region,'_similarity.pdf', sep='')
      pdf(file=plotname, width=9, height=12)
      par(mfrow=c(2,2))
      #x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
      ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      ecospat.plot.niche(z2,title="PCA-env - naturalised niche",name.axis1="PC1",name.axis2="PC2")
      
      
      
      #plot.niche.dyn(z1,z2,0.25,int=2, title=paste("Overlap in native and naturalised climatic niches"))  
      plot(ntv.overlap)
      plot(shape,main=paste("Native and naturalised distributions",region))
      points(native$x,native$y,col="blue",pch=".",cex=4)
      points(natur$x,natur$y,col="red",pch=".",cex=4)
      
      dev.off()
      
      #there are a lot of measurements here but for report here are the major conclusion:
      #level of expansion (summed density in expansion space: expansion
      #occupied gridcells with expansion: Occupied_expcells
      #magnitude of niche expansion (average length): exp.dis
      #direction of niche expansion: exp_rad
      #variance of niche expansion: CV
      
      scores.sp1$range<-"native"
      scores.sp2$range<-"natur"
      
      
      
      
      #test<-raster(paste("niche_shift_work/shift_PCArasters/",sp_name,"_PCAraster.tif",sep=""))
      
      
      # 
      # D.df <- c(sp_name,region,nrow(occ.sp1),nrow(occ.sp2),D_result,dyn.100$dynamic.index.w, nnddyn.100$nnd,
      #           native_center,natur_center,natur_poten_center,
      #           native_cinfo[c("bio5","bio6","bio12")],
      #           circ_model,circ_rstat,circ_pvalue,
      #           exp.mean.dir1,exp.median.dir1,exp.ten.dir1,exp.ninety.dir1,exp.max.dir1,
      #           exp.mean.dir2,exp.median.dir2,exp.ten.dir2,exp.ninety.dir2,exp.max.dir2,
      #           tmax_pro,tmin_pro,precip_pro,cold_pro,dry_pro,
      #           CV,dir1,dir2,conc_dir1,conc_dir2,circ_altmodel)
      # 
      # 
      # names(D.df) <- c('species name', 'region', 'native_occurences','naturalised_occurences',
      #                  'Dstat_overlap', names(dyn.100$dynamic.index.w), names(nnddyn.100$nnd),
      #                  'native_center_x','native_center_y','natur_center_x','natur_center_y',"natur_poten_x","natur_poten_y",
      #                  'native_center_tmax','native_center_tmin','native_center_precip',
      #                  'circular_model','circ_rstat','circ_pstat',
      #                  'mean_exp_dist1','median_exp_dist1','ten_exp_dist1','ninety_exp_dist1','max_exp_dist1',
      #                  'mean_exp_dist2','median_exp_dist2','ten_exp_dist2','ninety_exp_dist2','max_exp_dist2',
      #                  'proportion_tmax','proportion_tmin','proportion_precip','proportion_cold','proportion_dry',
      #                  'CV',"exp_dir1","exp_dir2","concentration_dir1","concentration_dir2","alt_circular_model")
      # 
      # 
      # print(paste("IntermediateOutput/plant_PCA_nichefilling/rstates_center3var/", sp_name, '_',region,'_nicheshift.R',sep=''))
      # save(D.df, file=paste("IntermediateOutput/plant_PCA_nichefilling/stabilityRegional/rstates_center3var/", sp_name, '_',region,'_nicheshift.R',sep=''))
      # print("finished species")
      # return(D.df)
      # 
      # 
      
      ##########
      
      
    }
  }
}





D.df<-unlist(yooo[[1]])
D.df

for (i in yooo){
  #print(i[1])
  if(i[1] != "" & i[1] != "not enough known occurrences to calculate niche"){D.df<-rbind(D.df,i)}
}

D.df<-D.df[-c(1),]



csvfile4 <- paste("IntermediateOutput/plant_PCA_nichefilling/3Varplant_D_FILLvalues","_",cor_sel,"_center.txt",sep="")
write.table(D.df, csvfile4, row.names = FALSE, sep = "\t")

D.df_convert$species.name<-paste(D.df_convert$row.names,D.df_convert$species.name)
D.df_convert$row.names<-NULL

write.table(D.df_convert, csvfile4, row.names = FALSE, sep = "\t")
D.df_convert<-read.table(csvfile4,row.names=NULL,header = T)
head(D.df_convert)


