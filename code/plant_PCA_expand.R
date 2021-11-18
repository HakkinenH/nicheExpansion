
##############################################
### DETECT SPECIES NICHE EXPANSION ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 
# This takes 1) a species list 2) occurrence data (desaggreagated AND non-analogous points removed), stored in RawData/GBIF and 3) compiled WorldClim data in a csv file
# In this script we project species' native and naturalised occurrences in PCA space to estimate their climatic niche.
# We find areas of niche expansion, and calculate several things about it (what areas, what proportion of niche expansion, what directions are expansion cells in?)
# We then save the output with as much detail as possible for later analysis
# Output is stored in: IntermediateOutput/plant_PCA_expand

# Based on Broenniman et al (2012) PCA method, adapted by Regan Early to measure niche NND and other dynamics
# Modified by HHakkinen 2017 to process species data. Points are read in, categorised, PCA space is created
# then save the outputs as raster for later summary, as well graphical output


### ###


rm(list=ls())

#set to location of repository
setwd("DIRECTORY_HERE")


#functions from the original broenniman paper, various versions of input
source("code/functions/niche_dynamic_functions_Nov11th2016HH.R")
source("code/functions/PCA_functions.R")
source("code/functions/miscFunctions.R")


library(dplyr)
library(biomod2)
library(ade4)
library(raster)
library(rworldmap)
library(ecospat)
library(adehabitatHR)
library(circular)
library(CircMLE)
library(rgeos)
library(sp)





###################################
#load files for processing
##################################
#CHECK FILE NAMES AND PATHS ARE CORRECT BEFORE RUNNING! 
#IF this section loads correctly then rest should run with no problems


#load shapefile with zones for range classification
shape <- shapefile("RawData/BiogeographicZones/biogeographic_zonesV2.shp")
#low-res map for plotting purposes
newmap <- getMap(resolution = "coarse")  # different resolutions available
#plot(shape,col=as.factor(shape@data$Name),axes=T)


#this is 3 climate variables
clim123_ref <- na.exclude(read.csv("RawData/WorldClim/bioclim_3Var_10min_region.csv"))
#this is a 4 variable version that includes precipitation seasonality
#clim123_ref <- na.exclude(read.csv("RawData/WorldClim/bioclim_4Var_10min_region.csv"))



#SELECT VARIABLES TO RUN THE ANALYSIS
#3 VAR IS TMAX, TMIN AND PRECIP: BIO5 BIO6 AND BIO12
#2 VAR is TMAX AND PRECIP: BIO5 AND BIO12
#4 VAR is TMAX. TMIN, PRECIP, PRECIP SEASON, BIO5, BIO6, BIO12, BIO15

varlis<-c("bio5", "bio6", "bio12")


#CHOOSE OUTPUT Folder (options are 2Var, 3Var, 4Var)
#set the name of the variable set. NOTE: this sets folder names and output files for later, so make sure it's memorable!
#we use this name to keep track of output and results
folOut<-"3Var"



clim123<-clim123_ref[,c("x","y", varlis)]
colnames(clim123)<-c("x","y","X1","X2","X3")


#if you want to look at what the PCA space will look like you can load these pre-made versions
#previously we've saved a copy of the global PCA (slightly circular but there it is)
#open and extract direction of each component vector
#if(folOut=="3Var"){load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecip")}
#if(folOut=="2Var"){load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TmaxPrecip")}
#if(folOut=="4Var"){load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecipPrSeas")}





#list of all species for processing. 
#this is the new summary list for informing analogue climate
plant_summary<-read.csv("IntermediateOutput/PCA_find_analogue_byRegion/plant_specieslist_analoguefiltered.csv",stringsAsFactors = F)


#blank dataframe that will be filled row by row for any species that fail. Useful for checking if anything fails
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







#################################
#DEFINE FUNCTIONS TO FIND EXPANSION IN PCA SPACE
#################################

#the actual function to process the data, it's bad form to have it this way round, but saves on scrolling
niche_cal<-function(sp_name,region){

  print(sp_name)
  print(region)
  

  #this is out current dataset
  #deaggregaated, and filtered so only analogue native/naturalised occurrence data is included
  #native (all regions) vs naturalised (the region we are currently looking at only)
  file_natur = paste("IntermediateOutput/PCA_find_analogue_byregion/",region,"/",sp_name,region,"_natur.csv",sep="")
  file_native = paste("IntermediateOutput/PCA_find_analogue_byregion/",region,"/",sp_name,region,"_native.csv",sep="")
  
  #read in all files and check formats and names are correct
  native<-read.csv(file_native)
  natur<-read.csv(file_natur)
  if(ncol(native)>2){native<-native[,1:2]}
  if(ncol(natur)>2){natur<-natur[,1:2]}
  
  native<-native[!is.na(native$x)&!is.na(native$y),]
  natur<-natur[!is.na(natur$x)&!is.na(natur$y),]
  
  colnames(native)<-c("x","y")
  colnames(natur)<-c("x","y")
  
  occ.sp_native<-native
  occ.sp_natur<-natur

  #global climate is treated as equally available, set climate
  clim1<-clim123_ref[,c("x","y", varlis)]
  clim2<-clim123_ref[,c("x","y", varlis)] 

  
  #clim12 is a combination of clim1 (native) and clim2 (naturalised)
  clim12 <- rbind(clim1, clim2)
  
  #filter out species with fewer than 5 occurences, PCA cannot be built on fewer.
  if (nrow(occ.sp_native)<=5 | nrow(occ.sp_natur)<=5){
    print ("not enough known occurrences to calculate niche")
  }
  
  else{
    
    #extract bioclim variables for species occurrences

    occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp_native,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=0.16666))
    occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp_natur,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=0.16666))
    
    
    
    if (nrow(occ.sp1)<5 | nrow(occ.sp2)<5){
      print ("Too many data point removals!")
      mergefail.df[nrow(mergefail.df)+1,1]<-c(sp_name)
      print ("not enough known occurrences to calculate niche")
      
    }
    else{
      
      
      ###
      ### ANALYSIS - selection of parameters 
      ###
      
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
      
      
      ###
      ###row weighting and grouping factors for ade4 functions 
      ###
      
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
      
      
      
      
      ###
      ### PCA-ENV 
      ###
      
      # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
    
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)

      
      #save global PCA for later use if you like
      #save(pca.cal,file="./IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecipPrSeas")
      #plot the loadings if needed
      ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
      arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))
      colnames(arr_tab)<-c("x","y")
      #table with points which give the direction as compared to the origin (0,0)
      #find the angle of direction for each PCA component
      arr_tab$angle<-apply(arr_tab,1,rad.ang.center,x2=c(0,0))
      
      #plot(0,0,xlim=c(-1,1),ylim=c(-1,1))
      #points(arr_tab$x,arr_tab$y)
      #arrows.circular(arr_tab[1:4,3],col="red")

      # predict the scores on the axes
      scores.clim12<- pca.cal$li[row.clim12,]
      scores.clim1<- pca.cal$li[row.clim1,]
      scores.clim2<- pca.cal$li[row.clim2,]
      scores.sp1<- pca.cal$li[row.sp1,]
      scores.sp2<- pca.cal$li[row.sp2,]
      
      
      #grid.clim.NNU is slightly different to the main plotting functions later (ecospat.grid.clim.dyn), but only in the output. Calculations are the same
      #return niche matrices both as a matrix and as a raster:
      z2<- grid.clim.NNU2(scores.clim12, scores.clim2, scores.sp2, R) ## Uses the PCA scores of the entire climate space, native climate space, and species distribution to calculate occurrence density on a grid. The entire climate space is only used to set the boundaries of the grid.
      z1<- grid.clim.NNU2(scores.clim12,scores.clim1,scores.sp1,R)

      
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
      
      ## Overlap between native and naturalised - all climate space
      z1.r100 <- z1.r > 0
      z2.r100 <- z2.r > 0
      overlap100 <- z1.r100 * z2.r100
      
      #in order to calculate what points lie where, what do I need?
      #outline of native occupied niche (100% analogue)
      ntv.distn.r.pa100<-ntv.distn.r.pa*overlap100
      
      #outline of naturalised occupied niche (100% analogue)
      natur.distn.r.pa100<-natur.distn.r.pa*overlap100
      
      #stab, exp and res in analogue space
      stab100<-stab*overlap100
      exp100<-exp*overlap100
      res100<-res*overlap100
      
      
      ###we have found all areas of overlapping niche space so now we can start to answer questions
      
      
      #########################1) do species expand? >0.1 expansion###################
      
      exp_stat<-dyn.100$dynamic.index.w["expansion"]
      
      exp_points<-scores.sp2[extract(exp100, scores.sp2)==1,]


      #is there any expansion?
      
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
      
      #an alternative method if more than 1 medoid
      #native_center<-as.data.frame(pam(scores.sp1, 1)$medoids)
      #natur_center<-as.data.frame(pam(scores.sp2, 1)$medoids)
      
      #what is the climate like in native center (center of expansion)
      native_clim_pca<-cbind(clim1,scores.clim1)
      natur_clim_pca<-cbind(clim2,scores.clim2)
      
      #link points in PCA space and their original lat/lon
      natur_occ<-cbind(occ.sp2,scores.sp2)
      natur_occ$exp<-0
      
      #specify columns needed
      col_ex<-c("Axis1","Axis2")

      #find real world coordinates that occur within expansion space
      expansion_coords<-extract(exp100, natur_occ[,col_ex])
      natur_occ[which(expansion_coords==1),"exp"]<-1


      #save these points for later
      csvfile4 <- paste("IntermediateOutput/plant_PCA_exp/",folOut,"/NaturPCAvalues/",sp_name,region,".csv",sep="")
      write.table(natur_occ, csvfile4, row.names = FALSE, sep = "\t")
      
      
      
      #find all naturalised climate cells that lie in native region and are therefore analogue
      natur_PCA<-raster(paste("IntermediateOutput/findRegionOutline/",folOut,"/",region,"PCA.tif",sep=""))
      
      
      #we also need info on the potential naturalised niche. Where is the center?
      #take naturalised climate and crop by native niche
      natur_poten <-natur_PCA*ntv.distn.r.pa100

      
      natur_poten_df<-as.data.frame(coordinates(natur_poten))
      natur_poten_df$layer<-as.data.frame(natur_poten$layer)$layer
      if(nrow(natur_poten_df[which(natur_poten_df$layer==0),])>0){natur_poten_df[which(natur_poten_df$layer==0),]<-NA}
      natur_poten_df<-natur_poten_df[!(is.na(natur_poten_df$layer)),]
      
      natur_poten_center<-COGravity(x=natur_poten_df$x, y=natur_poten_df$y, wt=natur_poten_df$layer)[c(1,3)]
      natur_poten_center<-data.frame(x=natur_poten_center[1],y=natur_poten_center[2])
      
      #now find information, what is it like at the center of the native niche? Get climate variable information
      native_cinfo<-near_point2(native_center,native_clim_pca)
      native_cinfo
      
      #prep output data
      circ_model<-NA;circ_rstat<-NA;circ_pvalue<-NA
      exp.mean.dir1<-NA;exp.median.dir1<-NA;exp.ten.dir1<-NA;exp.ninety.dir1<-NA;exp.max.dir1<-NA;
      exp.mean.dir2<-NA;exp.median.dir2<-NA;exp.ten.dir2<-NA;exp.ninety.dir2<-NA;exp.max.dir2<-NA;
      tmax_pro<-NA;tmin_pro<-NA;precip_pro<-NA;precipSea_pro<-NA;
      CV<-NA;dir1<-NA;dir2<-NA;circ_altmodel<-NA; conc_dir1<-NA;conc_dir2<-NA
      
      
      
      #only do the following if there is expansion
      if(nrow(exp_points)>0){

        ############## 2) is there strong expansion in one (or two) directions? Or is it in every direction? ##############
        #circ_mle comparing the central point of naturalised density with every gridcell of expansion (can I weight by density of gridcell?)
        #what model describes it best? Null is every direction at once by not very much
        #can we categorise it roughly? 1) no expansion 2) unidirectional expansion 3) bidirectional expansion 4)multidirectional exp
        

        #find all expansion cells and their associated density
        exp.distn.r<-exp100*z2$zz
        
        natur_PCA<-raster(paste("IntermediateOutput/findRegionOutline/",folOut,"/",region,"PCA.tif",sep=""))
        natur_PCA[natur_PCA==0]<-NA
        
        #rescale expansion for easier analysis (0-1000)
        exp.distn.r[exp.distn.r==0]<-NA
        #set threshold of expansion to discount
        thresh_exp<-quantile(exp.distn.r,probs=thresh_exp_per)

        exp.distn.r[exp.distn.r<thresh_exp]<-NA
        

        #correct for climate availability in the naturalised region
        exp.distn.r<-exp.distn.r/natur_PCA
        #rescale 0-1 to make it comparable for other species
        exp.distn.r<-exp.distn.r/cellStats(exp.distn.r,stat="max")

        #if we're not interested in compensating for climate availability then use these lines
        ### naturalised distribution
        #natur.distn<-z2[cor_sel][[1]]
        #natur.distn<-z2$z.cor
        #remove occurrences that fall under the threshold of occupied climate
        #natur.distn.r.pa<-z2$z.raw>thresh
        #natur.distn.r<-natur.distn*natur.distn.r.pa
        #natur.distn.dens <- cellStats(natur.distn.r, stat='sum')
        #exp.distn.r<-exp100*z2$natur.distn.r
        

        exp.distn.r2<-exp.distn.r*1000
        exp.distn.r2<-round(exp.distn.r2)
        exp.distn.r2[exp.distn.r2<1]<-NA
        
        
        #find distance from the central native density to each expansion point
        exp_df<-as.data.frame(coordinates(exp.distn.r2))
        exp_df$layer_orig<-as.data.frame(exp.distn.r$layer)$layer
        exp_df$layer<-as.data.frame(exp.distn.r2$layer)$layer


        exp_df<-exp_df[!is.na(exp_df$layer),]
        

        #calculate distance to each expansion cell
        exp_df$exp_dist<-unlist(apply(exp_df[,1:2], 1, euc.dist.center, x2=natur_poten_center))

        #calculate angle to each cell, 
        exp_df$exp_dir<-unlist(apply(exp_df[,1:2], 1, rad.ang.center, x2=natur_poten_center))

        
        #calculate circular variance http://resources.esri.com/help/9.3/arcgisdesktop/com/gp_toolref/spatial_statistics_tools/how_linear_directional_mean_spatial_statistics_works.htm
        CV<-1-(sqrt((sum(sin(exp_df$exp_dir))^2)+(sum(cos(exp_df$exp_dir))^2))/nrow(exp_df))
        
        
        #this is a nice summary table but in order to compensate for varying density of each expansion point
        #we are doing a weighting process. Each density point has been rescaled to have a minimum of one
        #we then multiply the angle by the weighted number
        
        exp_dir<-rep(exp_df$exp_dir,exp_df$layer)
        #exp_dir<-exp_df$exp_dir
        
        circ_res<-circ_mle(exp_dir)
        
        #what is the best model to describe expansion?
        
        circ_model<-circ_res$bestmodel
        circ_rstat<-circ_res$rt[1]
        circ_pvalue<-circ_res$rt[2]
        
        dAIC<-circ_res$results$deltaAIC
        
        #select names of other viable models
        dAIC2<-rownames(circ_res$results[which(dAIC>0&dAIC<10),])
        
        if(length(dAIC2)>0){circ_altmodel<-paste(dAIC2,collapse=";")}
        
        dir1<-circ_res$results[1,2]
        dir2<-circ_res$results[1,5]
        conc_dir1<-circ_res$results[1,3]
        conc_dir2<-circ_res$results[1,6]
        
        ntv_plot<-ntv.distn.r
        ntv_plot[ntv_plot==0]<-NA

        #save output for inspection
        plotname=paste('IntermediateOutput/plant_PCA_exp/',folOut,'/bucketDirection/',sp_name,'_',region,'_shift.jpg', sep='')
        png(file=plotname, width=900, height=700)
        
        par(mfrow=c(1,2))
        plot(exp100,legend=F,main=sp_name,col=c("grey","blue"))
        #plot(exp.distn.r2,add=T)
        plot(exp.distn.r,add=T)
        
        points(natur_poten_center,col="red")
        arrows(x1=exp_df$x,y1=exp_df$y,
               x0=natur_poten_center$x,y0=natur_poten_center$y,
               col="red",length=0.05)
        
        plot_circMLE(exp_dir, circ_res)

        mtext(circ_model)
        dev.off()
        
        #plot breakdown to look at
        plotname=paste('IntermediateOutput/plant_PCA_exp/',folOut,'/circleBreakdown/',sp_name,'_',region,'.jpg', sep='')
        png(file=plotname, width=1000, height=1000, res=75)
        
        par(mfrow=c(3,4))
        
        model.names<-row.names(circ_res$results)
        for(i in model.names){
          plot_circMLE(exp_dir, circ_res, model = i)
          mtext(i)
        }
        
        dev.off()
        
        
        

        ####################################
        #DIRECTION CUT DOWN AND SUMMARISE
        ###################################
        
        ###we have dir1 (and maybe dir2) which are directions of expansion
        #we have exp_df, find all points that are within 45
        ###what is the mean expansion distance in this direction?
        
        #made function for this because we're going to do this several times
        #find all points that lie in a 90 degree slice around our primary direction of expansion
        exp_df_dir1<-find_slice(exp_df,dir1)
        
        
        ####calculate mean,median,max and 5% 95% of expansion for points within this direction
        ####is there any reason to want to know overall mean? 
        #this is the length of the line
        exp.mean.dir1<-mean(exp_df_dir1$exp_dist)
        exp.median.dir1<-median(exp_df_dir1$exp_dist)
        exp.ten.dir1<-quantile(exp_df_dir1$exp_dist,0.1)
        exp.ninety.dir1<-quantile(exp_df_dir1$exp_dist,0.9)
        exp.max.dir1<-max(exp_df_dir1$exp_dist)
        if(exp.max.dir1==-Inf){
          #infinite means something has gone wrong (usually becuase we try and divide by 0)
          #make a note and investigate if needed
          print("oh noooooo")
          resu<-c(sp_name,"dir1",dir1,exp.max.dir1,nrow(exp_df_dir1))
          write.csv(resu,paste("./IntermediateOutput/plant_PCA_exp/",folOut,"/failure/",sp_name,"-",region,".csv",sep=""))

          }
        
        #if direction 2 exists in model, extract the same distance info
        if(!is.na(dir2)){
          
          #find all points that lie in a 90 degree slice around our secondary direction of expansion (if it exists)
          exp_df_dir2<-find_slice(exp_df,dir2)
        
          summary(exp_df_dir2)
          
          ####calculate mean,median,max and 5% 95% of expansion for points within this direction
          ####is there any reason to want to know overall mean? 
          #this is the length of the line
          exp.mean.dir2<-mean(exp_df_dir2$exp_dist)
          exp.median.dir2<-median(exp_df_dir2$exp_dist)
          exp.ten.dir2<-quantile(exp_df_dir2$exp_dist,0.1)
          exp.ninety.dir2<-quantile(exp_df_dir2$exp_dist,0.9)
          exp.max.dir2<-max(exp_df_dir2$exp_dist)
          
          if(exp.max.dir2==-Inf){
            
            print("oh noooooo")
            resu<-c(sp_name,"dir2",dir2,exp.max.dir2,nrow(exp_df_dir2))
            write.csv(resu,paste("./IntermediateOutput/plant_PCA_exp/",folOut,"/failure/",sp_name,region,".csv",sep=""))
            
          }
          
        }
        
        
        
        #we also want to know how much niche filling and expansion is in each climatic direction
        #for that we take a slice and take some summary stats
        #find all points that lie in a 90 degree slice around Precip vector (increasing)
        if ("bio12" %in% varlis){
          precip_dir<-arr_tab["bio12","angle"]
          precip_df<-find_slice(exp_df,precip_dir)
          #what proportion of exp density is in this slice?
          precip_pro<-sum(precip_df$layer_orig)/sum(exp_df$layer_orig)
        }
        
        #find all points that lie in a 90 degree slice around tmin vector (increasing)
        if ("bio6" %in% varlis){
          tmin_dir<-arr_tab["bio6","angle"]
          tmin_df<-find_slice(exp_df,tmin_dir)
          #what proportion of exp density is in this slice?
          tmin_pro<-sum(tmin_df$layer_orig)/sum(exp_df$layer_orig)
        }

        if ("bio5" %in% varlis){
          #find all points that lie in a 90 degree slice around tmax vector (increasing)
          tmax_dir<-arr_tab["bio5","angle"]
          tmax_df<-find_slice(exp_df,tmax_dir)
          #what proportion of exp density is in this slice?
          tmax_pro<-sum(tmax_df$layer_orig)/sum(exp_df$layer_orig)
        }
        
        if ("bio15" %in% varlis){
          #find all points that lie in a 90 degree slice around tmax vector (increasing)
          tmax_dir<-arr_tab["bio15","angle"]
          tmax_df<-find_slice(exp_df,tmax_dir)
          #what proportion of exp density is in this slice?
          tmax_pro<-sum(tmax_df$layer_orig)/sum(exp_df$layer_orig)
        }
        
        
        
        
        ###we extract the ordination values and densities to build into our summary table
        
        dir1_plot<-as.circular(dir1,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
        
        precip_dir_plot<-as.circular(precip_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
        tmin_dir_plot<-as.circular(tmin_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')
        tmax_dir_plot<-as.circular(tmax_dir,type='angles',units='radians',template='none',modulo='asis',zero=0,rotation='counter')

        
        plotname=paste('IntermediateOutput/plant_PCA_exp/',folOut,'/bucketDirection/',sp_name,'_',region,'_shiftpie.jpg', sep='')
        png(file=plotname, width=900, height=700)
        
        #let's plot
        par(mfrow=c(2,2))
        plot(exp.distn.r)
        #precip
        plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Precip")
        points(exp_df$x,exp_df$y,col="blue")
        
        points(natur_poten_center$x,natur_poten_center$y,col="red")
        
        arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
        
        arrows.circular(precip_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        arrows.circular(precip_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        arrows.circular(precip_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        
        points(precip_df$x,precip_df$y,col="red")
        text(2,2,paste("proportion:", round(precip_pro, digits=2)))
        
        #tmin
        plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Tmin")
        points(exp_df$x,exp_df$y,col="blue")

        points(natur_poten_center$x,natur_poten_center$y,col="red")
        arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)

        arrows.circular(tmin_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        arrows.circular(tmin_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        arrows.circular(tmin_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)

        points(tmin_df$x,tmin_df$y,col="red")
        text(2,2,paste("proportion:", round(tmin_pro, digits=2)))

        #tmax
        plot(ntv.overlap,xlim=c(-4,3),ylim=c(-4,4),main="Tmax")
        points(exp_df$x,exp_df$y,col="blue")
        
        points(natur_poten_center$x,natur_poten_center$y,col="red")
        arrows.circular(dir1_plot, shrink=1, length=0.1, col="blue", x0=natur_poten_center$x, y0=natur_poten_center$y)
        
        arrows.circular(tmax_dir_plot, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        arrows.circular(tmax_dir_plot+0.7853982, lty=2, shrink=3, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        arrows.circular(tmax_dir_plot-0.7853982, shrink=3, lty=2, length=0.1, x0=natur_poten_center$x, y0=natur_poten_center$y)
        
        points(tmax_df$x,tmax_df$y,col="red")
        text(2,2,paste("proportion:", round(tmax_pro, digits=2)))
        
        
        dev.off()
        

        exp.distn.r2_save<-exp.distn.r2>0

        #now save the final output
        writeRaster(exp.distn.r, paste("IntermediateOutput/plant_PCA_exp/",folOut,"/shift_PCArasters/",sp_name,'_',region,"_PCAraster.tif",sep=""), options= c("COMPRESSION=LZW","INTERLEAVE=PIXEL"), overwrite=TRUE)
        
        write.csv(exp_df,paste("IntermediateOutput/plant_PCA_exp/",folOut,"/DirData/",sp_name,'_',region,"_dir.csv",sep=""))
        
        
      }
      
      ######show plots to summarise these points
      plotname=paste('IntermediateOutput/plant_PCA_exp/',folOut,'/shift_PCArasters/',sp_name,'_',region,'_similarity.pdf', sep='')
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
      
      #there are a lot of measurements here but for report here are the major conclusions:
      #level of expansion (summed density in expansion space: expansion
      #occupied gridcells with expansion: Occupied_expcells
      #magnitude of niche expansion (average length): exp.dis
      #direction of niche expansion: exp_rad
      #variance of niche expansion: CV
      
      scores.sp1$range<-"native"
      scores.sp2$range<-"natur"


      D.df <- c(sp_name,region,nrow(occ.sp1),nrow(occ.sp2),D_result,dyn.100$dynamic.index.w, nnddyn.100$nnd,
                native_center,natur_center,natur_poten_center,
                native_cinfo[c(varlis)],
                circ_model,circ_rstat,circ_pvalue,
                exp.mean.dir1,exp.median.dir1,exp.ten.dir1,exp.ninety.dir1,exp.max.dir1,
                exp.mean.dir2,exp.median.dir2,exp.ten.dir2,exp.ninety.dir2,exp.max.dir2,
                tmax_pro,
                tmin_pro,
                precip_pro, 
                precipSea_pro,
                CV,dir1,dir2,conc_dir1,conc_dir2,circ_altmodel)
      
      if(folOut=="2Var"){
        ntv_inf<-c('native_center_tmax','native_center_precip')
        propinfo<-c("proportion_tmax", "proportion_precip") }
      if(folOut=="3Var"){
        ntv_inf<-c('native_center_tmax','native_center_tmin','native_center_precip')
        propinfo<-c("proportion_tmax", "proportion_tmin", "proportion_precip") }
      if(folOut=="4Var"){
        ntv_inf<-c('native_center_tmax','native_center_tmin','native_center_precip', 'native_center_precipSeas')
        propinfo<-c("proportion_tmax", "proportion_tmin", "proportion_precip", "proportion_precipSeas") } 
      
      
      names(D.df) <- c('species name', 'region', 'native_occurences','naturalised_occurences',
                       'Dstat_overlap', names(dyn.100$dynamic.index.w), names(nnddyn.100$nnd),
                       'native_center_x','native_center_y','natur_center_x','natur_center_y',"natur_poten_x","natur_poten_y",
                       ntv_inf,
                       'circular_model','circ_rstat','circ_pstat',
                       'mean_exp_dist1','median_exp_dist1','ten_exp_dist1','ninety_exp_dist1','max_exp_dist1',
                       'mean_exp_dist2','median_exp_dist2','ten_exp_dist2','ninety_exp_dist2','max_exp_dist2',
                       propinfo,
                       'CV',"exp_dir1","exp_dir2","concentration_dir1","concentration_dir2","alt_circular_model")
      
      
      print(paste("IntermediateOutput/plant_PCA_exp/",folOut,"/rstates_center/", sp_name, '_',region,'_nicheshift.R',sep=''))
      save(D.df, file=paste("IntermediateOutput/plant_PCA_exp/",folOut,"/rstates_center/", sp_name, '_',region,'_nicheshift.R',sep=''))
      print("finished species")
      return(D.df)
      
      
      
      ##########
      
      
    }
  }
}


#application of the function
yooo<-mapply(niche_cal, plant_summary$species_name,plant_summary$region)


#can run a single species for testing purposes
#sp_name<-"Cerastium fontanum"
#region<-"Neotropical"
#niche_cal(sp_name, region)



