rm(list=ls())


#set folder with datapoints to analyse
library(raster)
library(rgeos)
library(rworldmap)
library(stringr)
library(RColorBrewer)
library(ecospat)

#functions from the original broenniman paper, various versions of input
#source("niche.overlap.functions.R")
#source("occ.prep.functions.R")
source("niche analysis/niche_dynamic_functions_Nov11th2016HH.R")
#source("niche_dynamic_functions_Nov11th2016HH.R")

##########################################################
#Functions for data processing

#takes lat lon point, projects it and then finds the continent/region the point is in
clookup<-function(x){
  point<-x[1:2]
  if (is.na(point[1]) != TRUE && is.na(point[2]) != TRUE){
    sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(shape)))
    rest<<-over(sp2, shape)
  }else{
    #print ("failed")
    fail<<-fail+1
  }
} 



###################################
#load files for processing

#load shapefile with zones for range classification
shape <- shapefile("biogeographic_zonesV2.shp")
#low-res map for plotting purposes
newmap <- getMap(resolution = "coarse")  # different resolutions available
#plot(shape,col=as.factor(shape@data$Name),axes=T)
#load world raster as base for data to be overlaid, and for info on grid cell area (km squared)
world<-raster("Global_10min_grid.grd")

#global climate space. This is basically the above .grd file but with region assigned to each cell (in csv form)
#it takes a while to calculate this everytime, so I just saved it with region assigned (check it matches shape)
clim123_ref <- na.exclude(read.csv("ProducData/bioclimNPP_10min_region.csv"))
head(clim123_ref)
clim123<-clim123_ref[,-c(6,7)]

colnames(clim123)<-c("x","y","X1","X2","X3","X4")


#list of all species for processing
plant_summary<-read.csv("species lists/plant_data_summary0607.csv")
#remove species we ruled out for whatever reason
plant_list<-plant_summary[plant_summary$Verdict == "Yes",]


sp_list<-as.character(plant_list[1:nrow(plant_list),1])
#sp_name<-as.character(plant_list[25:30,1])  

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

x<-sp_list[1]
sp_list2<-sp_list[1:20]
sp_list2

yooo<-lapply(sp_list2, niche_cal)
yooo
#the actual function to process the data, it's bad form to have it this way round, but saves on scrolling
niche_cal<-function(x){
  sp_name<-x
  print(sp_name)
  
  #these are original raw files but we use the deaggregated points for speed
  #file_natur=paste('species lists/plant_dist_data_0802/naturalised/',sp_name,'_naturalised.csv', sep='')
  #file_native=paste('species lists/plant_dist_data_0802/native/',sp_name,'_native.csv', sep='')
  #file_unclass=paste('species lists/plant_dist_data_0802/',sp_name,'_unclass.csv', sep='')
  
  file_natur=paste('species lists/plant_dist_data_0802/desaggregated_data/',sp_name,'_naturalised.csv', sep='')
  file_native=paste('species lists/plant_dist_data_0802/desaggregated_data/',sp_name,'_native.csv', sep='')
  #file_natur=paste('desaggregated_data/',sp_name,'_naturalised.csv', sep='')
  #file_native=paste('desaggregated_data/',sp_name,'_native.csv', sep='')
  
  
  
  native<-read.csv(file_native)
  natur<-read.csv(file_natur)
  if(ncol(native)>2){native<-native[,1:2]}
  if(ncol(natur)>2){natur<-natur[,1:2]}
  native<-native[!is.na(native$lon)&!is.na(native$lat),]
  natur<-natur[!is.na(natur$lon)&!is.na(natur$lat),]
  
  colnames(native)<-c("x","y")
  colnames(natur)<-c("x","y")
  
  
  #we comment out because we're just pre-desaggregated data
  #occ.sp_native<-ecospat.occ.desaggregation(df=native,colxy=1:2,min.dist=0.16666,plot=F) 
  #occ.sp_natur<-ecospat.occ.desaggregation(df=natur,colxy=1:2,min.dist=0.16666,plot=F) 
  
  occ.sp_native<-native
  occ.sp_natur<-natur
  
  #we look up all native points, find the regions they occur in and treat that as "available" climate
  clookup(occ.sp_native)
  native$region<-rest$Name
  native_regions<-unique(rest$Name)
  native_regions<-native_regions[!is.na(native_regions)]

  clim_sel<-c(1:5,8)
  clim1<-clim123_ref[clim123_ref$region %in%native_regions  , clim_sel]  
  
  #do the same for naturalised points, find regions naturalised points occur and in use this to build clim2
  clookup(occ.sp_natur)
  natur$region<-rest$Name
  natur_regions<-unique(rest$Name)
  natur_regions<-natur_regions[!is.na(natur_regions)]
  
  clim2<-clim123_ref[clim123_ref$region %in%natur_regions  , clim_sel]  
  
  native2<-native
  natur2<-natur
  
  #clim12 is a combination of clim1 (native) and clim2 (naturalised)
  clim12 <- rbind(clim1, clim2)
  
  #filter out occurences with fewer than 5 occurences, PCA cannot be built on fewer.
  if (nrow(native2)<=5 | nrow(natur2)<=5){
    print ("not enough known occurrences to calculate niche")
  }
  
  else{
    
    #extract bioclim variables for species occurences
    
    occ.sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp_native,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=0.16666))
    occ.sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp_natur,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=0.16666))
    
    
    
    if (nrow(occ.sp1)<5 | nrow(occ.sp2)<5){
      print ("Too many data point removals!")
      mergefail.df[nrow(mergefail.df)+1,1]<-c(sp_name)
      
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
      Xvar<-c(3:ncol(clim12)) # This means use columns 3 to x of the climate files. If you need to adjust the variables used, see where climate is defined above
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
      
      
      # predict the scores on the axes
      scores.clim12<- pca.cal$li[row.clim12,]
      scores.clim1<- pca.cal$li[row.clim1,]
      scores.clim2<- pca.cal$li[row.clim2,]
      scores.sp1<- pca.cal$li[row.sp1,]
      scores.sp2<- pca.cal$li[row.sp2,]
      
      
      #grid.clim.NNU is slightly different to the main plotting functions later (ecospat.grid.clim.dyn), but only in the output. Calculatations are the same
      #return niche matrices both as a matrix and as a raster:
      #at some point I will update this due to adehabitat being outdated
      z2<- grid.clim.NNU(scores.clim12, scores.clim2, scores.sp2, R) ## Uses the PCA scores of the entire climate space, native climate space, and species distribution to calculate occurrence density on a grid. The entire climate space is only used to set the boundaries of the grid.
      z1<- grid.clim.NNU(scores.clim12,scores.clim1,scores.sp1,R)

      
      
      #test of equivalency taken out due to time delay
      #a<-ecospat.niche.equivalency.test(z1,z2,rep=1)# test of niche equivalency and similarity according to Warren et al. 2008
      #b<-ecospat.niche.similarity.test(z1,z2,rep=100) ## can take some time
      #b2<-ecospat.niche.similarity.test(z2,z1,rep=100) ## can take some time
      
      #take D result (overlap) just in case we want it
      D_result<-round(as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]),3)
      #plot equivalency values, not used currently			
      #x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
      #ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
      #ecospat.plot.niche(z2,title="PCA-env - naturalised niche",name.axis1="PC1",name.axis2="PC2")
      #ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
      #plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]),3)))
      #ecospat.plot.overlap.test(a,"D","Equivalency")
      #ecospat.plot.overlap.test(b,"D","Similarity 2->1")
      #ecospat.plot.overlap.test(b2,"D","Similarity 1->2")
      
      
      #set threshold for 'occupied' climate
      thresh=0.00
      
      #get all of the native distribution, this is scaled from 0:1
      ntv.distn <- z1$z.uncor
      #ntv.distn.r <- raster(ntv.distn, xmn=min(z1$x), xmx=max(z1$x), ymn=min(z1$y), ymx=max(z1$y)) ## rows are the y values
      
      #remove occurrences that fall under the threshold of occupied climate
      #z.raw is z/Z, so scaled to availability of climate. remove any points from z.uncor that fall under this threshold
      ntv.distn.r.pa<-z1$z.raw>thresh
      ntv.distn.r<-ntv.distn*ntv.distn.r.pa
      ntv.distn.dens <- cellStats(ntv.distn.r, stat='sum')
      
      
      ### naturalised distribution
      
      natur.distn<-z2$z.uncor
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
      
      # # Climate in naturalised range
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
      
      
      
      
      #########PREP MAP DATAFRAMES
      native_clim_pca<-cbind(clim1,scores.clim1)
      #over_9000<-(native_clim_pca[native_clim_pca$Axis2<=-10,])
      #plot(native_clim_pca$Axis1,native_clim_pca$Axis2,pch=".",cex=2,col="black")
      #points(over_9000$Axis1,over_9000$Axis2,pch=".",cex=2,col="red")
      natur_clim_pca<-cbind(clim2,scores.clim2)
      
      
      
      ######NICHE DYNAMICS IN NATURALISED RANGE
      #create 3 raster layers, climate that represents expansion in climate (may or may not be occupied in space)
      #stability, climate that is occupied in both ranges (though we may or may not have a record)
      #restriction, climate that should be occupied (according to native range), but is not in the naturalised
      natur_occ<-cbind(occ.sp2,scores.sp2)
      native_occ<-cbind(occ.sp1,scores.sp1)

      col_ex<-c(ncol(natur_clim_pca)-1,ncol(natur_clim_pca))
      
      head(natur_clim_pca)
      #find real world coordinates that occur within analogue space
      analogue_coords<-extract(overlap100, natur_occ[,col_ex])
      naturalised_analogue100<-natur_occ[which(analogue_coords==1),]
      naturalised_nonanalogue100<-natur_occ[which(analogue_coords==0),]
      
      #find real world coordinates that occur within shared climate space
      analogue_coords2<-extract(overlap100, native_occ[,col_ex])
      native_analogue100<-native_occ[which(analogue_coords2==1),]
      native_nonanalogue100<-native_occ[which(analogue_coords2==0),]
      
      natur_in<-nrow(naturalised_analogue100)
      natur_out<-nrow(natur_occ) - nrow(naturalised_analogue100)
      
      native_in<-nrow(native_analogue100)
      native_out<-nrow(native_occ) - nrow(native_analogue100)
      
 
      

      
      
      ##########################################
      #compare against old method, this is Regan's original method and function, is the easiest way to extract NND summary statistics
      # calculation of occurence density and test of niche equivalency and similarity 
      
      #z1_old<- ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1,R,th.sp=0)
      #z2_old<- ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2,R,th.sp=0)
      #remove occurences that do not at least inhabit the threshold level of available climate type
      #z1_old$z.uncor[(z1_old$z/z1_old$Z)<thresh]<-0
      #image(z1_old$z.uncor)
      #image(z1_old$z/z1_old$Z)
      
      #z2_old$z.uncor[(z2_old$z/z2_old$Z)<thresh]<-0
      #image(z2_old$z.uncor)
      #image(z2_old$z/z2_old$Z)
      
      
      dyn.75<-dynamic.index(z1,z2,thresh=25) ## dynamics in 75% most common  analagous climate space
      dyn.100<-dynamic.index(z1,z2,thresh=0) ## dynamics in all analagous climate space
      dyn.all<-dynamic.index(z1,z2,thresh='xxx') ## dynamics in all climate space
      
      nnddyn.75<-dynamic.index.nnd(z1,z2,thresh=25) ## dynamics in 75% most common  analagous climate space
      nnddyn.100<-dynamic.index.nnd(z1,z2,thresh=0) ## dynamics in all analagous climate space
      
      nnddyn.all<-dynamic.index.nnd(z1,z2,thresh='xxx') ## dynamics in all climate space
      
      ## Make pretty plot of the niche shift
      #pdf(file=paste(wd,sp,'_KS PCA-env 3vars.pdf',sep=''), width=8, height=8, pointsize=12) ## This will save the plot to a pdf file in the working directory specified above. Uncomment if you wish to use
      x11(width=180,height=200)

      a<-cellStats(z2$z.uncor,'max')
      b<-cellStats(z1$z.uncor,'max')
      if(a>0 & b>0){
        par(mfrow=c(3,2))
        
        #t1 <- paste("PCA-env - Native niche", sep='')
        #plot.niche(z1,title=t1,name.axis1="PC1",name.axis2="PC2")
        #t2 <- paste("PCA-env - Naturalised niche", sep='')
        #plot.niche(z2,title=t2,name.axis1="PC1",name.axis2="PC2")
        
        plot.contrib(pca.cal$co,pca.cal$eig)
        plot.niche.dyn(z1,z2,0.25,int=2, title=paste("niche overlap between native and naturalised"))
        #   fun.arrows(scores.sp1,scores.sp2,scores.clim1,scores.clim2) ## Don't use because just weird
        
        t.75 <- paste('75% CS: exp=',dyn.75[[2]][1],' stab=',dyn.75[[2]][2],' unf=',dyn.75[[2]][3], sep='')
        t.100 <- paste('100% CS: exp=',dyn.100[[2]][1],' stab=',dyn.100[[2]][2],' unf=',dyn.100[[2]][3], sep='')
        t.all <- paste('All: exp=',dyn.all[[2]][1],' stab=',dyn.all[[2]][2],' unf=',dyn.all[[2]][3], sep='')
        text(-0.5,1.9,t.75, cex=0.8)
        text(-0.5,0.9,t.100, cex=0.8)
        text(-0.5,-.1,t.all, cex=0.8)
        # 
        natur_plot<-natur.distn.r.pa
        natur_plot[natur_plot==0]<-NA
        plot(overlap100,legend=F,main=paste("Natur: In:",natur_in,"Out:",natur_out))
        plot(natur_plot,add=T,col=rainbow(1,alpha=0.7),legend=F)
        points(naturalised_nonanalogue100$Axis1,naturalised_nonanalogue100$Axis2,col="Blue",pch=".",cex=3)
        points(naturalised_analogue100$Axis1,naturalised_analogue100$Axis2,col="Yellow",pch=".",cex=3)
        
        
        native_plot<-ntv.distn.r.pa
        native_plot[native_plot==0]<-NA
        plot(overlap100,legend=F,main=paste("Native: In:",native_in,"Out:",native_out))
        plot(native_plot,add=T,col=rainbow(1,alpha=0.7),legend=F)
        points(native_nonanalogue100$Axis1,native_nonanalogue100$Axis2,col="Blue",pch=".",cex=3)
        points(native_analogue100$Axis1,native_analogue100$Axis2,col="Yellow",pch=".",cex=3)
        
        
        
        plot(shape,main="Native distribution")
        points(native_analogue100$x,native_analogue100$y,col="red",pch=".",cex=4)
        points(native_nonanalogue100$x,native_nonanalogue100$y,col="Blue",pch=".",cex=4)

        

        plot(shape,main="Naturalised distribution")
        points(naturalised_nonanalogue100$x,naturalised_nonanalogue100$y,col="red",pch=".",cex=4)
        points(naturalised_analogue100$x,naturalised_analogue100$y,col="Yellow",pch=".",cex=4)

        

        plotname=paste('niche_shift_work/analogue_disagg_21112018/Figures/',sp_name,'_PCA.png', sep='')
        dev.copy(png, filename=plotname, width=1000, height=1200,res=150)

        dev.off() # Uncomment if you wish to use the pdf function above
      }
      
      D.df<-c(sp_name,nrow(naturalised_analogue100),nrow(naturalised_nonanalogue100),
             nrow(native_analogue100),nrow(native_nonanalogue100))
      names(D.df) <- c('species_name', 'kept_natur','thrown_natur','kept_native','thrown_native')
      
      
      write.csv(native_analogue100,paste('niche_shift_work/analogue_disagg_21112018/',sp_name,'_native.csv', sep=''), row.names = FALSE)
      write.csv(naturalised_analogue100,paste('niche_shift_work/analogue_disagg_21112018/',sp_name,'_natur.csv', sep=''), row.names = FALSE)
      
      
      return(D.df)
      

      
      
      ##########
      
      
    }
  }
}

yooo=="not enough known occurrences to calculate niche"
D.df<-yooo[yooo!="not enough known occurrences to calculate niche"]
D.df

l_check<-lapply(D.df,length)
print(D.df[l_check!=5])
D.df<-D.df[l_check==5]

D.df2<-do.call(rbind, D.df)

write.csv(D.df2,"niche_shift_work/analogue_disagg_21112018/Figures/summary_loss.csv",row.names=F)


