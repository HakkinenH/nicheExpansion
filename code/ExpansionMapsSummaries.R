

##############################################
### EXPANSION MAPS AND SUMMARIES ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 
# This code produces a variety of different maps to summarise our data
# 1) It produces a histogram of how many species expand overall and by how much (it's not a map but it's short and goes with the maps in the paper)
# 2) It prints a map of the biogeographic zones
# 3) It summarises how many native/naturalised species occur in each region and subregion. (Useful way to visualise data)
# 4) It summarises how many species are undergoing niche expansion in each region and subregion

# The output from this file is used in Figure 2 and Appendix Figure S1.1

### ###

rm(list=ls())

#set to current repo
setwd("DIRECTORY_HERE")


library(raster)
library(rworldmap)
library(rgeos)
library(RColorBrewer)
library(ggplot2)
library(wesanderson)

source("./code/functions/miscFunctions.R")

###################################
#load files for processing
##################################

#load output from plant_pca_expansion and Rstate_compile.R
#niche expansion direction based on 3 variables
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")
niche_shift<-read.delim(csvfile4, sep="\t")


table_name<-"plant"

#load a species list
plant_summary<-read.csv("IntermediateOutput/PCA_find_analogue_byRegion/plant_specieslist_analoguefiltered.csv",stringsAsFactors = F)
sp_list<-unique(plant_summary$species_name)



#get a subregion map based on L4 administrative data
shape<-shapefile("./RawData/countryData/level4.shp")
#cut some larger landmasses out
shape<-shape[which(shape$LEVEL_4_NA!="Greenland"&shape$LEVEL_4_NA!="Antarctica"),]


#load biogeographic region info
region_shape <- shapefile("./RawData/BiogeographicZones/biogeographic_zonesV2.shp")



#################################################
#1) histogram of expansion
#################################################

pdf(file="./FinalOutput/ExpansionMapsSummaries/ExpSumm_hist.pdf",width=5,height=5)
hist(niche_shift$expansion,w=0.1,xlab="Proportion of Niche Expansion",main="",col="white")
abline(v=0.1,col="red",lwd=3)
dev.off()



#################################################
#2) plot a map of biogeographic zones
#################################################


#create a new aesthetic df. Colour by a.diff.sum
shp_df <- broom::tidy(region_shape, region = "Name")
table(shp_df$id)
#as.numeric(shp_df$id)


plotname=paste("./FinalOutput/Supplementary/ExpansionMapsSummaries/BiogeoZoneMap.pdf",sep="")
pdf(file=plotname, width=8, height=4)
map <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = (id)), colour = "black") +
  theme_void()+labs(fill = "Realm")
print(map)
dev.off()


#################################################
#3) 4) compile and plot native/naturalised/expansion maps
#################################################
region_shape$native_count<-0
region_shape$natur_count<-0
region_shape$exp_count<-0

region_shape$native_subcount<-0
region_shape$natur_subcount<-0
region_shape$exp_subcount<-0

shape$native_count<-0
shape$natur_count<-0
shape$exp_count<-0

shape$native_subcount<-0
shape$natur_subcount<-0
shape$exp_subcount<-0


nunavut_list<-c()


for (i in 1: length(sp_list)){
#for (i in 1: 20){
  sp_name<-sp_list[i]
  print(sp_name)
  
  sp_rows<-plant_summary[plant_summary$species_name==sp_name,]
  
  
  #reset counter to avoid double counts
  region_shape$native_subcount<-0
  region_shape$natur_subcount<-0
  region_shape$exp_subcount<-0
  
  shape$native_subcount<-0
  shape$natur_subcount<-0
  shape$exp_subcount<-0
  
  if(nrow(sp_rows)==0){break()}

  
  file_natur=paste(data_path,sp_name,'_naturalised.csv', sep='')
  file_native=paste(data_path,sp_name,'_native.csv', sep='')
  
  
  native<-read.csv(file_native)
  natur<-read.csv(file_natur)
  
  colnames(native)<-c("x","y")
  colnames(natur)<-c("x","y")
  
  native<-native[!is.na(native$x)&!is.na(native$y),]
  natur<-natur[!is.na(natur$x)&!is.na(natur$y),]
  
  
  if(nrow(native)>0){
    
    native_points<-SpatialPoints(native[,1:2])
    crs(native_points)<-crs(region_shape)
    
    
    #which continent/region has native species in it
    region_native<-unique(na.omit(over(native_points,region_shape)))
    region_shape$native_subcount[region_shape$Name%in%region_native$Name]<-region_shape$natur_subcount[region_shape$Name%in%region_native$Name]+1
    region_native_shape<-region_shape[region_shape$Name%in%region_native$Name,]
    
    #which subregion has native species in it
    subregion_native<-unique(na.omit(over(native_points,shape)))
    shape$native_subcount[shape$LEVEL_4_NA%in%subregion_native$LEVEL_4_NA]<-shape$natur_subcount[shape$LEVEL_4_NA%in%subregion_native$LEVEL_4_NA]+1
    subregion_native_shape<-shape[shape$LEVEL_4_NA%in%subregion_native$LEVEL_4_NA,]
    
  }
  
  if(nrow(natur)>0){
    natur_points<-SpatialPoints(natur[,1:2])
    crs(natur_points)<-crs(region_shape)
    
    #which continent/region has naturalised species in it
    region_natur<-unique(na.omit(over(natur_points,region_shape)))
    region_shape$natur_subcount[region_shape$Name%in%region_natur$Name]<-region_shape$natur_subcount[region_shape$Name%in%region_natur$Name]+1
    region_natur_shape<-region_shape[region_shape$Name%in%region_natur$Name,]
    
    
    #which subregion has naturalised species in it
    subregion_natur<-unique(na.omit(over(natur_points,shape)))
    shape$natur_subcount[shape$LEVEL_4_NA%in%subregion_natur$LEVEL_4_NA]<-shape$natur_subcount[shape$LEVEL_4_NA%in%subregion_natur$LEVEL_4_NA]+1
    subregion_natur_shape<-shape[shape$LEVEL_4_NA%in%subregion_natur$LEVEL_4_NA,]
    
    if("Nunavut" %in% subregion_natur$LEVEL_4_NA){print("in nunavut!");nunavut_list<-c(nunavut_list,sp_name)}
    
    
  }
  
  
  #part 2: now we check whether it expanded or not, and add a count for any occurrences that count as an expansion
  
  niche_sp<-niche_shift[which(niche_shift$species.name ==sp_name),]
  #limit to areas where there has been significant expansion
  niche_sp<-niche_sp[which(niche_sp$expansion>0.1),]
  
  if(nrow(niche_sp)>0){
    
    #now we look up where the species was expanding
    for(j in 1:nrow(niche_sp)){
      regn<-niche_sp[j,"region"]
      
      pcaval <- paste0("./IntermediateOutput/plant_PCA_exp/3Var/NaturPCAvalues/",sp_name,regn,".csv")
      expco<-read.delim(pcaval, sep = "\t")
      
      expco<-expco[expco$exp==1,]
      
      if(nrow(expco)>0){
        exp_points<-SpatialPoints(expco[,1:2])
        crs(exp_points)<-crs(region_shape)
        
        #which continent/region has native species in it
        region_exp<-unique(na.omit(over(exp_points,region_shape)))
        region_shape$exp_subcount[region_shape$Name%in%region_exp$Name]<-region_shape$exp_subcount[region_shape$Name%in%region_exp$Name]+1
        region_exp_shape<-region_shape[region_shape$Name%in%region_exp$Name,]
        
        #which subregion has native species in it
        subregion_exp<-unique(na.omit(over(native_points,shape)))
        shape$exp_subcount[shape$LEVEL_4_NA%in%subregion_exp$LEVEL_4_NA]<-shape$exp_subcount[shape$LEVEL_4_NA%in%subregion_exp$LEVEL_4_NA]+1
        subregion_exp_shape<-shape[shape$LEVEL_4_NA%in%subregion_exp$LEVEL_4_NA,]
      }
    }
    
    
  }
  
  
  
  #remove doublecounts
  shape$native_subcount[shape$native_subcount>1]<-1
  shape$natur_subcount[shape$natur_subcount>1]<-1
  shape$exp_subcount[shape$exp_subcount>1]<-1
  
  region_shape$native_subcount[region_shape$native_subcount>1]<-1
  region_shape$natur_subcount[region_shape$natur_subcount>1]<-1
  region_shape$exp_subcount[region_shape$exp_subcount>1]<-1
  
  shape$native_count<-shape$native_count+shape$native_subcount
  shape$natur_count<-shape$natur_count+shape$natur_subcount
  shape$exp_count<-shape$exp_count+shape$exp_subcount
  
  region_shape$native_count<-region_shape$native_count+region_shape$native_subcount
  region_shape$natur_count<-region_shape$natur_count+region_shape$natur_subcount
  region_shape$exp_count<-region_shape$exp_count+region_shape$subexp_count
  
  
}

shapefile(region_shape, "./IntermediateOutput/ExpansionMapsSummaries/composite_region_map.shp", overwrite=TRUE)
shapefile(shape, "./IntermediateOutput/ExpansionMapsSummaries/composite_subregion_map.shp", overwrite=TRUE)


#can skip it necessary to this point
#shape<-shapefile("./IntermediateOutput/ExpansionMapsSummaries/composite_subregion_map.shp")

a<-gTouches(shape,byid=T)
alist<-apply(a, 1, function(r) any(r %in% TRUE))
mainland<-shape[alist,]
islands<-shape[!alist,]
islands$LEVEL_4_N

mainland<-mainland[which(mainland$LEVEL_4_N!="Haiti"&mainland$LEVEL_4_N!="Dominican Republic"),]
mainlandpart2<-shape[which(shape$LEVEL_4_N=="Ireland"|shape$LEVEL_4_N=="Northern Ireland"|
                             shape$LEVEL_4_N=="Great Britain"|shape$LEVEL_4_N=="Madagascar"|
                             shape$LEVEL_4_N=="New Zealand South"|shape$LEVEL_4_N=="New Zealand North"|
                             shape$ISO_COD=="JP"),]
shape<-bind(mainland,mainlandpart2)



shape<-spTransform(shape,CRS="+proj=moll +units=m +ellps=WGS84")

names(shape)[17]<-"native_count"
names(shape)[18]<-"natur_count"
names(shape)[19]<-"exp_count"


#####################################
###NOW START PLOTTING MAPS
####################################

#set palette
pal <- (wes_palette("Zissou1", 100, type = "continuous"))


rng = range(c(shape$native_count, shape$natur_count,shape$exp_count))
rng



###NATIVE
#plot native species on world map
shp_df <- broom::tidy(shape, region = "native_count")

shp_df$id[shp_df$id==0]<-NA



plotname<-paste("./FinalOutput/ExpansionMapsSummaries/subregion_nativeMap.pdf",sep="")
pdf(file=plotname,width=9,height=4)

map <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = as.numeric(id)), colour = "black") +
  theme_void()+
  scale_fill_gradientn(colours = pal,na.value = "white", limits=c(1, ceiling(rng[2])))+
  labs(fill = "# Species")
print(map)

dev.off()



###NATURALISED
#plot naturalised species on world map
shp_df <- broom::tidy(shape, region = "natur_count")

shp_df$id[shp_df$id==0]<-NA
#shp_df$id<-log10(as.numeric(shp_df$id))

plotname<-paste("./FinalOutput/ExpansionMapsSummaries/subregion_naturMap.pdf",sep="")
pdf(file=plotname,width=9,height=4)

map <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = as.numeric(id)), colour = "black") +
  theme_void()+
  scale_fill_gradientn(colours = pal,na.value = "white", limits=c(1, ceiling(rng[2])))+
  labs(fill = "# Species")

map

dev.off()



###EXPANSION
#plot expanding species on world map

shp_df <- broom::tidy(shape, region = "exp_count")

shp_df$id[shp_df$id==0]<-NA
#shp_df$id<-log10(as.numeric(shp_df$id))


plotname<-paste("./FinalOutput/ExpansionMapsSummaries/subregion_expMap.pdf",sep="")
pdf(file=plotname,width=9,height=4)

map <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = as.numeric(id)), colour = "black") +
  theme_void()+
  scale_fill_gradientn(colours = pal,na.value = "white")+
  labs(fill = "# Species")

map

dev.off()






