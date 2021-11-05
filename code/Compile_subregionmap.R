


#step one

#read in occurrence data used in analysis

#layer over subregion map

#run a count


rm(list=ls())

setwd("U:/Data/niche_shift_work/")

setwd("/Users/henry_hakkinen/Documents/ExeterWork/Data/niche_shift_work")

library(raster)
library(rworldmap)
library(rgeos)
library(RColorBrewer)


table_name<-"plant"

####plant paths
raster_path<-"exp_unfill_rasters000_3clim/"
range_path<-"U:/Data/species lists/plant_dist_data_0802/estimated_range_rasters/"

#####plant paths
plant_summary<-read.csv("analogue_disagg3Var_04032019/plant_data_analogue05042019.csv",stringsAsFactors = F)
sp_list<-unique(plant_summary$species_name)



#remove species we ruled out for whatever reason
d_values<-read.csv("plant_D_values000_3clim.csv",stringsAsFactors = F)


#define where dat is located:
data_path<-"analogue_disagg3Var_04032019/"

#define output folder for species rasters
output_path<-"Figures/3VarEXP/"
output_figs<-"Figures/3VarEXP/"

shape<-shapefile("U:/Data/countrydata/level4.shp")
shape<-shapefile("/Users/henry_hakkinen/Documents/ExeterWork/Data/countrydata/level4.shp")


newmap<-getMap(resolution="coarse")
world<-raster("U:/Data/Global_10min_grid.grd")

#load region info
region_shape <- shapefile("U:/Data/biogeographic_zonesV2.shp")
region_shape <- shapefile("/Users/henry_hakkinen/Documents/ExeterWork/Data/biogeographic_zonesV2.shp")


#define islands vs nonislands
shape<-shape[which(shape$LEVEL_4_NA!="Greenland"&shape$LEVEL_4_NA!="Antarctica"),]




region_shape$native_count<-0
region_shape$natur_count<-0

region_shape$native_subcount<-0
region_shape$natur_subcount<-0

shape$native_count<-0
shape$natur_count<-0

shape$native_subcount<-0
shape$natur_subcount<-0

failure<-c()
nunavut_list<-c()

x11()
i
#for (i in 1: length(sp_list)){
for (i in 1: 20){
  sp_name<-sp_list[i]
  print(sp_name)
  
  sp_rows<-plant_summary[plant_summary$species_name==sp_name,]
  
  
  #reset counter to avoid double counts
  region_shape$native_subcount<-0
  region_shape$natur_subcount<-0
  
  shape$native_subcount<-0
  shape$natur_subcount<-0
  
  if(nrow(sp_rows)==0){break()}

  for(q in 1:nrow(sp_rows)){
    
    region<-sp_rows[q,2]
    
    file_natur=paste(data_path,region,"/",sp_name,region,'_natur.csv', sep='')
    file_native=paste(data_path,region,"/",sp_name,region,'_native.csv', sep='')
    
    
    native<-read.csv(file_native)
    natur<-read.csv(file_natur)

    colnames(native)<-c("x","y")
    colnames(natur)<-c("x","y")

    native<-native[!is.na(native$x)&!is.na(native$y),]
    natur<-natur[!is.na(natur$x)&!is.na(natur$y),]
    
    #if(nrow(sp_points)>0){
    
    native_points<-SpatialPoints(native)
    crs(native_points)<-crs(region_shape)
    natur_points<-SpatialPoints(natur)
    crs(natur_points)<-crs(region_shape)
    

    

    #which continent/region has native species in it
    region_native<-unique(na.omit(over(native_points,region_shape)))
    region_shape$native_subcount[region_shape$Name%in%region_native$Name]<-region_shape$natur_subcount[region_shape$Name%in%region_native$Name]+1
    region_native_shape<-region_shape[region_shape$Name%in%region_native$Name,]
    

    #which continent/region has naturalised species in it
    region_natur<-unique(na.omit(over(natur_points,region_shape)))
    region_shape$natur_subcount[region_shape$Name%in%region_natur$Name]<-region_shape$natur_subcount[region_shape$Name%in%region_natur$Name]+1
    region_natur_shape<-region_shape[region_shape$Name%in%region_natur$Name,]
    
    #which subregion has native species in it
    subregion_native<-unique(na.omit(over(native_points,shape)))
    shape$native_subcount[shape$LEVEL_4_NA%in%subregion_native$LEVEL_4_NA]<-shape$natur_subcount[shape$LEVEL_4_NA%in%subregion_native$LEVEL_4_NA]+1
    subregion_native_shape<-shape[shape$LEVEL_4_NA%in%subregion_native$LEVEL_4_NA,]
    
    
    #which subregion has naturalised species in it
    subregion_natur<-unique(na.omit(over(natur_points,shape)))
    shape$natur_subcount[shape$LEVEL_4_NA%in%subregion_natur$LEVEL_4_NA]<-shape$natur_subcount[shape$LEVEL_4_NA%in%subregion_natur$LEVEL_4_NA]+1
    subregion_natur_shape<-shape[shape$LEVEL_4_NA%in%subregion_natur$LEVEL_4_NA,]
    
    if("Nunavut" %in% subregion_natur$LEVEL_4_NA){print("in nunavut!");nunavut_list<-c(nunavut_list,sp_name)}
    
    
    par(mfrow=c(1,2))
    plot(newmap,main=paste(sp_name,": Native"))
    plot(region_native_shape,col="green",add=T)
    plot(subregion_native_shape,col="red",add=T)
    points(native$x,native$y,col="cyan",pch=".",cex=3)
    
    plot(newmap,main=paste(sp_name,": Naturalised",region))
    plot(region_natur_shape,col="green",add=T)
    plot(subregion_natur_shape,col="red",add=T)
    points(natur$x,natur$y,col="cyan",pch=".",cex=3)


    
  }
  
  #remove doublecounts
  shape$native_subcount[shape$native_subcount>1]<-1
  shape$natur_subcount[shape$natur_subcount>1]<-1
  
  region_shape$native_subcount[region_shape$native_subcount>1]<-1
  region_shape$natur_subcount[region_shape$natur_subcount>1]<-1
  
  shape$native_count<-shape$native_count+shape$native_subcount
  shape$natur_count<-shape$natur_count+shape$natur_subcount
  
  region_shape$native_count<-region_shape$native_count+region_shape$native_subcount
  region_shape$natur_count<-region_shape$natur_count+region_shape$natur_subcount

  
}

shapefile(region_shape, "composite_region_map.shp", overwrite=TRUE)
shapefile(shape, "composite_subregion_map.shp", overwrite=TRUE)



shape<-shapefile("composite_subregion_map.shp")

shape<-shape[which(shape$LEVEL_4_N!="Greenland"&shape$LEVEL_4_N!="Antarctica"),]
shape<-shape[which(shape$LEVEL_4_N!="Antarctica"),]



shape<-spTransform(shape,CRS="+proj=moll +units=m +ellps=WGS84")


x11()
names(shape)[17]<-"native_count"
names(shape)[18]<-"natur_count"



table_name<-"Plant"
rbPal <-  colorRampPalette(rev(brewer.pal(11,"Spectral")))


######plot expansion by region


region_shape$native_Col <- rbPal(10)[as.numeric(cut(region_shape$native_count,breaks = 10))]
region_shape$native_Col[region_shape$native_count==0]<-"#FFFFFF"

region_shape$natur_Col <- rbPal(10)[as.numeric(cut(region_shape$natur_count,breaks = 10))]
region_shape$natur_Col[region_shape$natur_count==0]<-"#FFFFFF"

#range of species from 0-110, cut at 10 interval
summary(shape$native_count)
summary(shape$natur_count)


#shape$native_Col <- rbPal(10)[as.numeric(cut(shape$native_count,breaks = 10))]
shape$native_Col <- rbPal(11)[as.numeric(cut(shape$native_count,breaks=seq(0,110,by=10)))]

shape$native_Col[shape$native_count==0]<-"#FFFFFF"

shape$natur_Col <- rbPal(11)[as.numeric(cut(shape$natur_count,breaks=seq(0,110,by=10)))]
shape$natur_Col[shape$natur_count==0]<-"#FFFFFF"



plotname<-paste("Figures/",table_name,"_subregion_rangesummary.pdf",sep="")

pdf(file=plotname,width=15,height=9)

plot(shape,col=shape$native_Col,main=paste(table_name,"Native"))
plot(shape,col=shape$natur_Col,main=paste(table_name,"Naturalised"))
dev.off()

plotname<-paste("Figures/subregion_legend.pdf",sep="")

pdf(file=plotname,width=2,height=5)
legend_image <- as.raster(matrix(rev(rbPal(10)), ncol=1))

plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,100,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()


#basically do the same with exp data
head(plant_summary)
sp_name<-"Abutilon theophrasti"
region<-"Nearctic"

test<-read.csv(paste("regional/DirData/",sp_name,'_',region,"_dir.csv",sep=""))
head(test)



