
rm(list=ls())


setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/IntermediateOutput/plant_PCA_exp/4Var/rstates_center/")

csvfile4<- paste("../../../Rstate_compile/4Var_plant_D_shiftvalues_zcor_center.txt",sep="")


yooo<-list.files(pattern="*_nicheshift.R")
fail<-list.files(path="flaws",pattern="*_nicheshift.R")

yooo<-yooo[!(yooo %in% fail)]
yooo

load(yooo[1])
D.df_total<-D.df
length(D.df_total)
D.df_total

for (n in yooo){
  print(n)
  load(n)
  i<-D.df

  #print(i[1])
  if(i[1] != "" & i[1] != "not enough known occurrences to calculate niche"){D.df_total<-rbind(D.df_total,i)}
}
length(i)

#Woodwardia radicans
D.df_total<-D.df_total[-c(1),]
names(D.df_total)
head(D.df_total)

#D.df<-do.call(rbind, lapply(seq_along(yooo), function(i){
#data.frame(CLUSTER=i, yooo[[i]])
#}))




write.table(D.df_total, csvfile4, row.names = F, sep = "\t")

D.df_convert<-read.table(csvfile4,row.names=NULL)

D.df_convert$species.name<-paste(D.df_convert$row.names,D.df_convert$species.name)
D.df_convert$row.names<-NULL

write.table(D.df_convert, csvfile4, row.names = F, sep = "\t")

shiftdf<-read.table(csvfile4,header=T)
tail(shiftdf)


library(raster)
#merge raster bricks
setwd("../shift_PCArasters")

remove_na<-function(x){
  #print(x)
  #print(class(x))
  #x$expansion[is.na(x$expansion)] <- 0
  #x$unfilling[is.na(x$unfilling) m] <- 0
  #x$natur_occ[is.na(x$natur_occ)] <- 0
  x[is.na(x)]<-0
  return(x)
}
brick_list<-list.files(pattern="*_PCAraster.tif")


raster_pile<-brick(brick_list[1])
raster_pile<-remove_na(raster_pile)

for (i in 2:length(brick_list)){
  print(i)
  raster_new<-brick(brick_list[i])
  raster_new<-remove_na(raster_new)
  
  raster_pile<-Reduce("+",list(raster_pile,raster_new))
  
  #close file, to avoid memory pileup
  removeTmpFiles(h=0)
}

names(raster_pile)
x11()

library(RColorBrewer)
cols <- rev(brewer.pal(11,"Spectral"))
raster_plot<-raster_pile
raster_plot[raster_plot==0]<-NA
plot(raster_plot$layer,col=cols)

csvfile5<- paste("../../../Rstate_compile/4Var_composite_rasterEXP_4clim_zcor.tif",sep="")
writeRaster(raster_pile, csvfile5, options= c("COMPRESSION=LZW","INTERLEAVE=PIXEL"), overwrite=TRUE)
raster_pile<-raster(csvfile5);names(raster_pile)<-"layer"





###OLD#####
###EXP SECTION
setwd("U:/Data/niche_shift_work/stabilityRegional/")
library(raster)
expRaster<-raster("composite_rasterSTAB_3clim_CENTER.tif")
cols <- brewer.pal(11,"Spectral")
raster_plot<-expRaster
raster_plot[raster_plot==0]<-NA
raster_plot<-raster_plot/cellStats(raster_plot,stat=max)

globalRasterR<-raster("../glob_clim_raster.tif")
globalRaster <-globalRasterR
globalRaster[globalRaster==0]<-NA
globalRaster<-globalRaster/cellStats(globalRaster,stat=max)

x11()
par(mfrow=c(2,2))
plot(globalRaster,col=cols)
plot(raster_plot,col=cols)

origin(raster_plot) <- origin(globalRaster)
extent(raster_plot) <- extent(globalRaster)
raster_plot[is.na(globalRaster)] <-NA
globalRaster[is.na(raster_plot)] <-NA

plot(globalRaster,col="green")
plot(raster_plot,col="red",add=T)
plot(globalRaster,col="green",add=T)

raster_int<-raster_plot/globalRaster
raster_int[raster_int>60]<-NA



plot(raster_int,col=cols)
plot(log10(raster_int),col=cols)



