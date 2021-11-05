rm(list=ls())
library(RColorBrewer)
library(CircStats)
library(scales)



setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/")


#plant_summary<-read.csv("plant_data_summary.csv")
#plant_summary_new<-read.csv("species lists/plant_dist_data_0802/plant_data_summary0607.csv")

plant_summary<-read.csv("IntermediateOutput/PCA_find_analogue_byRegion/plant_data_analogue05042019.csv",stringsAsFactors = F)


#niche expansion direction based on 3 variables
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")
shift_dyn2<-read.delim(csvfile4, sep="\t")


rad.ang <- function(x1, x2){
  dx = x1[1] - x2[1]
  dy = x1[2] - x2[2] 
  
  if(!is.numeric(dx)){dx<-as.numeric(dx)}
  if(!is.numeric(dy)){dy<-as.numeric(dy)}
  
  theta = atan2(dy,dx)
  
  if(theta<0){ theta<- theta+ 2*pi }
  return(theta)
} 
load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecip")
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))


arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))

rownames(arr_tab)<-c("TMax","TMin","Precip")


#make a second table to describe the 4 quarters of climate circle
quar_tab<-arr_tab[c(1,3),]

quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)

rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")


ampF<-2.2
ampF2<-1.5

plotname=paste("./IntermediateOutput/AdHoc_nicheshiftPlot/NicheVersusExp.png",sep="")
png(file=plotname, width=1080, height=360)

par(mfrow=c(1,3))
plot(shift_dyn2$native_center_tmax, shift_dyn2$expansion, 
     xlab="Central TMax of native niche", ylab="Proportion of expansion",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
plot(shift_dyn2$native_center_tmin, shift_dyn2$expansion, 
     xlab="Central TMin of native niche", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
plot(shift_dyn2$native_center_precip, shift_dyn2$expansion, 
     xlab="Central Precip of native niche", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
dev.off()


plotname=paste("./IntermediateOutput/AdHoc_nicheshiftPlot/NicheVersusDir.png",sep="")
png(file=plotname, width=1080, height=360)

par(mfrow=c(1,3))
plot(shift_dyn2$native_center_tmax, shift_dyn2$exp_dir1, 
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF,
     xlab="Central TMax of native niche", ylab="Direction of expansion")
points(shift_dyn2$native_center_tmax, shift_dyn2$exp_dir2,     
       col=alpha("black", 0.3), pch=16)

plot(shift_dyn2$native_center_tmin, shift_dyn2$exp_dir1, 
     xlab="Central TMin of native niche", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
points(shift_dyn2$native_center_tmin, shift_dyn2$exp_dir2,
       col=alpha("black", 0.3), pch=16)

plot(shift_dyn2$native_center_precip, shift_dyn2$exp_dir1, 
     xlab="Central Precip of native niche", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
points(shift_dyn2$native_center_precip, shift_dyn2$exp_dir2,
       col=alpha("black", 0.3), pch=16)

dev.off()






