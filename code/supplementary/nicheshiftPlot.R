##############################################
### nicheshiftPlot.R ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 

# as a supplementary analysis we want to compare conditions at the centre of the native niche and the proportion and direction of expansion

# in plot 1 we plot what the tmax, tmin and precip at the centre of the native niche is against the proportion of niche expansion
# do species from warm/cold/dry/wet places expand more?

# plot what the tmax, tmin and precip at the centre of the native niche is against the proportion of niche expansion
# do species from warm/cold/dry/wet places expand in certain directions?

# output from this is used in Appendix Figure S1.10

### ###


rm(list=ls())



###################################
#set paths and variables
##################################


#set to current repo
setwd("DIRECTORY_HERE")


library(RColorBrewer)
library(CircStats)
library(scales)
library(ecospat)

source("./code./functions/circ_plot_functions.R")


#load output from plant_pca_expansion and Rstate_compile.R
#niche expansion direction based on 3 variables
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")
shift_dyn2<-read.delim(csvfile4, sep="\t")


#read in information on the PCA and what the loadings/directions are. Should match the csvfile4 spec above
#this file is generated in "plant_PCA_expand.R"
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



###################################
#plot niche shift and expansion versus niche info
##################################

#set label size
ampF<-2.2
ampF2<-1.5


#plot what the tmax, tmin and precip of the native niche is against the proportion of niche expansion
#do species from warm/cold/dry/wet places expand more?
plotname=paste("./FinalOutput/supplementary/nicheshiftPlot/NicheVersusExp.png",sep="")
png(file=plotname, width=1080, height=360)

#adjust margin
#par(mar = c(bottom, left, top, right)

par(mfrow=c(1,3),
    mar = c(5.1, 4.4, 4.1, 2.1))
plot(shift_dyn2$native_center_tmax/10, shift_dyn2$expansion, 
     xlab="Central TMax of native niche (°C)", ylab="Proportion of expansion",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
plot(shift_dyn2$native_center_tmin/10, shift_dyn2$expansion, 
     xlab="Central TMin of native niche (°C)", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
plot(shift_dyn2$native_center_precip, shift_dyn2$expansion, 
     xlab="Central Precip of native niche (mm)", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
dev.off()



#plot what the tmax, tmin and precip of the native niche is against the proportion of niche expansion
#do species from warm/cold/dry/wet places expand in certain directions?
plotname=paste("./FinalOutput/supplementary/nicheshiftPlot/NicheVersusDir.png",sep="")
png(file=plotname, width=1080, height=360)


par(mfrow=c(1,3),
    mar = c(5.1, 4.4, 4.1, 2.1))
plot(shift_dyn2$native_center_tmax/10, shift_dyn2$exp_dir1, 
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF,
     xlab="Central TMax of native niche (°C)", ylab="Direction of expansion")
points(shift_dyn2$native_center_tmax/10, shift_dyn2$exp_dir2,     
       col=alpha("black", 0.3), pch=16)

plot(shift_dyn2$native_center_tmin/10, shift_dyn2$exp_dir1, 
     xlab="Central TMin of native niche (°C)", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
points(shift_dyn2$native_center_tmin/10, shift_dyn2$exp_dir2,
       col=alpha("black", 0.3), pch=16)

plot(shift_dyn2$native_center_precip, shift_dyn2$exp_dir1, 
     xlab="Central Precip of native niche (mm)", ylab="",
     col=alpha("black", 0.3), pch=16, cex.axis=ampF2, cex.lab=ampF)
points(shift_dyn2$native_center_precip, shift_dyn2$exp_dir2,
       col=alpha("black", 0.3), pch=16)

dev.off()






