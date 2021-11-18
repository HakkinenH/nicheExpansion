
##############################################
### PLOT NICHE FILLING HISTOGRAMS ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 
# In plant_PCA_nichefilling we calculated how much of their potential niche each species/population has colonised in their naturalised range
# E.g. we work out how much climate is available in the "wettest" quarter of climate, and how many cells are occupied in that quarter
# we can then work out what proportion of all occurrences are in the "wettest quarter" and also how much of the available "wet" climate is occupied
# We plot this in several ways, firstly we make some plots of where in their niche species have naturalised
# this output is used for Figure 4 in the manuscript

#we then do some subsidiary plots for the supplementary material
# we then remove any species/population that showed niche expansion, just to show it's not being driven by expanding species
# used in Appendix Figure S1.3 
# then we split niche filling  by region
# used in Appendix Appendix Figure S1.9

### ###


rm(list=ls())


#set to path of local repo
setwd("DIRECTORY_HERE")

library("ecospat")

source("./code./functions/circ_plot_functions.R")


###################################
#load files for processing
##################################

output<-"3VarEXP"

#niche filling direction based on 3 variables
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_Fillingvalues_zcor_center.txt",sep="")

shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")



#read in information on the PCA and what the loadings/directions are. Should match the csvfile4 spec above
#this file is generated in "plant_PCA_expand.R"
#for 3 var
load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecip")

ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))

colnames(arr_tab)<-c("x","y")
#table with points which give the direction as compared to the origin (0,0)
#find the angle of direction for each PCA component
arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))
rownames(arr_tab)<-c("TMax","TMin","Precip")



####



names(shiftdf)

##the old way to plot
test<-shiftdf[shiftdf$expansion>=0.1,]

x11(width=6,height=4)
par(mfrow=c(1,3))
hist(shiftdf$proportion_precip)
hist(test$proportion_precip,col="blue",add=T)

hist(shiftdf$proportion_tmax)
hist(test$proportion_tmax,col="blue",add=T)
hist(shiftdf$proportion_tmin)
hist(test$proportion_tmin,col="blue",add=T)

names(shiftdf)


##the new way, split into four windows (wet, hot, dry, cold)
pdf(file="./FinalOutput/PlotNichefillingHists/NicheFilling_hists.pdf", width=6, height=6)

par(mfrow=c(2,2))
hist(shiftdf$proportion_precip,main="Wet",xlab="Proportion of naturalised occurrence")
hist(shiftdf$proportion_tmax,main="Hot",xlab="Proportion of naturalised occurrence")

hist(shiftdf$proportion_dry,main="Dry",xlab="Proportion of naturalised occurrence")
hist(shiftdf$proportion_cold,main="Cold",xlab="Proportion of naturalised occurrence")

dev.off()



expdf<-shiftdf[which(shiftdf$expansion<=0.1),]
#for supp material i need when where only expanding species are included
pdf(file="./FinalOutput/Supplementary/PlotNichefillingHists/NicheFilling_hists_EXPONLY.pdf", width=6, height=6)

par(mfrow=c(2,2))
hist(expdf$proportion_precip,main="Wet",xlab="Proportion of naturalised occurrence")
hist(expdf$proportion_tmax,main="Hot",xlab="Proportion of naturalised occurrence")

hist(expdf$proportion_dry,main="Dry",xlab="Proportion of naturalised occurrence")
hist(expdf$proportion_cold,main="Cold",xlab="Proportion of naturalised occurrence")

dev.off()

median(shiftdf$proportion_precip,na.rm = T)
median(shiftdf$proportion_tmax,na.rm = T)
median(shiftdf$proportion_dry,na.rm = T)
median(shiftdf$proportion_cold,na.rm = T)


#loop all regions
reg_list<-unique(shiftdf$region)
reg_list<-sort(reg_list)


#work in sections
#1:6 then 7:12


for(sec in 1:2){
  if(sec==1){reg_now<-reg_list[1:6]}
  if(sec==2){reg_now<-reg_list[7:12]}

  
  
  pdf(file=paste0("./FinalOutput/Supplementary/PlotNichefillingHists/NicheFilling_regionalhists",sec,".pdf"), width=12, height=18)
  par(mfrow=c(6,4))
  for (i in reg_now){
    print(i)
    shift_sub<-shiftdf[shiftdf$region==i,]
    
    hist(shift_sub$proportion_precip,main=paste0(i,": Wet"),xlab="Proportion of naturalised occurrence",cex.lab=1.5,cex.main=1.5)
    hist(shift_sub$proportion_tmax,main=paste0(i,": Hot"),xlab="Proportion of naturalised occurrence",cex.lab=1.5,cex.main=1.5)
    
    hist(shift_sub$proportion_dry,main=paste0(i,": Dry"),xlab="Proportion of naturalised occurrence",cex.lab=1.5,cex.main=1.5)
    
    hist(shift_sub$proportion_cold,main=paste0(i,": Cold"),xlab="Proportion of naturalised occurrence",cex.lab=1.5,cex.main=1.5)
    
  }
  dev.off()  
}


#Oc = Dry
#Pan = Dry/Col
#Afrotropical = Cold
#Neotropical = Wet/Cold
#Australian = Cold Wet




