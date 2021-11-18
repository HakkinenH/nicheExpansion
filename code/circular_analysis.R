##############################################
### circular_analysis ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 
# After running plant_PCA_expand.R and Rstate_compile we should have a compiled dataframe with information on 
  #1) whether a species expanded (i.e. >0.1 proportion of their naturalised niche is outside of native niche)
  #2) what direction (if any) their niche expanded
#These are statistics for individual species and naturalised populations
#we now want to find trends across ALL populations and species
#This file takes all available information on niche expansions (or lack thereof) and produces some summary statistics, plots and output
#specifically it will:
  #1) run a circular parametric model to find out if  species expand more often in certain directions
  #2) is direction expansion predicted by a particular climatic direction?
  #3) get the residuals from the non-parametric model and find the circular R-squared
#you can optionally use the "PA_test" to account for pseudo-replication in the data (since species are typically repeated)

#The output from this script is used in Table 1, Table 2, Figure 3, Appendix Figure S1.4, Appendix Figure S1.5, Appendix Figure S1.6

### ###

rm(list=ls())



###################################
#set paths and variables
##################################
setwd("DIRECTORY_HERE")

library("ecospat")
library("NPCirc")
library("CircMLE")



source("./code./functions/circ_plot_functions.R")


#the following files are included in the repo and can be used as examples or to confirm out published analysis
#choose one for the main analysis

#POSSIBLE FILE LIST:
  # niche expansion direction based on 3 variables
  # 3Var_plant_D_shiftvalues_zcor_center.txt

  # niche expansion direction based on 2 variables (precipitation and maximum temperature)
  # 2Var_plant_D_shiftvalues_zcor_center.txt

  # niche expansion direction based on 4 variables (precipitation, precipitation seasonality and minimum temperature and maximum temperature)
  # 4Var_plant_D_shiftvalues_zcor_center.txt

  # direction of unfilling (where species succeed/fail to fail their niche)
  # 3Var0304plant_D_shiftvalues_zcor_center.txt




#select which output file you would like to analyse (published version is niche expansion direction based on 3 variables)
csvfile4 <- paste0("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt")

#choose the name of the output file
output<-"3VarEXP"

PA_test<-"FALSE"


#read in the main file with results
shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")

#check it looks right
head(shiftdf)



#read in information on the PCA and what the loadings/directions are. Should match the csvfile4 spec above
#this file is generated in "plant_PCA_expand.R"
#for 2var
#load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TmaxPrecip")
#for 3 var
load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecip")
#for 4 var
#load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecipPrSeas")






###################################
#BASIC EXPLORATION AND CHECKS
##################################

#check the names of the datafile
names(shiftdf)



#x11()
#basic plot to show whether expansion/filling happens most in wettest/highest tmin/highest tmax portion of climate niche
par(mfrow=c(1,3))
hist(shiftdf$proportion_precip)
hist(shiftdf$proportion_tmin)
hist(shiftdf$proportion_tmax)



#cut out invalid options, remove any species with no valid expansion
shiftdf<-shiftdf[which(shiftdf$max_exp_dist1!=-Inf | is.na(shiftdf$max_exp_dist1)),]
shiftdf<-shiftdf[which(shiftdf$max_exp_dist2!=-Inf | is.na(shiftdf$max_exp_dist2)),]

#we occasionally get expansion where we don't have observations because of the way kernel density is calculated. We cut these out
shiftdf<-shiftdf[!is.na(shiftdf$exp_dir1)|!is.na(shiftdf$exp_dir2),]

#how many (valid) species/populations do we have that expand significantly? (i.e. over 10% of their naturalised density lies outside of the native niche?)
sum(shiftdf$expansion>0.1)


#only keep these species/populations, we will now analyse species that expand only
shiftdf<-shiftdf[shiftdf$expansion>0.1,]
#split by region
table(shiftdf$region)
#split by circular model
table(shiftdf$circular_model)

#build PCA circle and take information on directions and climate
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))
colnames(arr_tab)<-c("x","y")


#table with points which give the direction as compared to the origin (0,0)
#find the angle of direction for each PCA component
#variable depending on how many variables are being put in
if (grepl("2Var", output)){
  arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))
  arr_tab
  rownames(arr_tab)<-c("TMax","Precip")
  #make a second table to describe the 4 quarters of climate circle
  quar_tab<-arr_tab[c(1,2),]
  
  quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
  quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)
  
  rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")
  
  
}else if (grepl("3Var", output)){
  
  arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))
  arr_tab
  rownames(arr_tab)<-c("TMax","TMin","Precip")
  
  
  #make a second table to describe the 4 quarters of climate circle
  quar_tab<-arr_tab[c(1,3),]
  
  quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
  quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)
  
  rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")
  
  
}else{
  arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))
  arr_tab
  rownames(arr_tab)<-c("TMax","TMin","Precip", "PrecipSeason")
  
  
  #make a second table to describe the 4 quarters of climate circle
  quar_tab<-arr_tab[c(1,3),]
  
  quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
  quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)
  
  rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")
  
  
}



###################################
#PREPARE CIRCULAR DATA
##################################



#subselect relevant columns from data frame
shiftdf_1<-shiftdf[,c("median_exp_dist1",
                      "exp_dir1",
                      "region",
                      "native_center_tmax",
                      "native_center_tmin",
                      "native_center_precip",
                      "species.name","expansion")]
shiftdf_2<-shiftdf[,c("median_exp_dist2",
                      "exp_dir2",
                      "region",
                      "native_center_tmax",
                      "native_center_tmin",
                      "native_center_precip",
                      "species.name","expansion")]

names(shiftdf_2)<-names(shiftdf_1)

#bind rows together so all information can be compared together
shiftdf_1<-rbind(shiftdf_1,shiftdf_2)

#make sure all rows have valid data
shiftdf_1<-shiftdf_1[complete.cases(shiftdf_1),]
#circular models don't constrain results to 0-2pi, they will include e.g. 5pi. Adjust so everything is on same scale
shiftdf_1[shiftdf_1[,2]>2*pi,2] <- shiftdf_1[shiftdf_1[,2]>2*pi,2] - 2*pi


#if we want to select species at random to account for pseudo-replication
#randomly takes one record per species name, discards the rest
if(PA_test){

  shiftdf_1R <- shiftdf_1[0,]
  
  sp_un<-unique(shiftdf_1$species.name)
  
  for(i in 1:length(sp_un)){
    sp_now<-sp_un[i]
    
    sp_shift<-shiftdf_1[shiftdf_1$species.name == sp_now,]
    sp_keep<-sp_shift[sample(1:nrow(sp_shift), 1),]
    
    shiftdf_1R<-rbind(shiftdf_1R, sp_keep)
  }
  shiftdf_1<-shiftdf_1R
}


#take the direction of expansion column and make it circular format
dir<-as.circular(shiftdf_1[,2], units="radians",rotation='counter')

#set which column has distance information in it
dist<-shiftdf_1[,1]




###############################################
###circular model, do species expand more often in certain directions?
################################################


#run a parametric circular model on direction of expansion
circ_res<-circ_mle(dir)

#make a circular histogram of direction of expansion
plotname=paste("./FinalOutput/circular_analysis/",output,"/DirectionRose.pdf",sep="")
pdf(file=plotname, width=5, height=5)

plot_circMLE.custom(dir, circ_res,shrink=1.4)
axis.circular(at=circular(seq(0, 2*pi-pi/2, pi/2)), 
              labels=c("0",expression(pi/2),expression(pi),expression(3*pi/2)))
s.corcircle(arr_tab[,c(1,2)],grid=F,add=T)

dev.off()



#how far is mean angles of expansion from main climate directions?
#if the result is close to 0 this means the average direction of expansion is in that direction
#if the result is close to 2pi then it is the OPPOSITE direction


#take mean direction
dir1<-circ_res$results[1,2]
dir2<-circ_res$results[1,5]

find_CLIMangle(dir1,dir2,"Wetter")
find_CLIMangle(dir1,dir2,"Drier")
find_CLIMangle(dir1,dir2,"Warmer")
find_CLIMangle(dir1,dir2,"Colder")

find_CLIMangle(dir1,NA,"Wetter")[1]
find_CLIMangle(dir1,NA,"Drier")[1]
find_CLIMangle(dir1,NA,"Warmer")[1]
find_CLIMangle(dir1,NA,"Colder")[1]

find_CLIMangle(dir2,NA,"Wetter")[1]
find_CLIMangle(dir2,NA,"Drier")[1]
find_CLIMangle(dir2,NA,"Warmer")[1]
find_CLIMangle(dir2,NA,"Colder")[1]



#when you run circ_mle(), several alternative models are compared, if multiple models are viable we can investigate
#in the example data (3Var_plant_D_shiftvalues_zcor_center.txt), the second model is also viable
dir1<-circ_res$results[2,2]
dir2<-circ_res$results[2,5]

find_CLIMangle(dir1,dir2,"Wetter")
find_CLIMangle(dir1,dir2,"Drier")
find_CLIMangle(dir1,dir2,"Warmer")
find_CLIMangle(dir1,dir2,"Colder")

#third model is also viable
dir1<-circ_res$results[3,2]
dir2<-circ_res$results[3,5]

find_CLIMangle(dir1,dir2,"Wetter")
find_CLIMangle(dir1,dir2,"Drier")
find_CLIMangle(dir1,dir2,"Warmer")
find_CLIMangle(dir1,dir2,"Colder")



#what is the best model to describe expansion?
circ_model<-circ_res$bestmodel
circ_rstat<-circ_res$rt[1]
circ_pvalue<-circ_res$rt[2]
circ_res

#save some key statistics on the best circular model
write.csv(circ_res$results,paste0("./FinalOutput/circular_analysis/",output,"/circ_summary.csv"))


#OPTIONAL: alternative approach with curve fitting
#provides a different way to analyse and visualise results
library(mclust)

mod4 <- densityMclust(shiftdf_1$exp_dir1)

summary(mod4)

plotname=paste0("./FinalOutput/circular_analysis/",output,"/CurveFit_Direction.pdf")
pdf(file=plotname, width=9, height=7)

plot(mod4, what = "density", data = shiftdf_1$exp_dir1, breaks = 20,xlab="Direction of niche shift")
abline(v=arr_tab[3,3])
abline(v=arr_tab[3,3]+pi)
dev.off()



###############################################
###circular model, is the distance of expansion predicted by direction?
################################################



#run a non-parametric circular regression, with kernel smoothing
#we use 2 alternative smoothing methods, Nadaraya-Watson and Local-Linear
#these 2 methods should give similar results
estNW <- kern.reg.circ.lin(dir, dist, method="NW")
estLL <- kern.reg.circ.lin(dir, dist, method="LL")


#make a plot of this regression
plotname=paste("./FinalOutput/circular_analysis/",output,"/DirectionMagnitude_reg.pdf",sep="")
pdf(file=plotname, width=4, height=4)


res<-plot(estLL, plot.type="circle", points.plot=TRUE, show.radial.grid=F, grid.col="black",
          labels=c("0",expression(pi/2),expression(pi),expression(3*pi/2)),
          label.pos=seq(0,7*pi/2,by=pi/2), zero=0, clockwise=FALSE,cex=0.6)

lines(estLL, plot.type="circle", plot.info=res, line.col=2)
#s.corcircle(arr_tab[,c(1,2)]*4,grid=F,add=T)


dev.off()

#plot differently, can help with visualising
res<-plot(estNW, plot.type="line", points.plot=TRUE)



#we want to get some information on which direction has maximum expansion, and what direction that's in

#get data on the output of the regression. x is angle, y is distance
estNW2<-as.data.frame(estNW[c("x","y")])
#find where the maximum distance is
estNWmax<-estNW2[which.max(estNW2$y),]
#where are the 90th percentile points?
estNW2<-estNW2[which(estNW2$y>quantile(estNW$y, 0.9)),]

#what direction is this maximum distance of expansion in?
#find_CLIMangle(estNWmax$x,NA,"Precip")
#find_CLIMangle(estNWmax$x,NA,"TMin")
#find_CLIMangle(estNWmax$x,NA,"TMax")

find_CLIMangle(estNWmax$x,NA,"Wetter")
find_CLIMangle(estNWmax$x,NA,"Drier")
find_CLIMangle(estNWmax$x,NA,"Warmer")
find_CLIMangle(estNWmax$x,NA,"Colder")


#plot some additional results based on this maximum distance
#plot where on the circle the area of maxmimum distance was, and where the angle of precipitation is
plotname=paste("./FinalOutput/circular_analysis/",output,"/DirectionMagnitude_Suppmaxplot.pdf",sep="")
pdf(file=plotname, width=4, height=4)


res<-plot(estLL, plot.type="circle", points.plot=TRUE,
          #labels=c("N","NE","E","SE","S","SO","O","NO"),
          label.pos=seq(0,7*pi/4,by=pi/4), zero=0, clockwise=FALSE,cex=0.6)
points(estNW2$x,estNW2$y,col="red",cex=0.75)
points(estNWmax$x,estNWmax$y,col="blue",cex=1.2)
lines(estLL, plot.type="circle", plot.info=res, line.col=1)
arrows.circular(arr_tab["Precip","angle"], length=0.1)
#s.corcircle(arr_tab[,c(1,2)]*4,grid=F,add=T)

dev.off()


###########################################
##### get the residuals from non-parametric regression and find the circular R-squared
###########################################


row.names(shiftdf_1)<-1:nrow(shiftdf_1)

#Step 1: find maximum residuals when gradient is set to 0 at mean (intercept)
restotal<-euc.dist.center(shiftdf_1[,c(2,1)], data.frame(shiftdf_1[,2],mean(shiftdf_1$median_exp_dist1)))

#step 2: find our residuals right now
#get fitted values
fitted<-data.frame(x=as.numeric(estNW$x),y=as.numeric(estNW$y))

#do some rounding, we need to find the nearest fitted points to each true data point
shift_val<-shiftdf_1[,c(2,1)]
shift_val$xround2<-round(shift_val$exp_dir1, digits=2)
shift_val$xround1<-round(shift_val$exp_dir1, digits=1)

fitted$x_round2<-round(fitted$x,digits=2)
fitted$x_round1<-round(fitted$x,digits=1)

sum(shift_val$xround %in% fitted$x_round)
newdf<-data.frame(x2=shift_val$xround2,x1=shift_val$xround1,y=shift_val$median_exp_dist1,
                  x_pred2=shift_val$xround2,x_pred1=shift_val$xround1)
newdf$y_pred<-NA
head(newdf)

for(i in 1:nrow(newdf)){
  xnow<-newdf[i,4]
  #print(i)
  #print(newdf[i,])
  
  y_fit<-fitted[which(fitted$x_round2==xnow),2] 
  #print("attempt 1:")
  #print(fitted[which(fitted$x_round2==xnow),] )
  
  if(length(y_fit)==0){
    xnow<-newdf[i,5]
    y_fit<-fitted[which(fitted$x_round1==xnow),2] 
    #print("attempt 2:")
    #print(fitted[which(fitted$x_round1==xnow),] )
  }
  
  #bad line of code, correct later!
  newdf[i,"y_pred"]<-y_fit[1]
}

#sum square of residuals with new model
resfitted<-euc.dist.center(newdf[,c(1,3)],newdf[,c(4,6)])
plot(newdf[,1],newdf[,3])
points(newdf[,4],newdf[,6],col="red")

#get the final pseudo R-squared
1 - (resfitted/restotal)




