##############################################
### circular_analysis_Regional ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 

#This is a modified version of "circular_analysis" which separates circular information by region, instead of running it all together

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


#The output from this script is used in Table 1, Table 2, Appendix Figure S1.7, Appendix Figure S1.8


### ###

rm(list=ls())



###################################
#set paths and variables
##################################

#set to current location of repo
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

# direction of unfilling (where species fail to fail their niche)
# unfillingRegional/3Var0304plant_D_shiftvalues_zcor_center.txt

# direction of stability (where species fail to fail their niche)
# stabilityRegional/3Var0304plant_D_shiftvalues_zcor_center.txt



#select which output file you would like to analyse (published version is niche expansion direction based on 3 variables)
csvfile4 <- paste0("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt")

#choose the name of the output file
output<-"3VarEXP"



shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")
#View(shiftdf)



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


#cut out invalid options, remove any species with no valid expansion
shiftdf<-shiftdf[which(shiftdf$max_exp_dist1!=-Inf | is.na(shiftdf$max_exp_dist1)),]
shiftdf<-shiftdf[which(shiftdf$max_exp_dist2!=-Inf | is.na(shiftdf$max_exp_dist2)),]

#we occasionally get expansion where we don't have observations because of the way kernel density is calculated. We cut these out
shiftdf<-shiftdf[!is.na(shiftdf$exp_dir1)|!is.na(shiftdf$exp_dir2),]

#only keep these species/populations, we will now analyse species that expand only
shiftdf<-shiftdf[shiftdf$expansion>0.1,]



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


#make a second table to describe the 4 quarters of climate circle
quar_tab<-arr_tab[c(1,3),]

quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)

rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")



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

#take the direction of expansion column and make it circular format
dir<-as.circular(shiftdf_1[,2], units="radians",rotation='counter')

#set which column has distance information in it
dist<-shiftdf_1[,1]





##############################################
#circular anova: is expansion direction dependent on region?
##############################################


#make a list of all regions to process
reg_list<-sort(unique(shiftdf_1$region))



#prepare output table for parametric model, to be filled per region
par_table<-data.frame(region=character(),circ_model=numeric(),
                      circ_params=numeric(),q1=numeric(),k1=numeric(),
                      q2=numeric(),k2=numeric(), R=numeric(), 
                      q1_nearang=numeric(),q1_climdir=numeric(), q2_nearang=numeric(),q2_climdir=numeric(),
                      AICweight=numeric(),stringsAsFactors = F)
#prepare output table for non-parametric model, to be filled per region
nonpar_table<-data.frame(region=character(),smoothBW=numeric(),smoothLL=numeric(),
                         Rsquared=numeric(), medDistanc=numeric(), maxDistanc=numeric(), MaxDistDir = numeric(),
                         WetAng=numeric(),DryAng=numeric(),WarmAng=numeric(), ColdAng=numeric(),
                         stringsAsFactors = F)

#loop through each region
for (i in 1: length(reg_list)){

  reg<-reg_list[i]
  print(as.character(reg))
  
  #subset to only this region
  shiftdf_reg<-shiftdf_1[shiftdf_1$region==reg,]

  
  dir<-shiftdf_reg$exp_dir1
  
  dir<-as.circular(dir, units="radians",rotation='counter')
  
  dist<-shiftdf_reg$median_exp_dist1
  
  
  #############part one: run parametric model###############
  circ_res<-circ_mle(dir)
  
  #make a circular histogram of direction of expansion
  plotname=paste("FinalOutput/supplementary/circular_analysis_Regional/",output,"/DirectionRose_",reg,".pdf",sep="")
  pdf(file=plotname, width=5, height=5)

  plot_circMLE.custom(dir, circ_res,shrink=1.4)
  s.corcircle(arr_tab[,c(1,2)],grid=F,add=T)
  title(paste(reg,"n= ",length(dir)))
  
  dev.off()
  
  
  #how far is mean angles of expansion from main climate directions?
  #if the result is close to 0 this means the average direction of expansion is in that direction
  #if the result is close to 2pi then it is the OPPOSITE direction
  dir1<-circ_res$results[1,2]
  dir2<-circ_res$results[1,5]

  
  wetang1<-find_CLIMangle(dir1,NA,"Wetter")[1]
  dryang1<-find_CLIMangle(dir1,NA,"Drier")[1]
  warmang1<-find_CLIMangle(dir1,NA,"Warmer")[1]
  coldang1<-find_CLIMangle(dir1,NA,"Colder")[1]
  dir1_summ<-cbind(wetang1,dryang1,warmang1,coldang1)
  colnames(dir1_summ)<-c("wetter","drier","warmer","colder")
  
  dir1_clim<-colnames(dir1_summ)[which.min(dir1_summ)]
  dir1_min<-min(dir1_summ)
  
  #if model has secondary direction then include this too
  if(!is.na(dir2)){
    wetang2<-find_CLIMangle(dir2,NA,"Wetter")[1]
    dryang2<-find_CLIMangle(dir2,NA,"Drier")[1]
    warmang2<-find_CLIMangle(dir2,NA,"Warmer")[1]
    coldang2<-find_CLIMangle(dir2,NA,"Colder")[1]
    dir2_summ<-cbind(wetang2,dryang2,warmang2,coldang2)
    colnames(dir2_summ)<-c("wetter","drier","warmer","colder")
    
    dir2_clim<-colnames(dir2_summ)[which.min(dir2_summ)]
    dir2_min<-min(dir2_summ)
    
    
  }else{dir2_clim<-NA;dir2_min<-NA}

  
  #compile resuts of the model and put into a data row
  circ_model<-circ_res$bestmodel
  circ_rstat<-circ_res$rt[1]
  circ_pvalue<-circ_res$rt[2]
  
  circ_params<-circ_res$results[1,1]
  circ_AICweight<-round(circ_res$results[1,19],digits=2)
  
  circ_q1<-round(circ_res$results[1,2],digits=2)
  circ_k1<-round(circ_res$results[1,3],digits=2)
  circ_R<-round(circ_res$results[1,4],digits=2)
  circ_q2<-round(circ_res$results[1,5],digits=2)
  circ_k2<-round(circ_res$results[1,6],digits=2)
  
  
  #compile results and add to result table
  res1<-c(as.character(reg),circ_model,circ_params,
          circ_q1,circ_k1,circ_q2,circ_k2,circ_R,
          dir1_min, dir1_clim,dir2_min, dir2_clim,circ_AICweight)
  par_table[nrow(par_table)+1,]<-res1
  

  
  
  ####################part 2: non-parametric regression#################
  
  #run a non-parametric circular regression, with kernel smoothing
  #we use 2 alternative smoothing methods, Nadaraya-Watson and Local-Linear
  #these 2 methods should give similar results
  estNW <- kern.reg.circ.lin(dir, dist, method="NW")
  estLL <- kern.reg.circ.lin(dir, dist, method="LL")
  
  
  
  
  ###########################################
  ##### get the residuals from non-parametric regression and find the circular R-squared
  ###########################################

  
  row.names(shiftdf_reg)<-1:nrow(shiftdf_reg)
  
  #Step 1: find maximum residuals when gradient is set to 0 at mean (intercept)
  restotal<-euc.dist.center(shiftdf_reg[,c(2,1)], data.frame(shiftdf_reg[,2],mean(shiftdf_reg$median_exp_dist1)))
  
  #step 2: find our residuals right now
  #get fitted values
  fitted<-data.frame(x=as.numeric(estNW$x),y=as.numeric(estNW$y))
  

  
  #do some rounding, we need to find the nearest fitted points to each true data point
  shift_val<-shiftdf_reg[,c(2,1)]
  shift_val$xround2<-round(shift_val$exp_dir1, digits=2)
  shift_val$xround1<-round(shift_val$exp_dir1, digits=1)
  
  fitted$x_round2<-round(fitted$x,digits=2)
  fitted$x_round1<-round(fitted$x,digits=1)
  
  sum(shift_val$xround %in% fitted$x_round)
  newdf<-data.frame(x2=shift_val$xround2,x1=shift_val$xround1,y=shift_val$median_exp_dist1,
                    x_pred2=shift_val$xround2,x_pred1=shift_val$xround1)
  newdf$y_pred<-NA

  
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
  rsq<-1 - (resfitted/restotal)
  
  
  #make a plot of our non-parametric regression
  plotname=paste("FinalOutput/supplementary/circular_analysis_Regional/",output,"/DirectionMagnitude_",reg,".pdf",sep="")
  pdf(file=plotname, width=4, height=4)
  

  res<-plot(estNW, plot.type="circle", points.plot=TRUE,main=paste(reg,"n:",length(dist)),
            #labels=c("N","NE","E","SE","S","SO","O","NO"),
            label.pos=seq(0,7*pi/4,by=pi/4), zero=0, clockwise=FALSE,cex=0.6)
  lines(estLL, plot.type="circle", plot.info=res, line.col=1)
  
  dev.off()


  
  #we want to get some information on which direction has maximum expansion, and what direction that's in
  #get median distance
  medDist<-summary(dist)["Median"]
  #get data on the output of the regression. x is angle, y is distance
  estNW2<-as.data.frame(estNW[c("x","y")])
  #find where the maximum distance is
  estNWmax<-estNW2[which.max(estNW2$y),]
  #where are the 90th percentile points?
  estNW2<-estNW2[which(estNW2$y>quantile(estNW2$y, 0.9)),]
  
  #what direction is this maximum distance of expansion in?
  wetang<-find_CLIMangle(estNWmax$x,NA,"Wetter")
  dryang<-find_CLIMangle(estNWmax$x,NA,"Drier")
  warmang<-find_CLIMangle(estNWmax$x,NA,"Warmer")
  coldang<-find_CLIMangle(estNWmax$x,NA,"Colder")
  
  
  #plot some additional results based on this maximum distance
  #plot where on the circle the area of maxmimum distance was, and where the angle of precipitation is
  plotname=paste("FinalOutput/supplementary/circular_analysis_Regional/",output,"/DirectionMagnitude_Suppmaxplot_",reg,".pdf",sep="")
  pdf(file=plotname, width=4, height=4)
  
  
  res<-plot(estNW, plot.type="circle", points.plot=TRUE,
            #labels=c("N","NE","E","SE","S","SO","O","NO"),
            label.pos=seq(0,7*pi/4,by=pi/4), zero=0, clockwise=FALSE,cex=0.6)
  points(estNW2$x,estNW2$y,col="red",cex=0.75)
  points(estNWmax$x,estNWmax$y,col="blue",cex=1.2)
  lines(estNW, plot.type="circle", plot.info=res, line.col=1)
  arrows.circular(arr_tab["Precip","angle"], length=0.1)
  arrows.circular(arr_tab["Tmin","angle"], length=0.1)
  arrows.circular(arr_tab["TMax","angle"], length=0.1)
  #s.corcircle(arr_tab[,c(1,2)]*4,grid=F,add=T)
  
  
  dev.off()
  
  
  
  
  res2<-c(as.character(reg),round(estNW$bw,digits=2),round(estNW$bw,digits=2),round(rsq,digits=2),
          round(medDist,2), round(estNWmax$y,2),round(estNWmax$x,2), wetang[1], dryang[1], warmang[1], coldang[1])
  
  
  nonpar_table[nrow(nonpar_table)+1,]<-res2
  
  
}


#now save output tables for later
write_path<-paste("FinalOutput/supplementary/circular_analysis_Regional/",output,"/ParStat_regional.csv",sep="")
write.csv(par_table,file=write_path,row.names = F)

write_path<-paste("FinalOutput/supplementary/circular_analysis_Regional/",output,"/NonParStat_regional.csv",sep="")
write.csv(nonpar_table,file=write_path,row.names = F)






