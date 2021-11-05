rm(list=ls())

setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/")



#niche expansion direction based on 4 variables (including NPP)
#csvfile4 <- paste("NPPplant_D_shiftvalues_zuncor_center.txt",sep="")

#niche expansion direction based on 3 variables
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")
#csvfile4 <- paste("regional/3Var0306plant_D_shiftvalues_zcor_center.txt",sep="")

#csvfile4 <- paste("regionalCUT/3Var0306plant_D_shiftvalues_zuncor_center.txt",sep="")

#direction of unfilling (where species fail to fail their niche)
#csvfile4 <- paste("unfillingRegional/3Var0304plant_D_shiftvalues_zcor_center.txt",sep="")

#direction of stability (where species fail to fail their niche)
#csvfile4 <- paste("stabilityRegional/3Var0304plant_D_shiftvalues_zcor_center.txt",sep="")


#COMPARE AGAINST THE 2VAR DATASET
#csvfile4<- paste("./IntermediateOutput/Rstate_compile/2Var_plant_D_shiftvalues_zcor_center.txt",sep="")

#COMPARE AGAINST THE 4VAR DATASET
#csvfile4<- paste("./IntermediateOutput/Rstate_compile/4Var_plant_D_shiftvalues_zcor_center.txt",sep="")


output<-"3VarEXP"
PA_test<-"FALSE"


shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")

head(shiftdf)






names(shiftdf)
#x11()
par(mfrow=c(1,3))
hist(shiftdf$proportion_precip)
hist(shiftdf$proportion_tmin)
hist(shiftdf$proportion_tmax)



#cut out invalid options
shiftdf<-shiftdf[which(shiftdf$max_exp_dist1!=-Inf | is.na(shiftdf$max_exp_dist1)),]
shiftdf<-shiftdf[which(shiftdf$max_exp_dist2!=-Inf | is.na(shiftdf$max_exp_dist2)),]

#we occasionally get expansion because of the kernal where we don't have observations. We cut these out
shiftdf<-shiftdf[!is.na(shiftdf$exp_dir1)|!is.na(shiftdf$exp_dir2),]


sum(shiftdf$expansion>0.1)

shiftdf<-shiftdf[shiftdf$expansion>0.1,]
table(shiftdf$region)
table(shiftdf$circular_model)

rad.ang <- function(x1, x2){
  dx = x1[1] - x2[1]
  dy = x1[2] - x2[2] 
  
  if(!is.numeric(dx)){dx<-as.numeric(dx)}
  if(!is.numeric(dy)){dy<-as.numeric(dy)}
  
  theta = atan2(dy,dx)
  
  if(theta<0){ theta<- theta+ 2*pi }
  return(theta)
} 

find_CLIMangle<-function(dir1,dir2,clim){
  if(dir1>2*pi){dir1<-dir1-2*pi}
  if(dir1>2*pi){dir1<-dir1-2*pi}
  
  if(!is.na(dir2)){
    if(dir2>2*pi){dir2<-dir2-2*pi}
    if(dir2>2*pi){dir2<-dir2-2*pi}
  }
  
  diff1<-abs(quar_tab[clim,3]-dir1)
  diff2<-abs(quar_tab[clim,3]-dir2)
  
  
  diff1_1<-abs(dir1-quar_tab[clim,3])
  diff2_2<-abs(dir2-quar_tab[clim,3])
  
  
  res_rad<-round(min(c(diff1,diff2,diff1_1,diff2_2),na.rm=T),digits=2)
  res_deg<-round((res_rad * (180/pi)),digits=2)
  return(c(res_rad,res_deg))
  
} 



library("ecospat")

source("./code./functions/circ_plot_functions.R")

#for 2var
#load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TmaxPrecip")
#for 3 var
#load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecip")
#for 4 var
load("IntermediateOutput/plant_PCA_exp/PCA_contrib_TminTmaxPrecipPrSeas")


ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))

colnames(arr_tab)<-c("x","y")
#table with points which give the direction as compared to the origin (0,0)
#find the angle of direction for each PCA component
if (grepl("2Var", output)){
  arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))
  arr_tab
  rownames(arr_tab)<-c("TMax","Precip")
  #make a second table to describe the 4 quarters of climate circle
  quar_tab<-arr_tab[c(1,2),]
  
  quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
  quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)
  
  rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")
  
  
}else if (grepl("4Var", output)){
  
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



#####DATA PREP####

library("NPCirc")
library("CircMLE")



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
shiftdf_1<-rbind(shiftdf_1,shiftdf_2)

shiftdf_1<-shiftdf_1[complete.cases(shiftdf_1),]
shiftdf_1[shiftdf_1[,2]>2*pi,2] <- shiftdf_1[shiftdf_1[,2]>2*pi,2] - 2*pi


#if we want to select species at random to account for pseudo-replication
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

dir<-as.circular(shiftdf_1[,2], units="radians",rotation='counter')

dist<-shiftdf_1[,1]


#################

###############################################
###circular model, do species expand more often in certain directions?
################################################


plotname=paste("./IntermediateOutput/circular_analysis/",output,"/DirectionRose.pdf",sep="")
pdf(file=plotname, width=5, height=5)

circ_res<-circ_mle(dir)
plot_circMLE.custom(dir, circ_res,shrink=1.4)
axis.circular(at=circular(seq(0, 2*pi-pi/2, pi/2)), 
              labels=c("0",expression(pi/2),expression(pi),expression(3*pi/2)))
s.corcircle(arr_tab[,c(1,2)],grid=F,add=T)

dev.off()



#how far is angle from precipitation angle?

circ_res$results$lamda

circ_res

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



#second model is also viable
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


write.csv(circ_res$results,paste0("./IntermediateOutput/circular_analysis/",output,"/circ_summary.csv"))


#alternative approach with curve fitting
library(mclust)

mod4 <- densityMclust(shiftdf_1$exp_dir1)

summary(mod4)

plotname=paste0("./IntermediateOutput/circular_analysis/",output,"/CurveFit_Direction.pdf")
pdf(file=plotname, width=9, height=7)

plot(mod4, what = "density", data = shiftdf_1$exp_dir1, breaks = 20,xlab="Direction of niche shift")
abline(v=arr_tab[3,3])
abline(v=arr_tab[3,3]+pi)
dev.off()



###############################################
###circular model, is direction expansion predicted by (as above)
################################################

dir<-as.circular(shiftdf_1[,2], units="radians",rotation='counter')
dist<-shiftdf_1[,1]


estNW <- kern.reg.circ.lin(dir, dist, method="NW")
estLL <- kern.reg.circ.lin(dir, dist, method="LL")


plotname=paste("./IntermediateOutput/circular_analysis/",output,"/DirectionMagnitude_reg.pdf",sep="")
pdf(file=plotname, width=4, height=4)


res<-plot(estLL, plot.type="circle", points.plot=TRUE, show.radial.grid=F, grid.col="black",
          labels=c("0",expression(pi/2),expression(pi),expression(3*pi/2)),
          label.pos=seq(0,7*pi/2,by=pi/2), zero=0, clockwise=FALSE,cex=0.6)

lines(estLL, plot.type="circle", plot.info=res, line.col=2)
#s.corcircle(arr_tab[,c(1,2)]*4,grid=F,add=T)


dev.off()
res<-plot(estNW, plot.type="line", points.plot=TRUE)

#we want to get some information on which direction has maximum expansion, and what direction that's in

medDist<-summary(dist)["Median"]
estNW2<-as.data.frame(estNW[c("x","y")])
estNWmax<-estNW2[which.max(estNW2$y),]
estNW2<-estNW2[which(estNW2$y>quantile(estNW$y, 0.9)),]

find_CLIMangle(estNWmax$x,NA,"Precip")
find_CLIMangle(estNWmax$x,NA,"TMin")
find_CLIMangle(estNWmax$x,NA,"TMax")

find_CLIMangle(estNWmax$x,NA,"Wetter")
find_CLIMangle(estNWmax$x,NA,"Drier")
find_CLIMangle(estNWmax$x,NA,"Warmer")
find_CLIMangle(estNWmax$x,NA,"Colder")


plotname=paste("./IntermediateOutput/circular_analysis/",output,"/DirectionMagnitude_Suppmaxplot.pdf",sep="")
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


###########
#####if interested get the residuals from this and find the circular R-squared
###########


euc.dist.center <- function(x1, x2) (sum((x1 - x2) ^ 2))

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

1 - (resfitted/restotal)

#R-Squared 0.3947819




#adhoc
#yes! Species expand more towards higher precipitation

#compare
shift3var<-read.table(paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep=""),header=T,row.names = NULL,sep="\t")
dim(shift3var) #1883 for 3 var
dim(shiftdf) #1883 for 2 var
length(unique(shift3var$species.name)) #606 species
length(unique(shiftdf$species.name)) #606 species
sum(shift3var$expansion>0.1) #885 events expand
sum(shiftdf$expansion>0.1) # 703 events expand

shift3var<-shift3var[shift3var$expansion>0.1,]; var3list<-unique(shift3var$species.name)
shiftdf<-shiftdf[shiftdf$expansion>0.1,]; var2list<-unique(shiftdf$species.name)

sum(var3list %in% var2list) #351 species in 3var are in 2var (out of 404)
sum(var2list %in% var3list) #351 species in 2var are in 3var (out of 384)

spl<-unique(shift3var$species.name)

spDF<-data.frame(sp_name=spl[1], numNat=sum(shift3var$species.name==spl[1]))


for (i in 1:length(spl)){
  print(spl[i]);print(i)
  spDF[nrow(spDF)+1,]<-c(spl[i], sum(shift3var$species.name==spl[i]))
}

spDF$numNat<-as.numeric(spDF$numNat)
summary(spDF)


