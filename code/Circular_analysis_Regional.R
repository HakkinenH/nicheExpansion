
rm(list=ls())

setwd("/Users/henry_hakkinen/Documents/Research/ExeterWork/Data/niche_shift_work")

setwd("U:/Data/niche_shift_work/")


#niche expansion direction based on 4 variables (including NPP)
#csvfile4 <- paste("NPPplant_D_shiftvalues_zuncor_center.txt",sep="")

#niche expansion direction based on 3 variables
csvfile4 <- paste("regional/3Var0304plant_D_shiftvalues_zcor_center.txt",sep="")

#direction of unfilling (where species fail to fail their niche)
#csvfile4 <- paste("unfillingRegional/3Var0304plant_D_shiftvalues_zcor_center.txt",sep="")

#direction of stability (where species fail to fail their niche)
#csvfile4 <- paste("stabilityRegional/3Var0304plant_D_shiftvalues_zcor_center.txt",sep="")


output<-"3VarEXP"

shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")
#View(shiftdf)



shiftdf<-shiftdf[which(shiftdf$max_exp_dist1!=-Inf | is.na(shiftdf$max_exp_dist1)),]
shiftdf<-shiftdf[which(shiftdf$max_exp_dist2!=-Inf | is.na(shiftdf$max_exp_dist2)),]

shiftdf<-shiftdf[shiftdf$expansion>0.1,]

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

source("circ_plot_functions.R")
#load("U:/Data/niche_shift_work/PCA_contrib")
load("/Users/henry_hakkinen/Documents/Research/ExeterWork/Data/niche_shift_work/PCA_contrib")


ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
arr_tab<-pca.cal$co[, 1:2]/max(abs(pca.cal$co[, 1:2]))

colnames(arr_tab)<-c("x","y")
#table with points which give the direction as compared to the origin (0,0)
#find the angle of direction for each PCA component
arr_tab$angle<-apply(arr_tab,1,rad.ang,x2=c(0,0))
rownames(arr_tab)<-c("Tmin","TMax","Precip")

#make a second table to describe the 4 quarters of climate circle
quar_tab<-arr_tab[c(1,3),]

quar_tab[3,]<-c(NA,NA,quar_tab[1,"angle"]-pi)
quar_tab[4,]<-c(NA,NA,quar_tab[2,"angle"]+pi)

rownames(quar_tab)<-c("Warmer","Wetter","Colder","Drier")

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

dir<-as.circular(shiftdf_1[,2], units="radians",rotation='counter')

dist<-shiftdf_1[,1]





##############################################
#circular anova: is expansion direction dependent on region?
##############################################



reg_list<-sort(unique(shiftdf_1$region))

reg_list




par_table<-data.frame(region=character(),circ_model=numeric(),
                      circ_params=numeric(),q1=numeric(),k1=numeric(),
                      q2=numeric(),k2=numeric(), R=numeric(), 
                      q1_nearang=numeric(),q1_climdir=numeric(), q2_nearang=numeric(),q2_climdir=numeric(),
                      AICweight=numeric(),stringsAsFactors = F)
nonpar_table<-data.frame(region=character(),smoothBW=numeric(),smoothLL=numeric(),
                         Rsquared=numeric(), medDistanc=numeric(), maxDistanc=numeric(), MaxDistDir = numeric(),
                         WetAng=numeric(),DryAng=numeric(),WarmAng=numeric(), ColdAng=numeric(),
                         stringsAsFactors = F)


i<-1
for (i in 1: length(reg_list)){
  
  
  reg<-reg_list[i]
  print(as.character(reg))
  
  shiftdf_reg<-shiftdf_1[shiftdf_1$region==reg,]

  
  dir<-shiftdf_reg$exp_dir1
  
  dir<-as.circular(dir, units="radians",rotation='counter')
  
  dist<-shiftdf_reg$median_exp_dist1
  
  
  #############part one: run parametric regression###############
  circ_res<-circ_mle(dir)
  
  plotname=paste("Figures/",output,"/Regional/DirectionRose_",reg,".pdf",sep="")
  pdf(file=plotname, width=5, height=5)
  


  plot_circMLE.custom(dir, circ_res,shrink=1.4)
  s.corcircle(arr_tab[,c(1,2)],grid=F,add=T)
  title(paste(reg,"n= ",length(dir)))
  
  dev.off()
  
  #how far is angle from precipitation angle?
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
  
  
  res1<-c(as.character(reg),circ_model,circ_params,
          circ_q1,circ_k1,circ_q2,circ_k2,circ_R,
          dir1_min, dir1_clim,dir2_min, dir2_clim,circ_AICweight)
  par_table[nrow(par_table)+1,]<-res1
  

  
  
  ####################part 2: non-parametric regression#################
  estNW <- kern.reg.circ.lin(dir, dist, method="NW")
  estLL <- kern.reg.circ.lin(dir, dist, method="LL")
  
  
  
  
  ##all of the below is to calculate a pseudo R-squared
  
  
  euc.dist.center <- function(x1, x2) (sum((x1 - x2) ^ 2))
  
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
  
  
  
  plotname=paste("Figures/",output,"/Regional/DirectionMagnitude_",reg,".pdf",sep="")
  pdf(file=plotname, width=4, height=4)
  

  res<-plot(estNW, plot.type="circle", points.plot=TRUE,main=paste(reg,"n:",length(dist)),
            #labels=c("N","NE","E","SE","S","SO","O","NO"),
            label.pos=seq(0,7*pi/4,by=pi/4), zero=0, clockwise=FALSE,cex=0.6)
  lines(estLL, plot.type="circle", plot.info=res, line.col=1)
  
  dev.off()


  
  
  medDist<-summary(dist)["Median"]
  estNW2<-as.data.frame(estNW[c("x","y")])
  estNWmax<-estNW2[which.max(estNW2$y),]
  estNW2<-estNW2[which(estNW2$y>quantile(estNW2$y, 0.9)),]
  
  
  wetang<-find_CLIMangle(estNWmax$x,NA,"Wetter")
  dryang<-find_CLIMangle(estNWmax$x,NA,"Drier")
  warmang<-find_CLIMangle(estNWmax$x,NA,"Warmer")
  coldang<-find_CLIMangle(estNWmax$x,NA,"Colder")
  
  
  plotname=paste("Figures/",output,"/Regional/Supp/",as.character(reg),"DirectionMagnitude_Suppmaxplot.pdf",sep="")
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


write_path<-paste("Figures/",output,"/ParStat_regional.csv",sep="")
write.csv(par_table,file=write_path,row.names = F)

write_path<-paste("Figures/",output,"/NonParStat_regional.csv",sep="")
write.csv(nonpar_table,file=write_path,row.names = F)


read.csv(write_path)

###################################################
##########do the same but with native region. Where do species go?
#####################################################

native_info<-read.csv("primary_native_info.csv")
head(shiftdf_1)
head(native_info)

shiftdf_11<-merge(shiftdf_1, native_info[,c(1,7)], by.x="species.name",by.y="sp_name",all.x=T)


reg_list<-sort(unique(shiftdf_11$primary_native))


par_table<-data.frame(region=character(),circ_model=numeric(),
                      circ_params=numeric(),q1=numeric(),k1=numeric(),
                      q2=numeric(),k2=numeric(), precip_rad=numeric(),precip_deg=numeric(),
                      AICweight=numeric(),stringsAsFactors = F)
nonpar_table<-data.frame(region=character(),smoothBW=numeric(),smoothLL=numeric(),
                         Rsquared=numeric(), stringsAsFactors = F)

for (i in 1: length(reg_list)){
  
  
  reg<-reg_list[i]
  print(as.character(reg))
  
  shiftdf_reg<-shiftdf_11[shiftdf_11$primary_native==reg,]

  
  dir<-shiftdf_reg$exp_dir1
  
  dir<-as.circular(dir, units="radians",rotation='counter')
  
  dist<-shiftdf_reg$median_exp_dist1
  
  
  #############part one: run parametric regression###############
  circ_res<-circ_mle(dir)
  
  plotname=paste("Figures/",output,"/Regional/NATIVEDirectionRose_",reg,".pdf",sep="")
  pdf(file=plotname, width=8, height=8)
  
  
  plot_circMLE.custom(dir, circ_res,shrink=1.2)
  s.corcircle(arr_tab[,c(1,2)],grid=F,add=T)
  title(paste(reg,"n= ",length(dir)))
  
  dev.off()
  
  #how far is angle from precipitation angle?
  dir1<-circ_res$results[1,2]
  
  if(dir1>2*pi){dir1<-dir1-2*pi}
  if(dir1>2*pi){dir1<-dir1-2*pi}
  
  dir2<-circ_res$results[1,5]
  
  if(!is.na(dir2)){
    if(dir2>2*pi){dir2<-dir2-2*pi}
    if(dir2>2*pi){dir2<-dir2-2*pi}
  }
  
  diff1<-abs(arr_tab[3,3]-dir1)
  diff2<-abs(arr_tab[3,3]-dir2)
  
  
  diff1_1<-abs(dir1-arr_tab[3,3])
  diff2_2<-abs(dir2-arr_tab[3,3])
  
  
  precip_rad<-round(min(c(diff1,diff2,diff1_1,diff2_2),na.rm=T),digits=2)
  precip_deg<-round((precip_rad * (180/pi)),digits=2)
  
  
  
  circ_model<-circ_res$bestmodel
  circ_rstat<-circ_res$rt[1]
  circ_pvalue<-circ_res$rt[2]
  
  circ_params<-circ_res$results[1,1]
  circ_AICweight<-round(circ_res$results[1,19],digits=2)
  
  circ_q1<-round(circ_res$results[1,2],digits=2)
  circ_k1<-round(circ_res$results[1,3],digits=2)
  circ_q2<-round(circ_res$results[1,5],digits=2)
  circ_k2<-round(circ_res$results[1,6],digits=2)
  
  
  res1<-c(as.character(reg),circ_model,circ_params,
          circ_q1,circ_k1,circ_q2,circ_k2,
          precip_rad,precip_deg,circ_AICweight)
  par_table[nrow(par_table)+1,]<-res1
  
  
  ####################part 2: non-parametric regression#################
  estNW <- kern.reg.circ.lin(dir, dist, method="NW")
  estLL <- kern.reg.circ.lin(dir, dist, method="LL")
  
  
  
  
  ##all of the below is to calculate a pseudo R-squared
  
  
  euc.dist.center <- function(x1, x2) (sum((x1 - x2) ^ 2))
  
  row.names(shiftdf_reg)<-1:nrow(shiftdf_reg)
  
  
  #Step 1: find maximum residuals when gradient is set to 0 at mean (intercept)
  restotal<-euc.dist.center(shiftdf_reg[,c("exp_dir1","median_exp_dist1")], data.frame(shiftdf_reg$exp_dir1,mean(shiftdf_reg$median_exp_dist1)))

  
  
  #step 2: find our residuals right now
  #get fitted values
  fitted<-data.frame(x=as.numeric(estNW$x),y=as.numeric(estNW$y))
  
  
  
  #do some rounding, we need to find the nearest fitted points to each true data point
  

  shift_val<-shiftdf_reg[,c("exp_dir1","median_exp_dist1")]

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
  
  plotname=paste("Figures/",output,"/Regional/NATIVEDirectionMagnitude_",reg,".pdf",sep="")
  pdf(file=plotname, width=8, height=8)
  
  
  res<-plot(estNW, plot.type="circle", points.plot=TRUE,main=paste(reg,"n:",length(dist)),
            #labels=c("N","NE","E","SE","S","SO","O","NO"),
            label.pos=seq(0,7*pi/4,by=pi/4), zero=0, clockwise=FALSE,cex=0.6)
  lines(estLL, plot.type="circle", plot.info=res, line.col=2)
  
  
  plot(newdf[,1],newdf[,3])
  points(newdf[,4],newdf[,6],col="red")
  
  dev.off()
  
  
  
  rsq<-1 - (resfitted/restotal)
  
  
  
  
  res2<-c(as.character(reg),round(estNW$bw,digits=2),round(estLL$bw,digits=2),round(rsq,digits=2))
  nonpar_table[nrow(nonpar_table)+1,]<-res2
  
  
}


write_path<-paste("Figures/",output,"/NATIVEParStat_regional.csv",sep="")
write.csv(par_table,file=write_path,row.names = F)

write_path<-paste("Figures/",output,"/NATIVENonParStat_regional.csv",sep="")
write.csv(nonpar_table,file=write_path,row.names = F)



