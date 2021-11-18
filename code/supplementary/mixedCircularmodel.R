##############################################
### mixedCircularmodel ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 

# the file circular_analysis.R analyses direction and magnitude trends in expansion across all populations
# But unfortunately it has pseudo-replication, it has multiple species repeating that colonise new regions 
#(e.g. A. theophrasti is naturalised in N. America and Europe, these are two rows in the data but are not independent)

#We account for this in circular_analysis.R by randomly drawing single species accounts and running again
#this isn't ideal either because it cuts down the data
#ideally we should run a mixed circular model to account for random effects
#this random effects circular model can't do many of the clever things circular_analysis.R can do BUT it can check how much REs influence our patterns of expansion

#the output of this file are used in the results section of the main paper as a subsidiary analysis


### ###

rm(list=ls())


###################################
#set paths and variables
##################################

#set to current repo
setwd("DIRECTORY_HERE")

library(bpnreg)

#niche expansion direction based on 3 variables. load about from Rstate_compile (from plant_PCA_expand)
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")

#load expansion based on 2 variables if you prefer
#csvfile4<- paste("./IntermediateOutput/Rstate_compile/2Var_plant_D_shiftvalues_zcor_center.txt",sep="")


#read in file
shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")



###################################
#BASIC EXPLORATION AND CHECKS
##################################

#cut out invalid options, remove any species with no valid expansion
shiftdf<-shiftdf[which(shiftdf$max_exp_dist1!=-Inf | is.na(shiftdf$max_exp_dist1)),]
shiftdf<-shiftdf[which(shiftdf$max_exp_dist2!=-Inf | is.na(shiftdf$max_exp_dist2)),]

#we occasionally get expansion because of the kernel where we don't have observations. We cut these out
shiftdf<-shiftdf[!is.na(shiftdf$exp_dir1)|!is.na(shiftdf$exp_dir2),]


sum(shiftdf$expansion>0.1)

shiftdf<-shiftdf[shiftdf$expansion>0.1,]
#get a df with direction, region and species

shiftdf_1<-shiftdf[,c("exp_dir1",
                      "region",
                      "native_center_tmax",
                      "native_center_tmin",
                      "native_center_precip",
                      "species.name","expansion")]
shiftdf_2<-shiftdf[,c("exp_dir2",
                      "region",
                      "native_center_tmax",
                      "native_center_tmin",
                      "native_center_precip",
                      "species.name","expansion")]

names(shiftdf_2)<-names(shiftdf_1)
shiftdf_1<-rbind(shiftdf_1,shiftdf_2)

shiftdf_1<-shiftdf_1[complete.cases(shiftdf_1),]
shiftdf_1[shiftdf_1[,1]>2*pi,1] <- shiftdf_1[shiftdf_1[,1]>2*pi,1] - 2*pi

length(unique(shiftdf_1$species.name))


levels(shiftdf_1$region)<-1:length(unique(shiftdf_1$region))


#set region to factor and assign level number
shiftdf_1$region_f<-NA

shiftdf_1$region_f[which(shiftdf_1$region=="Afrotropical")]<-1
shiftdf_1$region_f[which(shiftdf_1$region=="Australian")]<-2
shiftdf_1$region_f[which(shiftdf_1$region=="Madagascan")]<-3
shiftdf_1$region_f[which(shiftdf_1$region=="Nearctic")]<-4
shiftdf_1$region_f[which(shiftdf_1$region=="Neotropical")]<-5
shiftdf_1$region_f[which(shiftdf_1$region=="Oceanian")]<-6
shiftdf_1$region_f[which(shiftdf_1$region=="Oriental")]<-7
shiftdf_1$region_f[which(shiftdf_1$region=="Palearctic_East")]<-8
shiftdf_1$region_f[which(shiftdf_1$region=="Palearctic_West")]<-9
shiftdf_1$region_f[which(shiftdf_1$region=="Panamanian")]<-10
shiftdf_1$region_f[which(shiftdf_1$region=="Saharo-Arabian")]<-11
shiftdf_1$region_f[which(shiftdf_1$region=="Sino-Japanese")]<-12
shiftdf_1$region_f<-as.factor(shiftdf_1$region_f)


#set species to factor and assign level number
shiftdf_1$species_f<-NA

sp_list<-unique(shiftdf_1$species.name)

for(i in 1:length(sp_list)){
  print(sp_list[i])
  shiftdf_1$species_f[which(shiftdf_1$species.name==sp_list[i])]<-i
}
shiftdf_1$species_f<-as.numeric(shiftdf_1$species_f)


###############################################
###Fit a Bayesian circular mixed-effects model
################################################

fit.Maps <- bpnme(pred.I = exp_dir1  ~ 1 + (1|species_f),
                  data = shiftdf_1,
                  its = 1000, burn = 50, n.lag = 3)

#check results
print(fit.Maps)

#check traceplot, make sure it's a "fuzzy caterpillar"
traceplot(fit.Maps, parameter="beta1")

#check linear coefficients
coef_lin(fit.Maps)
#check circular coefficients
coef_circ(fit.Maps)

#check random effects
coef_ran(fit.Maps)
#check circular random effects.  posterior summaries of the circular or linear random effect variances, overall influence of random effects
coef_ran(fit.Maps, type = "circular")
#check linear random effects
coef_ran(fit.Maps, type = "linear")
#check model statistics
fit(fit.Maps)






