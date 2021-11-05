library(bpnreg)


bp<-bpnme(Error.rad ~ Maze + Trial.type + (1|Subject), Maps, its = 1000, burn=1)

x11()
traceplot(bp)
BFc(bp)


library(bpnreg)
fit.Maps <- bpnme(pred.I = Error.rad ~ Maze + Trial.type + L.c + (1|Subject),
                  data = Maps,
                  its = 100, burn = 10, n.lag = 3)

print(fit.Maps)

traceplot(fit.Maps, parameter="beta2")
BFc(fit.Maps, hypothesis = "Maze1 < Trial.type1")
coef_lin(fit.Maps)
coef_circ(fit.Maps)

coef_ran(fit.Maps)
coef_ran(fit.Maps, type = "circular")
fit(fit.Maps)



rm(list=ls())

setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/")


#niche expansion direction based on 3 variables
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")

#COMPARE AGAINST THE 2VAR DATASET
#csvfile4<- paste("./IntermediateOutput/Rstate_compile/2Var_plant_D_shiftvalues_zcor_center.txt",sep="")


#compare
shiftdf<-read.table(paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep=""),header=T,row.names = NULL,sep="\t")


#cut out invalid options
shiftdf<-shiftdf[which(shiftdf$max_exp_dist1!=-Inf | is.na(shiftdf$max_exp_dist1)),]
shiftdf<-shiftdf[which(shiftdf$max_exp_dist2!=-Inf | is.na(shiftdf$max_exp_dist2)),]

#we occasionally get expansion because of the kernal where we don't have observations. We cut these out
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

shiftdf_1$species_f<-NA

sp_list<-unique(shiftdf_1$species.name)

for(i in 1:length(sp_list)){
  print(sp_list[i])
  shiftdf_1$species_f[which(shiftdf_1$species.name==sp_list[i])]<-i
}
shiftdf_1$species_f<-as.numeric(shiftdf_1$species_f)

fit.Maps <- bpnme(pred.I = exp_dir1  ~ 1 + (1|species_f),
                  data = shiftdf_1,
                  its = 1000, burn = 50, n.lag = 3)


print(fit.Maps)

traceplot(fit.Maps, parameter="beta1")
#BFc(fit.Maps, hypothesis = "Maze1 < Trial.type1")
coef_lin(fit.Maps)
coef_circ(fit.Maps)

coef_ran(fit.Maps)
coef_ran(fit.Maps, type = "circular")
coef_ran(fit.Maps, type = "linear")
fit(fit.Maps)






