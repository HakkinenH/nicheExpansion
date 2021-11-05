
#correlate how much available climate there is in the native niche, and how much is actually occupied



setwd("C:/Users/Henry/Documents/Research/RepoCode/nicheExpansion/")

yooo<-list.files(path="./IntermediateOutput/plant_PCA_nichefilling/findFilling",pattern="*.csv")

occRaw<-read.csv(paste0("./IntermediateOutput/plant_PCA_nichefilling/findFilling/",yooo[1]))

for (i in 2:length(yooo)){
  print(i)
  occRaw_nw<-read.csv(paste0("./IntermediateOutput/plant_PCA_nichefilling/findFilling/", yooo[i]))
  occRaw<-rbind(occRaw, occRaw_nw)
}


write.csv(occRaw,paste("IntermediateOutput/plant_PCA_nichefilling/3Varplant_D_FILLvalues_z.uncor_center.csv",sep=""))


cold<-occRaw[,c("sp_name", "region", "coldOccdens", "coldPotdens", "coldPro")]
names(cold)<-c("sp_name", "region", "Occdens", "Potdens", "Pro")
cold$quar<-"cold"

hot<-occRaw[,c("sp_name", "region", "tmaxOccdens", "tmaxPotdens", "tmaxPro")]
names(hot)<-c("sp_name", "region", "Occdens", "Potdens", "Pro")
hot$quar<-"hot"

dry<-occRaw[,c("sp_name", "region", "dryOccdens", "dryPotdens", "dryPro")]
names(dry)<-c("sp_name", "region", "Occdens", "Potdens", "Pro")
dry$quar<-"dry"

wet<-occRaw[,c("sp_name", "region", "precipOccdens", "precipPotdens", "precipPro")]
names(wet)<-c("sp_name", "region", "Occdens", "Potdens", "Pro")
wet$quar<-"wet"



occFull<-rbind(cold, hot, dry, wet)
head(occFull)

dim(occFull)
sum(occFull$Occdens==0)
sum(occFull$Potdens==0)
sum(occFull$Pro==0, na.rm=T)


occFull$logOccdens<-log10(occFull$Occdens+min(occFull$Occdens[occFull$Occdens>0], na.rm=T)/10)
occFull$logPotdens<-log10(occFull$Potdens+min(occFull$Potdens[occFull$Potdens>0], na.rm=T)/10)
occFull$logPro<-log10(occFull$Pro+min(occFull$Pro[occFull$Pro>0], na.rm=T)/10)

hist(occFull$logOccdens)
hist(occFull$logPotdens)
hist(occFull$Pro)


boxplot(logOccdens  ~ quar, data=occFull)
boxplot(logPotdens  ~ quar, data=occFull, xlab = "Type of climate available", ylab="Total climate available (summed density logged)")
boxplot(logPro  ~ quar, data=occFull)



#basic comparison, is filling predicted by the quarter
occFull1<-occFull[occFull$Pro>0,]

anovRes<-glm(logPro ~ quar, data=occFull1)
plot(anovRes)

boxplot(logPro  ~ quar, data=occFull1)
summary(anovRes)

#mixed linear regression, is the slope of niche filling predicted by the direction?
library(lme4)
library(scales)

occFull2<-occFull[occFull$Occdens>0,]
plot(occFull2$logPotdens, occFull2$logOccdens, cex=0.4)


model2<-lmer(logOccdens ~ logPotdens*quar + (1|region) + (1|sp_name), data=occFull2, REML=F)

model3<-lmer(logOccdens ~ logPotdens + (1|region) + (1|sp_name) , data=occFull2, REML=F)
anova(model3, model2)

plot(model2)

summary(model2)
fixef(model2)
summary(ranef(model2)$sp_name)

boxplot(logOccdens~region, data=occFull2)
hotDF<-data.frame(logPotdens=seq(min(occFull2$logPotdens), max(occFull2$logPotdens),by=0.01) , 
                   quar="hot", 
                   region="Australian", 
                   sp_name="Abutilon theophrasti")

hotDF$predVals<-predict(model2, hotDF)

coldDF<-data.frame(logPotdens=seq(min(occFull2$logPotdens), max(occFull2$logPotdens),by=0.01) , 
                  quar="cold", 
                  region="Australian", 
                  sp_name="Abutilon theophrasti")

coldDF$predVals<-predict(model2, coldDF)

wetDF<-data.frame(logPotdens=seq(min(occFull2$logPotdens), max(occFull2$logPotdens),by=0.01) , 
                  quar="wet", 
                  region="Australian", 
                  sp_name="Abutilon theophrasti")

wetDF$predVals<-predict(model2, wetDF)


dryDF<-data.frame(logPotdens=seq(min(occFull2$logPotdens), max(occFull2$logPotdens),by=0.01) , 
                  quar="dry", 
                  region="Australian", 
                  sp_name="Abutilon theophrasti")

dryDF$predVals<-predict(model2, dryDF, se.fit = T)


plot(occFull2$logPotdens[occFull2$quar=="hot"], occFull2$logOccdens[occFull2$quar=="hot"], cex=0.5, col="red",
     xlab = "Total climate availabile (summed density logged)", ylab="Total climate occuped (summed density logged)")
points(occFull2$logPotdens[occFull2$quar=="cold"], occFull2$logOccdens[occFull2$quar=="cold"], cex=0.5, col="blue")
points(occFull2$logPotdens[occFull2$quar=="wet"], occFull2$logOccdens[occFull2$quar=="wet"], cex=0.5, col="green")
points(occFull2$logPotdens[occFull2$quar=="dry"], occFull2$logOccdens[occFull2$quar=="dry"], cex=0.5, col="orange")

plot(occFull2$logPotdens, occFull2$logOccdens, cex=0.5, col=alpha("black", 0.3),
     xlab = "Total climate availabile (summed density logged)", ylab="Total climate occuped (summed density logged)")

points(hotDF$logPotdens, hotDF$predVals, cex=0.2, col="red", pch=4)
points(coldDF$logPotdens, coldDF$predVals, cex=0.2, col="blue", pch=4)
points(wetDF$logPotdens, wetDF$predVals, cex=0.2, col="green", pch=4)
points(dryDF$logPotdens, dryDF$predVals, cex=0.2, col="orange", pch=4)






#load
csvfile4 <- paste("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt",sep="")
shiftdf<-read.table(csvfile4,header=T,row.names = NULL,sep="\t")

names(shiftdf)


