

##############################################
### correlateFillingAvail ###
##############################################



### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 

#Does the amount of available climate affect our results across climatic quarters?

#we find differences in expansion proportion and niche filling across the 4 quarters of climate (dry, cold, hot, wet)
#but are these because there is MORE wet climate available (for example)?

#To investigate we made a supplementary analysis. This file runs two models
#this script takes output from plant_PCA_nichefilling.R

#ANOVA/kruskal-wallis: is there more available cliamte in one quarter? Is the proportion of niche filling higher in any given quarter?
#mixed linear regression: Is the relationship between climate availability and niche filling DIFFERENT across different quarters of climate
#ie. correlate how much available climate there is in the naturalised niche, how much is occupied


#The output from this script is used in Appendix Figure S1.11


### ###

rm(list=ls())


###################################
#set paths and variables
##################################

#set to current repo
setwd("DIRECTORY_HERE")

library(lme4)
library(scales)

#get a list of all the output from plant_PCA_nichefilling.R
yooo<-list.files(path="./IntermediateOutput/plant_PCA_nichefilling/findFilling",pattern="*.csv")


###################################
#LOAD FILES AND PROCESS FOR ANALYSIS
##################################


#load first file as a template
occRaw<-read.csv(paste0("./IntermediateOutput/plant_PCA_nichefilling/findFilling/",yooo[1]))

#load each file in turn and add to result dataframe
for (i in 2:length(yooo)){
  print(i)
  occRaw_nw<-read.csv(paste0("./IntermediateOutput/plant_PCA_nichefilling/findFilling/", yooo[i]))
  occRaw<-rbind(occRaw, occRaw_nw)
}


#write out the compiled file in case we need it?
#write.csv(occRaw,paste("IntermediateOutput/correlateFillingAvail/3Varplant_D_FILLvalues_zcor_center.csv",sep=""))


#seperate by quarters of climate niche, the coldest, hottest, driest and wettest quarters
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


#bind them together to make a single table in new format
occFull<-rbind(cold, hot, dry, wet)
head(occFull)

#what is size of table?
dim(occFull)
#how many have no valid occurrence density (in their climatic quarter)
sum(occFull$Occdens==0)
#how many have NO available climate (in their climatic quarter)
sum(occFull$Potdens==0)
#how many have occupied NONE of the available climate (in their climatic quarter)
sum(occFull$Pro==0, na.rm=T)

#log for normality, but add a very minor number to avoid zeroes
occFull$logOccdens<-log10(occFull$Occdens+min(occFull$Occdens[occFull$Occdens>0], na.rm=T)/10)
occFull$logPotdens<-log10(occFull$Potdens+min(occFull$Potdens[occFull$Potdens>0], na.rm=T)/10)
occFull$logPro<-log10(occFull$Pro+min(occFull$Pro[occFull$Pro>0], na.rm=T)/10)

#check histograms, should look vaguely normal
hist(occFull$logOccdens)
hist(occFull$logPotdens)
hist(occFull$Pro)



###########################################
#MODEL 1: is there a statistically significant amount of niche filling in different quarters?
##########################################

#basic boxplots, are there any obvious differences?
boxplot(logOccdens  ~ quar, data=occFull)

#save this one for later
plotname=paste("FinalOutput/supplementary/correlateFillingAvail/availableClimateBoxplot.pdf",sep="")
pdf(file=plotname, width=6, height=5)
boxplot(logPotdens  ~ quar, data=occFull, xlab = "Type of climate available", ylab="Total climate available (summed density logged)")
dev.off()
boxplot(logPro  ~ quar, data=occFull)


#basic comparison, is filling predicted by the quarter

#fudge a little and cut out zeroes (will come back to this)
occFull1<-occFull[occFull$Pro>0,]

#does proportion of niche filling vary by climatic quarter, run ANOVA
anovRes<-glm(logPro ~ quar, data=occFull1)
#check residuals
plot(anovRes)
#they look ok

boxplot(logPro  ~ quar, data=occFull1)
summary(anovRes)

#significant differences, cold niche filling is less than the others
#wet niche filling is highest. Niche filling in dry and hot are in the middle


#PART 2: do some quarters simply have more climate available?


#fudge a little and cut out zeroes (will come back to this)
occFull1<-occFull[occFull$Pro>0,]

#does proportion of niche filling vary by climatic quarter, run ANOVA
anovRes2<-glm(logPotdens ~ quar, data=occFull1)
#check residuals
plot(anovRes2)
#they're not great

boxplot(logPotdens  ~ quar, data=occFull1)
summary(anovRes2)
#wet has slightly less climate available (significantly so based on p<0.05), all other differences are very minor and non-significant

#as the residuals weren't great, run a non-parametric equivalent
kruskal.test(logPotdens ~ quar, data = occFull1)

#significant, run pair-wise
pairwise.wilcox.test(occFull1$logPotdens, occFull1$quar,
                     p.adjust.method = "BH")

#wet is significantly different (slightly lower) than the rest, all others are non-significant
#So availability is not the same

tapply(occFull1$logPotdens, occFull1$quar, summary)
#but differences are very small overall and p value is likely driven by large sample size



###########################################
#MODEL 2: mixed linear regression, is the slope of niche filling predicted by the direction?
##########################################

#we would expect niche filling to get higher as total climate available and occurrence density get higher
# a significant slope is fine, but  we expect this pattern to be the same across all 4 climatic directions 
#(i.e. availability, occupation rate and relationship are roughly the same)
#to do this we correlate climate availability, occupancy, and check for differences between climatic quarters
#we also account for random effects of region and species, since these are replicated across the dataset


#remove species with no occurrences
occFull2<-occFull[occFull$Occdens>0,]
plot(occFull2$logPotdens, occFull2$logOccdens, cex=0.4)

#start with complex model, add interaction
model2<-lmer(logOccdens ~ logPotdens*quar + (1|region) + (1|sp_name), data=occFull2, REML=F)

#step-wise removal. cut out the interaction, doe it make a difference?
model3<-lmer(logOccdens ~ logPotdens + (1|region) + (1|sp_name) , data=occFull2, REML=F)

#it does! Keep the complex model
anova(model3, model2)

#check residuals
plot(model2)
#not perfect, but looks ok

#check summaries
summary(model2)
fixef(model2)
summary(ranef(model2)$sp_name)

#how much variation across regions?
boxplot(logOccdens~region, data=occFull2)


#some differences in intercept (probably due to sparse data)
#but when climate is available estimates are very similar, minor differences in slope


#make some example data so we can plot model output
#there's not that much variation across regions. There is quite a bit over species, so selected one near the median
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


#plot in colour if needed
#plot(occFull2$logPotdens[occFull2$quar=="hot"], occFull2$logOccdens[occFull2$quar=="hot"], cex=0.5, col="red",
     #xlab = "Total climate availabile (summed density logged)", ylab="Total climate occuped (summed density logged)")
#points(occFull2$logPotdens[occFull2$quar=="cold"], occFull2$logOccdens[occFull2$quar=="cold"], cex=0.5, col="blue")
#points(occFull2$logPotdens[occFull2$quar=="wet"], occFull2$logOccdens[occFull2$quar=="wet"], cex=0.5, col="green")
#points(occFull2$logPotdens[occFull2$quar=="dry"], occFull2$logOccdens[occFull2$quar=="dry"], cex=0.5, col="orange")


#plot in BW with lines, save output
plotname=paste("FinalOutput/supplementary/correlateFillingAvail/availFillregression.pdf",sep="")
pdf(file=plotname, width=6, height=5)

plot(occFull2$logPotdens, occFull2$logOccdens, cex=0.5, col=alpha("black", 0.3),
     xlab = "Total climate availabile (summed density logged)", ylab="Total climate occuped (summed density logged)")

points(hotDF$logPotdens, hotDF$predVals, cex=0.2, col="red", pch=4)
points(coldDF$logPotdens, coldDF$predVals, cex=0.2, col="blue", pch=4)
points(wetDF$logPotdens, wetDF$predVals, cex=0.2, col="green", pch=4)
points(dryDF$logPotdens, dryDF$predVals, cex=0.2, col="orange", pch=4)

dev.off()






