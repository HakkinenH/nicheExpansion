
##############################################
### CompileSuppTables.R ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 

# There are a number of tables and maps to compile for the supplementary material
# Each is quite minor, so I compiled them into a single file

#this file can be run at any time without running anything else,
#however you need (in rawdata):

#1) a species list
#2) a list of sources for native/naturalised (provided is our list)
#3) a map of biogeographic regions

#In this file it will 
#1) make a neater species list table for the supplementary material
#2) make a list of citations of all sources used to categorise data as native/naturalised from SQL table
# NOTE THAT THIS IS NOT THAT HELPFUL BY ITSELF. The citations are incomplete. To see the full citations see "supp_citations_Dict.xlsx" in rawdata/speciesLists
#3) Make a neater map of the biogeographic zones, with colour and legends

#Output from this code is used in Appendix Table S1.1; Appendix Table S1.2;


### ###



rm(list=ls())



###################################
#set paths and variables
##################################


#set to current repo
setwd("DIRECTORY_HERE")


library("raster")
library("ggplot2")

#this is the final species list, use throughout for analysis
plant_summary<-read.csv("./RawData/speciesList/plant_data_summary0607.csv", stringsAsFactors=F)

#needed for part 2:
#Plant_nativesources_SQLtable.csv and natur equivalent are taken from our SQL table with out full species list and info
#this is just a dump into a csv

#read in the native table info
native_sql<-read.csv("./RawData/speciesList/Plant_nativesources_SQLtable.csv",stringsAsFactors = F,sep=";")
#read in the naturalised table info
natur_sql<-read.csv("./RawData/speciesList/Plant_natursources_SQLtable.csv",stringsAsFactors = F,sep=";")



###################################
#make a neater table for the supplementary material
##################################


plant_list<-unique(plant_summary$Species_name)
length(plant_list)

table_height<-ceiling(length(plant_list)/3)

species_mat<-matrix(nrow=table_height,ncol=3)

species_mat[,1] <- plant_list[1:table_height]
species_mat[,2] <- plant_list[(table_height+1):(table_height*2)]
species_mat[1:(length(plant_list)-2*table_height),3] <- plant_list[((table_height*2+1):length(plant_list))]


write.csv(species_mat,"./FinalOutput/Supplementary/CompileSuppTablesMaps/FullSPlist.csv")




###################################
#make a list of citations of all sources used 
#to categorise data as native/naturalised from SQL table
##################################


#cut blanks
native_sql<-native_sql[native_sql$country!="",]
natur_sql<-natur_sql[natur_sql$country!="",]


#set status
native_sql$status<-"native"
natur_sql$status<-"natur"

#cut down to only those in the final species list
native_sql<-native_sql[native_sql$scientific_name %in% plant_list,]
natur_sql<-natur_sql[natur_sql$scientific_name %in% plant_list,]


#plant_list[!(plant_list %in% natur_sql$scientific_name)]
#if there are some conflicting entries, we need to sort these out
#if in doubt use the manual entry

native_sql$check<-paste(native_sql$scientific_name,native_sql$country,sep="-")
natur_sql$check<-paste(natur_sql$scientific_name,natur_sql$country,sep="-")

sum(native_sql$check %in% natur_sql$check)
sum(natur_sql$check %in% native_sql$check)


native_sql_unique<-native_sql[!(native_sql$check %in% natur_sql$check),]
natur_sql_unique<-natur_sql[!(natur_sql$check %in% native_sql$check),]
dim(native_sql_unique)
dim(natur_sql_unique)

native_sql_unique<-native_sql_unique[order(native_sql_unique$source),]

natur_sql_unique<-natur_sql_unique[order(natur_sql_unique$citation),]

#neaten and combine tables
native_sql_shared<-native_sql[native_sql$check %in% natur_sql$check,]
natur_sql_shared<-natur_sql[natur_sql$check %in% native_sql$check,]
colnames(natur_sql_shared)[1]<-"source";colnames(natur_sql_shared)[4]<-"citation"


shared_info<-rbind(native_sql_shared[,c("source","citation","country","scientific_name","status")],
                  natur_sql_shared[,c("source","citation","country","scientific_name","status")])
shared_info$source<-paste0(shared_info$source, " (", shared_info$citation, ")")
shared_info$citation<-NULL

#save the output
#NOTE that the output is not that useful as is, as citations are not properly formatted
#Use "./RawData/speciesList/supp_citations_Dict.xlsx" for full citations
write.csv(shared_info,"./FinalOutput/Supplementary/CompileSuppTablesMaps/NativeNaturSources.csv")





