
####################################
######RANGE FILLING CHAPTER#######
####################################

rm(list=ls())



#TABLE 1
#read in unfilling data
setwd("U:/Data/trait_data/PLANT TRAITS")
setwd("/Users/henry_hakkinen/Documents/ExeterWork/Data/trait_data/PLANT TRAITS")
unfill_plants<-read.csv("unfilling_region_kmsummary.csv",stringsAsFactors = F)
length(unique(unfill_plants$sp_name))

head(unfill_plants)
sp_list<-unique(unfill_plants$sp_name)
table_height<-ceiling(length(sp_list)/3)
species_mat<-matrix(nrow=table_height,ncol=3)

species_mat[,1] <- sp_list[1:table_height]
species_mat[,2] <- sp_list[(table_height+1):(table_height*2)]
species_mat[1:(length(sp_list)-2*table_height),3] <- sp_list[((table_height*2+1):length(sp_list))]
dim(species_mat)

write.csv(species_mat,"/Users/henry_hakkinen/Documents/ExeterWork/writeup/UnfillingChapter/FullPlantSPlist.csv",row.names = F)

#species list for birds
setwd("U:/Data/trait_data/BIRD TRAITS")
setwd("/Users/henry_hakkinen/Documents/ExeterWork/Data/trait_data/BIRD TRAITS")
unfill_birds<-read.csv("unfilling_region_kmsummary.csv")

length(unique(unfill_birds$sp_name))

sp_list<-unique(unfill_birds$sp_name)
table_height<-ceiling(length(sp_list)/3)
species_mat<-matrix(nrow=table_height,ncol=3)

species_mat[,1] <- sp_list[1:table_height]
species_mat[,2] <- sp_list[(table_height+1):(table_height*2)]
species_mat[1:(length(sp_list)-2*table_height),3] <- sp_list[((table_height*2+1):length(sp_list))]
dim(species_mat)

write.csv(species_mat,"/Users/henry_hakkinen/Documents/ExeterWork/writeup/UnfillingChapter/FullBirdSPlist.csv",row.names = F)


#species list for mammals
setwd("U:/Data/trait_data/MAMMAL TRAITS")
setwd("/Users/henry_hakkinen/Documents/ExeterWork/Data/trait_data/MAMMAL TRAITS")
table_name<-"Mammal"
unfill_mammals<-read.csv("unfilling_region_kmsummary.csv")

length(unique(unfill_mammals$sp_name))

sp_list<-unique(unfill_mammals$sp_name)
table_height<-ceiling(length(sp_list)/3)
species_mat<-matrix(nrow=table_height,ncol=3)

species_mat[,1] <- sp_list[1:table_height]
species_mat[,2] <- sp_list[(table_height+1):(table_height*2)]
species_mat[1:(length(sp_list)-2*table_height),3] <- sp_list[((table_height*2+1):length(sp_list))]
dim(species_mat)

write.csv(species_mat,"/Users/henry_hakkinen/Documents/ExeterWork/writeup/UnfillingChapter/FullMammalsSPlist.csv",row.names = F)


#Table 2
#####information on species and where they were introduced

setwd("U:/Data/trait_data/beta_reg")
setwd("/Users/henry_hakkinen/Documents/ExeterWork/Data/trait_data/beta_reg")


plant_regdf<-read.csv("plant_3var1catdataforwriteup.csv",stringsAsFactors = F)
plant_regdf<-plant_regdf[,c("sp_name","xcat1")]
names(plant_regdf)<-c("Species Name","Realm")
mammal_regdf<-read.csv("mammal_3var1catdataforwriteup.csv",stringsAsFactors = F)
head(mammal_regdf)
mammal_regdf<-mammal_regdf[,c("sp_name","xcat1")]
names(mammal_regdf)<-c("Species Name","Realm")
bird_regdf<-read.csv("bird_3var1catdataforwriteup.csv",stringsAsFactors = F)
bird_regdf<-bird_regdf[,c("sp_name","xcat1")]
names(bird_regdf)<-c("Species Name","Realm")

write.csv(plant_regdf,"/Users/henry_hakkinen/Documents/ExeterWork/writeup/UnfillingChapter/Plant_sp_region.csv")
write.csv(mammal_regdf,"/Users/henry_hakkinen/Documents/ExeterWork/writeup/UnfillingChapter/Mammal_sp_region.csv")
write.csv(bird_regdf,"/Users/henry_hakkinen/Documents/ExeterWork/writeup/UnfillingChapter/Bird_sp_region.csv")





######################################
#####EXPANSION CHAPTER####################
#########################################


#28/05/2019
#Make supplementary material info

rm(list=ls())

setwd("U:/Data/niche_shift_work")

#this is the final species list, 3 climate variables and filtered for analogue climate
plant_summary<-read.csv("analogue_disagg3Var_04032019/plant_data_analogue05042019.csv",stringsAsFactors = F)

plant_list<-unique(plant_summary$species_name)
length(plant_list)

table_height<-ceiling(length(plant_list)/3)

species_mat<-matrix(nrow=table_height,ncol=3)

species_mat[,1] <- plant_list[1:table_height]
species_mat[,2] <- plant_list[(table_height+1):(table_height*2)]
species_mat[1:(length(plant_list)-2*table_height),3] <- plant_list[((table_height*2+1):length(plant_list))]



write.csv(species_mat,"U:/writeup/ExpansionChapter/FullSPlist.csv")



#make a list of citations from SQL table


#read in the native table
native_sql<-read.csv("U:/writeup/ExpansionChapter/Plant_native_SQLtable.csv",stringsAsFactors = F,sep=";")
natur_sql<-read.csv("U:/writeup/ExpansionChapter/Plant_natur_SQLtable.csv",stringsAsFactors = F,sep=";")
#read in the naturalised table

native_sql<-native_sql[native_sql$country!="",]
natur_sql<-natur_sql[natur_sql$country!="",]



native_sql$status<-"native"
natur_sql$status<-"natur"

native_sql<-native_sql[native_sql$scientific_name %in% plant_list,]
natur_sql<-natur_sql[natur_sql$scientific_name %in% plant_list,]




#plant_list[!(plant_list %in% natur_sql$scientific_name)]
#there are some conflicting entries, we need to sort these out
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
#write.csv(native_sql_unique,"U:/writeup/ExpansionChapter/native_citations_part1.csv")

natur_sql_unique<-natur_sql_unique[order(natur_sql_unique$citation),]
#write.csv(natur_sql_unique,"U:/writeup/ExpansionChapter/natur_citations_part1.csv")


native_sql_shared<-native_sql[native_sql$check %in% natur_sql$check,]
natur_sql_shared<-natur_sql[natur_sql$check %in% native_sql$check,]
colnames(natur_sql_shared)[1]<-"source"
head(native_sql_shared)
head(natur_sql_shared)

head(native_sql_shared)


shared_info<-rbind(native_sql_shared[,c("source","citation","country","scientific_name","status")],
                  natur_sql_shared[,c("source","citation","country","scientific_name","status")])
write.csv(shared_info,"U:/writeup/ExpansionChapter/shared_statusinfo.csv")


shared_info<-read.csv("U:/writeup/ExpansionChapter/shared_statusinfo.csv",stringsAsFactors = F)
sp_shared<-unique(shared_info$scientific_name)

save_df<-shared_info[0,]
fail_df<-shared_info[0,]

i<-4
#for(i in 1:length(sp_shared)){
for(i in 1:10){
  sp_name<-sp_shared[i]
  print(as.character(sp_name))
  
  sp_rows<-shared_info[shared_info$scientific_name==sp_name,]
  
  
  
  search<-grepl("manual",sp_rows$source)
  
  if(sum(search>0)){
    save_row<-sp_rows[grepl("manual",sp_rows$source),]
    fail_row<-sp_rows[ -c(which(grepl("manual",sp_rows$source))),]
    
    
    save_df[nrow(save_df+1),]<-save_row
    fail_df[nrow(fail_df+1),]<-fail_row
    
  }else{
    fail_df<-sp_rows
    
  }

}




#make plot of biogeographic regions
library("raster")
library("ggplot2")
shape <- shapefile("U:/Data/biogeographic_zonesV2.shp")


plot(shape,col=as.factor(shape$Name))
#create a new aesthetic df. Colour by a.diff.sum
shp_df <- broom::tidy(shape, region = "Name")
table(shp_df$id)
as.numeric(shp_df$id)

x11()

map <- ggplot() + geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = (id)), colour = "black") +
  theme_void()+labs(fill = "Realm")
map


