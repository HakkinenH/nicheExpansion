

#lookup function to find country and info that a point lies in
clookup<-function(x){
  point<-x[1:2]
  if (is.na(point[1]) != TRUE && is.na(point[2]) != TRUE){
    sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(shape)))
    rest<<-over(sp2, shape)
  }else{
    #print ("failed")
    fail<<-fail+1
  }
} 

#lookup function to find country and info that a point lies in
clookup2<-function(x, shape){
  point<-x[1:2]
  if (is.na(point[1]) != TRUE && is.na(point[2]) != TRUE){
    sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(shape)))
    rest<<-over(sp2, shape)
  }else{
    #print ("failed")
    fail<<-fail+1
  }
} 


#function to remove NA values from raster bricks
remove_na<-function(x){
  #print(x)
  #print(class(x))
  #x$expansion[is.na(x$expansion)] <- 0
  #x$unfilling[is.na(x$unfilling) m] <- 0
  #x$natur_occ[is.na(x$natur_occ)] <- 0
  x[is.na(x)]<-0
  return(x)
}

