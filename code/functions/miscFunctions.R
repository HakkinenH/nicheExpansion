

#lookup function
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