#miscellaneous functions used throughout repo for PCA purposes

euc.dist.center <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

euc.dist.edge <- function(x){
  x1<-x[1:2]
  x2<-x[3:4]
  sqrt(sum((x1 - x2) ^ 2))
  
} 

#find the angle between two points and return a number of times

rad.ang.center <- function(x1, x2){
  dx = x1[1] - x2[1]
  dy = x1[2] - x2[2] 
  
  if(!is.numeric(dx)){dx<-as.numeric(dx)}
  if(!is.numeric(dy)){dy<-as.numeric(dy)}
  
  theta = atan2(dy,dx)
  
  if(theta<0){ theta<- theta+ 2*pi }
  return(theta)
} 

rad.ang.edge <- function(x){
  x1<-x[1:2]
  x2<-x[3:4]
  
  dx = x1[1] - x2[1]
  dy = x1[2] - x2[2] 
  
  if(!is.numeric(dx)){dx<-as.numeric(dx)}
  if(!is.numeric(dy)){dy<-as.numeric(dy)}
  
  theta = atan2(dy,dx)
  
  if(theta<0){ theta<- theta+ 2*pi }
  return(theta)
} 


#takes lat lon point, projects it and then finds the continent/region the point is in
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




near_point<-function(x,stab_p){
  x<-cbind(x[1],x[2])
  point<-SpatialPoints(x)
  
  if(class(stab_p)!="SpatialPoints"){stab_p<-SpatialPoints(stab_p)}
  
  dis_df<-gDistance(point,stab_p,byid = T)
  stab_p2<-as.data.frame(stab_p)
  near_p<-stab_p2[which.min(dis_df),1:2]
  near_p<-cbind(min(dis_df),near_p)
  return(near_p)
}

#running near point takes ages, so we filter via another function (don't murder the workflow)
near_point2<-function(x,df){
  search_df<-paste(round(df[,"Axis1"],digits=1),round(df[,"Axis2"],digits=1),sep="_")
  point<-paste(round(x[,1],digits=1),round(x[,2],digits=1),sep="_")
  
  find_list<-which(search_df %in% point)
  new_search<-df[find_list,c("Axis1","Axis2")]
  
  #more formal search now that we have a sublist
  row_near<-near_point(x,new_search)
  res_near<-df[rownames(row_near),]
  
  return(res_near)
}




find_slice<-function(data,dir1){
  
  ####find 45degrees either side of this angle
  #there will always be at least one direction
  #45 degrees in radians is 0.7853982

  
  dir1_high<-dir1+0.7853982
  dir1_low<-dir1-0.7853982
  

  
  if(dir1_low<0){dir1_low<-dir1_low+2*pi}
  
  
  #there is an annoying habit of circres of producing directions of over 2*pi, correct for this.
  
  if(dir1_high>(2*pi)){dir1_high<-dir1_high-(2*pi)}
  if(dir1_low>(2*pi)){dir1_low<-dir1_low-(2*pi)}
  
  
  
  #do it again!
  if(dir1_high>(2*pi)){dir1_high<-dir1_high-(2*pi)}
  if(dir1_low>(2*pi)){dir1_low<-dir1_low-(2*pi)}
  
  #do it again!
  if(dir1_high>(2*pi)){dir1_high<-dir1_high-(2*pi)}
  if(dir1_low>(2*pi)){dir1_low<-dir1_low-(2*pi)}
  
  if(dir1>(2*pi)){print("FUUUCK")}

  
  ####take only points that lie within this angle
  exp_df_dir1<-data[which(data$exp_dir>dir1_low & data$exp_dir<dir1_high),]
  
  #build a special scenario if our angle of error crosses 0.
  if((dir1_high>0 &dir1_high<0.5*pi)  & (dir1_low<2*pi & dir1_low>1.5*pi)){
    exp_df_dir1<-data[which(data$exp_dir>dir1_low | data$exp_dir<dir1_high),]
  }
  
  
  return(exp_df_dir1)
}


#wt.mean
wt.mean <- function(x,wt) {
  s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
  return( sum(wt * x)/sum(wt) ) #return the mean
}

#' @rdname wt.mean
#' @export
wt.var <- function(x,wt) {
  s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
  xbar = wt.mean(x,wt) #get the weighted mean
  return( sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))) ) #return the variance
} 

#' @rdname wt.mean
#' @export
wt.sd <- function(x,wt) { 
  return( sqrt(wt.var(x,wt)) ) #return the standard deviation
} 


#Centre of Gravity or Mass calculations for spatial data

COGravity <- function(x,y=NULL,z=NULL,wt=NULL) {
  #check if raster from sp or raster package and convert if necessary
  if (any(class(x) %in% 'RasterLayer')) x = asc.from.raster(x)
  if (any(class(x) == 'SpatialGridDataFrame')) x = asc.from.sp(x)
  #check if x is vector or matrix
  if (is.vector(x)) { #if the data is a vector...do calculations
    if (is.null(wt)) {	#if no weighting supplied, calculate means & standard deviations
      out = c(COGx=mean(x,na.rm=TRUE),COGx.sd=sd(x,na.rm=TRUE))
      if (!is.null(y)) out = c(out,COGy=mean(y,na.rm=TRUE),COGy.sd=sd(y,na.rm=TRUE))
      if (!is.null(z)) out = c(out,COGz=mean(z,na.rm=TRUE),COGz.sd=sd(z,na.rm=TRUE))
    } else { #if weighting supplied, calculate weighted means and variances to get COG
      out = c(COGx=wt.mean(x,wt),COGx.sd=wt.sd(x,wt))
      if (!is.null(y)) out = c(out,COGy=wt.mean(y,wt),COGy.sd=wt.sd(y,wt))
      if (!is.null(z)) out = c(out,COGz=wt.mean(z,wt),COGz.sd=wt.sd(z,wt))
    }
  } else if (any(class(x) == 'asc')) { #if x is of class 'asc'
    if (is.null(wt)) { #if wt is null then assume that values in x are the weights
      pos = as.data.frame(which(is.finite(x),arr.ind=TRUE))
      pos$x = getXYcoords(x)$x[pos$row]
      pos$y = getXYcoords(x)$y[pos$col]
      pos$wt = x[cbind(pos$row,pos$col)]
      out = c(COGx=wt.mean(pos$x,pos$wt),COGx.sd=wt.sd(pos$x,pos$wt),COGy=wt.mean(pos$y,pos$wt),COGy.sd=wt.sd(pos$y,pos$wt))
    } else { #if wt is supplied, it must be of the same dim as x and then the values of x are assumed to be your z
      if (!all(dim(x)==dim(wt))) stop('the grids for x & weights must be of the same dimensions')
      pos = as.data.frame(which(is.finite(x),arr.ind=TRUE))
      pos$x = getXYcoords(x)$x[pos$row]
      pos$y = getXYcoords(x)$y[pos$col]
      pos$z = x[cbind(pos$row,pos$col)]
      pos$wt = wt[cbind(pos$row,pos$col)]
      out = c(COGx=wt.mean(pos$x,pos$wt),COGx.sd=wt.sd(pos$x,pos$wt),COGy=wt.mean(pos$y,pos$wt),COGy.sd=wt.sd(pos$y,pos$wt),COGz=wt.mean(pos$z,pos$wt),COGz.sd=wt.sd(pos$z,pos$wt))
    }
    
  }
  # return the output	
  return(out)
}
