
#various circular calculations used throughout repo
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

plot_circMLE.custom<-function (data, table, model, bins, shrink, col, lwd, lty) {
  if (is.null(data) | missing(data)) 
    stop("Please provide input data vector")
  if (is.null(table) | missing(table)) 
    stop("Please provide the output list from the \"circ_mle\" function")
  if (missing(model)) 
    model = rownames(table$results)[1]
  if (length(model) != 1) 
    stop("Only 1 model can be specified")
  if (!all(model %in% rownames(table$results))) 
    stop("Model not specified correctly.")
  if (missing(bins)) 
    bins = 18
  else bins = bins
  if (length(bins) != 1 | !is.numeric(bins)) 
    stop("Must specify the number of bins as a numeric vector of length 1")
  if (missing(shrink)) 
    shrink = 1.5
  else shrink = shrink
  if (length(shrink) != 1 | !is.numeric(shrink)) 
    stop("Must specify \"shrink\" as a numeric vector of length 1")
  if (missing(col)) 
    col = c("grey", "red", "black", "black")
  else col = col
  if (missing(lwd)) 
    lwd = c(2, 2, 2)
  else lwd = lwd
  if (missing(lty)) 
    lty = c("solid", "dashed", "dashed")
  else lty = lty
  data = check_data(data)
  params = circularp(data)

  
  rose.diag(data, bins = bins, col = col[1], shrink = shrink, axes=FALSE)
  #arrows.circular(mean.circular(data, na.rm = T), y = rho.circular(data, 
                                                                   #na.rm = T), col = col[2], lwd = lwd[1], lty = lty[1])
  model.vector = table$results[which(rownames(table$results) == 
                                       model), 2:6]
  q1 = as.circular(model.vector[1], control.circular = params)
  q2 = as.circular(model.vector[4], control.circular = params)
  k1 = as.numeric(model.vector[2])
  k2 = as.numeric(model.vector[5])
  l = as.numeric(model.vector[3])
  plot.function.circular(function(x) dmixedvonmises(x, q1, 
                                                    q2, k1, k2, l), add = T, lwd = lwd[2], shrink = shrink, 
                         lty = lty[2], col = col[3], zero = params$zero, rotation = params$rotation)
  if (any(c("M2A", "M2B", "M2C") == model)) {
    arrows.circular(q1, col = col[4], lwd = lwd[3], lty = lty[3])
  }
  if (any(c("M3A", "M3B", "M4A", "M4B", "M5A", "M5B") == model)) {
    arrows.circular(circular(q1, zero = 0, rotation = "counter"), 
                    col = col[4], lwd = lwd[3], lty = lty[3])
    arrows.circular(circular(q2, zero = 0, rotation = "counter"), 
                    col = col[4], lwd = lwd[3], lty = lty[3])
  }
}


plot.circustom<-function(x, plot.type=c("circle", "line"), points.plot=FALSE, rp.type="p", type="l",
                                   line.col=1, points.col="grey", points.pch=1, xlim=NULL, ylim=NULL, radial.lim=NULL, xlab=NULL, ylab=NULL, 
                                   labels=NULL, label.pos=NULL, units=NULL, zero=NULL, clockwise=NULL, main=NULL){
  
  xcircularp <- attr(x$x, "circularp")
  ycircularp <- attr(x$y, "circularp")
  if (is.null(xcircularp) && is.null(ycircularp)) 
    stop("the component 'x' and/or the component 'y' of the object must be of class circular")
  plot.type <- match.arg(plot.type)
  if (is.circular(x$datax) && !is.circular(x$datay)){
    if (is.null(units)) units <- xcircularp$units
    template <- xcircularp$template
    x$x <- conversion.circular(x$x, units = "radians", modulo = "2pi")
    x$datax <- conversion.circular(x$datax, units = "radians", modulo = "2pi")
    attr(x$x, "class") <- attr(x$x, "circularp") <- NULL  
    attr(x$datax, "class") <- attr(x$datax, "circularp") <- NULL
    if (plot.type=="line" & units == "degrees") {
      x$x <- x$x/pi * 180
      x$datax <- x$datax/pi * 180
      if (is.null(xlab)) xlab <- "degrees"
    }
    if (plot.type=="line" & units == "hours") {
      x$x <- x$x/pi * 12
      x$datax <- x$datax/pi * 12
      if (is.null(xlab)) xlab <- "hours"
    }
    if (is.null(xlab)) xlab <- "radians"
    if (is.null(ylab)) ylab <- ""
    
  } else if (is.circular(x$datax) && is.circular(x$datay)){
    if (plot.type=="circle") plot.type <- "torus"
    template <- xcircularp$template
    if (is.null(units)) units <- xcircularp$units
    x$x <- conversion.circular(x$x, units = "radians", modulo = "2pi")
    x$datax <- conversion.circular(x$datax, units = "radians", modulo = "2pi")
    attr(x$x, "class") <- attr(x$x, "circularp") <- NULL  
    attr(x$datax, "class") <- attr(x$datax, "circularp") <- NULL
    
    template <- ycircularp$template
    x$y <- conversion.circular(x$y, units = "radians", modulo = "2pi")
    x$datay <- conversion.circular(x$datay, units = "radians", modulo = "2pi")
    attr(x$y, "class") <- attr(x$y, "circularp") <- NULL  
    attr(x$datay, "class") <- attr(x$datay, "circularp") <- NULL
    
    x$datax[x$datax>pi]<-x$datax[x$datax>pi]-2*pi
    x$datay[x$datay>pi]<-x$datay[x$datay>pi]-2*pi
    x$x[x$x>pi]<-x$x[x$x>pi]-2*pi
    x$y[x$y>pi]<-x$y[x$y>pi]-2*pi

    if (plot.type=="line" & units == "degrees") {
      x$x <- x$x/pi * 180
      x$datax <- x$datax/pi * 180
      x$y <- x$y/pi * 180
      x$datay <- x$datay/pi * 180
      if (is.null(xlab)) xlab <- "degrees"
      if (is.null(ylab)) ylab <- "degrees"
    }
    if (plot.type=="line" & units == "hours") {
      x$x <- x$x/pi * 12
      x$datax <- x$datax/pi * 12
      x$y <- x$y/pi * 12
      x$datay <- x$datay/pi * 12
      if (is.null(xlab)) xlab <- "hours"
      if (is.null(ylab)) ylab <- "hours"
    }
    if (is.null(xlab)) xlab <- "radians"
    if (is.null(ylab)) ylab <- "radians"
    
  } else if (!is.circular(x$datax) && is.circular(x$datay)){
    if (plot.type=="circle") plot.type <- "cylinder"
    if (is.null(units)) units <- ycircularp$units
    template <- ycircularp$template
    x$y <- conversion.circular(x$y, units = "radians", modulo = "2pi")
    x$datay <- conversion.circular(x$datay, units = "radians", modulo = "2pi")
    attr(x$y, "class") <- attr(x$y, "circularp") <- NULL  
    attr(x$datay, "class") <- attr(x$datay, "circularp") <- NULL
    
    if (plot.type=="line" & units == "degrees") {
      x$y <- x$y/pi * 180
      x$datay <- x$datay/pi * 180
      if (is.null(ylab)) ylab <- "degrees"
    }
    if (plot.type=="line" & units == "hours") {
      x$y <- x$y/pi * 12
      x$datay <- x$datay/pi * 12
      if (is.null(ylab)) ylab <- "hours"
    }
    if (is.null(xlab)) xlab <- ""
    if (is.null(ylab)) ylab <- "radians"
  }	
  if (is.null(main)) main <- deparse(x$call)
  
  if (plot.type == "line") {
    if (is.null(xlim)) xlim <- range(c(x$x, x$datax))
    if (is.null(ylim)) ylim <- range(c(x$y, x$datay))
    xorder <- order(x$x)
    x$x <- x$x[xorder]
    x$y <- x$y[xorder]
    plot.default(x$x, x$y, type = type, col=line.col, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, main=main)
    if (points.plot) points(x$datax, x$datay, col=points.col, pch=points.pch)
  } else {
    if (plot.type=="torus" | plot.type=="circle"){
      if (is.null(labels) || is.null(label.pos)){
        if (is.null(units)){
          if (xcircularp$units=="degrees"){
            labels <- c("0","45","90","135","180","225","270","315")
          }else if (xcircularp$units=="radians"){
            labels <- c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),
                        expression(5*pi/4),expression(3*pi/2),expression(7*pi/4))
          }else if (units=="hours"){
            labels <- c("0h","3h","6h","9h","12h","15h","18h","21h")
          }
        } else if (units=="degrees"){
          labels <- c("0","45","90","135","180","225","270","315")
        } else if (units=="radians"){
          labels <- c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),
                      expression(5*pi/4),expression(3*pi/2),expression(7*pi/4))
        } else if (units=="hours"){
          labels <- c("0h","3h","6h","9h","12h","15h","18h","21h")
        }
        label.pos <- rad.pos <- seq(0,7*pi/4,by=pi/4)
      }
    }
    if (plot.type=="torus"){
      torus <- parametric3d(fx = function(u,v) (1+0.25*cos(v))*cos(u),
                            fy = function(u,v) (1+0.25*cos(v))*sin(u),
                            fz = function(u,v) 0.25*sin(v),
                            u = seq(0,2*pi,length.out=50),
                            v = seq(0,2*pi,length.out=50), 
                            engine = "none", color="paleturquoise", alpha=0.5)
      drawScene.rgl(torus)
      xx <-cos(x$x)*(1+0.25*cos(x$y))
      yy <- sin(x$x)*(1+0.25*cos(x$y))
      zz <- 0.25*sin(x$y)
      lines3d(xx, yy, zz, col=line.col)
      #texts3d(1.25*cos(label.pos), 1.25*sin(label.pos), 0, texts=labels)
      labels
      
      if (points.plot) {
        xx <- cos(x$datax)*(1+0.25*cos(x$datay))
        yy <- sin(x$datax)*(1+0.25*cos(x$datay))
        zz <- 0.25*sin(x$datay)
        points3d(xx, yy, zz, col=points.col)
      }
      
    } else if (plot.type=="cylinder"){
      R<- diff(range(x$datax))/8
      cylinder <- parametric3d(fx = function(u,v) v,
                               fy = function(u,v) R*cos(u),
                               fz = function(u,v) R*sin(u),
                               u = seq(0,2*pi,length.out=50),
                               v = seq(min(x$datax),max(x$datax),length.out=50),
                               engine = "none", color="paleturquoise", alpha=0.5)
      drawScene.rgl(cylinder)	
      xx <- x$x
      yy <- R*cos(x$y)
      zz <- R*sin(x$y)
      lines3d(xx, yy, zz, col=line.col)
      if (is.null(units)){
        if (xcircularp$units=="degrees"){
          labels <- c("0","90","180","270")
        }else if (xcircularp$units=="radians"){
          labels <- c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2))
        }else if (units=="hours"){
          labels <- c("0h","6h","12h","18h")
        }
      } else if (units=="degrees"){
        labels <- c("0","90","180","270")
      } else if (units=="radians"){
        labels <- c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2))
      } else if (units=="hours"){
        labels <- c("0h","6h","12h","18h")
      }
      label.pos <- rad.pos <- seq(0,3*pi/2,by=pi/2)
      yy <- R*cos(label.pos)
      zz <- R*sin(label.pos)
      #texts3d(min(x$datax), yy, zz, texts=labels)
      
      
      #texts3d(max(x$datax), yy, zz, texts=labels)
      xlim <- range(dist)
      xx <- pretty(xlim)
      xx <- xx[xx>=min(dist) & xx<=max(dist)]
      yy <- R
      zz <- R
      #texts3d(xx, yy, zz, texts=xx)
      if (points.plot) {
        xx <- x$datax
        yy <- R*cos(x$datay)
        zz <- R*sin(x$datay)
        points3d(xx, yy, zz, col=points.col)
      }
    }else{
      if (is.null(zero)){
        if (template == "geographics" | template == "clock24") zero <- pi/2
        else zero <- xcircularp$zero
      }
      if (is.null(clockwise)){
        if (template == "geographics" | template == "clock24") clockwise <- TRUE
        else clockwise <- ifelse(xcircularp$rotation=="counter",FALSE,TRUE)
      }
      if (is.null(radial.lim)) radial.lim <- range(c(x$datay,x$y))
      options (warn=-1)
      radial.plot(x$y, x$x, rp.type=rp.type, line.col=line.col, labels=labels, label.pos=label.pos,
                  start=zero, clockwise=clockwise, radial.lim=radial.lim, main=main) 
      if (points.plot) {
        radial.plot(x$datay, x$datax, rp.type="s", start=zero, clockwise=clockwise, radial.lim=radial.lim, 
                    point.col=points.col, point.symbols=points.pch, add=TRUE)
      }
      return(invisible(list(zero=zero, clockwise=clockwise, radial.lim=radial.lim)))
    }
  }
}

euc.dist.center <- function(x1, x2) (sum((x1 - x2) ^ 2));


