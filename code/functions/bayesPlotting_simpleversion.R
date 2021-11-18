
#run the OccExpan_BetaBayes file first and then run this

#provide model run
library("plyr")

pars = model_run$BUGSoutput$sims.list



####produce your output tables
par_list<-names(pars)[4:(cont_n+4)]
par_list


####SUMMARY1:################
#overall effect of parameters on filling, ignoring random/categorical effects
par_summary1<-data.frame(matrix(nrow=(cont_n+1),ncol=5))
colnames(par_summary1)<-c("parameter","mean_estimate","95%_estimate","5%_estimate","consistent")


i<-2
for(i in 1:length(par_list)){
  par_name<-par_list[i]
  
  if(cat_n>0){
    
    if(par_name=="alpha"){par_name<-"mu_alpha"}
    
    if(grepl("beta", par_name)){par_name<-paste("mu_",par_name,sep="")}
  }
  
  
  
  #sub<-pars[[par_name]]
  #par_mean<-mean(sub)
  #par_05<-quantile(sub, probs=0.05)
  #par_95<-quantile(sub, probs=0.95)
  
  par_summary1[i,1]<-par_name
  par_summary1[i,2]<-summary_tab1[par_name,1]
  par_summary1[i,3]<-summary_tab1[par_name,7]
  par_summary1[i,4]<-summary_tab1[par_name,3]
  
  if(sign(par_summary1[i,3])==sign(par_summary1[i,4])){
    par_summary1[i,5]<-"TRUE"
  }
  
}

x_list_add<-c("intercept",x_list)
par_summary1$parameter_name<-x_list_add


par_summary1

write.csv(par_summary1,paste(output,table_name,"_bayesresults_mean.csv",sep=""),row.names = F)



####SUMMARY 2######################
#separate by category (whatever they are, and give useful summary tables)

if(cat_n>0){
  for(i in 1:cat_n){
    print(i)
    cat_list<-levels(unfilling_df[,paste("xcat",i,sep="")])
    
    n_nrow<-length(par_list)*length(cat_list)
    par_summary2<-data.frame(matrix(nrow=n_nrow,ncol=7))
    colnames(par_summary2)<-c("factor","parameter","mean_estimate","95%_estimate","5%_estimate","consistent","parameter_name")
    
    par_summary2[1:length(cat_list),"parameter"]<-"alpha"
    par_summary2[1:length(cat_list),"factor"]<-cat_list
    par_summary2[1:length(cat_list),"parameter_name"]<-x_list_add[1]
    
    
    
    for(q in 2:length(par_list)){
      par_now<-par_list[q]
      ind1<-length(cat_list)*(q-1)+1
      ind2<-length(cat_list)*(q)
      
      par_summary2[c(ind1:ind2),"parameter"]<-par_now
      par_summary2[c(ind1:ind2),"factor"]<-cat_list
      par_summary2[c(ind1:ind2),"parameter_name"]<-x_list_add[q]
      
    }
    
    #empty framework for results, now fill it up!
    for(n in 1:length(par_list)){
      par_name<-par_list[n]
      
      sub<-pars[[par_name]]
      #we have to be a bit clever to make this script generic
      #if we're working in dimension 2 (columns of categories), we need to rbind all the 3rd dimensions
      #if we're working in dimension 3 (separate arrays of categories), we need to cbind all the 2nd dimension elements
      if(i==1){par_summ<-apply(sub, 2, function(x) rbind(x))}
      if(i==2){par_summ<-apply(sub, 3, function(x) cbind(x))}
      #if you have more than 2 categories, more code needed!
      
      par_mean<-c()
      par_05<-c()
      par_95<-c()
      
      for(q in 1:length(cat_list)){
        sub2<-par_summ[,q]
        
        par_mean<-c(par_mean,mean(sub2))
        par_05<-c(par_05,quantile(sub2, probs=0.05))
        par_95<-c(par_95,quantile(sub2, probs=0.95))
        
      }
      par_summary2[which(par_summary2$parameter==par_name),3]<-par_mean
      par_summary2[which(par_summary2$parameter==par_name),4]<-par_95
      par_summary2[which(par_summary2$parameter==par_name),5]<-par_05
      
    }
    
    for(p in 1:nrow(par_summary2)){
      if(sign(par_summary2[p,3]) == sign(par_summary2[p,4]) &
         sign(par_summary2[p,3]) == sign(par_summary2[p,5]) &
         sign(par_summary2[p,4]) == sign(par_summary2[p,5])){
        par_summary2[p,6] <- "TRUE"
      }
      
    }
    
    output_name<-paste(output,table_name,"_bayesresults_",cat_names[i],".csv",sep="")
    write.csv(par_summary2,output_name,row.names = F)
    assign(paste("rem_summary",i,sep=""), par_summary2)
    
  }
}



####produce your output plots

param_table<-data.frame(matrix(ncol=3,nrow=cont_n))
colnames(param_table) <- c("param_name","par_ind","mod_ind")

param_table$param_name <- x_list
param_table$par_ind <- which(grepl("xvar",colnames(unfilling_df)))
param_table$mod_ind <- par_list[-1]


param_table

#shrink the palette if necessary
pal_n<-pal[names(pal)%in% levels(unfilling_df$xcat1)]

pal_df<-as.data.frame(pal_n)
pal_df$xcat1<-rownames(pal_df)
pal_df

#assign those colours to the unfilling table
#col_pal<-join(unfilling_df,pal_df,by="xcat1")
#col_pal<-as.character(col_pal$pal_n)


##################plot type 1: plot mean effects ignoring hierachical effects#######
par_list
par_mean<-par_list

if (cat_n>0){
  par_mean[grepl("beta",par_mean)]<-paste("mu_",par_mean[grepl("beta",par_mean)],sep="")
}

param_table
n1<-1
for(n1 in 1:nrow(param_table)){
  
  
  #select 
  par_ind<-param_table[n1,2]
  mod_ind<-param_table[n1,3]
  mod_ind<-paste("mu_",mod_ind,sep="")
  
  param_name<-param_table[n1,1]
  
  n_sims = model_run$BUGSoutput$n.sims
  
  #assign a median to all other parameters to hold steady while building predictor plots
  for(n2 in 1:nrow(param_table)){
    par_vec<-unfilling_df[,param_table[n2,"par_ind"]]
    assign(paste("x_grid", n2, sep = ""), median(par_vec))
  }
  
  #build a range of x values for line of best fit
  x_int = pretty(unfilling_df[,par_ind], n = 100)
  assign(paste("x_grid", n1, sep = ""), x_int) 
  
  x_plot<-pretty(unfilling_df[,par_ind], n = 100)
  
  grid_list<-c(ls(pattern="x_grid"))
  
  grid_df<-vector("list",length(grid_list))
  for(gr in 1:length(grid_list)){
    now_gr<-grid_list[gr]
    ndf<-get(now_gr)
    grid_df[gr]<-list(ndf)
  }
  
  
  pars = model_run$BUGSoutput$sims.list
  #skip fit and fitnew parameters
  pars<-pars[4:16]
  
  now_par<-grid_df[[n1]]
  
  
  
  pred = matrix(NA, ncol = length(x_plot), nrow = n_sims)

  for(j in 1:n_sims) {
    if(cat_n==0){sub_t<-pars$alpha[j,]}
    if(cat_n==1){sub_t<-pars$mu_alpha[j,]}
    if(cat_n==2){sub_t<-pars$mu_alpha[j,,]}
    
    for (k in 2:length(par_mean)){
      
      curr_par<-par_mean[k]

      if(cat_n==0){
        sub_t<-sub_t + pars[[curr_par]][j,] * unlist(grid_df[k-1])
      }
      pars[[curr_par]]
      if(cat_n==1){
        sub_t<-sub_t + pars[[curr_par]][j,] * unlist(grid_df[k-1])
      }
      if(cat_n==2){
        sub_t<-sub_t + pars[[curr_par]][j,,] * unlist(grid_df[k-1])
      }
    }
    pred[j,] = inv.logit(sub_t)
  }

  
  pdf(file=paste(output,table_name,"_",param_name,"_mean.pdf",sep=""),width=4,height=4)
  
  #par(xpd = T, mar = par()$mar + c(0,0,0,8))
  if (cat_n==0){
    plot(unfilling_df[,par_ind], unfilling_df$prop_filling, col="DarkBlue",
         ylab = 'Proportion niche expansion', xlab = param_name,
         xlim=c(min(unfilling_df[,par_ind]),max(unfilling_df[,par_ind])),ylim=c(0,1), pch=19,cex=0.8)
    #legend(2,1.04, sort(as.character(dis_names)), col = pal_n[order(dis_names)],lty = 1,lwd=2,bty="n")
    
  }else{
    plot(unfilling_df[,par_ind], unfilling_df$prop_filling, col=col_pal,
         ylab = 'Proportion niche expansion', xlab = param_name,
         xlim=c(min(unfilling_df[,par_ind]),max(unfilling_df[,par_ind])),ylim=c(0,1), pch=19,cex=0.8)
    #legend(2,1.04, sort(as.character(dis_names)), col = pal_n[order(dis_names)],lty = 1,lwd=2,bty="n")
  }
  

  
  int_test<-as.data.frame(HPDinterval(as.mcmc(pred), prob=0.50))
  polygon(c(now_par,rev(now_par)),c(int_test$lower,rev(int_test$upper)), col="#0000ff60", border=NA)
  
  int_test<-as.data.frame(HPDinterval(as.mcmc(pred), prob=0.95))
  polygon(c(now_par,rev(now_par)),c(int_test$lower,rev(int_test$upper)), col="#0000ff60", border=NA)
  
  pred_summary = apply(pred,2,'quantile',probs = c(0.05, 0.5, 0.95))
  line_r<-par_summary1[which(par_summary1[,"parameter_name"] == param_name),"consistent"]
  
  
  line_type<-2
  if(!is.na(line_r)){line_type<-1}
  #points(curr_dat[,par_ind], curr_dat$prop_filling, col = pal[i],pch=".",cex=5) 
  

  lines(x_plot , pred_summary[1,], col = "black",lty=3,lwd=2)  
  lines(x_plot , pred_summary[2,], col = "black",lty=line_type,lwd=2)  
  lines(x_plot , pred_summary[3,], col = "black",lty=3,lwd=2)  
  
  dev.off()
}

