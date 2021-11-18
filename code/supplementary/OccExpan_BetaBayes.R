
##############################################
### OccExpan_BetaBayes.R ###
##############################################

### META ###
# HHakkinen
# Complete Date: 01/07/2021
# University of Exeter
# Code repo used to support:
#   "Plant naturalisations are constrained by temperature but released by precipitation"
# 

#This is a simpler, generic version of several bayesian scripts I made.
#it runs a hierarchical bayesian model using a beta distribution to explain the proportion of species range filling/expansion

# 1) input a file. The file expects the y response column to be called 'prop_filling' though this can be changed
# 2) Specify the continuous parameters, give the column names and give the label name (what would you like the parameter to look like on plots)
# the continuous parameters will be normalised to have a mean of 0 and an SD of 2, BUT do any data transformations BEFORE you do this.
# 3) Specify the categorical parameters. As above, specify which columns are categorical predictors
# 4) is variance constant over x (nVar), or variable? I've included a few options for several eventualities. Specified as "var_type"
# 5) there are also some plotting parameters, these should work with any range of info, but can be modified!

#the script will assemble a model dataframe based on this info. The actual models are placed in 'JAGSBetaModels'
#the script selects an appropriate model with the information you've provided, BUT I haven't covered all eventualities
#for instance with a dataframe with 1 categorical predictor and 3 continuous predictors, we run "3Beta1Cat_nVar.txt"
#but I've never had to run a model with 6 continuous predictors and 2 categories, so it will throw an error and you should make the model yourself.
#in addition I've set some reasonable priors for the data I have, but ALWAYS CHECK YOUR PRIORS.

#This script has been modified specifically to examine "Number of naturalised grid cells" versus "Proportion of niche expansion"
#We want to make sure expansion isn't strongly driven by number of naturalised occurrences, and that other factors are at play
#This model is 1 continuous predictor, and 1 beta-distribution proportional response
#It saves various diagnostic plots in IntermediateOutput and summaries of the main output in FinalOutput
#Output is used in Appendix Figure S1.2. The underlying model spec can be found in Appendix Table 1.3 and in code/supplementary/JAGSBetaModels

rm(list=ls())

library(ggplot2)
library(RColorBrewer)
library(R2jags)
library(boot)
library(loo)



#set to current repo
setwd("DIRECTORY_HERE")


#give output folder
output<-"./FinalOutput/supplementary/OccExpan_BetaBayes/"

#tax group
table_name<-"Plant"

annot<-""



################################
####the following lines are used to correlate number of occurrences to filling/expansion
#they are not used for normal regression
################################
sp_info<-read.delim("./IntermediateOutput/Rstate_compile/3Var_plant_D_shiftvalues_zcor_center.txt", sep="\t")


unfilling<-sp_info[,c("species.name","naturalised_occurences","expansion")]
names(unfilling)[1:2]<-c("sp_name","nat_occ")

#set the y response
unfilling$prop_filling<-unfilling$expansion
unfilling$nat_occ<-log10(unfilling$nat_occ)

################################


#is data variance constant? Options are nVar (normal) and expVar (variance increase exponentially over x1 (first predictive parameter))
var_type<-"nVar"

summary(unfilling$prop_filling)

unfilling$prop_filling[unfilling$prop_filling==1]<-0.999
unfilling$prop_filling[unfilling$prop_filling==0]<-0.001

#list all continuous parameter columns!
unfilling$xvar1<-unfilling$nat_occ  

#list all categorical parameter columns
#unfilling$xcat1<-unfilling$region
#unfilling$xcat1<-factor(unfilling$xcat1)




###give some plotting parameters!
#What labels do you want for each parameter? This will go on tables and on plots

x_name1<-"Number of Naturalised Gridcells"

#cat_name1<-"Region"


#what is the colour palette? how many levels do you need?
#region levels are as follows (Oceanian and Madagascan are removed as they are never used)
region_list<-c("Afrotropical","Australian","Nearctic","Neotropical","Oriental","Palearctic-East",
"Palearctic-West","Panamanian","Saharo-Arabian","Sino-Japanese")


pal<-brewer.pal(length(region_list),"Spectral")
names(pal)<-region_list
colScale <- scale_colour_manual(name = "region_list",values = pal)


#pal<-brewer.pal(length(unique(unfilling_df$xcat1)),"Spectral")
#names(pal)<-unique(unfilling_df$xcat1)





####everything under here should be automated#########

#assemble arbitrary length lists
x_temp<-c(ls(pattern="x_name"))
x_list <- unlist(lapply(x_temp, get))

cat_temp<-c(ls(pattern="cat_name"))
cat_names <- unlist(lapply(cat_temp, get))

#stitch a new data frame together
cont_index<-which(grepl("xvar",colnames(unfilling)))
cat_index<-which(grepl("xcat",colnames(unfilling)))

#changing continuous variables to have mean of 0 and SD of 1
unfilling[,cont_index]<-scale(unfilling[,cont_index], center=T, scale=T)[,1:length(cont_index)]


#create new dataframe and trim to complete cases only
unfilling_df<-data.frame(unfilling$sp_name,unfilling$prop_filling,unfilling[,cont_index],unfilling[,cat_index])


colnames(unfilling_df)<-colnames(unfilling)[c(1,6,cont_index,cat_index)]
colnames(unfilling_df)[c(1,2)]<-c("sp_name","prop_filling")


unfilling_df<-unfilling_df[complete.cases(unfilling_df),]


#in case levels have been removed then relevel the cat
#unfilling_df$xcat1<-droplevels(unfilling_df$xcat1)
#unfilling_df$xcat1<-factor(unfilling_df$xcat1)


#if any factor levels have been removed, script will break, remove empty factor levels
new_catref<-which(grepl("xcat",colnames(unfilling_df)))
if(length(new_catref)>1){
  unfilling_df[,new_catref]<-lapply(unfilling_df[,new_catref],factor)
}else{
  unfilling_df[,new_catref]<-factor(unfilling_df[,new_catref])
}

#x11()
par(mfrow=c(2,3))
for(i in 1:ncol(unfilling_df)){
  if(is.numeric(unfilling_df[,i])){
    hist(unfilling_df[,i],breaks=20,main=colnames(unfilling_df)[i])
  }
}



#########select model############
#unfilling_df<-unfilling_df[,colnames(unfilling_df)!="xvar4"]

#how many categorical parameters are there?
cat_n<-sum(grepl("xcat",colnames(unfilling_df)))
#how many continuous predictive parameters are there?
cont_n<-sum(grepl("xvar",colnames(unfilling_df)))

model_file<-paste("./code/supplementary/JAGSBetaModels/",cont_n,"Beta",cat_n,"Cat_",var_type,".txt",sep="")
model_code = readChar(model_file, file.info(model_file)$size)

model_file

#####REMEMBER TO CHECK YOUR PRIORS!
#in this case we can get away with it as all parameters have similar null models (slope is 0, intercept is inv.logit(-1))
#phi prior (dispersion) is set to uniform distibution 0-10

##create parameter list for model. This is done in several stages because of the weird way R deals with lists
baselist<-list(N = nrow(unfilling_df), 
           y = unfilling_df$prop_filling)

#get all continuous variables into list
cont_col<-unfilling_df[,which(grepl("xvar",colnames(unfilling_df)))]
if(class(cont_col)!="data.frame"){cont_col<-list(cont_col)}
names(cont_col)<-paste("x",1:length(cont_col),sep="")


#get all categorical variables
cat_col<-unfilling_df[,which(grepl("xcat",colnames(unfilling_df)))]
if(class(cat_col)!="data.frame"){cat_col<-list(cat_col)}
if(length(cat_col)>0){names(cat_col)<-paste("cat",1:length(cat_col),sep="")}


#specify number of levels of categories
cat_s<-lapply(cat_col,function(x){length(unique(x))})
if(length(cat_col)>0){names(cat_s)<-paste("N_cat",1:length(cat_col),sep="")}

#define phi variables
if (var_type=="nVar"){phi_param<-"phi"}
if (var_type=="expVar"){phi_param<-c("phi_alpha","phi_beta")}

model_data<-c(baselist,cont_col,cat_col,cat_s)


# Choose the parameters to watch
model_parameters =  c(model_parameters =  c("alpha",
                                             paste("beta",1:cont_n,sep=""),
                                            phi_param,
                                            "mu",
                                            "mu_alpha",
                                            paste("mu_beta",1:cont_n,sep=""),
                                            "sigma_alpha",
                                            paste("sigma_beta",1:cont_n,sep=""),
                                            "y",
                                            "y_pred","LogLik", "Fit","FitNew"))
model_parameters

#set an exception for if cont_n=1 and cat_n=0
if (cont_n==1 & cat_n==0){
  
  model_parameters =  c(model_parameters =  c("alpha",
                                              "beta1",
                                              phi_param,
                                              "mu",
                                              "y",
                                              "y_pred","LogLik", "Fit","FitNew"))
  
}




# Run the model
model_run = jags(data = model_data, n.iter=20000,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code))


print(model_run)


rhat_check<-model_run$BUGSoutput$summary[,"Rhat"]

rm1<-which(grepl("mu", names(rhat_check)))
rm2<-which(grepl("y_pred", names(rhat_check)))
rm3<-which(grepl("y", names(rhat_check)))
rm4<-which(grepl("LogLik", names(rhat_check)))
rm5<-which(grepl("Fit", names(rhat_check)))
rm6<-which(grepl("FitNew", names(rhat_check)))
rm_all<-c(rm1,rm2,rm3,rm4,rm5,rm6)

rhat_check<-rhat_check[-rm_all]



#check for convergence
conver<-"fine"
if(max(rhat_check)>=1.1){print("WARNING: CONVERGENCE PROBLEMS");conver<-"PROBLEMS!"}


#get some summary statistics
DIC_res<-model_run$BUGSoutput$DIC
pd_res<-model_run$BUGSoutput$pD

loglik_full <- model_run$BUGSoutput$sims.list$LogLik
waic_t <- waic(loglik_full)$estimates["waic",]

#can also calculate loo results but not currently used
loo_t<-loo(loglik_full)
#plot(loo_t)


#check 


#0 basically means total lack of fit, how well does our beta distribution approximate the observed data?
#ideally it should be about 0.5, since residuals should be random around the beta distribution (above and below equally)


####PRODUCE some diagnostic plots
  #for each value of x in range we predicted a y value with our equation


#how well did it fit with out actual data? Take the mean estimate for all x range values and plot
#against actual y values. It fits worse and worse as phi increases


#calculate pearson residuals
mu<-model_run$BUGSoutput$sims.list$mu


linear_mu<-logit(mu)
Resid <- -1*sweep(mu,2,unfilling_df$prop_filling,'-')/sqrt(mu*(1-mu))
#store our linear predictors and our link transformed response
cor_df<-data.frame(prop=logit(unfilling_df$prop_filling),linear_mu=apply(linear_mu,2,mean))


#can calculate a pseudo r-squared (squared correlation of linear predictor and link-transformed response)
#adapted from betareg package
rsq <- function (x, y) cor(x, y) ^ 2
cor_test<-cor.test(cor_df$prop,cor_df$linear_mu, method="pearson")
pseudoR<-rsq(cor_df$prop,cor_df$linear_mu)


summary_tab1<-model_run$BUGSoutput$summary

#check fit with residuals
out <- model_run$BUGSoutput
sumsq_compare<-mean(out$sims.list$FitNew > out$sims.list$Fit)





#also do a post posterior check, let's see how accurate our predicted values were
y_pred<-model_run$BUGSoutput$sims.list$y_pred
y<-unfilling_df$prop_filling
y_pred_quant = apply(y_pred, 2, quantile, probs = c(0.25, 0.5, 0.75))



res_dif = (logit(y) - logit(y_pred_quant))^2
rms<-sqrt(mean(res_dif,na.rm=T))



pdf(file=paste("./IntermediateOutput/OccExpan_BetaBayes/Diagnostics/Diagnostics_",table_name,".pdf",sep=""),width=9,height=9)
par(mfrow=c(2,2))


#plot residuals
plot(apply(Resid,2,mean)~apply(linear_mu,2,mean),main="residuals of regression")
abline(a=0,b=0,col="red")

#plot linear predictor versus link tranformed response (used for pearson's residuals)
plot(apply(linear_mu,2,mean),logit(unfilling_df$prop_filling))
abline(lm(logit(unfilling_df$prop_filling)~apply(linear_mu,2,mean)))
legend("topleft",paste("pseudo R-squared:", round(pseudoR,digits=3)))

# Create a plot of the true values against the posterior predicted values
plot((y_pred_quant[2,]), (y), pch = 19,ylab="True y",xlab="predicted y",ylim=c(0,1),cex=0.6)
for(i in 1:ncol(y_pred_quant)) {
  lines(c(y_pred_quant[1,i], y_pred_quant[3,i]), c(y[i], y[i]))
}
abline(a = 0, b = 1,lty=2,col="red")

#plot fitted over true y
#plot(logit(y_pred_quant[2,]), logit(y), pch = 19,ylab="logit(True y)",xlab="logit(predicted y)",cex=0.6)
plot(out$sims.list$FitNew,out$sims.list$Fit,xlab="Fitted",ylab="True Values")
mtext(paste("RMS:",round(rms,2)),3)
#for(i in 1:ncol(y_pred_quant)) {
  #lines(c(logit(y[i]), logit(y[i])), c(logit(y_pred_quant[1,i]), logit(y_pred_quant[3,i])))
#}
abline(a = 0, b = 1)


dev.off()



reg_names<-unique(unfilling_df$xcat1)
pars = model_run$BUGSoutput$sims.list

diag_par<-names(pars)
diag_par<-diag_par[diag_par!="y_pred"]

#create some histograms of parameters to look at, do they deviate from 0?
pdf(file=paste("./IntermediateOutput/OccExpan_BetaBayes/Diagnostics/parameter_hists_",table_name,".pdf",sep=""),width=9,height=6)

for (i in 4:(cont_n+4)){
  cur_name<-diag_par[i]
  cur_par<-pars[[i]]
  print(cur_name)
  par(mfrow=c(3,3))
  for(k in 1:length(reg_names)) {
    if(cat_n==1){hist(cur_par[,k],breaks = 30, main = reg_names[k], xlim = range(cur_par), xlab = cur_name)  }
    if(cat_n==2){hist(cur_par[,k,],breaks = 30, main = reg_names[k], xlim = range(cur_par), xlab = cur_name)  }
  }
}

for(i in (cont_n+5):length(diag_par)){
  cur_name<-diag_par[i]
  cur_par<-pars[[i]]
  print(cur_name)
  par(mfrow=c(1,1))
  hist(cur_par,breaks = 30, main = cur_name, xlim = range(cur_par), xlab = cur_name)  
}

dev.off()

#check our priors haven't screwed up!
#values here should match priors in the model
#wher did priors end up?
pdf(file=paste("./IntermediateOutput/OccExpan_BetaBayes/Diagnostics/check_priors_",table_name,".pdf",sep=""),width=12,height=12)

par(mfrow=c(4,4))
dens = density(pars$alpha)
curve(dnorm(x, mean = -1, sd = 2), -5,5, xlab = 'alpha', ylab = '', ylim = range(dens$y))
lines(dens, col='red')
dens = density(pars$beta1)
curve(dnorm(x, mean = 0, sd = 1), -5, 5, xlab = 'beta1', ylab = '', ylim = range(dens$y))
lines(dens, col='red')
# dens = density(pars$beta2)
# curve(dnorm(x, mean = 0, sd = 1), -5, 5, xlab = 'beta2', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$beta3)
# curve(dnorm(x, mean = 0, sd = 1), -5, 5, xlab = 'beta3', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$beta4)
# curve(dnorm(x, mean = 0, sd = 1), -5, 5, xlab = 'beta4', ylab = '', ylim = range(dens$y))
#lines(dens, col='red')
dens = density(pars$phi)
curve(dunif(x, 0, 10), -10, 20, xlab = 'phi', ylab = '', ylim = range(dens$y))
lines(dens, col='red')
# dens = density(pars$phi_alpha)
# curve(dunif(x, 0, 10), -10, 20, xlab = 'phi_alpha', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$phi_beta)
# curve(dnorm(x, mean = 0, sd = 2), -10, 20, xlab = 'phi_beta', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')


# library("extraDistr")
# #used to be dunif but half cauchy is more conservative (tries to minimise differences)
# dens = density(pars$mu_alpha)
# curve(dnorm(x, mean=0,sd=5), -5, 5, xlab = 'mu_alpha', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$mu_beta1)
# curve(dnorm(x, mean=0,sd=5), -5, 5, xlab = 'mu_beta1', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$mu_beta2)
# curve(dnorm(x, mean=0,sd=5), -5, 5, xlab = 'mu_beta2', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$mu_beta3)
# curve(dnorm(x, mean=0,sd=5), -5, 5, xlab = 'mu_beta3', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$mu_beta4)
# curve(dnorm(x, mean=0,sd=5), -5, 5, xlab = 'mu_beta4', ylab = '', ylim = range(dens$y))
# lines(dens, col='red')
# dens = density(pars$sigma_beta1)
# curve(dnorm(x, mean=0,sd=5), 0, 5, xlab = 'sigma_beta', ylab = '', ylim = range(dens$y),xlim=c(0,10))
# lines(dens, col='red')
# dens = density(pars$sigma_beta2)
# curve(dnorm(x, mean=0,sd=5), 0, 5, xlab = 'sigma_beta2', ylab = '', ylim = range(dens$y),xlim=c(0,10))
# lines(dens, col='red')
# dens = density(pars$sigma_beta3)
# curve(dnorm(x, mean=0,sd=5), 0, 5, xlab = 'sigma_beta3', ylab = '', ylim = range(dens$y),xlim=c(0,10))
# lines(dens, col='red')
# dens = density(pars$sigma_beta4)
# curve(dnorm(x, mean=0,sd=5), 0, 5, xlab = 'sigma_beta3', ylab = '', ylim = range(dens$y),xlim=c(0,10))
# lines(dens, col='red')

dev.off()

source("./code/functions/bayesPlotting_simpleversion.R")


cat("Convergence:",conver,"\n",
    "DIC:", DIC_res,"\n",
    "pD",pd_res,"\n",
    "WAIC",waic_t[1],"\n",
    "pseudo R-Squared",pseudoR,"\n",
    "Residual Fit",sumsq_compare,"\n",
    "RMS: ", rms, "\n",
    "sample size", nrow(unfilling_df),
    file=paste(output,table_name,annot,"summary_stats.csv",sep=""))


#double check for convergence problems!
cat("Convergence:",conver,"\n",
    "DIC:", DIC_res,"\n",
    "pD",pd_res,"\n",
    "WAIC",waic_t[1],"\n",
    "pseudo R-Squared",pseudoR,"\n",
    "Residual Fit",sumsq_compare,"\n",
    "RMS: ", rms, "\n",
    "sample size", nrow(unfilling_df))

loo_t

summary(loo_t$diagnostics$pareto_k)
length(loo_t$diagnostics$pareto_k)

#anything over 0.7 is not a good fit and should be a concern.
which(loo_t$diagnostics$pareto_k>0.7)
unfilling_df[which(loo_t$diagnostics$pareto_k>0.7),]



