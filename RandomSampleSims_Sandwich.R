#Simulation code for "Practical  considerations for sandwich variance estimation in two-stage regression settings"
#Random sample - Logistic Regression

library(MASS)
library(mvtnorm)
library(survival)
library(sandwich)
library(xtable)
library(dplyr)
library(numDeriv)
library(survey)


#Functions for this program:
#Function for assessing performance (mean)
Data_output_singlecov_mean <- function(betaest,betase,betatrue){
  
  ASE<-mean(betase)
  ESE<-sd(betaest)
  #Median percent bias corrected = (estimated-target)/target
  bias<-(betaest-betatrue)/betatrue
  mean_percent_bias<-mean(bias)*100
  
  #Coverage Probability Calculation
  beta_1_l95<-betaest-qnorm(0.975)*(betase)
  beta_1_r95<-betaest+qnorm(0.975)*(betase)
  
  beta_1_ci<-as.data.frame(cbind(beta_1_l95,beta_1_r95))
  
  # check the proportion of intervals containing the parameter
  CP1<-mean(apply(beta_1_ci, 1, findInterval, x = betatrue) == 1)
  
  results<-cbind(mean_percent_bias,ASE,ESE,CP1)
  
  return(results)
}

#Function for assessing performance (median)
Data_output_singlecov_median <- function(betaest,betase,betatrue){
  
  ASE<-median(betase)
  ESE<-sd(betaest)
  MAD<-mad(betaest)
  
  #Median percent bias corrected = (estimated-target)/target
  bias<-(betaest-betatrue)/betatrue
  median_percent_bias<-median(bias)*100
  
  #Coverage Probability Calculation
  beta_1_l95<-betaest-qnorm(0.975)*(betase)
  beta_1_r95<-betaest+qnorm(0.975)*(betase)
  
  beta_1_ci<-as.data.frame(cbind(beta_1_l95,beta_1_r95))
  
  # check the proportion of intervals containing the parameter
  CP1<-mean(apply(beta_1_ci, 1, findInterval, x = betatrue) == 1)
  
  results<-cbind(median_percent_bias,ASE,MAD,CP1)
  
  return(results)
}

#Introduce expit function
expit <-function(xb){return(exp(xb)/(1+exp(xb)))}

add <- function(x) Reduce("+", x)

#Functions for A and B matrices for stacked estimating equations approach to sandwich
#Derivatives are calculated exactly for regression calibration and logistic regression outcome model
A_mat_RC_Logistic<-function(Xstar,Z,XHat,Beta_0,Beta_X,Beta_Z,Y,V){
  p<-expit(Beta_0+Beta_X*XHat+Beta_Z*Z)
  A11<- -V
  A21<- A12<- -Xstar*V
  A31<- A13<- -Z*V
  A41<- -Beta_X*p*(1-p)
  A51<- XHat*(Beta_X*p^2-Beta_X*p)+Y-p
  A61<- -Beta_X*Z*p*(1-p)
  A22<- -(Xstar^2)*V
  A32<- A23<- -Xstar*Z*V
  A42<- -Beta_X*Xstar*p*(1-p)
  A52<-  XHat*(Beta_X*Xstar*p^2-Beta_X*Xstar*p)+Xstar*(Y-p)
  A62<- -Beta_X*Z*Xstar*p*(1-p)
  A33<- -(Z^2)*V
  A43<- -Beta_X*Z*p*(1-p)
  A53<-  XHat*(Beta_X*Z*p^2-Beta_X*Z*p)+Z*(Y-p)
  A63<- -Beta_X*(Z^2)*p*(1-p)
  A14<-A24<-A34<-A15<-A25<-A35<-A16<-A26<-A36<- 0
  A44<- -p*(1-p)
  A54<- A45<- -XHat*p*(1-p)
  A64<-A46<- -Z*p*(1-p)
  A55<- -(XHat^2)*p*(1-p)
  A65<- A56 <- -XHat*Z*p*(1-p)
  A66<- -(Z^2)*p*(1-p)
  
  mymat<-matrix(c(A11,A21,A31,A41,A51,A61,
                  A12,A22,A32,A42,A52,A62,
                  A13,A23,A33,A43,A53,A63,
                  A14,A24,A34,A44,A54,A64,
                  A15,A25,A35,A45,A55,A65,
                  A16,A26,A36,A46,A56,A66), nrow=6, ncol=6)
  
  return(mymat)
}

B_mat_RC_Logistic<-function(Xstar,Z,XHat,Beta_0,Beta_X,Beta_Z,Y,V,alpha_0,alpha_X,alpha_Z,Xstarstar){
  p<-expit(Beta_0+Beta_X*XHat+Beta_Z*Z)
  eq1<-V*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq2<-V*Xstar*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq3<-V*Z*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq4<-Y-p
  eq5<-(Y-p)*XHat
  eq6<-(Y-p)*Z
  
  mymat<-matrix(c(eq1,eq2,eq3,eq4,eq5,eq6), nrow=6, ncol=1)
  
  mymatT<-mymat%*%t(mymat)
  return(mymatT)
  
}

RC_Logistic_Est_Eq<-function(parametervec){ #Xstar,Z,XHat,Beta_0,Beta_X,Beta_Z,Y,V,alpha_0,alpha_X,alpha_Z,Xstarstar
  alpha_0<-parametervec[1]
  alpha_X<-parametervec[2]
  alpha_Z<-parametervec[3]
  Beta_0<-parametervec[4]
  Beta_X<-parametervec[5]
  Beta_Z<-parametervec[6]
  
  XHat<-alpha_0+alpha_X*Xstar + alpha_Z*Z
  
  p<-expit(Beta_0+Beta_X*XHat+Beta_Z*Z)
  eq1<-V*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq2<-V*Xstar*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq3<-V*Z*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq4<-Y-p
  eq5<-(Y-p)*XHat
  eq6<-(Y-p)*Z
  
  #mymat<-matrix(c(eq1,eq2,eq3,eq4,eq5,eq6), nrow=6, ncol=1)
  mymat<-c(eq1,eq2,eq3,eq4,eq5,eq6)
  
  return(mymat)
  
}

######################################################################
################Begin data generating mechanism#######################
######################################################################

#set seed for reproducibility
set.seed(2903)

#Main study sample size 
N<-1000

#subset sample size 
n_sub<-450

#Set a reasonable value for our true beta_1, true beta_2, and true beta_3
beta_0<-log(0.2) #-1.609
beta_1<-log(1.5) #0.405
beta_2<-log(1.2) #0.182
nbeta<-3

#Set number of simulations AND number of bootstrap replicates
NSIM<-1000
NBOOT<-500

#Create empty object to store simulation data
myresults<-NULL

#Create a vector of means for simulated covariates 
mu <- rep(0,nbeta-1)

correlation<-0.3
ME_var<-0.5
#Create the covariance matrix for the simulated data
Sigma1<-matrix(c(1,correlation,correlation,1), nrow=nbeta-1, ncol=nbeta-1)


for(iter in 1:NSIM){
  
  #Generate x's from a mvrnorm distribution
  x_data <- (mvrnorm(n=N, mu=mu, Sigma=Sigma1))
  y <- ifelse(runif(N)<expit(beta_0+beta_1*x_data[,1]+beta_2*x_data[,2]),1,0)
  prop.table(table(y))
  mydata<-data.frame(x_data,y)
  colnames(mydata)<-c("x","z","y")
  mydata$ID<-1:nrow(mydata)
  
  ### Measurement error model
  #Create error-prone x subject to both systematic + random error as well as x with just classical error
  mydata$xstar<- .37*mydata$x+.15*mydata$z + rnorm(N,0,sqrt(ME_var))  #R-squared about .4
  
  #Create measure with classical error
  mydata$xstarstar<- mydata$x + rnorm(N,0,sqrt(.2))
  
  #Subset the data
  subset_index<-sample(mydata$ID,n_sub,replace=FALSE)
  mydata$valid_indicator<-as.numeric(mydata$ID%in%subset_index)
  mydata$xstarstar<-ifelse(mydata$valid_indicator==0, NA, mydata$xstarstar) 

  #Random sample design
  sampdesign <- svydesign(id=~1, data=mydata)
  sampdesign_sub<-subset(sampdesign,valid_indicator==1)
  cfit<-svyglm(xstarstar~xstar+z,design=sampdesign,family=gaussian(),subset=valid_indicator==1)

  #update samp design
  sampdesign <- update(sampdesign,xhat =predict(cfit,newdata=sampdesign$variables) )
  #sampdesign$variables
  rcfit<-svyglm(y ~ xhat+z, design=sampdesign,family=quasibinomial)
  tfit<-svyglm(y ~ x+z, design=sampdesign,family=quasibinomial)
  nfit<-svyglm(y ~ xstar+z, design=sampdesign,family=quasibinomial)

  #Now obtain A and B matrices for sandwich
  Amat_All<-Bmat_All<-NULL
  for(j in c(1:N)) {
    Amat_All[[j]]<-A_mat_RC_Logistic(Xstar=sampdesign$variables$xstar[j],Z=sampdesign$variables$z[j],XHat=sampdesign$variables$xhat[j],Beta_0=rcfit$coefficients[1],
                                     Beta_X=rcfit$coefficients[2],Beta_Z=rcfit$coefficients[3],
                                     Y=sampdesign$variables$y[j],V=sampdesign$variables$valid_indicator[j])
    Bmat_All[[j]]<-B_mat_RC_Logistic(Xstar=sampdesign$variables$xstar[j],
                                     Z=sampdesign$variables$z[j],
                                     XHat=sampdesign$variables$xhat[j],
                                     Beta_0=rcfit$coefficients[1],
                                     Beta_X=rcfit$coefficients[2],
                                     Beta_Z=rcfit$coefficients[3],
                                     Y=sampdesign$variables$y[j],
                                     V=sampdesign$variables$valid_indicator[j],
                                     alpha_0=cfit$coefficients[1],
                                     alpha_X=cfit$coefficients[2],
                                     alpha_Z=cfit$coefficients[3],
                                     Xstarstar=ifelse(is.na(sampdesign$variables$xstarstar[j]),
                                                      0,sampdesign$variables$xstarstar[j]))
    
  }
  length(Bmat_All)
  sampdesign$variables$xstarstar[j]
  mydata$xstarstar[j]
  
  summary(sampdesign$variables$xstarstar)
  summary(mydata$xstarstar)
  
  
  AmatSUM<-add(Amat_All)/N
  BmatSUM<-add(Bmat_All)/N
  A<-solve(AmatSUM)
  B<-BmatSUM
  sandwichvar<-(A)%*%B%*%(t((A)))/N

  #Bootstrap instead!
  bans<-NULL
  for(b in 1:NBOOT){
    ### Need to bootstrap stratified on membership to calibration subset
    bID1<-sample(1:n_sub,n_sub,replace=T)
    bID2<-sample((n_sub+1):N,(N-n_sub),replace=T)
    bID<-c(bID1,bID2)
    bdata<-mydata[bID,]
    bcfit<-lm(xstarstar~xstar+z,data=bdata)
    bdata$xhat<-predict(bcfit,newdata=bdata)
    brcfit<-glm(y ~ xhat+z ,data=bdata,family=quasibinomial())
    bans<-rbind(bans,coef(brcfit))
  }
  
  #Save
  myresults<-rbind(myresults,c(summary(cfit)$coefficients[,1],summary(cfit)$coefficients[,2],
                               summary(tfit)$coefficients[,1],summary(tfit)$coefficients[,2],
                               summary(nfit)$coefficients[,1],summary(nfit)$coefficients[,2],
                               summary(rcfit)$coefficients[,1],summary(rcfit)$coefficients[,2],
                               sd(bans[,1]),sd(bans[,2]),sd(bans[,3]),
                               sqrt(diag(sandwichvar))))
  
  
  print(iter)
} ### End of simulation loop
myresults<-data.frame(myresults)
names(myresults)<-c("a0","aX","aZ","se_a0","se_aX","se_aZ","b0_True","bX_True","bZ_True","se_b0_True","se_bX_True","se_bZ_True",
                    "b0_Naive","bX_Naive","bZ_Naive","se_b0_Naive","se_bX_Naive","se_bZ_Naive",
                    "b0_RC","bX_RC","bZ_RC","se_b0_RC","se_bX_RC","se_bZ_RC","se_b0_boot","se_bX_boot","se_bZ_boot",
                    "se_a0_sand","se_aX_sand","se_aZ_sand","se_b0_sand","se_bX_sand","se_bZ_sand")

head(myresults)

#### Simulation results
mytable<-rbind(Data_output_singlecov_median(myresults$bX_True,myresults$se_bX_True,beta_1),
               Data_output_singlecov_median(myresults$bX_Naive,myresults$se_bX_Naive,beta_1),
               Data_output_singlecov_median(myresults$bX_RC,myresults$se_bX_RC,beta_1),
               Data_output_singlecov_median(myresults$bX_RC,myresults$se_bX_sand,beta_1),
               Data_output_singlecov_median(myresults$bX_RC,myresults$se_bX_boot,beta_1))
xtable(mytable,digits=3)


