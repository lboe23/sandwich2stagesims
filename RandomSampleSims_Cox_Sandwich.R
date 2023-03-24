#Simulation code for "Practical  considerations for sandwich variance estimation in two-stage regression settings"
#Random sample - Cox Proportional Hazards Model Regression

library(MASS)
library(mvtnorm)
library(survival)
library(sandwich)
library(survey)
library(xtable)
library(dplyr)
library(numDeriv)


#Functions for this program:
My.Cox.Residuals<-function(alphas){
  #code
  type<-"score"
  
  otype <- type
  n <- length(PH_rcfit_svy$residuals)
  rr <- PH_rcfit_svy$residuals
  y <- PH_rcfit_svy$y
  x <- PH_rcfit_svy[['x']]  # avoid matching PH_rcfit_svy$xlevels
  x2 <- PH_rcfit_svy[['x']]  # avoid matching PH_rcfit_svy$xlevels
  
  vv <- drop(PH_rcfit_svy$naive.var)
  
  if (is.null(vv)) vv <- drop(PH_rcfit_svy$var)
  weights <- PH_rcfit_svy$weights
  if (is.null(weights)) weights <- rep(1,n)
  strat <- PH_rcfit_svy$strata
  method <- PH_rcfit_svy$method
  
  # I need Y, and perhaps the X matrix (and strata)
  Terms <- PH_rcfit_svy$terms
  if (!inherits(Terms, 'terms'))
    stop("invalid terms component of object")
  strats <- attr(Terms, "specials")$strata
  
  
  
  
  if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
    temp <- survival:::coxph.getdata(PH_rcfit_svy, y=TRUE, x=TRUE, stratax=TRUE)
    y <- temp$y
    x <- temp$x
    
    myvarscox<-all.vars(formula(PH_rcfit_svy)[[3]])[-1]
    zvars<-x[,(myvarscox)]
    
    svydata<-PH_rcfit_svy$survey.design$variables
    xmatstage1<-as.matrix(cbind(rep(1,N),svydata[,all.vars(formula(cfit)[[3]])]))
    dim(xmatstage1)
    myxhat<-as.vector(xmatstage1%*%alphas)
    xmatstage2 <- as.matrix(cbind(myxhat,zvars))
    
    if (length(strats)) strat <- temp$strata
  }
  ny <- ncol(y)
  status <- y[,ny,drop=TRUE]
  
  
  
  nstrat <- as.numeric(strat)
  nvar <- ncol(xmatstage2)
  
  if (is.null(strat)) {
    ord <- order(y[,ny-1], -status)
    newstrat <- rep(0,n)
  }
  newstrat[n] <- 1
  
  # sort the data
  x <- x[ord,]
  xmatstage2 <- xmatstage2[ord,]
  
  y <- y[ord,]
  linpred<-exp(xmatstage2%*%coef(PH_rcfit_svy))
  mymeans<-apply(xmatstage2,2,mean)
  expmean<-as.vector(exp(mymeans%*%coef(PH_rcfit_svy)))
  score<-as.vector(linpred/expmean)
  
  
  
  #Score
  storage.mode(y) <- storage.mode(x) <- storage.mode(xmatstage2) <- "double"
  storage.mode(newstrat) <- "integer"
  storage.mode(score) <- storage.mode(weights) <- "double"
  #if (ny==2) {
  resid <-  .Call(survival:::Ccoxscore2, 
                  y, 
                  xmatstage2, 
                  newstrat,
                  score,
                  weights[ord],
                  as.integer(method=='efron'))
  
  summary(resid)
  #} else {
  #  resid<- .Call(survival:::Cagscore2,
  #                y, 
  #                x, 
  #                newstrat,
  #                score,
  #                weights[ord],
  #                as.integer(method=='efron'))
  #}
  
  if (nvar >1) {
    rr <- matrix(0, n, nvar)
    rr[ord,] <- resid
    dimnames(rr) <- list(names(PH_rcfit_svy$residuals), 
                         names(PH_rcfit_svy$coefficients))
  }else {
    rr[ord] <- resid
  }
  
  
  # rtest<-resid[50,]
  return(rr)
}

My.Cox.Residuals.1person<-function(alphas){
  #code
  type<-"score"
  
  otype <- type
  n <- length(PH_rcfit_svy$residuals)
  rr <- PH_rcfit_svy$residuals
  y <- PH_rcfit_svy$y
  x <- PH_rcfit_svy[['x']]  # avoid matching PH_rcfit_svy$xlevels
  x2 <- PH_rcfit_svy[['x']]  # avoid matching PH_rcfit_svy$xlevels
  
  vv <- drop(PH_rcfit_svy$naive.var)
  
  if (is.null(vv)) vv <- drop(PH_rcfit_svy$var)
  weights <- PH_rcfit_svy$weights
  if (is.null(weights)) weights <- rep(1,n)
  strat <- PH_rcfit_svy$strata
  method <- PH_rcfit_svy$method
  
  # I need Y, and perhaps the X matrix (and strata)
  Terms <- PH_rcfit_svy$terms
  if (!inherits(Terms, 'terms'))
    stop("invalid terms component of object")
  strats <- attr(Terms, "specials")$strata
  
  
  
  
  if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
    temp <- survival:::coxph.getdata(PH_rcfit_svy, y=TRUE, x=TRUE, stratax=TRUE)
    y <- temp$y
    x <- temp$x
    
    myvarscox<-all.vars(formula(PH_rcfit_svy)[[3]])[-1]
    zvars<-x[,(myvarscox)]
    
    svydata<-PH_rcfit_svy$survey.design$variables
    xmatstage1<-as.matrix(cbind(rep(1,N),svydata[,all.vars(formula(cfit)[[3]])]))
    dim(xmatstage1)
    myxhat<-as.vector(xmatstage1%*%alphas)
    xmatstage2 <- as.matrix(cbind(myxhat,zvars))
    
    if (length(strats)) strat <- temp$strata
  }
  ny <- ncol(y)
  status <- y[,ny,drop=TRUE]
  
  
  
  nstrat <- as.numeric(strat)
  nvar <- ncol(xmatstage2)
  
  if (is.null(strat)) {
    ord <- order(y[,ny-1], -status)
    newstrat <- rep(0,n)
  }
  newstrat[n] <- 1
  
  # sort the data
  x <- x[ord,]
  xmatstage2 <- xmatstage2[ord,]
  
  y <- y[ord,]
  linpred<-exp(xmatstage2%*%coef(PH_rcfit_svy))
  mymeans<-apply(xmatstage2,2,mean)
  expmean<-as.vector(exp(mymeans%*%coef(PH_rcfit_svy)))
  score<-as.vector(linpred/expmean)
  
  
  
  #Score
  storage.mode(y) <- storage.mode(x) <- storage.mode(xmatstage2) <- "double"
  storage.mode(newstrat) <- "integer"
  storage.mode(score) <- storage.mode(weights) <- "double"
  #if (ny==2) {
  resid <-  .Call(survival:::Ccoxscore2, 
                  y, 
                  xmatstage2, 
                  newstrat,
                  score,
                  weights[ord],
                  as.integer(method=='efron'))
  
  summary(resid)
  #} else {
  #  resid<- .Call(survival:::Cagscore2,
  #                y, 
  #                x, 
  #                newstrat,
  #                score,
  #                weights[ord],
  #                as.integer(method=='efron'))
  #}
  
  if (nvar >1) {
    rr <- matrix(0, n, nvar)
    rr[ord,] <- resid
    dimnames(rr) <- list(names(PH_rcfit_svy$residuals), 
                         names(PH_rcfit_svy$coefficients))
  }else {
    rr[ord] <- resid
  }
  
  
  # rtest<-resid[50,]
  return(rr[2678,])
}

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


######################################################################
################Begin data generating mechanism#######################
######################################################################

#set seed for reproducibility
set.seed(2903)

#Main study sample size - 16000 like HCHS/SOL
N<-1000

#subset sample size - 450 like HCHS/SOL
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

  #Create Cox outcome
  #Begin time to event: choose beta carotene for highest R2
  lambda_sim <- 0.10 * exp(c(beta_1*x_data[,1]+beta_2*x_data[,2]))
  Event_time <- rexp(N, lambda_sim)
  
  #censoring time
  finaltime<-2
  
  #Create times censored at final time and indicator
  timetoevent<-ifelse(Event_time<finaltime,Event_time,finaltime)
  event_indicator<-ifelse(timetoevent<finaltime,1,0)
  eventrate<-prop.table((table(event_indicator)))
  
  #Finalize data frame
  mydata<-data.frame(x_data,y,timetoevent,event_indicator)
  
  colnames(mydata)<-c("x","z","y","Time","delta")
  mydata$ID<-1:nrow(mydata)
  
  ### Measurement error model
  #Create error-prone x subject to both systematic + random error as well as x with just classical error
  mydata$xstar<- .37*mydata$x+.15*mydata$z + rnorm(N,0,sqrt(ME_var))  #R-squared about .4
  
  #Create measure with classical error
  mydata$xstarstar<- mydata$x + rnorm(N,0,sqrt(.2))
  
  #Now, since value with classical ME can only be available on a subset due to cost, make available for subset only
  mydata$xstarstar<-ifelse(mydata$ID > n_sub, NA, mydata$xstarstar) 
  mydata$valid_indicator<-ifelse(is.na(mydata$xstarstar),0,1)
  
  #before fitting any models, make sure we're ordered with all missing gold standard last 
  mydata<-arrange(mydata, desc(is.na(xstarstar))) 
  #Create survey design
  sampdesign <- svydesign(id=~1, data=mydata)
  sampdesign_sub<-subset(sampdesign,valid_indicator==1)
  cfit<-svyglm(xstarstar~xstar+z,design=sampdesign,family=gaussian(),subset=valid_indicator==1)

  mydata$xhat<-predict(cfit,newdata=mydata)
  
  #update samp design
  sampdesign <- update(sampdesign,xhat =predict(cfit,newdata=sampdesign$variables) )

  
  #Now let's run Cox time-to-event models
  PH_rcfit_svy<-svycoxph(Surv(Time, delta) ~ xhat+z, design=sampdesign)
  PH_tfit_svy<-svycoxph(Surv(Time, delta) ~ x+z, design=sampdesign)
  PH_nfit_svy<-svycoxph(Surv(Time, delta) ~ xstar+z, design=sampdesign)
  
  #Everything below here is the stuff we want to include in our sandwich variance function
  
  #j unknown alphas, k unknown betas
  j_dim<-length(coef(cfit))
  k_dim<-length(coef(PH_rcfit_svy))
  
  #Now we can begin constructing meat of sandwich -- use existing functions
  #Create the estfun for the stage 1 model by first constructing a N by j_dim matrix of zeroes
  #  and insert estfun contributions for those who are in stage 1 model
  #  this way, all contributions for those NOT in stage 1 model are zero
  # with estimating equation contributions for those in the stage 1 model
  xstarstar<-sampdesign$variables$xstarstar
  
  estfun.stage1<-matrix(0,nrow=N,ncol=j_dim)
  is.calibration <- !is.na(xstarstar)
  estfun.stage1[is.calibration,] <- survey:::estfuns.lm(cfit)
  
  #estfun for stage 2 is straightfroward -- just pull from estfun 
  estfun.stage2<-as.matrix(estfun(PH_rcfit_svy))
  
  #Combine estfuns for both stage 1 and stage 2 models 
  estfun.all<-cbind(estfun.stage1,estfun.stage2)
  
  #Function to add matrices
  add <- function(x) Reduce("+", x)
  
  #upper right is zeroes 
  A_upperright<- matrix(0,nrow=j_dim,ncol=k_dim)
  
  #Bottom right is Hessian for stage 2 model
  A_bottomright<--solve(PH_rcfit_svy$inv.info)/N #bottom right hand

  #Upper left is Hessian for stage 1 model 
  A_upperleft<--solve(cfit$naive.cov)/N #bottom right hand

  #Save alphas
  alphas<-coef(cfit)
  
  #Create cox model jacobian
  jacobiancoxall<-jacobian(My.Cox.Residuals,coef(cfit))
  
  mat1jacobiancoxall<-jacobiancoxall[c(1:N),]
  mat2jacobiancoxall<-jacobiancoxall[c((N+1):dim(jacobiancoxall)[1]),]

  JacobianList.stage2<-NULL
  for(j in c(1:N)) {
    JacobianList.stage2[[j]]<-as.matrix(rbind(mat1jacobiancoxall[j,],mat2jacobiancoxall[j,]))
  }
  
  A_bottomleft<-add(JacobianList.stage2)/N
  
  AAll<-rbind(cbind(A_upperleft,A_upperright),cbind(A_bottomleft,A_bottomright))

  A.inv<-solve(AAll)

  #Compute influence functions
  infl<- as.matrix(estfun.all)%*%t(A.inv)/N
  
  #Compute the sandwich 2 ways for SRS -- directly in sandwich from, or using vcov(svytotal())
 # sand1<-(A.inv)%*%UUt%*%(t(A.inv))/N
  sand2<-vcov(svytotal(infl,PH_rcfit_svy$survey.design))
  
  sqrt(diag(sand2))
  
  
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
     brcfit<-coxph(Surv(Time, delta) ~ xhat+z, data=bdata)
     bans<-rbind(bans,coef(brcfit))
  }
  
  
  #Save
  
  myresults<-rbind(myresults,c(coef(PH_rcfit_svy),coef(PH_tfit_svy),coef(PH_nfit_svy),
                               sqrt(diag(vcov(PH_rcfit_svy))),sqrt(diag(vcov(PH_tfit_svy))),
                               sqrt(diag(vcov(PH_nfit_svy))),  
                               sqrt(diag(sand2)),sd(bans[,1]),sd(bans[,2])))
  
  
  print(iter)
} ### End of simulation loop



myresults<-data.frame(myresults)
names(myresults)<-c("bX_RC","bZ_RC",
                    "bX_True","bZ_True",
                    "bX_Naive","bZ_Naive",
                    "se_bX_RC","se_bZ_RC",
                    "se_bX_True","se_bZ_True",
                    "se_bX_Naive","se_bZ_Naive",
                    "se_a0_sand","se_aX_sand","se_aZ_sand","se_bX_sand","se_bZ_sand",
                    "se_bX_boot","se_bZ_boot")

head(myresults)

#### Simulation results
mytable<-rbind(Data_output_singlecov_median(myresults$bX_True,myresults$se_bX_True,beta_1),
               Data_output_singlecov_median(myresults$bX_Naive,myresults$se_bX_Naive,beta_1),
               Data_output_singlecov_median(myresults$bX_RC,myresults$se_bX_RC,beta_1),
               Data_output_singlecov_median(myresults$bX_RC,myresults$se_bX_sand,beta_1),
               Data_output_singlecov_median(myresults$bX_RC,myresults$se_bX_boot,beta_1))

mytablenew<-cbind(mytable[,-c(2,4)],mytable[,c(2,4)])
xtable(mytablenew,digits=3)





