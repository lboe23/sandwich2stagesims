#Simulation code for "Practical  considerations for sandwich variance estimation in two-stage regression settings"
#Complex Survey Data Simulation - Logistic Regression

library(sandwich)
library(survey)
library(Hmisc)
library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(car)
library(data.table)
library(microbenchmark)
library(splitstackshape)


#Function for assessing performance
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



B_mat_RC_Linear<-function(Xstar,Z,XHat,Beta_0,Beta_X,Beta_Z,Y,V,alpha_0,alpha_X,alpha_Z,Xstarstar){
  eq1<-V*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq2<-V*Xstar*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq3<-V*Z*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq4<-1*(Y-Beta_0-Beta_X-XHat-Beta_Z*Z)
  eq5<-XHat*(Y-Beta_0-Beta_X-XHat-Beta_Z*Z)
  eq6<-Z*(Y-Beta_0-Beta_X-XHat-Beta_Z*Z)
  
  mymat<-matrix(c(eq1,eq2,eq3,eq4,eq5,eq6), nrow=6, ncol=1)
  
  mymatT<-mymat%*%t(mymat)
  return(mymatT)
  
}

EstFun_RC_Logistic<-function(Xstar,Z,XHat,Beta_0,Beta_X,Beta_Z,Y,V,alpha_0,alpha_X,alpha_Z,Xstarstar){
  p<-expit(Beta_0+Beta_X*XHat+Beta_Z*Z)
  eq1<-V*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq2<-V*Xstar*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq3<-V*Z*(Xstarstar-alpha_0-alpha_X*Xstar-alpha_Z*Z)
  eq4<-Y-p
  eq5<-(Y-p)*XHat
  eq6<-(Y-p)*Z
  
  estfun<-c(eq1,eq2,eq3,eq4,eq5,eq6)
  
  
  return(estfun)
  
}

RC_Logistic_Est_Eq<-function(parametervec){ #Xstar,Z,XHat,Beta_0,Beta_X,Beta_Z,Y,V,alpha_0,alpha_X,alpha_Z,Xstarstar
  alphas<-parametervec[1:N_alpha]
  betas<-parametervec[(N_alpha+1):length(parametervec)]
  AlphaMat<-as.matrix(cbind(1,Xstar,Z))
  
  XHat<-AlphaMat%*%alphas 
  BetaMat<-as.matrix(cbind(1,XHat,Z))
  
  p<-expit(BetaMat%*%betas)
  eq1<-V*(Xstarstar-XHat)
  eq2<-V*Xstar*(Xstarstar-XHat)
  eq3<-V*Z*(Xstarstar-XHat)
  eq4<-Y-p
  eq5<-(Y-p)*XHat
  eq6<-(Y-p)*Z
  
  mymat<-unlist(c(eq1,eq2,eq3,eq4,eq5,eq6))
  
  return(mymat)
  
}

######################################################################
################Load in Complex Survey Data#######################
######################################################################

sampname = 'SampleData'
S=1000 #Number of samples drawn
idx=1:1000

#set seed for reproducibility
set.seed(2903)

beta_1<-log(1.5)


options(survey.adjust.domain.lonely=TRUE)
options(survey.lonely.psu="adjust")


#Most important step: Create path for reading in samples (simulated from "GenerateData_ComplexSurveySims.R")
mycorr<-0.3;MEvar<-0.25;meanN<-10000
popname = paste0('TargetPopulationData',MEvar,'_corr',mycorr)
sampname = 'SampleData'

files = paste0(popname,'_N',meanN,'/',sampname,"_",str_pad(1:S,nchar(S),pad=0),".RData")

#Number of betas
nbeta<-3

#Set number of simulations AND number of bootstrap replicates
NBOOT<-25 

boot_list = list()


#Create some empty objects to store simulation data
myresults<-NULL
for(i in idx){ #i=1
  boot_list[[i]] = list()
  
  load(file=files[i])
  samp <- as.data.table(samp)
  N<-nrow(samp)
  samp$ID<-c(1:N)
  prop.table(table(samp$y))
  
  #create SOLNAS indicator, V
  samp$V<-ifelse(samp$solnas==T,1,0)
  samp$xstarstar<-ifelse(samp$V==0,NA,samp$xstarstar)
  
  summary(samp$xstarstar)
  #new strata = 8
  samp$newstrat<- interaction(samp$strat, samp$V)
  
  #create sampling design
  sampdesign <- svydesign(id=~BGid, strata=~newstrat, weights=~bghhsub_s2, data=samp,nest=TRUE) #let solnas design be subset 
  sampdesignsubset<-subset(sampdesign,V==1)
  
  #create RANDOM sampling design
  sampdesignrand <- svydesign(id=~1, data=samp) #let solnas design be subset 
  
  #Fit calibration model
  cfit<-svyglm(xstarstar~xstar+X2,design=sampdesignsubset,family = stats::gaussian())

  #predict xhat and update samp design
  sampdesign <- update(sampdesign,xhat = predict(cfit,newdata=sampdesign$variables))
  N_alpha<-length(cfit$coefficients)
  
  sampdesign<-subset(sampdesign,V==1)
  
  #Now get regression coefficients for true, naive and RC model
  tfit<-svyglm(y ~ X1+X2 ,design=sampdesign,family=quasibinomial())
  nfit<-svyglm(y ~ xstar+X2 ,design=sampdesign,family=quasibinomial())
  rcfit<-svyglm(y ~ xhat+X2 ,design=sampdesign,family=quasibinomial(),influence=TRUE)

  N<-nobs(rcfit)
  #Try to get correct 
  U.stage2<-survey:::estfuns.glm(rcfit)/weights(rcfit$survey.design,"sampling")

  #U.stage2<-survey:::estfuns.glm(rcfit)/rcfit$prior.weights
  summary(rcfit$prior.weights)
  mybread<- -solve((rcfit$naive.cov)/mean(rcfit$prior.weights))/N
  mybreadinv<-solve(mybread)
  
  bread(rcfit)
  t(mybreadinv)
  
  infl.stage2<- U.stage2%*%t(mybreadinv)/N
  summary(svytotal(infl.stage2,  sampdesign))
  summary((rcfit$prior.weights))
  summary(weights(rcfit$survey.design,"sampling"))
  
  #infl<- as.matrix(estfun.all)%*%t(A.inv)/N
  
  v.stage2<-vcov(svytotal(infl.stage2,  sampdesign))
  sqrt(diag(v.stage2))
  SE(rcfit)
  su

  workingweightstest <- as.vector(((rcfit$family$variance(mu.j))^2)/rcfit$family$variance(mu.j))
  

  Ainv<-summary(rcfit)$cov.unscaled
  estfuntest<-survey:::estfuns.glm(rcfit)
  estfuntest<- residuals(rcfit, "working") * rcfit$weights * model.matrix(rcfit)
  testagain<-svyrecvar(estfuntest%*%Ainv,rcfit$survey.design$cluster,rcfit$survey.design$strata,rcfit$survey.design$fpc,postStrata=rcfit$survey.design$postStrata)
  sqrt(diag(testagain))
  
  
  
  Amat_All<-Bmat_All<-NULL
  for(j in c(1:N)) {
    Amat_All[[j]]<-A_mat_RC_Logistic(Xstar=sampdesign$variables$xstar[j],Z=sampdesign$variables$X2[j],XHat=sampdesign$variables$xhat[j],Beta_0=rcfit$coefficients[1],
                                     Beta_X=rcfit$coefficients[2],Beta_Z=rcfit$coefficients[3],
                                     Y=sampdesign$variables$y[j],V=sampdesign$variables$V[j])
    Bmat_All[[j]]<-B_mat_RC_Logistic(Xstar=sampdesign$variables$xstar[j],Z=sampdesign$variables$X2[j],XHat=sampdesign$variables$xhat[j],Beta_0=rcfit$coefficients[1],
                                     Beta_X=rcfit$coefficients[2],Beta_Z=rcfit$coefficients[3],
                                     Y=sampdesign$variables$y[j],V=sampdesign$variables$V[j],alpha_0=cfit$coefficients[1],
                                     alpha_X=cfit$coefficients[2],alpha_Z=cfit$coefficients[3],Xstarstar=ifelse(is.na(sampdesign$variables$xstarstar[j]),0,sampdesign$variables$xstarstar[j]))
    
  }
  
  
  AmatSUM<-add(Amat_All)/N
  BmatSUM<-add(Bmat_All)/N
  A<-(solve(AmatSUM))
  Asvy<-t(solve(AmatSUM))
  
  B<-BmatSUM
  sandwichvar<-(A)%*%B%*%(t((A)))/N
  sqrt(diag(sandwichvar))
  
  parametervec<-c(cfit$coefficients,rcfit$coefficients)
  
  sandwichvar<-(A)%*%B%*%(t((A)))/N
  
  allestfun<-NULL
  for(j in c(1:N)) {
    Xstar<-sampdesign$variables$xstar[j]
    Z<-sampdesign$variables$X2[j]
    Y<-sampdesign$variables$y[j]
    V<-sampdesign$variables$V[j]
    Xstarstar<-ifelse(is.na(sampdesign$variables$xstarstar[j]),0,sampdesign$variables$xstarstar[j])
    estfunj<-RC_Logistic_Est_Eq(parametervec)
    allestfun<-rbind(allestfun,estfunj)
  }
  
  

  
  U<-as.matrix(allestfun)
  
  summary(allestfun)
  summary(U.stage2)
  summary(survey:::estfuns.glm(rcfit)/(weights(rcfit$survey.design,"sampling")))
  summary(survey:::estfuns.glm(rcfit)/(rcfit$prior.weights))
  
#Influence functions
  infl<- U%*%Asvy/N
  v<-vcov(svytotal(infl,  sampdesign))
  
  #Save models
  stage1.model.new<-cfit
  stage2.model.new<-rcfit
  ###Now let's obtain the pieces for matrix A 
  #upper-right quadrant (all elements are zero), j times k
  A.upperright<- matrix(0,nrow=j_dim,ncol=k_dim)
  
  #Here is the upper left computed using the Hessian from stage 1 
  A.upperleft<- -solve(stage1.model.new$naive.cov)/N
  
  #Bottom-right quadrant is just Hessian for stage 2 model
  A.bottomright<- -solve(stage2.model.new$naive.cov)/N
  
  #Save alphas - coefficients from stage 1 model
  alphas.stage1<-coef(stage1.model.new)
  
  #Now, for stage 2 model, we write an estfun function as a function of alphas, holding beta constant 
  stage2.alphas.estfuns.glm<-function (alphas) {
    MyX_Naive<-model.matrix(stage2.model.naive)
    MyX_RC<-model.matrix(stage2.model.new)
    y.j=stage2.model.new$y
    myxhat<-as.vector(MyX_Naive%*%alphas)
    xmat.j <- as.matrix(cbind(MyX_RC[,1],myxhat,MyX_RC[,c(3:ncol(MyX_RC))]))
    mu.j<-as.vector(stage2.model.new$family$linkinv(xmat.j%*%coef(stage2.model.new)))
    r.j<-(y.j - mu.j)/(stage2.model.new$family$variance(mu.j))
    workingweights.j <- as.vector(((stage2.model.new$family$variance(mu.j))^2)/stage2.model.new$family$variance(mu.j))
    
    myestfun<-r.j * workingweights.j * xmat.j
    
    return(myestfun)
  }
  
  
  #Now, use Jacobian function from numderiv to obtain derivatives of stage 2 estimating equation wrt alpha
  #Then split it up, because it computes large matrix with each person's contributions for each parameter
  JacobianAll.stage2<-jacobian(func=stage2.alphas.estfuns.glm, x=alphas.stage1)
  JacobianList.stage2<-lapply(split(JacobianAll.stage2,rep(c(1:N),each=1)),matrix,nrow=k_dim)
  
  #Add these matrices for subjects 1...N then divide by N
  A.bottomleft<-add(JacobianList.stage2)/N
  
  #Now we can combine all four quadrants of our A matrix and obtain AAll
  AAll<-rbind(cbind(A.upperleft,A.upperright),cbind(A.bottomleft,A.bottomright))
  
  #Invert this matrix; needed for sandwich 
  A.inv<-solve(AAll)
  
  for(b in 1:NBOOT){
    # Fit calibration equation
    sampdesign_solnas_boot <- sampdesign_solnas[sample(1:nrow(sampdesign_solnas),nrow(sampdesign_solnas),replace = T),]
    
    #Fit calibration model -- boot
    cfit_boot<-svyglm(xstarstar~xstar+X2,design=sampdesign_solnas_boot,family = stats::gaussian())
    
    #predict xhat and update samp design
    sampdesign <- update(sampdesign,xhat_boot = predict(cfit_boot,newdata=samp))
    
    #Now get regression coefficients for true, naive and RC model
    
    rcfit_boot<-svyglm(y ~ xhat_boot+X2 ,design=sampdesign,family=quasibinomial())
    betas_boot<-summary((rcfit_boot))$coefficients[,1]
    ses_boot<- summary((rcfit_boot))$coefficients[,2]
    #bans<-rbind(bans,c(betas_boot,ses_boot))
    
    boot_list[[i]][[b]] <- data.frame(Sim=i,MI=b,Coeff=names(rcfit_boot$coefficients),
                                      Est=as.numeric(rcfit_boot$coefficients),SE=as.numeric(SE(rcfit_boot)))
    print(b)
  }
  
  
  #Save
  myresults<-rbind(myresults,c(summary(cfit)$coefficients[,1],summary(cfit)$coefficients[,2],
                               summary(tfit)$coefficients[,1],summary(tfit)$coefficients[,2],
                               summary(nfit)$coefficients[,1],summary(nfit)$coefficients[,2],
                               summary(rcfit)$coefficients[,1],summary(rcfit)$coefficients[,2],
                               sqrt(diag(v)),
                               sqrt(diag(sandwichvar))))
  
  boot_list[[i]] = do.call(rbind,boot_list[[i]])
  
  print(i)
} ### End of simulation loop

myresults<-data.frame(myresults)

names(myresults)<-c("a0","aX","aZ","se_a0","se_aX","se_aZ","b0_True","bX_True","bZ_True","se_b0_True","se_bX_True","se_bZ_True",
                    "b0_Naive","bX_Naive","bZ_Naive","se_b0_Naive","se_bX_Naive","se_bZ_Naive",
                    "b0_RC","bX_RC","bZ_RC","se_b0_RC","se_bX_RC","se_bZ_RC",
                    "se_a0_survey","se_aX_survey","se_aZ_survey","se_b0_survey","se_bX_survey","se_bZ_survey",
                    "se_a0_sand","se_aX_sand","se_aZ_sand","se_b0_sand","se_bX_sand","se_bZ_sand")
