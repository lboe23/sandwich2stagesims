# sandwich2stagesims
Code for simulation studies in paper, "Practical Considerations for sandwich variance estimation in two-stage regression settings", Boe, Lumley, and Shaw

Code files include:

(1) R code that runs simulations for data from a simple random sample to assess the performance of the sandwich variance estimator compared to the boostrap, when the outcome (stage 2) model is a logistic regression model: RandomSampleSims_Sandwich.R

(2) R code that runs simulations for data from a simple random sample to assess the performance of the sandwich variance estimator compared to the boostrap, when the outcome (stage 2) model is a Cox proportional hazards regression model: RandomSampleSims_Cox_Sandwich.R

(3) R code modified from Baldoni et al. 2021 that simulates a superpopulation and draws 1000 samples from it using a three-stage, stratified sampling scheme. Our code combines the files generateTargetPopulation.R, generateTargetPopulationData.R, and sampleTargetPopulation.R from https://github.com/plbaldoni/HCHSsim/tree/master/R into one code file and modifies the code in order to simulate our covariate and event times of interest and adjust the sample size. Our resulting code file is called GenerateData_ComplexSurveySims.R, which produces 1000 simulated data sets to be used for the simulation in file (4).

(4) R code that runs simulations for complex survey design data to assess the performance of the sandwich variance estimator compared to the multiple-imputation based variance estimator, when the outcome (stage 2) model is a logistic regression model: ComplexSurveySims_Sandwich.R
