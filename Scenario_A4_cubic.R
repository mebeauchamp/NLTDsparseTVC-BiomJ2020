#######################################################################################
##                                                                                   ##
## Simulations for scenario A4, with non-linear (NL) and time-dependent (TD) effects ##
## estimated with cubic splines                                                      ##
##                                                                                   ##
#######################################################################################

# Last update: December 24, 2019

### Simulations settings ###

rm(list=ls())
gc()

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")
SBPdata <- read.csv("SBPdata.csv")

set.seed(11781)
n.sims <- 300

N <- 1000
incidence.rate <- 0.6
fupt <- 296

# True NL function
nlx5<-function(x){
  return((x+0.5)^2/6)
}

# True TD function
tdx5<-function(x){
  g<-1.5*exp(-((x-80)/80)^{2})+0.2
  return(g)
}

SBPdata$sbp.ideal.nltd <- tdx5(SBPdata$Stop)*nlx5(SBPdata$sbp.ideal)

### Declaration fof objects to be saved ###

# Ideal case, true and full model without TEL
coefA4iM8cTD <- list()
knotsA4iM8cTD <- list()
coefA4iM8cNL <- list()
knotsA4iM8cNL <- list()
devA4iM8c <- rep(NA, n.sims)

# Baseline only, true and full model without TEL
coefA4bM8cTD <- list()
knotsA4bM8cTD <- list()
coefA4bM8cNL <- list()
knotsA4bM8cNL <- list()
devA4bM8c <- rep(NA, n.sims)

# Sparse, true and full model without TEL
coefA4sM8cTD <- list()
knotsA4sM8cTD <- list()
coefA4sM8cNL <- list()
knotsA4sM8cNL <- list()
devA4sM8c <- rep(NA, n.sims)

# Sparse, true and full model with TEL
coefA4sM4cTD <- list()
knotsA4sM4cTD <- list()
coefA4sM4cNL <- list()
knotsA4sM4cNL <- list()
coefA4sM4cTEL <- list()
knotsA4sM4cTEL <- list()
devA4sM4c <- rep(NA, n.sims)

### Simulation loop ###

options(warn=1)
library(PermAlgo)

for ( i in 1:n.sims){
  cat("i=",i,'\n')

  eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/fupt))
  censortimes <- rep(fupt, times=N)

  datA4 <- permalgorithm(numSubjects=N, maxTime=fupt, 
                         Xmat=as.matrix(SBPdata[,c('sbp.ideal.nltd','sbp.ideal','sbp.bl','sbp.sparse','tel')]), 
                         XmatNames=c('sbp.ideal.nltd','sbp.ideal','sbp.bl','sbp.sparse','tel'),
                         eventRandom=eventtimes, censorRandom=censortimes,
                         betas=c(1,0,0,0,0))

  ### Ideal case ###
  
  # True and full model without TEL
  modA4iM8c <- last_prog(data=datA4, Type=c("Start","Stop","Event"),
                         variables=c("sbp.ideal"), TD=c(1), NL=c(1), m=1, p=3,
                         knots=-999)
  coefA4iM8cNL[[i]] <- modA4iM8c$coefficients_splines_NL
  knotsA4iM8cNL[[i]] <- modA4iM8c$knots_covariates
  coefA4iM8cTD[[i]] <- modA4iM8c$coefficients_splines_TD
  knotsA4iM8cTD[[i]] <- modA4iM8c$knots_time
  devA4iM8c[i]<- -2*modA4iM8c$Partial_Log_Likelihood

  rm(modA4iM8c); gc()
  
  ### Baseline only ###
  
  # True and full model without TEL
  modA4bM8c <- last_prog(data=datA4, Type=c("Start","Stop","Event"),
                         variables=c("sbp.bl"), TD=c(1), NL=c(1), m=1, p=3,
                         knots=-999)
  coefA4bM8cNL[[i]] <- modA4bM8c$coefficients_splines_NL
  knotsA4bM8cNL[[i]] <- modA4bM8c$knots_covariates
  coefA4bM8cTD[[i]] <- modA4bM8c$coefficients_splines_TD
  knotsA4bM8cTD[[i]] <- modA4bM8c$knots_time
  devA4bM8c[i]<- -2*modA4bM8c$Partial_Log_Likelihood

  rm(modA4bM8c); gc()
  
  ### Sparse ###

  # True and full model without TEL
  modA4sM8c <- last_prog(data=datA4, Type=c("Start","Stop","Event"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), m=1, p=3,
                         knots=-999)
  coefA4sM8cNL[[i]] <- modA4sM8c$coefficients_splines_NL
  knotsA4sM8cNL[[i]] <- modA4sM8c$knots_covariates
  coefA4sM8cTD[[i]] <- modA4sM8c$coefficients_splines_TD
  knotsA4sM8cTD[[i]] <- modA4sM8c$knots_time
  devA4sM8c[i]<- -2*modA4sM8c$Partial_Log_Likelihood
  
  # True and full model with TEL
  modA4sM4c <- last_prog_TEL(data=datA4, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), TEL=c(1), m=1, p=3,
                         knots=-999)
  coefA4sM4cTEL[[i]] <- modA4sM4c$coefficients_splines_TEL
  knotsA4sM4cTEL[[i]] <- modA4sM4c$knots_TEL
  coefA4sM4cTD[[i]] <- modA4sM4c$coefficients_splines_TD
  knotsA4sM4cTD[[i]] <- modA4sM4c$knots_time
  coefA4sM4cNL[[i]] <- modA4sM4c$coefficients_splines_NL
  knotsA4sM4cNL[[i]] <- modA4sM4c$knots_covariates
  devA4sM4c[i] <- -2*modA4sM4c$Partial_Log_Likelihood

  rm(modA4sM8c, modA4sM4c, datA4, eventtimes, censortimes); gc()  
}
gc()

save(n.sims, N, incidence.rate, fupt, nlx5, tdx5,
  coefA4iM8cTD, knotsA4iM8cTD, coefA4iM8cNL, knotsA4iM8cNL, devA4iM8c, 
  coefA4bM8cTD, knotsA4bM8cTD, coefA4bM8cNL, knotsA4bM8cNL, devA4bM8c, 
  coefA4sM8cTD, knotsA4sM8cTD, coefA4sM8cNL, knotsA4sM8cNL, devA4sM8c, 
  coefA4sM4cTD, knotsA4sM4cTD, coefA4sM4cNL, knotsA4sM4cNL, coefA4sM4cTEL, knotsA4sM4cTEL, devA4sM4c, 
  file="Scenario_A4_cubic.RData")

