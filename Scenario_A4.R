#######################################################################################
##                                                                                   ##
## Simulations for scenario A4, with non-linear (NL) and time-dependent (TD) effects ##
##                                                                                   ##
#######################################################################################

# Last update: December 24, 2019

## NOTE: It includes investigations on the 

### Simulations settings ###

rm(list=ls())
gc()

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")
SBPdata <- read.csv("SBPdata.csv")

set.seed(5378111)
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
coefA4iM8TD <- list()
knotsA4iM8TD <- list()
coefA4iM8NL <- list()
knotsA4iM8NL <- list()
devA4iM8 <- rep(NA, n.sims)

# Baseline only, true and full model without TEL
coefA4bM8TD <- list()
knotsA4bM8TD <- list()
coefA4bM8NL <- list()
knotsA4bM8NL <- list()
devA4bM8 <- rep(NA, n.sims)

# Sparse, true and full model without TEL
coefA4sM8TD <- list()
knotsA4sM8TD <- list()
coefA4sM8NL <- list()
knotsA4sM8NL <- list()
devA4sM8 <- rep(NA, n.sims)

# Sparse, true and full model with TEL
coefA4sM4TD <- list()
knotsA4sM4TD <- list()
coefA4sM4NL <- list()
knotsA4sM4NL <- list()
coefA4sM4TEL <- list()
knotsA4sM4TEL <- list()
devA4sM4 <- rep(NA, n.sims)

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
  modA4iM8 <- last_prog(data=datA4, Type=c("Start","Stop","Event"),
                        variables=c("sbp.ideal"), TD=c(1), NL=c(1), m=1, p=2,
                        knots=-999)
  coefA4iM8NL[[i]] <- modA4iM8$coefficients_splines_NL
  knotsA4iM8NL[[i]] <- modA4iM8$knots_covariates
  coefA4iM8TD[[i]] <- modA4iM8$coefficients_splines_TD
  knotsA4iM8TD[[i]] <- modA4iM8$knots_time
  devA4iM8[i]<- -2*modA4iM8$Partial_Log_Likelihood
   
  rm(modA4iM8); gc()
   
  ### Baseline only ###
   
  # True and full model without TEL
  modA4bM8 <- last_prog(data=datA4, Type=c("Start","Stop","Event"),
                        variables=c("sbp.bl"), TD=c(1), NL=c(1), m=1, p=2,
                        knots=-999)
  coefA4bM8NL[[i]] <- modA4bM8$coefficients_splines_NL
  knotsA4bM8NL[[i]] <- modA4bM8$knots_covariates
  coefA4bM8TD[[i]] <- modA4bM8$coefficients_splines_TD
  knotsA4bM8TD[[i]] <- modA4bM8$knots_time
  devA4bM8[i]<- -2*modA4bM8$Partial_Log_Likelihood
   
  rm(modA4bM8); gc()
   
  ### Sparse ###
   
  # True and full model without TEL
  modA4sM8 <- last_prog(data=datA4, Type=c("Start","Stop","Event"),
                        variables=c("sbp.sparse"), TD=c(1), NL=c(1), m=1, p=2,
                        knots=-999)
  coefA4sM8NL[[i]] <- modA4sM8$coefficients_splines_NL
  knotsA4sM8NL[[i]] <- modA4sM8$knots_covariates
  coefA4sM8TD[[i]] <- modA4sM8$coefficients_splines_TD
  knotsA4sM8TD[[i]] <- modA4sM8$knots_time
  devA4sM8[i]<- -2*modA4sM8$Partial_Log_Likelihood
  
  # True and full model with TEL
  modA4sM4 <- last_prog_TEL(data=datA4, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA4sM4TEL[[i]] <- modA4sM4$coefficients_splines_TEL
  knotsA4sM4TEL[[i]] <- modA4sM4$knots_TEL
  coefA4sM4TD[[i]] <- modA4sM4$coefficients_splines_TD
  knotsA4sM4TD[[i]] <- modA4sM4$knots_time
  coefA4sM4NL[[i]] <- modA4sM4$coefficients_splines_NL
  knotsA4sM4NL[[i]] <- modA4sM4$knots_covariates
  devA4sM4[i] <- -2*modA4sM4$Partial_Log_Likelihood
  
  rm(modA4sM8, modA4sM4, datA4, eventtimes, censortimes); gc()
}
gc()

save(n.sims, N, incidence.rate, fupt, nlx5, tdx5,
           coefA4iM8TD, knotsA4iM8TD, coefA4iM8NL, knotsA4iM8NL, devA4iM8, 
           coefA4bM8TD, knotsA4bM8TD, coefA4bM8NL, knotsA4bM8NL, devA4bM8, 
           coefA4sM8TD, knotsA4sM8TD, coefA4sM8NL, knotsA4sM8NL, devA4sM8, 
           coefA4sM4TD, knotsA4sM4TD, coefA4sM4NL, knotsA4sM4NL, coefA4sM4TEL, knotsA4sM4TEL, devA4sM4, 
           file="Scenario_A4.RData")
