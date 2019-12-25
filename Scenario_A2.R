#####################################################################
##                                                                 ##
## Simulations for scenario A2, with non-linear (NL) and PH effect ##
##                                                                 ##
#####################################################################

# Last update: December 24, 2019

### Simulations settings ###

rm(list=ls())
gc()

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")
SBPdata <- read.csv("SBPdata.csv")

set.seed(123235)
n.sims <- 300

N <- 1000
incidence.rate <- 0.6
fupt <- 296

# True NL function
nlx5<-function(x){
  return((x+0.5)^2/6)
}

SBPdata$sbp.ideal.nl <- nlx5(SBPdata$sbp.ideal)

### Declaration fof objects to be saved ###

# Ideal case, true model without TEL
coefA2iM6NL <- list()
knotsA2iM6NL <- list()
devA2iM6 <- rep(NA, n.sims)

# Ideal case, full model without TEL
coefA2iM8TD <- list()
knotsA2iM8TD <- list()
coefA2iM8NL <- list()
knotsA2iM8NL <- list()
devA2iM8 <- rep(NA, n.sims)

# Baseline only, true model without TEL
coefA2bM6NL <- list()
knotsA2bM6NL <- list()
devA2bM6 <- rep(NA, n.sims)

# Baseline only, full model without TEL
coefA2bM8TD <- list()
knotsA2bM8TD <- list()
coefA2bM8NL <- list()
knotsA2bM8NL <- list()
devA2bM8 <- rep(NA, n.sims)

# Sparse, true model without TEL
coefA2sM6NL <- list()
knotsA2sM6NL <- list()
devA2sM6 <- rep(NA, n.sims)

# Sparse, true model with TEL
coefA2sM2NL <- list()
knotsA2sM2NL <- list()
coefA2sM2TEL <- list()
knotsA2sM2TEL <- list()
devA2sM2 <- rep(NA, n.sims)

# Sparse, full model without TEL
coefA2sM8TD <- list()
knotsA2sM8TD <- list()
coefA2sM8NL <- list()
knotsA2sM8NL <- list()
devA2sM8 <- rep(NA, n.sims)

# Sparse, full model with TEL
coefA2sM4TD <- list()
knotsA2sM4TD <- list()
coefA2sM4NL <- list()
knotsA2sM4NL <- list()
coefA2sM4TEL <- list()
knotsA2sM4TEL <- list()
devA2sM4 <- rep(NA, n.sims)

### Simualtion loop ###

options(warn=1)
library(PermAlgo)

for ( i in 1:n.sims){
  cat("i=",i,'\n')

  eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/fupt))
  censortimes <- rep(fupt, times=N)

  datA2 <- permalgorithm(numSubjects=N, maxTime=fupt, 
                         Xmat=as.matrix(SBPdata[,c('sbp.ideal.nl','sbp.ideal','sbp.bl','sbp.sparse','tel')]), 
                         XmatNames=c('sbp.ideal.nl','sbp.ideal','sbp.bl','sbp.sparse','tel'),
                         eventRandom=eventtimes, censorRandom=censortimes,
                         betas=c(1,0,0,0,0))

  ### Ideal case ###
  
  # True model without TEL
  modA2iM6 <- last_prog(data=datA2, Type=c("Start","Stop","Event"),
                         variables=c("sbp.ideal"), TD=c(0), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA2iM6NL[[i]] <- modA2iM6$coefficients_splines_NL
  knotsA2iM6NL[[i]] <- modA2iM6$knots_covariates
  devA2iM6[i]<- -2*modA2iM6$Partial_Log_Likelihood

  # Full model without TEL
  modA2iM8 <- last_prog(data=datA2, Type=c("Start","Stop","Event"),
                         variables=c("sbp.ideal"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA2iM8NL[[i]] <- modA2iM8$coefficients_splines_NL
  knotsA2iM8NL[[i]] <- modA2iM8$knots_covariates
  coefA2iM8TD[[i]] <- modA2iM8$coefficients_splines_TD
  knotsA2iM8TD[[i]] <- modA2iM8$knots_time
  devA2iM8[i]<- -2*modA2iM8$Partial_Log_Likelihood

  rm(modA2iM6, modA2iM8); gc()
  
  ### Baseline only ###
  
  # True model without TEL
  modA2bM6 <- last_prog(data=datA2, Type=c("Start","Stop","Event"),
                         variables=c("sbp.bl"), TD=c(0), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA2bM6NL[[i]] <- modA2bM6$coefficients_splines_NL
  knotsA2bM6NL[[i]] <- modA2bM6$knots_covariates
  devA2bM6[i]<- -2*modA2bM6$Partial_Log_Likelihood  

  # Full model without TEL
  modA2bM8 <- last_prog(data=datA2, Type=c("Start","Stop","Event"),
                         variables=c("sbp.bl"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA2bM8NL[[i]] <- modA2bM8$coefficients_splines_NL
  knotsA2bM8NL[[i]] <- modA2bM8$knots_covariates
  coefA2bM8TD[[i]] <- modA2bM8$coefficients_splines_TD
  knotsA2bM8TD[[i]] <- modA2bM8$knots_time
  devA2bM8[i]<- -2*modA2bM8$Partial_Log_Likelihood

  rm(modA2bM6, modA2bM8); gc()
  
  ### Sparse ###
  
  # True model without TEL
  modA2sM6 <- last_prog(data=datA2, Type=c("Start","Stop","Event"),
                         variables=c("sbp.sparse"), TD=c(0), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA2sM6NL[[i]] <- modA2sM6$coefficients_splines_NL
  knotsA2sM6NL[[i]] <- modA2sM6$knots_covariates
  devA2sM6[i]<- -2*modA2sM6$Partial_Log_Likelihood  
  
  # True model with TEL
  modA2sM2 <- last_prog_TEL(data=datA2, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(0), NL=c(1), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA2sM2NL[[i]] <- modA2sM2$coefficients_splines_NL
  knotsA2sM2NL[[i]] <- modA2sM2$knots_covariates
  coefA2sM2TEL[[i]] <- modA2sM2$coefficients_splines_TEL
  knotsA2sM2TEL[[i]] <- modA2sM2$knots_TEL
  devA2sM2[i] <- -2*modA2sM2$Partial_Log_Likelihood
  
  # Full model without TEL
  modA2sM8 <- last_prog(data=datA2, Type=c("Start","Stop","Event"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA2sM8NL[[i]] <- modA2sM8$coefficients_splines_NL
  knotsA2sM8NL[[i]] <- modA2sM8$knots_covariates
  coefA2sM8TD[[i]] <- modA2sM8$coefficients_splines_TD
  knotsA2sM8TD[[i]] <- modA2sM8$knots_time
  devA2sM8[i]<- -2*modA2sM8$Partial_Log_Likelihood
  
  # Full model with TEL
  modA2sM4 <- last_prog_TEL(data=datA2, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA2sM4TEL[[i]] <- modA2sM4$coefficients_splines_TEL
  knotsA2sM4TEL[[i]] <- modA2sM4$knots_TEL
  coefA2sM4TD[[i]] <- modA2sM4$coefficients_splines_TD
  knotsA2sM4TD[[i]] <- modA2sM4$knots_time
  coefA2sM4NL[[i]] <- modA2sM4$coefficients_splines_NL
  knotsA2sM4NL[[i]] <- modA2sM4$knots_covariates
  devA2sM4[i] <- -2*modA2sM4$Partial_Log_Likelihood

  rm(modA2sM6, modA2sM2, modA2sM8, modA2sM4, datA2, eventtimes, censortimes); gc()  
}
gc()

save(n.sims, N, incidence.rate, fupt, nlx5, 
  coefA2iM6NL, knotsA2iM6NL, devA2iM6, 
  coefA2iM8TD, knotsA2iM8TD, coefA2iM8NL, knotsA2iM8NL, devA2iM8, 
  coefA2bM6NL, knotsA2bM6NL, devA2bM6, 
  coefA2bM8TD, knotsA2bM8TD, coefA2bM8NL, knotsA2bM8NL, devA2bM8, 
  coefA2sM6NL, knotsA2sM6NL, devA2sM6, 
  coefA2sM2NL, knotsA2sM2NL, coefA2sM2TEL, knotsA2sM2TEL, devA2sM2, 
  coefA2sM8TD, knotsA2sM8TD, coefA2sM8NL, knotsA2sM8NL, devA2sM8,
  coefA2sM4TD, knotsA2sM4TD, coefA2sM4NL, knotsA2sM4NL, coefA2sM4TEL, knotsA2sM4TEL, devA2sM4, 
  file="Scenario_A2.RData")

