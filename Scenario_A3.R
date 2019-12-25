#############################################################################
##                                                                         ##
## Simulations for scenario A3, with time-dependent (TD) and linear effect ##
##                                                                         ##
#############################################################################

# Last update: December 24, 2019

### Simulations settings ###

rm(list=ls())
gc()

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")
SBPdata <- read.csv("SBPdata.csv")

set.seed(999935)
n.sims <- 300

N <- 1000
incidence.rate <- 0.6
fupt <- 296

# True TD function
tdx5<-function(x){
  g<-1.5*exp(-((x-80)/80)^{2})+0.2
  return(g)
}

SBPdata$sbp.ideal.td <- tdx5(SBPdata$Stop)*SBPdata$sbp.ideal

### Declaration fof objects to be saved ###

# Ideal case, true model without TEL
coefA3iM7TD <- list()
knotsA3iM7TD <- list()
devA3iM7 <- rep(NA, n.sims)

# Ideal case, full model without TEL
coefA3iM8TD <- list()
knotsA3iM8TD <- list()
coefA3iM8NL <- list()
knotsA3iM8NL <- list()
devA3iM8 <- rep(NA, n.sims)

# Baseline only, true model without TEL
coefA3bM7TD <- list()
knotsA3bM7TD <- list()
devA3bM7 <- rep(NA, n.sims)

# Baseline only, full model without TEL
coefA3bM8TD <- list()
knotsA3bM8TD <- list()
coefA3bM8NL <- list()
knotsA3bM8NL <- list()
devA3bM8 <- rep(NA, n.sims)

# Sparse, true model without TEL
coefA3sM7TD <- list()
knotsA3sM7TD <- list()
devA3sM7 <- rep(NA, n.sims)

# Sparse, true model with TEL
coefA3sM3TD <- list()
knotsA3sM3TD <- list()
coefA3sM3TEL <- list()
knotsA3sM3TEL <- list()
devA3sM3 <- rep(NA, n.sims)

# Sparse, full model without TEL
coefA3sM8TD <- list()
knotsA3sM8TD <- list()
coefA3sM8NL <- list()
knotsA3sM8NL <- list()
devA3sM8 <- rep(NA, n.sims)

# Sparse, full model with TEL
coefA3sM4TD <- list()
knotsA3sM4TD <- list()
coefA3sM4NL <- list()
knotsA3sM4NL <- list()
coefA3sM4TEL <- list()
knotsA3sM4TEL <- list()
devA3sM4 <- rep(NA, n.sims)

### Simulation loop ###

options(warn=1)
library(PermAlgo)

for ( i in 1:n.sims){
  cat("i=",i,'\n')

  eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/fupt))
  censortimes <- rep(fupt, times=N)

  datA3 <- permalgorithm(numSubjects=N, maxTime=fupt, 
                         Xmat=as.matrix(SBPdata[,c('sbp.ideal.td','sbp.ideal','sbp.bl','sbp.sparse','tel')]), 
                         XmatNames=c('sbp.ideal.td','sbp.ideal','sbp.bl','sbp.sparse','tel'),
                         eventRandom=eventtimes, censorRandom=censortimes,
                         betas=c(1,0,0,0,0))

  ### Ideal case ###
  
  # True model without TEL
  modA3iM7 <- last_prog(data=datA3, Type=c("Start","Stop","Event"),
                         variables=c("sbp.ideal"), TD=c(1), NL=c(0), m=1, p=2,
                         knots=-999)
  coefA3iM7TD[[i]] <- modA3iM7$coefficients_splines_TD
  knotsA3iM7TD[[i]] <- modA3iM7$knots_time
  devA3iM7[i]<- -2*modA3iM7$Partial_Log_Likelihood

  # Full model without TEL
  modA3iM8 <- last_prog(data=datA3, Type=c("Start","Stop","Event"),
                         variables=c("sbp.ideal"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA3iM8NL[[i]] <- modA3iM8$coefficients_splines_NL
  knotsA3iM8NL[[i]] <- modA3iM8$knots_covariates
  coefA3iM8TD[[i]] <- modA3iM8$coefficients_splines_TD
  knotsA3iM8TD[[i]] <- modA3iM8$knots_time
  devA3iM8[i]<- -2*modA3iM8$Partial_Log_Likelihood

  rm(modA3iM7, modA3iM8); gc()
  
  ### Baseline only ###
  
  # True model without TEL
  modA3bM7 <- last_prog(data=datA3, Type=c("Start","Stop","Event"),
                         variables=c("sbp.bl"), TD=c(1), NL=c(0), m=1, p=2,
                         knots=-999)
  coefA3bM7TD[[i]] <- modA3bM7$coefficients_splines_TD
  knotsA3bM7TD[[i]] <- modA3bM7$knots_time
  devA3bM7[i]<- -2*modA3bM7$Partial_Log_Likelihood  

  # Full model without TEL
  modA3bM8 <- last_prog(data=datA3, Type=c("Start","Stop","Event"),
                         variables=c("sbp.bl"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA3bM8NL[[i]] <- modA3bM8$coefficients_splines_NL
  knotsA3bM8NL[[i]] <- modA3bM8$knots_covariates
  coefA3bM8TD[[i]] <- modA3bM8$coefficients_splines_TD
  knotsA3bM8TD[[i]] <- modA3bM8$knots_time
  devA3bM8[i]<- -2*modA3bM8$Partial_Log_Likelihood

  rm(modA3bM7, modA3bM8); gc()
  
  ### Sparse ###
  
  # True model without TEL
  modA3sM7 <- last_prog(data=datA3, Type=c("Start","Stop","Event"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(0), m=1, p=2,
                         knots=-999)
  coefA3sM7TD[[i]] <- modA3sM7$coefficients_splines_TD
  knotsA3sM7TD[[i]] <- modA3sM7$knots_time
  devA3sM7[i]<- -2*modA3sM7$Partial_Log_Likelihood  
  
  # True model with TEL
  modA3sM3 <- last_prog_TEL(data=datA3, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(0), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA3sM3TD[[i]] <- modA3sM3$coefficients_splines_TD
  knotsA3sM3TD[[i]] <- modA3sM3$knots_time
  coefA3sM3TEL[[i]] <- modA3sM3$coefficients_splines_TEL
  knotsA3sM3TEL[[i]] <- modA3sM3$knots_TEL
  devA3sM3[i] <- -2*modA3sM3$Partial_Log_Likelihood
  
  # Full model without TEL
  modA3sM8 <- last_prog(data=datA3, Type=c("Start","Stop","Event"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA3sM8NL[[i]] <- modA3sM8$coefficients_splines_NL
  knotsA3sM8NL[[i]] <- modA3sM8$knots_covariates
  coefA3sM8TD[[i]] <- modA3sM8$coefficients_splines_TD
  knotsA3sM8TD[[i]] <- modA3sM8$knots_time
  devA3sM8[i]<- -2*modA3sM8$Partial_Log_Likelihood
  
  # Full model with TEL
  modA3sM4 <- last_prog_TEL(data=datA3, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA3sM4TEL[[i]] <- modA3sM4$coefficients_splines_TEL
  knotsA3sM4TEL[[i]] <- modA3sM4$knots_TEL
  coefA3sM4TD[[i]] <- modA3sM4$coefficients_splines_TD
  knotsA3sM4TD[[i]] <- modA3sM4$knots_time
  coefA3sM4NL[[i]] <- modA3sM4$coefficients_splines_NL
  knotsA3sM4NL[[i]] <- modA3sM4$knots_covariates
  devA3sM4[i] <- -2*modA3sM4$Partial_Log_Likelihood

  rm(modA3sM7, modA3sM3, modA3sM8, modA3sM4, datA3, eventtimes, censortimes); gc()  
}
gc()

save(n.sims, N, incidence.rate, fupt, tdx5,
  coefA3iM7TD, knotsA3iM7TD, devA3iM7, 
  coefA3iM8TD, knotsA3iM8TD, coefA3iM8NL, knotsA3iM8NL, devA3iM8, 
  coefA3bM7TD, knotsA3bM7TD, devA3bM7, 
  coefA3bM8TD, knotsA3bM8TD, coefA3bM8NL, knotsA3bM8NL, devA3bM8, 
  coefA3sM7TD, knotsA3sM7TD, devA3sM7, 
  coefA3sM3TD, knotsA3sM3TD, coefA3sM3TEL, knotsA3sM3TEL, devA3sM3, 
  coefA3sM8TD, knotsA3sM8TD, coefA3sM8NL, knotsA3sM8NL, devA3sM8, 
  coefA3sM4TD, knotsA3sM4TD, coefA3sM4NL, knotsA3sM4NL, coefA3sM4TEL, knotsA3sM4TEL, devA3sM4,
  file="Scenario_A3.RData")

