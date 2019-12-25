############################################################
##                                                        ##
## Simulations for scenario A1, with linear and PH effect ##
##                                                        ##
############################################################

# Last update: December 24, 2019

### Simulations settings ###

rm(list=ls())
gc()

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")
SBPdata <- read.csv("SBPdata.csv")

set.seed(703444)
n.sims <- 300

N <- 1000
incidence.rate <- 0.6
fupt <- 296

### Declaration of objects to be saved ###

# Ideal case, true model without TEL
coefA1iM5 <- rep(NA, n.sims)
devAA1iM5 <- rep(NA, n.sims)

# Ideal case, full model without TEL
coefA1iM8TD <- list()
knotsA1iM8TD <- list()
coefA1iM8NL <- list()
knotsA1iM8NL <- list()
devA1iM8 <- rep(NA, n.sims)

# Baseline only, true model without TEL
coefA1bM5 <- rep(NA, n.sims)
devAA1bM5 <- rep(NA, n.sims)

# Baseline only, full model without TEL
coefA1bM8TD <- list()
knotsA1bM8TD <- list()
coefA1bM8NL <- list()
knotsA1bM8NL <- list()
devA1bM8 <- rep(NA, n.sims)

# Sparse, true model without TEL
coefA1sM5 <- rep(NA, n.sims)
devAA1sM5 <- rep(NA, n.sims)

# Sparse, true model with TEL
coefA1sM1TEL <- list()
knotsA1sM1TEL <- list()
devA1sM1 <- rep(NA, n.sims)

# Sparse, full model without TEL
coefA1sM8TD <- list()
knotsA1sM8TD <- list()
coefA1sM8NL <- list()
knotsA1sM8NL <- list()
devA1sM8 <- rep(NA, n.sims)

# Sparse, full model with TEL
coefA1sM4TD <- list()
knotsA1sM4TD <- list()
coefA1sM4NL <- list()
knotsA1sM4NL <- list()
coefA1sM4TEL <- list()
knotsA1sM4TEL <- list()
devA1sM4 <- rep(NA, n.sims)

### Simualtion loop ###

options(warn=1)
library(PermAlgo)

for (i in 1:n.sims){
  cat("i=",i,'\n')

  eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/fupt))
  censortimes <- rep(fupt, times=N)

  datA1 <- permalgorithm(numSubjects=N, maxTime=fupt, 
                         Xmat=as.matrix(SBPdata[,c('sbp.ideal','sbp.bl','sbp.sparse','tel')]), 
                         XmatNames=c('sbp.ideal','sbp.bl','sbp.sparse','tel'),
                         eventRandom=eventtimes, censorRandom=censortimes,
                         betas=c(log(1.5),0,0,0))

  ### Ideal case ###
  
  # True model without TEL
  modA1iM5 <- coxph(Surv(Start,Stop, Event) ~ sbp.ideal, data=datA1)
  coefA1iM5[i] <- summary(modA1iM5)$coef[1]
  devAA1iM5[i] <- -2*summary(modA1iM5)$loglik[2] 

  # Full model without TEL
  modA1iM8 <- last_prog(data=datA1, Type=c("Start","Stop","Event"),
                         variables=c("sbp.ideal"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA1iM8NL[[i]] <- modA1iM8$coefficients_splines_NL
  knotsA1iM8NL[[i]] <- modA1iM8$knots_covariates
  coefA1iM8TD[[i]] <- modA1iM8$coefficients_splines_TD
  knotsA1iM8TD[[i]] <- modA1iM8$knots_time
  devA1iM8[i]<- -2*modA1iM8$Partial_Log_Likelihood

  rm(modA1iM5, modA1iM8); gc()
  
  ### Baseline only ###
  
  # True model without TEL
  modA1bM5 <- coxph(Surv(Start,Stop, Event) ~ sbp.bl, data=datA1)
  coefA1bM5[i] <- summary(modA1bM5)$coef[1]
  devAA1bM5[i] <- -2*summary(modA1bM5)$loglik[2]   

  # Full model without TEL
  modA1bM8 <- last_prog(data=datA1, Type=c("Start","Stop","Event"),
                         variables=c("sbp.bl"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA1bM8NL[[i]] <- modA1bM8$coefficients_splines_NL
  knotsA1bM8NL[[i]] <- modA1bM8$knots_covariates
  coefA1bM8TD[[i]] <- modA1bM8$coefficients_splines_TD
  knotsA1bM8TD[[i]] <- modA1bM8$knots_time
  devA1bM8[i]<- -2*modA1bM8$Partial_Log_Likelihood

  rm(modA1bM5, modA1bM8); gc()
  
  ### Sparse ###
  
  # True model without TEL
  modA1sM5 <- coxph(Surv(Start,Stop, Event) ~ sbp.sparse, data=datA1)
  coefA1sM5[i] <- summary(modA1sM5)$coef[1]
  devAA1sM5[i] <- -2*summary(modA1sM5)$loglik[2]    
  
  # True model with TEL
  modA1sM1 <- last_prog_TEL(data=datA1, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(0), NL=c(0), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA1sM1TEL[[i]] <- modA1sM1$coefficients_splines_TEL
  knotsA1sM1TEL[[i]] <- modA1sM1$knots_TEL
  devA1sM1[i] <- -2*modA1sM1$Partial_Log_Likelihood
  
  # Full model without TEL
  modA1sM8 <- last_prog(data=datA1, Type=c("Start","Stop","Event"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), m=1, p=2,
                         knots=-999)
  coefA1sM8NL[[i]] <- modA1sM8$coefficients_splines_NL
  knotsA1sM8NL[[i]] <- modA1sM8$knots_covariates
  coefA1sM8TD[[i]] <- modA1sM8$coefficients_splines_TD
  knotsA1sM8TD[[i]] <- modA1sM8$knots_time
  devA1sM8[i]<- -2*modA1sM8$Partial_Log_Likelihood
  
  # Full model with TEL
  modA1sM4 <- last_prog_TEL(data=datA1, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(1), NL=c(1), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA1sM4TEL[[i]] <- modA1sM4$coefficients_splines_TEL
  knotsA1sM4TEL[[i]] <- modA1sM4$knots_TEL
  coefA1sM4TD[[i]] <- modA1sM4$coefficients_splines_TD
  knotsA1sM4TD[[i]] <- modA1sM4$knots_time
  coefA1sM4NL[[i]] <- modA1sM4$coefficients_splines_NL
  knotsA1sM4NL[[i]] <- modA1sM4$knots_covariates
  devA1sM4[i] <- -2*modA1sM4$Partial_Log_Likelihood

  rm(modA1sM5, modA1sM1, modA1sM8, modA1sM4, datA1, eventtimes, censortimes); gc()  
}
gc()

save(n.sims, N, incidence.rate, fupt,
  coefA1iM5, devAA1iM5, 
  coefA1iM8TD, knotsA1iM8TD, coefA1iM8NL, knotsA1iM8NL, devA1iM8, 
  coefA1bM5, devAA1bM5, coefA1bM8TD, knotsA1bM8TD, coefA1bM8NL, knotsA1bM8NL, devA1bM8, 
  coefA1sM5, devAA1sM5, 
  coefA1sM1TEL, knotsA1sM1TEL, devA1sM1, 
  coefA1sM8TD, knotsA1sM8TD, coefA1sM8NL, knotsA1sM8NL, devA1sM8, 
  coefA1sM4TD, knotsA1sM4TD, coefA1sM4NL, knotsA1sM4NL, coefA1sM4TEL, knotsA1sM4TEL, devA1sM4, 
  file="Scenario_A1.RData")

