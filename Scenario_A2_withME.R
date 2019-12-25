############################################################################################################
##                                                                                                        ##
## Simulations for scenario A2, with non-linear (NL) and PH effect, when measurment errors (ME) is added. ##
## Application of SIMEX at the data analysis stage.                                                       ##
## For sparse data only.                                                                                  ##
##                                                                                                        ##
############################################################################################################

# Last update: December 24, 2019

### Simulation settings ###

rm(list=ls())
gc()

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")
SBPdata <- read.csv("SBPdata.csv")

# Choose one scenario for standard deviation (SD) for measuremtent error (ME), by putting 1 line or  
# the otherin comment 
  # sd.ME is the SD of Normal distribution from which ME is generated 
  # sd.ME.SIMEX is the SD of ME assumed when applying SIMEX, on purpose underestimated
#sd.ME <- 0.5; sd.ME.SIMEX <- 0.4; set.seed(477210)   # Results reported in Figure 5
sd.ME <- 0.25; sd.ME.SIMEX <- 0.2; set.seed(5989960)  # Results reported in Figure S5

n.sims <- 300

N <- 1000
incidence.rate <- 0.6
fupt <- 296

n.sims.byrho <- 50
rho <- c(0.25, 0.5, 1, 1.5, 2, 2.5, 3)

# True NL function
nlx5<-function(x){
  return((x+0.5)^2/6)
}

SBPdata$sbp.ideal.nl <- nlx5(SBPdata$sbp.ideal)

ref.value.NL <- 0

### Declaration of objects to be saved ###

mean.sbp.noME <- rep(NA,n.sims)
sd.sbp.noME <- rep(NA,n.sims)
mean.sbp.ME <- rep(NA,n.sims)
sd.sbp.ME <- rep(NA, n.sims)

# Sparse, true model with TEL, without ME
coefA2sM2NL.noME <- list()
knotsA2sM2NL.noME <- list()
coefA2sM2TEL.noME <- list()
knotsA2sM2TEL.noME <- list()
devA2sM2.noME <- rep(NA, n.sims)

# Sparse, true model with TEL, with ME
coefA2sM2NL.ME <- list()
knotsA2sM2NL.ME <- list()
coefA2sM2TEL.ME <- list()
knotsA2sM2TEL.ME <- list()
devA2sM2.ME <- rep(NA, n.sims)

# Sparse, SIMEX
coefA2sM2.SIMEXcox <- list()
coefA2sM2.SIMEXlin <- list()
coefA2sM2.SIMEXlin.m1 <- rep(NA, n.sims)

### Simulation loop ###

options(warn=1)
library(PermAlgo)

for ( i in 1:n.sims){
  cat("i=",i,'\n')

  ### Data generation ###
    
  eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/fupt))
  censortimes <- rep(fupt, times=N)

  datA2 <- permalgorithm(numSubjects=N, maxTime=fupt, 
                         Xmat=as.matrix(SBPdata[,c('sbp.ideal.nl','sbp.sparse','tel')]), 
                         XmatNames=c('sbp.ideal.nl','sbp.sparse','tel'),
                         eventRandom=eventtimes, censorRandom=censortimes,
                         betas=c(1,0,0))
  mean.sbp.noME[i] <- mean(datA2$sbp.sparse)
  sd.sbp.noME[i] <- sd(datA2$sbp.sparse)
  
  ### Add ME ~ N(0, sd.ME) ###
  
  temp <- datA2[datA2$tel==0, c('Id','Fup','Stop','sbp.sparse','tel')]
  temp$sbp.sparse.ME <- temp$sbp.sparse + rnorm(nrow(temp), 0, sd.ME)
  rep.vec <- c(diff(temp$Stop), 0)     
  rep.vec[rep.vec<=0] <- temp$Fup[rep.vec<=0] - temp$Stop[rep.vec<=0] + 1
  
  datA2$sbp.sparse.ME <- rep(temp$sbp.sparse.ME, times=rep.vec)
  mean.sbp.ME[i] <- mean(datA2$sbp.sparse.ME)
  sd.sbp.ME[i] <- sd(datA2$sbp.sparse.ME)   
  
  ### Sparse ###
  
  # True model with TEL on data without ME
  modA2sM2.noME <- last_prog_TEL(data=datA2, Type=c("Start","Stop","Event","tel"),
                         variables=c("sbp.sparse"), TD=c(0), NL=c(1), TEL=c(1), m=1, p=2,
                         knots=-999)
  coefA2sM2NL.noME[[i]] <- modA2sM2.noME$coefficients_splines_NL
  knotsA2sM2NL.noME[[i]] <- modA2sM2.noME$knots_covariates
  coefA2sM2TEL.noME[[i]] <- modA2sM2.noME$coefficients_splines_TEL
  knotsA2sM2TEL.noME[[i]] <- modA2sM2.noME$knots_TEL
  devA2sM2.noME[i] <- -2*modA2sM2.noME$Partial_Log_Likelihood

  # True model with TEL on data with ME
  modA2sM2.ME <- last_prog_TEL(data=datA2, Type=c("Start","Stop","Event","tel"),
                            variables=c("sbp.sparse.ME"), TD=c(0), NL=c(1), TEL=c(1), m=1, p=2,
                            knots=-999)
  coefA2sM2NL.ME[[i]] <- modA2sM2.ME$coefficients_splines_NL
  knotsA2sM2NL.ME[[i]] <- modA2sM2.ME$knots_covariates
  coefA2sM2TEL.ME[[i]] <- modA2sM2.ME$coefficients_splines_TEL
  knotsA2sM2TEL.ME[[i]] <- modA2sM2.ME$knots_TEL
  devA2sM2.ME[i] <- -2*modA2sM2.ME$Partial_Log_Likelihood
  
  ### Application of SIMEX ###
  
  nlest.modA2sM2.ME.fct <- function(x){
    NLestim.sim(x, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) - 
      NLestim.sim(ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])
  }
  
  coefA2sM2.SIMEXcox[[i]] <- matrix(NA, nrow=n.sims.byrho, ncol=length(rho))
  temp.ME <- datA2[datA2$tel==0, c('Id','Fup','Stop','sbp.sparse.ME','tel')]
  
  for (j in 1:length(rho)){
    cat("   j=",j,"\n")
    for (k in 1:50){
      temp.ME$sbp.sparse.ME.SIMEX <- temp.ME$sbp.sparse.ME + rnorm(nrow(temp.ME), 0, (sd.ME.SIMEX*rho[j]))
      
      datA2$sbp.sparse.ME.SIMEX.NL <- nlest.modA2sM2.ME.fct(rep(temp.ME$sbp.sparse.ME.SIMEX, times=rep.vec))
      coefA2sM2.SIMEXcox[[i]][k,j] <- coxph(Surv(Start, Stop, Event) ~ sbp.sparse.ME.SIMEX.NL, data=datA2)$coef
    }
    gc()
  } 

  rm(modA2sM2.noME, modA2sM2.ME, temp, temp.ME, nlest.modA2sM2.ME.fct, datA2, eventtimes, censortimes); gc()  
  
}
gc()

### SIMEX extrapolation ###

x.rho <- rep(rho, each=n.sims.byrho)
x.rho2 <- x.rho*x.rho
for (i in 1:n.sims){
  coefA2sM2.SIMEXlin[[i]] <- lm(c(coefA2sM2.SIMEXcox[[i]]) ~ x.rho)$coef
  coefA2sM2.SIMEXlin.m1[i] <- coefA2sM2.SIMEXlin[[i]][1] + -1*coefA2sM2.SIMEXlin[[i]][2]
}

save(n.sims, N, incidence.rate, fupt, nlx5,
  sd.ME, sd.ME.SIMEX, n.sims.byrho, rho,
  mean.sbp.noME, sd.sbp.noME, mean.sbp.ME, sd.sbp.ME, 
  coefA2sM2NL.noME, knotsA2sM2NL.noME, coefA2sM2TEL.noME, knotsA2sM2TEL.noME, devA2sM2.noME, 
  coefA2sM2NL.ME, knotsA2sM2NL.ME, coefA2sM2TEL.ME, knotsA2sM2TEL.ME, devA2sM2.ME,
  coefA2sM2.SIMEXcox, coefA2sM2.SIMEXlin, coefA2sM2.SIMEXlin.m1,
  file=paste0("Scenario_A2_withME_sdME", sub("[:.:]", "", sd.ME), ".RData"))

