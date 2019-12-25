############################################################################################
##                                                                                        ##
## Code for Figures 1-5 in the manuscript and Figures S1-S5 in the Supporting Information ##
##                                                                                        ##
############################################################################################

# Last update: December 24, 2019

##### Preparation for all figures #####

# Set the current folder
setwd("C://path//to//your//directory")

source("Functions.R")

# X axis for time
axist <- seq(0, 275)

# X axis for X(u) values
axisx <- seq(-3, 4, length=100)

# X axis for TEL values
axistel <- seq(0, 42, length=100)

# Reference value for NL effect plots
ref.value.NL <- 0

# Load of simulation results from scenarios A1, A4, and A4 for cubic splines
load("Scenario_A1.RData")
load("Scenario_A2.RData")
load("Scenario_A3.RData")
load("Scenario_A4.RData")
load("Scenario_A4_cubic.RData")


########################################################################################################
## Figure 1 
## Ideal case (measurements updated at each time unit), NL and TD estimates from full model
## Scenarios A1 and A4

#pdf("Figure 1.pdf", width=9, height=9)
png("Figure 1.png", units="in", width=9, height=9, res=600)

par(mfrow=c(2,2))

### Panel (a): TD estimates from full model in A1

# For scaling the NL effect 
delta <- rep(NA, n.sims)
for (i in 1:n.sims){
  delta[i] <- NLestim.sim(x=2, m=1, p=2, coefA1iM8NL[[i]], knotsA1iM8NL[[i]]) - 
                NLestim.sim(x=-2, m=1, p=2, coefA1iM8NL[[i]], knotsA1iM8NL[[i]])
}
plot(axist, rep(log(1.5), times=length(axist)), xaxt="n", xlab="Follow-up time (years)", ylab="log(HR)", 
     type="l", ylim=c(-0.5,2))
axis(1, at=c(0,60,120,180,240), labels=c("0","10","20","30","40"))
title("(a) TD estimates in scenario A1")
for (i in 1:n.sims){
  lines(axist, TDestim.sim(axist, m=1, p=2, coefA1iM8TD[[i]], knotsA1iM8TD[[i]]) * (delta[i]/4), 
        col="grey", xlab="time", ylab="log(HR)", type="l")
}
A1iM8TD <- matrix(NA, nrow=n.sims, ncol=length(axist))
for(i in 1:n.sims){
  A1iM8TD[i,] <- TDestim.sim(axist, m=1, p=2, coefA1iM8TD[[i]], knotsA1iM8TD[[i]]) * (delta[i]/4)
}
meanA1iM8TD <- apply(A1iM8TD,2,mean)
lines(axist, meanA1iM8TD, col="white", lwd=5)
lines(axist, rep(log(1.5), times=length(axist)), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate","true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (b): NL estimates from full model in A1

plot(axisx, axisx*log(1.5), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-3,4))
title("(b) NL estimates in scenario A1")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA1iM8NL[[i]], knotsA1iM8NL[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA1iM8NL[[i]], knotsA1iM8NL[[i]])), 
        col="grey")
}
A1iM8NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A1iM8NL[i,] <- NLestim.sim(axisx, m=1, p=2, coefA1iM8NL[[i]], knotsA1iM8NL[[i]])-
    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA1iM8NL[[i]], knotsA1iM8NL[[i]])
}
meanA1iM8NL <- apply(A1iM8NL, 2, mean)
lines(axisx, meanA1iM8NL, col="white", lwd=5)
lines(axisx, axisx*log(1.5), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (c): TD estimates from full model in A4

plot(axist, tdx5(axist), xaxt="n", xlab="Follow-up time (years)", ylab="log(HR)", type="l", ylim=c(-0.5,2))
axis(1, at=c(0,60,120,180,240), labels=c("0","10","20","30","40"))
title("(c) TD estimates in scenario A4")
for (i in 1:n.sims){
  lines(axist, TDestim.sim(axist, m=1, p=2, coefA4iM8TD[[i]], knotsA4iM8TD[[i]]), col="grey", 
        xlab="time", ylab="log(HR)", type="l")
}
A4iM8TD <- matrix(NA, nrow=n.sims, ncol=length(axist))
for(i in 1:n.sims){
  A4iM8TD[i,] <- TDestim.sim(axist, m=1, p=2, coefA4iM8TD[[i]], knotsA4iM8TD[[i]])
}
meanA4iM8TD <- apply(A4iM8TD, 2, mean)
lines(axist, meanA4iM8TD, col="white", lwd=5)
lines(axist, tdx5(axist), col="black", lwd=5)
legend('bottomleft', inset=0.05, c("mean estimate", "individual estimate","true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (d): NL estimates from full model in A4

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-3,4))
title("(d) NL estimates in scenario A4")
for (i in 1:n.sims){
  lines(axisx, NLestim.sim(axisx, m=1, p=2, coefA4iM8NL[[i]], knotsA4iM8NL[[i]]) - 
          NLestim.sim(x=ref.value.NL, m=1, p=2, coefA4iM8NL[[i]], knotsA4iM8NL[[i]]), 
        col="grey")
}
A4iM8NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A4iM8NL[i,] <- NLestim.sim(axisx, m=1, p=2, coefA4iM8NL[[i]], knotsA4iM8NL[[i]]) -
    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA4iM8NL[[i]], knotsA4iM8NL[[i]])
}
meanA4iM8NL <- apply(A4iM8NL, 2, mean)
lines(axisx, meanA4iM8NL, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('bottomleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE,
       bg="gray90")

dev.off()


########################################################################################################
## Figure 2 
## Baseline only measurement, TD estimates from full model
## Scenarios A1, A4 

#pdf("Figure 2.pdf", width=11, height=5.5)
png("Figure 2.png", units="in", width=11, height=5.5, res=600)

par(mfrow=c(1,2))

### Panel (a): TD estimates from full model in A1

# For scaling the NL effect 
delta <- rep(NA, n.sims)
for (i in 1:n.sims){
  delta[i] <- NLestim.sim(x=2, m=1, p=2, coefA1bM8NL[[i]], knotsA1bM8NL[[i]]) - 
              NLestim.sim(x=-2, m=1, p=2, coefA1bM8NL[[i]], knotsA1bM8NL[[i]])
}
plot(axist, rep(log(1.5), times=length(axist)), xaxt="n", xlab="Follow-up time (years)", ylab="log(HR)", 
     type="l", ylim=c(-0.5,2))
axis(1, at=c(0,60,120,180,240), labels=c("0","10","20","30","40"))
title("(a) TD estimates in scenario A1")
for (i in 1:n.sims){
  lines(axist, TDestim.sim(axist, m=1, p=2, coefA1bM8TD[[i]], knotsA1bM8TD[[i]]) * (delta[i]/4), 
        col="grey")
}
A1bM8TD <- matrix(NA, nrow=n.sims, ncol=length(axist))
for(i in 1:n.sims){
  A1bM8TD[i,] <- TDestim.sim(axist, m=1, p=2, coefA1bM8TD[[i]], knotsA1bM8TD[[i]]) * (delta[i]/4)
}
meanA1bM8TD <- apply(A1bM8TD, 2, mean)
lines(axist, meanA1bM8TD, col="white", lwd=5)
lines(axist, rep(log(1.5), times=length(axist)), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate","true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (b): TD estimates from full model in A4

plot(axist, tdx5(axist), xaxt="n", xlab="Follow-up time (years)", ylab="log(HR)", type="l", ylim=c(-0.5,2))
axis(1, at=c(0,60,120,180,240), labels=c("0","10","20","30","40"))
title("(b) TD estimates in scenario A4")
for (i in 1:n.sims){
  lines(axist, TDestim.sim(axist, m=1, p=2, coefA4bM8TD[[i]], knotsA4bM8TD[[i]]), col="grey", 
        xlab="time", ylab="log(HR)", type="l")
}
A4bM8TD <- matrix(NA, nrow=n.sims, ncol=length(axist))
for(i in 1:n.sims){
  A4bM8TD[i,] <- TDestim.sim(axist, m=1, p=2, coefA4bM8TD[[i]], knotsA4bM8TD[[i]])
}
meanA4bM8TD <- apply(A4bM8TD, 2, mean)
lines(axist, meanA4bM8TD, col="white", lwd=5)
lines(axist, tdx5(axist), col="black", lwd=5)
legend('topright', inset=0.05, c("mean estimate", "individual estimate","true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

dev.off()


########################################################################################################
## Figure 3 
## Sparse measurements, NL estimates at fixed times u from full model
## Scenario A1

### Function to create panels (b)-(d) and (f)-(h)

fig.NLmTD.A1 <- function(panel1, panel2, u){
  
  # NL estimates at time u from full model without TEL model in A1
  
  plot(axisx, axisx*log(1.5), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
       type="l", ylim=c(-1.5,2.5))
  title(paste0("(", panel1, ") NL estimates at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA1sM8NL[[i]], knotsA1sM8NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA1sM8NL[[i]], knotsA1sM8NL[[i]])) *
                  TDestim.sim(x=u, m=1, p=2, coefA1sM8TD[[i]], knotsA1sM8TD[[i]]),
          col="grey")
  }
  A1sM8NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A1sM8NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA1sM8NL[[i]], knotsA1sM8NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2, coefA1sM8NL[[i]], knotsA1sM8NL[[i]])) *
                    TDestim.sim(x=u, m=1, p=2, coefA1sM8TD[[i]], knotsA1sM8TD[[i]])
  }
  meanA1sM8NL<-apply(A1sM8NL, 2, mean)
  lines(axisx, meanA1sM8NL, col="white", lwd=5)
  lines(axisx, axisx*log(1.5), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
    
  # NL estimates at time u from full model with TEL model in A1
  
  plot(axisx, axisx*log(1.5), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
       type="l", ylim=c(-1.5,2.5))
  title(paste0("(", panel2, ") NL estimates TEL-corrected at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA1sM4NL[[i]], knotsA1sM4NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA1sM4NL[[i]], knotsA1sM4NL[[i]])) * 
                  coefA1sM4TEL[[i]][1,] *
                  TDestim.sim(x=u, m=1, p=2, coefA1sM4TD[[i]], knotsA1sM4TD[[i]]), 
          col="grey")
  }
  A1sM4NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A1sM4NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA1sM4NL[[i]], knotsA1sM4NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA1sM4NL[[i]], knotsA1sM4NL[[i]])) * 
                    coefA1sM4TEL[[i]][1,] *
                    TDestim.sim(x=u, m=1, p=2, coefA1sM4TD[[i]], knotsA1sM4TD[[i]])
  }
  meanA1sM4NL <- apply(A1sM4NL, 2, mean)
  lines(axisx, meanA1sM4NL, col="white", lwd=5)
  lines(axisx, axisx*log(1.5), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
}

#pdf("Figure 3.pdf", width=13, height=6.5)
png("Figure 3.png", units="in", width=13, height=6.5, res=600)

par(mfcol=c(2,4))

### Panel (a): TD effect from full model in A1

# Scaling the NL effect 
delta <- rep(NA, n.sims)
for (i in 1:n.sims){
  delta[i] <- NLestim.sim(x=2, m=1, p=2, coefA1sM8NL[[i]], knotsA1sM8NL[[i]]) - 
              NLestim.sim(x=-2, m=1, p=2, coefA1sM8NL[[i]], knotsA1sM8NL[[i]])
}
plot(axist, rep(log(1.5), times=length(axist)), xlab="Follow-up time (years)", ylab="log(HR)", 
     type="l", ylim=c(-0.5,2))
title("(a) TD estimates")
for (i in 1:n.sims){
  lines(axist, TDestim.sim(axist, m=1, p=2, coefA1sM8TD[[i]], knotsA1sM8TD[[i]]) * (delta[i]/4), 
        col="grey")
}
A1sM8TD <- matrix(NA, nrow=n.sims, ncol=length(axist))
for(i in 1:n.sims){
  A1sM8TD[i,] <- TDestim.sim(axist, m=1, p=2, coefA1sM8TD[[i]], knotsA1sM8TD[[i]]) * (delta[i]/4)
}
meanA1sM8TD <- apply(A1sM8TD, 2, mean)
lines(axist, meanA1sM8TD, col="white", lwd=5)
lines(axist, rep(log(1.5), times=length(axist)), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate","true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (e): scaled TD effect from full model with TEL in A1

# For scaling the NL effect 
delta <- rep(NA, n.sims)
for (i in 1:n.sims){
  delta[i] <- NLestim.sim(x=2, m=1, p=2, coefA1sM4NL[[i]], knotsA1sM4NL[[i]]) - 
              NLestim.sim(x=-2, m=1, p=2, coefA1sM4NL[[i]], knotsA1sM4NL[[i]])
}
plot(axist, rep(log(1.5), times=length(axist)), xlab="Follow-up time (years)", ylab="log(HR)", 
     type="l", ylim=c(-0.5,2))
title("(e) TD estimates TEL-corrected")
for (i in 1:n.sims){
  lines(axist, TDestim.sim(axist, m=1, p=2, coefA1sM4TD[[i]], knotsA1sM4TD[[i]]) * (delta[i]/4) *
          coefA1sM4TEL[[i]][1,], 
        col="grey")
}
A1sM4TD <- matrix(NA, nrow=n.sims, ncol=length(axist))
for(i in 1:n.sims){
  A1sM4TD[i,] <- TDestim.sim(axist, m=1, p=2, coefA1sM4TD[[i]], knotsA1sM4TD[[i]]) * (delta[i]/4) *
    coefA1sM4TEL[[i]][1,]
}
meanA1sM4TD <- apply(A1sM4TD, 2, mean)
lines(axist, meanA1sM4TD, col="white", lwd=5)
lines(axist, rep(log(1.5), times=length(axist)), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate","true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panels (b)-(d) and (f)-(h)

fig.NLmTD.A1(panel1="b", panel2="f", u=30)  # u=5 years
fig.NLmTD.A1(panel1="c", panel2="g", u=90)  # u=15 years
fig.NLmTD.A1(panel1="d", panel2="h", u=210) # u=35 years

dev.off()


########################################################################################################
## Figure 4 
## Sparse measurements, NL estimates at fixed times u from full model
## Scenario A4

### Function to create panels (a)-(h)

fig.NLmTD.A4 <- function(panel1, panel2, u){
  
  # NL estimates at time u from full model without TEL model in A4
  
  plot(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), xlab=expression(italic(X) * '(' * italic(u) * ')'), 
       ylab="log(HR)", type="l", ylim=c(-1,4.5))
  title(paste0("(", panel1, ") NL estimates at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA4sM8NL[[i]], knotsA4sM8NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA4sM8NL[[i]], knotsA4sM8NL[[i]])) *
                  TDestim.sim(x=u, m=1, p=2, coefA4sM8TD[[i]], knotsA4sM8TD[[i]]),
        col="grey")
  }
  A4sM8NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A4sM8NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA4sM8NL[[i]], knotsA4sM8NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA4sM8NL[[i]], knotsA4sM8NL[[i]])) *
                    TDestim.sim(x=u, m=1, p=2, coefA4sM8TD[[i]], knotsA4sM8TD[[i]])
  }
  meanA4sM8NL<-apply(A4sM8NL, 2, mean)
  lines(axisx, meanA4sM8NL, col="white", lwd=5)
  lines(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

  # NL estimates at time u from full model with TEL model in A4
  
  plot(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), xlab=expression(italic(X) * '(' * italic(u) * ')'), 
       ylab="log(HR)", type="l", ylim=c(-1,4.5))
  title(paste0("(", panel2, ") NL estimates TEL-corrected at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA4sM4NL[[i]], knotsA4sM4NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA4sM4NL[[i]], knotsA4sM4NL[[i]])) * 
                  coefA4sM4TEL[[i]][1,] *
                  TDestim.sim(x=u, m=1, p=2, coefA4sM4TD[[i]], knotsA4sM4TD[[i]]), 
        col="grey")
  }
  A4sM4NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A4sM4NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA4sM4NL[[i]], knotsA4sM4NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA4sM4NL[[i]], knotsA4sM4NL[[i]])) * 
                    coefA4sM4TEL[[i]][1,] *
                    TDestim.sim(x=u, m=1, p=2, coefA4sM4TD[[i]], knotsA4sM4TD[[i]])
  }
  meanA4sM4NL <- apply(A4sM4NL, 2, mean)
  lines(axisx, meanA4sM4NL, col="white", lwd=5)
  lines(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
}

#pdf("Figure 4.pdf", width=13, height=6.5)
png("Figure 4.png", units="in", width=13, height=6.5, res=600)

par(mfcol=c(2,4))

fig.NLmTD.A4(panel1="a", panel2="e", u=30)  # u=5 years
fig.NLmTD.A4(panel1="b", panel2="f", u=90)  # u=15 years
fig.NLmTD.A4(panel1="c", panel2="g", u=150) # u=25 years
fig.NLmTD.A4(panel1="d", panel2="h", u=210) # u=35 years

dev.off()


########################################################################################################
## Figure 5 
## Sparse measurements, with measurement errors ~ N(0, sd=0.5), NL estimates from true model, 
##   with SIMEX correction
## Scenario A2

load("Scenario_A2_withME_sdME05.RData")

#pdf(file="Figure 5.pdf", height=9, width=9)
png(file="Figure 5.png", units="in", height=9, width=9, res=600)

par(mfrow=c(2,2))

### Panel (a): NL estimates from true model with TEL model in A2, without ME added

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-1,4))
title("(a) NL estimates TEL-corrected when\n no measurement errors were added")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]])) * 
                coefA2sM2TEL.noME[[i]][1,], 
        col="grey")
}
A2sM2NL.noME <- matrix(NA, nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A2sM2NL.noME[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]])-
                         NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]])) * 
                       coefA2sM2TEL.noME[[i]][1,]
}
meanA2sM2NL.noME <- apply(A2sM2NL.noME, 2, mean)
lines(axisx, meanA2sM2NL.noME, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (b): NL estimates from full with TEL model in A2, with ME added

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-1,4))
title("(b) NL estimates TEL-corrected when\n measurement errors were added")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                  coefA2sM2TEL.ME[[i]][1,], 
        col="grey")
}
A2sM2NL.ME <- matrix(NA, nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A2sM2NL.ME[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) -
                       NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                      coefA2sM2TEL.ME[[i]][1,]
}
meanA2sM2NL.ME <- apply(A2sM2NL.ME, 2, mean)
lines(axisx, meanA2sM2NL.ME, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (c): SIMEX extrapolation in A2 with ME added

mean.beta.byrho <- matrix(NA, nrow=n.sims, ncol=length(rho))
for (i in 1:n.sims){
  mean.beta.byrho[i,] <- apply(coefA2sM2.SIMEXcox[[i]], 2, mean)
}
plot(c(-1,3), c(0.2, 1.4), type="n", xlab=expression(rho), ylab='')
title("(c) SIMEX extrapolation when\n measurement errors were added")
mtext(expression(hat(eta)), side=2, line=2.5, cex=0.8)
for (i in 1:n.sims){
  points(rho, mean.beta.byrho[i,], col="grey")
}
points(rho, apply(mean.beta.byrho, 2, mean), pch=4, lty=2, lwd=2)
points(-1, mean(coefA2sM2.SIMEXlin.m1[1:n.sims]), pch=4, lty=2, lwd=2)

### Panel (d): NL estimates from true with TEL model in A2, with ME added and SIMEX correction

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), 
     ylab="log(HR)", type="l", ylim=c(-1,4))
title("(d) NL estimates TEL- and SIMEX-corrected\n when measurement errors were added")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                coefA2sM2TEL.ME[[i]][1,] * coefA2sM2.SIMEXlin.m1[i], 
        col="grey")
}
A2sM2NL.ME.SIMEX <- matrix(NA, nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A2sM2NL.ME.SIMEX[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) -
                             NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                            coefA2sM2TEL.ME[[i]][1,] * coefA2sM2.SIMEXlin.m1[i]
}
meanA2sM2NL.ME.SIMEX <- apply(A2sM2NL.ME.SIMEX, 2, mean)
lines(axisx, meanA2sM2NL.ME.SIMEX, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

dev.off()

rm(sd.ME, sd.ME.SIMEX, n.sims.byrho, rho,
   mean.sbp.noME, sd.sbp.noME, mean.sbp.ME, sd.sbp.ME, 
   coefA2sM2NL.noME, knotsA2sM2NL.noME, coefA2sM2TEL.noME, knotsA2sM2TEL.noME, devA2sM2.noME, 
   coefA2sM2NL.ME, knotsA2sM2NL.ME, coefA2sM2TEL.ME, knotsA2sM2TEL.ME, devA2sM2.ME,
   coefA2sM2.SIMEXcox, coefA2sM2.SIMEXlin, coefA2sM2.SIMEXlin.m1,
   mean.beta.byrho, A2sM2NL.ME, A2sM2NL.ME.SIMEX, A2sM2NL.noME,
   meanA2sM2NL.ME, meanA2sM2NL.ME.SIMEX, meanA2sM2NL.noME)


########################################################################################################
## Figure S1 
## Sparse measurements, TEL estimates from full model
## Scenarios A1, A2, A3 and A4

#pdf("Figure S1.pdf", width=9, height=9)
png("Figure S1.png", units="in", width=9, height=9, res=600)

par(mfrow=c(2,2))

### Panel (a): TEL effect from full with TEL model in A1

plot(c(0,42), c(0,2), xaxt="n", xlab="Time elapsed since last observation (years)", ylab="log(HR)", type="n")
axis(1, at=c(0,12,24,36), labels=c("0","2","4","6"))
title("(a) TEL estimates in scenario A1")
for (i in 1:n.sims){
  lines(axistel, TELestim.sim(axistel, m=1, p=2, coefA1sM4TEL[[i]], knotsA1sM4TEL[[i]]), col="grey")
}
A1sM4TEL <- matrix(NA,nrow=n.sims,ncol=length(axisx))
for(i in 1:n.sims){
  A1sM4TEL[i,] <- TELestim.sim(axistel, m=1, p=2, coefA1sM4TEL[[i]], knotsA1sM4TEL[[i]])
}
meanA1sM4TEL<-apply(A1sM4TEL, 2, mean)
lines(axistel, meanA1sM4TEL, col="white", lwd=5)
legend('top', inset=0.05, c("mean estimate", "individual estimate"), col=c("white","grey"),
       lty=c(1,1), lwd=c(2,2), merge=TRUE,
       bg="gray90")

### Panel (b): TEL effect from full with TEL model in A2

plot(c(0,42), c(0,2), xaxt="n", xlab="Time elapsed since last observation (years)", ylab="log(HR)", type="n")
axis(1, at=c(0,12,24,36), labels=c("0","2","4","6"))
title("(b) TEL estimates in scenario A2")
for (i in 1:n.sims){
  lines(axistel, TELestim.sim(axistel, m=1, p=2, coefA2sM4TEL[[i]], knotsA2sM4TEL[[i]]), col="grey")
}
A2sM4TEL <- matrix(NA,nrow=n.sims,ncol=length(axisx))
for(i in 1:n.sims){
  A2sM4TEL[i,] <- TELestim.sim(axistel, m=1, p=2, coefA2sM4TEL[[i]], knotsA2sM4TEL[[i]])
}
meanA2sM4TEL<-apply(A2sM4TEL, 2, mean)
lines(axistel, meanA2sM4TEL, col="white", lwd=5)
legend('top', inset=0.05, c("mean estimate", "individual estimate"), col=c("white","grey"),
       lty=c(1,1), lwd=c(2,2), merge=TRUE,
       bg="gray90")

### Panel (c): TEL effect from full with TEL model in A3

plot(c(0,42), c(0,2), xaxt="n", xlab="Time elapsed since last observation (years)", ylab="log(HR)", type="n")
axis(1, at=c(0,12,24,36), labels=c("0","2","4","6"))
title("(c) TEL estimates in scenario A3")
for (i in 1:n.sims){
  lines(axistel, TELestim.sim(axistel, m=1, p=2, coefA3sM4TEL[[i]], knotsA3sM4TEL[[i]]), col="grey")
}
A3sM4TEL <- matrix(NA,nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A3sM4TEL[i,] <- TELestim.sim(axistel, m=1, p=2, coefA3sM4TEL[[i]], knotsA3sM4TEL[[i]])
}
meanA3sM4TEL<-apply(A3sM4TEL, 2, mean)
lines(axistel, meanA3sM4TEL, col="white", lwd=5)
legend('top', inset=0.05, c("mean estimate", "individual estimate"), col=c("white","grey"),
       lty=c(1,1), lwd=c(2,2), merge=TRUE, bg="gray90")

### Panel (d): TEL effect from full with TEL model in A4

plot(c(0,42), c(0,2), xaxt="n", xlab="Time elapsed since last observation (years)", ylab="log(HR)", type="n")
axis(1, at=c(0,12,24,36), labels=c("0","2","4","6"))
title("(d) TEL estimates in scenario A4")
for (i in 1:n.sims){
  lines(axistel, TELestim.sim(axistel, m=1, p=2, coefA4sM4TEL[[i]], knotsA4sM4TEL[[i]]), col="grey")
}
A4sM4TEL <- matrix(NA,nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A4sM4TEL[i,] <- TELestim.sim(axistel, m=1, p=2, coefA4sM4TEL[[i]], knotsA4sM4TEL[[i]])
}
meanA4sM4TEL<-apply(A4sM4TEL, 2, mean)
lines(axistel, meanA4sM4TEL, col="white", lwd=5)
legend('top', inset=0.05, c("mean estimate", "individual estimate"), col=c("white","grey"),
       lty=c(1,1), lwd=c(2,2), merge=TRUE, bg="gray90")

dev.off()


########################################################################################################
## Figure S2 
## Sparse measurements, NL estimates at fixed times u from full model
## Scenario A2

### Function to create panels (a)-(h)

fig.NLmTD.A2 <- function(panel1, panel2, u){
  
  # NL estimates at time u from full model without TEL model in A2
  
  plot(axisx, (nlx5(axisx)-nlx5(0)), xlab=expression(italic(X) * '(' * italic(u) * ')'), 
       ylab="log(HR)", type="l", ylim=c(-1,4.5))
  title(paste0("(", panel1, ") NL estimates at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM8NL[[i]], knotsA2sM8NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM8NL[[i]], knotsA2sM8NL[[i]])) *
                  TDestim.sim(u, m=1, p=2, coefA2sM8TD[[i]], knotsA2sM8TD[[i]]),
          col="grey")
  }
  A2sM8NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A2sM8NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM8NL[[i]], knotsA2sM8NL[[i]])-
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA2sM8NL[[i]], knotsA2sM8NL[[i]])) *
                    TDestim.sim(x=u, m=1, p=2, coefA2sM8TD[[i]], knotsA2sM8TD[[i]])
  }
  meanA2sM8NL<-apply(A2sM8NL, 2, mean)
  lines(axisx, meanA2sM8NL, col="white", lwd=5)
  lines(axisx, (nlx5(axisx)-nlx5(0)), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
  
  # NL estimates at time u from full model with TEL model in A2
  
  plot(axisx, (nlx5(axisx)-nlx5(0)), xlab=expression(italic(X) * '(' * italic(u) * ')'), 
       ylab="log(HR)", type="l", ylim=c(-1,4.5))
  title(paste0("(", panel2, ") NL estimates TEL-corrected at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM4NL[[i]], knotsA2sM4NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM4NL[[i]], knotsA2sM4NL[[i]])) *
                  coefA2sM4TEL[[i]][1,] * 
                  TDestim.sim(x=u, m=1, p=2, coefA2sM4TD[[i]], knotsA2sM4TD[[i]]), 
          col="grey")
  }
  A2sM4NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A2sM4NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM4NL[[i]], knotsA2sM4NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA2sM4NL[[i]], knotsA2sM4NL[[i]])) * 
                      coefA2sM4TEL[[i]][1,] *
                      TDestim.sim(x=u, m=1, p=2, coefA2sM4TD[[i]], knotsA2sM4TD[[i]])
  }
  meanA2sM4NL <- apply(A2sM4NL, 2, mean)
  lines(axisx, meanA2sM4NL, col="white", lwd=5)
  lines(axisx, (nlx5(axisx)-nlx5(0)), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
}

#pdf("Figure S2.pdf", width=13, height=6.5)
png("Figure S2.png", units="in", width=13, height=6.5, res=600)

par(mfcol=c(2,4))

fig.NLmTD.A2(panel1="a", panel2="e", u=30)  # u=5 years
fig.NLmTD.A2(panel1="b", panel2="f", u=90)  # u=15 years
fig.NLmTD.A2(panel1="c", panel2="g", u=150) # u=25 years
fig.NLmTD.A2(panel1="d", panel2="h", u=210) # u=35 years

dev.off()


########################################################################################################
## Figure S3 
## Sparse measurements, NL estimates at fixed times u from full model
## Scenario A3

### Function to create panels (a)-(h)

fig.NLmTD.A3 <- function(panel1, panel2, u){
  
  # NL estimates at time u from full model without TEL model in A3
  
  plot(axisx, tdx5(u)*axisx, xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
       type="l", ylim=c(-5,7))
  title(paste0("(", panel1, ") NL estimates at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA3sM8NL[[i]], knotsA3sM8NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA3sM8NL[[i]], knotsA3sM8NL[[i]])) *
                  TDestim.sim(x=u, m=1, p=2, coefA3sM8TD[[i]], knotsA3sM8TD[[i]]),
          col="grey")
  }
  A3sM8NL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A3sM8NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA3sM8NL[[i]], knotsA3sM8NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA3sM8NL[[i]], knotsA3sM8NL[[i]])) *
                    TDestim.sim(x=u, m=1, p=2, coefA3sM8TD[[i]], knotsA3sM8TD[[i]])
  }
  meanA3sM8NL<-apply(A3sM8NL, 2, mean)
  lines(axisx, meanA3sM8NL, col="white", lwd=5)
  lines(axisx, tdx5(u)*axisx, col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"),
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
  
  # NL estimates at time u from full model with TEL model in A3
  
  plot(axisx, tdx5(u)*axisx, xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)",
       type="l", ylim=c(-5,7))
  title(paste0("(", panel2, ") NL estimates TEL-corrected at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA3sM4NL[[i]], knotsA3sM4NL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=2, coefA3sM4NL[[i]], knotsA3sM4NL[[i]])) * 
                  coefA3sM4TEL[[i]][1,] *
                  TDestim.sim(x=u, m=1, p=2, coefA3sM4TD[[i]], knotsA3sM4TD[[i]]), 
          col="grey")
  }
  A3sM4NL <- matrix(NA,nrow=n.sims,ncol=length(axisx))
  for(i in 1:n.sims){
    A3sM4NL[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA3sM4NL[[i]], knotsA3sM4NL[[i]]) -
                      NLestim.sim(x=ref.value.NL, m=1, p=2,coefA3sM4NL[[i]], knotsA3sM4NL[[i]])) *
                    coefA3sM4TEL[[i]][1,] *
                    TDestim.sim(x=u, m=1, p=2, coefA3sM4TD[[i]], knotsA3sM4TD[[i]])
  }
  meanA3sM4NL <- apply(A3sM4NL, 2, mean)
  lines(axisx, meanA3sM4NL, col="white", lwd=5)
  lines(axisx, tdx5(u)*axisx, col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
}

#pdf("Figure S3.pdf", width=13, height=6.5)
png("Figure S3.png", units="in", width=13, height=6.5, res=600)

par(mfcol=c(2,4))

fig.NLmTD.A3(panel1="a", panel2="e", u=30)
fig.NLmTD.A3(panel1="b", panel2="f", u=90)
fig.NLmTD.A3(panel1="c", panel2="g", u=150)
fig.NLmTD.A3(panel1="d", panel2="h", u=210)

dev.off()


########################################################################################################
## Figure S4 
## Cubic splines, Sparse measurements, NL estimates at fixed times u from full model
## Scenario A4

### Function to create panels (a)-(h)

fig.NLmTD.A4c <- function(panel1, panel2, u){
  
  # NL estimates at time u from full model without TEL model in A4
  
  plot(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), xlab=expression(italic(X) * '(' * italic(u) * ')'), 
       ylab="log(HR)", type="l", ylim=c(-1,4.5))
  title(paste0("(", panel1, ") NL estimates at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=3, coefA4sM8cNL[[i]], knotsA4sM8cNL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=3, coefA4sM8cNL[[i]], knotsA4sM8cNL[[i]])) *
                  TDestim.sim(x=u, m=1, p=3, coefA4sM8cTD[[i]], knotsA4sM8cTD[[i]]),
          col="grey")
  }
  A4sM8cNL <- matrix(NA, nrow=n.sims, ncol=length(axisx))
  for(i in 1:n.sims){
    A4sM8cNL[i,] <- (NLestim.sim(axisx, m=1, p=3, coefA4sM8cNL[[i]], knotsA4sM8cNL[[i]]) -
                       NLestim.sim(x=ref.value.NL, m=1, p=3,coefA4sM8cNL[[i]], knotsA4sM8cNL[[i]])) *
                     TDestim.sim(x=u, m=1, p=3, coefA4sM8cTD[[i]], knotsA4sM8cTD[[i]])
  }
  meanA4sM8cNL <- apply(A4sM8cNL, 2, mean)
  lines(axisx, meanA4sM8cNL, col="white", lwd=5)
  lines(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
  
  # NL estimates at time u from full model with TEL model in A4
  
  plot(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
       type="l", ylim=c(-1,4.5))
  title(paste0("(", panel2, ") NL estimates TEL-corrected at ", u/6, " years"))
  for (i in 1:n.sims){
    lines(axisx, (NLestim.sim(axisx, m=1, p=3, coefA4sM4cNL[[i]], knotsA4sM4cNL[[i]]) - 
                    NLestim.sim(x=ref.value.NL, m=1, p=3, coefA4sM4cNL[[i]], knotsA4sM4cNL[[i]])) * 
                  coefA4sM4cTEL[[i]][1,] *
                  TDestim.sim(x=u, m=1, p=3, coefA4sM4cTD[[i]], knotsA4sM4cTD[[i]]), 
          col="grey")
  }
  A4sM4cNL <- matrix(NA,nrow=n.sims,ncol=length(axisx))
  for(i in 1:n.sims){
    A4sM4cNL[i,] <- (NLestim.sim(axisx, m=1, p=3, coefA4sM4cNL[[i]], knotsA4sM4cNL[[i]]) -
                       NLestim.sim(x=ref.value.NL, m=1, p=3, coefA4sM4cNL[[i]], knotsA4sM4cNL[[i]])) * 
                     coefA4sM4cTEL[[i]][1,] *
                     TDestim.sim(x=u, m=1, p=3, coefA4sM4cTD[[i]], knotsA4sM4cTD[[i]])
  }
  meanA4sM4cNL <- apply(A4sM4cNL, 2, mean)
  lines(axisx, meanA4sM4cNL, col="white", lwd=5)
  lines(axisx, (nlx5(axisx)-nlx5(0))*tdx5(u), col="black", lwd=5)
  legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
         col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")
}

#pdf("Figure S4.pdf", width=13, height=6.5)
png("Figure S4.png", units="in", width=13, height=6.5, res=600)

par(mfcol=c(2,4))

fig.NLmTD.A4c(panel1="a", panel2="e", u=30)  # u=5 years
fig.NLmTD.A4c(panel1="b", panel2="f", u=90)  # u=15 years
fig.NLmTD.A4c(panel1="c", panel2="g", u=150) # u=25 years
fig.NLmTD.A4c(panel1="d", panel2="h", u=210) # u=35 years

dev.off()


########################################################################################################
## Figure S5 
## Sparse measurements, with measurement error ~ N(0, sd=0.25), NL estimates from true model and 
##   SIMEX correction
## Scenario A2

load("Scenario_A2_withME_sdME025.RData")

#pdf(file="Figure S5.pdf", height=9, width=9)
png(file="Figure S5.png", units="in", height=9, width=9, res=600)

par(mfrow=c(2,2))

### Panel (a): NL estimates from true model with TEL model in A2, without ME added

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-1,4))
title("(a) NL estimates TEL-corrected when\n no measurement errors weres added")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]])) *
                coefA2sM2TEL.noME[[i]][1,], col="grey")
}
A2sM2NL.noME <- matrix(NA, nrow=n.sims, ncol=length(axisx))
for(i in 1:n.sims){
  A2sM2NL.noME[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]]) -
                         NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.noME[[i]], knotsA2sM2NL.noME[[i]])) *
                       coefA2sM2TEL.noME[[i]][1,]
}
meanA2sM2NL.noME <- apply(A2sM2NL.noME, 2, mean)
lines(axisx, meanA2sM2NL.noME, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"),
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (b): NL estimates from full with TEL model in A2, with ME added

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-1,4))
title("(b) NL estimates TEL-corrected when\n measurement errors were added")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                coefA2sM2TEL.ME[[i]][1,], col="grey")
}
A2sM2NL.ME <- matrix(NA,nrow=n.sims,ncol=length(axisx))
for(i in 1:n.sims){
  A2sM2NL.ME[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) -
                       NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                     coefA2sM2TEL.ME[[i]][1,]
}
meanA2sM2NL.ME <- apply(A2sM2NL.ME, 2, mean)
lines(axisx, meanA2sM2NL.ME, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

### Panel (c): SIMEX extrapolation in A2 with ME added

mean.beta.byrho <- matrix(NA, nrow=n.sims, ncol=length(rho))
for (i in 1:n.sims){
  mean.beta.byrho[i,] <- apply(coefA2sM2.SIMEXcox[[i]], 2, mean)
}
plot(c(-1,3), c(0.2, 1.4), type="n", xlab=expression(rho), ylab='')
title("(c) SIMEX extrapolation when\n measurement errors were added")
mtext(expression(hat(eta)), side=2, line=2.5, cex=0.8)
for (i in 1:n.sims){
  points(rho, mean.beta.byrho[i,], col="grey")
}
points(rho, apply(mean.beta.byrho, 2, mean), pch=4, lty=2, lwd=2)
points(-1, mean(coefA2sM2.SIMEXlin.m1[1:n.sims]), pch=4, lty=2, lwd=2)

### Panel (d): NL estimates from true with TEL model in A2, with ME added and SIMEX correction

plot(axisx, nlx5(axisx)-nlx5(0), xlab=expression(italic(X) * '(' * italic(u) * ')'), ylab="log(HR)", 
     type="l", ylim=c(-1,4))
title("(d) NL estimates TEL- and SIMEX-corrected\n when measurement errors were added")
for (i in 1:n.sims){
  lines(axisx, (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) - 
                  NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                coefA2sM2TEL.ME[[i]][1,] * coefA2sM2.SIMEXlin.m1[i], col="grey")
}
A2sM2NL.ME.SIMEX <- matrix(NA,nrow=n.sims,ncol=length(axisx))
for(i in 1:n.sims){
  A2sM2NL.ME.SIMEX[i,] <- (NLestim.sim(axisx, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]]) -
                             NLestim.sim(x=ref.value.NL, m=1, p=2, coefA2sM2NL.ME[[i]], knotsA2sM2NL.ME[[i]])) *
                           coefA2sM2TEL.ME[[i]][1,] * coefA2sM2.SIMEXlin.m1[i]
}
meanA2sM2NL.ME.SIMEX <- apply(A2sM2NL.ME.SIMEX, 2, mean)
lines(axisx, meanA2sM2NL.ME.SIMEX, col="white", lwd=5)
lines(axisx, nlx5(axisx)-nlx5(0), col="black", lwd=5)
legend('topleft', inset=0.05, c("mean estimate", "individual estimate", "true effect"), 
       col=c("white","grey","black"), lty=c(1,1,1), lwd=c(2,2,2), merge=TRUE, bg="gray90")

dev.off()

rm(sd.ME, sd.ME.SIMEX, n.sims.byrho, rho,
   mean.sbp.noME, sd.sbp.noME, mean.sbp.ME, sd.sbp.ME, 
   coefA2sM2NL.noME, knotsA2sM2NL.noME, coefA2sM2TEL.noME, knotsA2sM2TEL.noME, devA2sM2.noME, 
   coefA2sM2NL.ME, knotsA2sM2NL.ME, coefA2sM2TEL.ME, knotsA2sM2TEL.ME, devA2sM2.ME,
   coefA2sM2.SIMEXcox, coefA2sM2.SIMEXlin, coefA2sM2.SIMEXlin.m1,
   mean.beta.byrho, A2sM2NL.ME, A2sM2NL.ME.SIMEX, A2sM2NL.noME,
   meanA2sM2NL.ME, meanA2sM2NL.ME.SIMEX, meanA2sM2NL.noME)
