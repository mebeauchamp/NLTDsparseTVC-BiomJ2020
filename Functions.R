#############################################################################
##                                                                         ##
## Functions to estimate non-linear (NL), time-dependent (NL), and time    ##
## elapsed since last observation (TEL) effects of time-varying covariates ##
## in survival analysis, i.e. models (1) and (2) of the manuscript         ##
##                                                                         ##
#############################################################################

# Last update: December 25, 2019

### Main functions:

  ## last_prog_TEL(data, Type, variables, TD, NL, TEL, m, p, knots=-999) ##

    # Function estimating model (2) of the manuscript
    # Arguments:
    #   data: data frame containing the data. The first column most be the variable identifying the subjects
    #   Type: a vector including the variable indicating the start and stop of time intervals, the event, and the 
    #         TEL variables in data, e.g. c("Start","Stop","Event","tel")
    #   variables: a vector containing all independent variables in the model
    #   TD: a vector indicating for each variable if the TD effect is modeled (0/1)
    #   NL: a vector indicating for each variable if the NL effect is modeled (0/1).
    #       Can be 1 only for continuous variables.
    #   TEL: a vector indicating for each variable if the TEL effect is modeled (0/1)
    #   m: number of interior knots (the same for all TD, NL, and TEL effects)
    #   p: degree of the splines (the same for all TD, NL, and TEL effects)
    #   knots: position of interior knots. -999 indicates that the knots are automatically allocated.
    
  ## last_prog(data, Type, variables, TD, NL, m, p, knots=-999) ##

    # Function estimating model (1) of the manuscript
    # Arguments:
    #   data: data frame containing the data. The first column most be the variable identifying the subjects
    #   Type: a vector including the variable indicating the start and stop of time intervals, and the event
    #         variable in data, e.g. c("Start","Stop","Event")
    #   variables: a vector containing all independent variables in the model
    #   TD: a vector indicating for each variable if the TD effect is modeled (0/1)
    #   NL: a vector indicating for each variable if the NL effect is modeled (0/1).
    #       Can be 1 only for continuous variables.
    #   m: number of interior knots (the same for all TD, NL, and TEL effects)
    #   p: degree of the splines (the same for all TD, NL, and TEL effects)
    #   knots: position of interior knots. -999 indicates that the knots are automatically allocated.

### Functions to calculate the estimated TD, NL, and TEL effect, respectively, for specific values of 
###  one variable from the simulation results:

  ## TDestim.sim(x, m, p, coefs, knots)       ##
  ## NLestim.sim(x, model.NL, coefs, knots)   ##
  ## TELestim.sim(x, model.TEL, coefs, knots) ##

    # NOTE: the three functions are identical but seperate names were kept to facilitate understanding of the code

    # Arguments:
    #   x: value or vector of values for which the TD, NL, or TEL effect is estimated
    #   m: number of interior knots used to estimate the effect
    #   p: degree of the splines used to estimate the effect  
    #   coefs: estimated splines coefficients for TD, NL, or TEL effect of the variable  
    #   knots: knots used to estimate the effect 

### Internal functions:

  ## spli(x, j, p, knots) ## 

   # Function calulating the spline value for x, for one spline in the basis
   # Arguments:
   #   x: value
   #   j: number of the spline
   #   p: degree of splines
   #   knots: position of interior knots

  ## DvlpMatrix(data, listeT, ncol, TypeNEW) ##

   # Function to rearrange the data to recreate the risk sets. It deletes rows not in any risk sets to improve efficiency.
   # Arguments:
   #  data: the data
   #  listeT: list of unique event times in data
   #  ncol: number of columns in data
   #  TypeNEW: position of variables in argument Type in data


library(survival)
library(splines)


last_prog_TEL <- function (data, Type, variables, TD, NL, TEL, m, p, knots=-999){
   
  ### IF THERE ARE NO TEL EFFECTS ###
  if (sum(TEL)==0){ 
   
    last_prog(data=data, Type, variables, TD, NL, m, p, knots=-999)
    
  ### IF THERE ARE TEL EFFECTS ###
  } else {    

    i1 <- sum((NL+TD+TEL)==0)                 # number of non-nl non-td non-tel variables
    i2 <- sum(((NL==1) & (TD==0) & (TEL==0))) # number of NL non-td non-tel variables
    i3 <- sum(((NL==0) & (TD==1) & (TEL==0))) # number of non-nl TD non-tel variables
    i4 <- sum((NL+TD)==2 & (TEL==0))          # number of NL TD non-tel variables
    i5 <- sum(((NL==0) & (TD==0) & (TEL==1))) # number of non-nl non-td TEL variables
    i6 <- sum(((NL==1) & (TD==0) & (TEL==1))) # number of NL non-td TEL variables
    i7 <- sum(((NL==0) & (TD==1) & (TEL==1))) # number of non-nl TD TEL variables
    i8 <- sum(((NL==1) & (TD==1) & (TEL==1))) # number of NL TD TEL variables
    V <- length(variables) # number of variables
    
    nTEL <- sum(TEL) # number of variables having a TEL effect (and possibly other effects)
    
    variablesNEW <- match(variables, names(data)) # positions in data of items in argument variables 
    TypeNEW <- match(Type, names(data))           # positions in data of items in argument Type           

    ## Store interior and exterior knots in (2*(p+1)+m) columns,
     # for variables (rows 1:V), time (row V+1), and TEL effects (rows (V+2):(V+3+nTEL))

    knotsNEW <- matrix(nrow = V + 1+nTEL, ncol = p+1 + m + p+1) 
    listeprobaquantiles <- seq(1, m)/(m+1)     # quantiles for interior knots for variables and time
    
    listeprobaquantilesTEL <- rep(NA, times=m) # quantiles imposed for interior knots for TEL variables
    if (m>3 & nTEL>0) stop("Current version cannot accomodate m > 3 for TEL effects. Use m <= 3.")
    if (m==1){
      listeprobaquantilesTEL <- 0.25
    } else if (m==2){
      listeprobaquantilesTEL <- c(0.25, 0.5)
    } else if (m==3){
      listeprobaquantilesTEL <- c(0.25, 0.5, 0.75)
    }

    # If no user-defined knots
    if (is.matrix(knots) == FALSE){   
      if (is.numeric(knots)==TRUE & knots==-999){ # use default knots 

        # Knots for variables:
          # p+1 exterior knots at min(variables[i])
          # interior knots at quantiles
          # p+1 exterior knots equally spaced between max(variable) and max(variable)+p
        for (i in 1:V) {
          knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p+1), 
                             quantile(data[, variables[i]], probs=listeprobaquantiles), 
                             seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1)) 
        }
        
        # Knots for time:
          # p+1 exterior knots at 0
          # interior knots at quantiles of event times
          # p+1 exterior knots equally spaced between max(time) and max(time)+p
        knotsNEW[V+1, ] <- c(rep(0, p+1), 
                            quantile(data[data[, TypeNEW[3]]==1, TypeNEW[2]], probs=listeprobaquantiles),  
                            seq(max(data[data[,TypeNEW[3]]==1, TypeNEW[2]]), 
                                max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1)) 
 
        # Knots for TEL effects:
          # p+1 exterior knots at 0
          # interior knots at imposed quantiles
          # p+1 exterior knots equally spaced between max(TEL variable) and max(TEL variable)+p 
        for(i in 1:nTEL){
          knotsNEW[(V+1+i),] <- c(rep(0, p+1), 
                                quantile(data[, TypeNEW[3+i]], probs=listeprobaquantilesTEL), 
                                seq(max(data[, TypeNEW[3+i]]), max(data[,TypeNEW[3+i]])+p, 1)) 
        }
      }                                                                          
    }

    # If (some) user-defined interior knots are specified
    if (is.matrix(knots)==TRUE){  
      if (dim(knots)[1] != (length(variables)+2) | dim(knots)[2] != m) 
        stop("Error Message: argument knots should be a matrix of dimension (length(variables)+2) by m")

      # Knots for variables
      for (i in 1:V){
        if (is.na(knots[i,1])==TRUE){  # if user-defined knots for that variable is NA, use default
          knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p+1),   
                             quantile(data[, variables[i]], probs=listeprobaquantiles), 
                             seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
        } else {
          knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p+1), knots[i,], 
                             seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
        }
      }

      # Knots for time      
      if (is.na(knots[(V+1),1])==TRUE){  # if user-defined knots for time is NA, use default
        knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                               quantile(data[data[, TypeNEW[3]]==1, TypeNEW[2]], 
                                        probs=listeprobaquantiles), 
                               seq(max(data[data[,TypeNEW[3]]==1, TypeNEW[2]]), 
                                   max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1))
      } else {
        knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                               knots[(V+1),], 
                               seq(max(data[data[,TypeNEW[3]]==1, TypeNEW[2]]), 
                                   max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1))
      }
      
      # Knots for TEL effects
      if (is.na(knots[(V+2),1]) == TRUE){  # if user-defined knots for TEL effect is NA, use default  
        knotsNEW[(V+2), ] <- c(rep(0, p + 1), 
                               quantile(data[, TypeNEW[4]], probs=listeprobaquantilesTEL), 
                               seq(max(data[, TypeNEW[4]]), max(data[, TypeNEW[4]])+p, 1))
      } else {
        knotsNEW[(V+2), ] <- c(rep(0, p+1), 
                               knots[(V+2),], 
                               seq(max(data[, TypeNEW[4]]), max(data[, TypeNEW[4]])+p, 1))
      }
    }

    ## Prepare data (QWR) for estimation
    
    data <- as.matrix(data)
    QWR <- data  

    # Add splines for NL effects in QWR    
    nbNL <- sum(NL) # number of NL effects 
    if (nbNL != 0){ # generate spline values for each NL effect, with constraint a=0 for first spline
      for (i in 1:nbNL){
        QWR <- cbind(QWR, splineDesign(knotsNEW[seq(1,V, 1)[NL==1][i],], 
                                       x=QWR[, variablesNEW[NL==1][i]], ord=p+1)[,-1])
      }
    }

    # Add splines for TD effects in QWR (splines for time)
    QWR <- cbind(QWR, splineDesign(knotsNEW[V+1,], x=QWR[, TypeNEW[2]], ord=p+1, outer.ok=TRUE))
    
    # Add splines for TEL effects in QWR      
    for(j in 1:nTEL){
      QWR <- cbind(QWR, splineDesign(knotsNEW[(V+1+j),], x=QWR[, TypeNEW[3+j]], ord=p+1))
    }

    ## Construct equation (tt) of Cox model (part to be estimated at all steps of all iterations)
    
    tt <- paste("modX <- coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,",TypeNEW[2], "],QWR[,", TypeNEW[3], "]) ~ ",
                sep="")
    Nn <- 0  # number of parameters added to model up to now
 
    # Add parameters for variables without NL, TD and TEL effects   
    nbpara <- sum((NL==0 & TD==0 & TEL==0 )) 
    if (nbpara != 0){ 
      for (k in 1:nbpara){
        tt <- paste(tt, "QWR[,", variablesNEW[NL==0 & TD==0 & TEL==0][k], "]+", sep="")
        Nn <- Nn + 1
      }
    }

    # Add parameters for variables with NL effect only 
    nbonlyNL <- sum((NL==1 & TD==0 & TEL==0)) 
    if (nbonlyNL != 0){
      onlyNL <- match(variablesNEW[NL==1 & TD==0 & TEL==0], variablesNEW[NL==1]) 
      for (k in 1:nbonlyNL){
        covp <- paste("QWR[,", (dim(data)[2] + (onlyNL[k]-1)*(m+p)+1):(dim(data)[2] + onlyNL[k]*(m+p)), "]",
                      sep="")
        tt <- paste(tt, paste(c(covp,""), collapse="+"))
        Nn <- Nn+m+p
      }
    }

    # Add parameters for variables with TD effect only 
    nbonlyTD <- sum((NL==0 & TD==1 & TEL==0)) 
    if (nbonlyTD != 0){
      for (k in 1:nbonlyTD){
        flag <- dim(QWR)[2]+1
        QWR <- cbind(QWR, QWR[, variablesNEW[NL==0 & TD==1 & TEL==0][k]] *                  # variable
                       QWR[, (dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2] +(m+p)*(nbNL+1)+1)]) # splines for time
        covp <- paste("QWR[,", flag: dim(QWR)[2], "]", sep="")
        tt <- paste(tt, paste(c(covp,""), collapse="+"))  
        Nn <- Nn +m+p+1
      }
    }

    # Add parameters for variables with TEL effect only
    nbonlyTEL <- sum((NL==0 & TD==0 & TEL==1)) 
    if (nbonlyTEL != 0){
      onlyTEL <- match(variablesNEW[NL==0 & TD==0 & TEL==1], variablesNEW[TEL==1]) 
      for (k in 1:nbonlyTEL){
        flag <- dim(QWR)[2]+1
        QWR <- cbind(QWR, QWR[, variablesNEW[NL==0 & TD==0 & TEL==1][k]] *           # variable
                     QWR[, (dim(data)[2]+(m+p)*(nbNL+1)+1+(onlyTEL[k]-1)*(m+p+1)+1): # splines for TEL variable
                           (dim(data)[2] +(m+p)*(nbNL+1)+1+onlyTEL[k]*(1+m+p))])      
        covp <- paste("QWR[,", flag: dim(QWR)[2], "]", sep="")
        tt <- paste(tt, paste(c(covp,""), collapse="+"))  
        Nn <- Nn+1+m+p
      }
    }

    ## Number of variables with each combination of 2 or 3 effects and their position in argument 
      # "variables"among variables with a NL effect or TD effect
    
    nbNLTDTEL <- sum((NL==1 & TD==1 & TEL==1)) # nbr variables with NL, TD and TEL effects
    pNLTDTEL <- match(variablesNEW[NL==1 & TD==1 & TEL==1], variablesNEW[NL==1])  
    
    nbNLTD <- sum((NL==1 & TD==1 & TEL==0)) # nbr variables with NL and TD no TEL effects
    pNLTD <- match(variablesNEW[NL==1 & TD==1 & TEL==0], variablesNEW[NL==1])  
    
    nbNLTEL <- sum((NL==1 & TD==0 & TEL==1)) # nbr variables with NL and TEL no TD effects
    pNLTEL <- match(variablesNEW[NL==1 & TD==0 & TEL==1], variablesNEW[NL==1]) 
    
    nbTDTEL <- sum((NL==0 & TD==1 & TEL==1)) # nbr variables with TD and TEL no NL effects
    pTDTEL <- match(variablesNEW[NL==0 & TD==1 & TEL==1], variablesNEW[TD==1])  

    # Position in argument "variables" of variables with TEL effects
    posTELNLTD <- match(variablesNEW[NL==1 & TD==1 & TEL==1], variablesNEW[TEL==1])
    
    posTELTD <- match(variablesNEW[NL==0 & TD==1 & TEL==1], variablesNEW[TEL==1])  
    
    posTELNL <- match(variablesNEW[NL==1 & TD==0 & TEL==1], variablesNEW[TEL==1])  

    tt2 <- tt  # tt accounts for variables with no or only 1 effect
    vrais <- c() # Likelihood
    
    
    ### If some variables have up to 3 effects to be estimated ###

    if (nbNLTDTEL != 0){ 
 
      ## Iteration 1, step 1 (estimate NL effects) ##
      
      # Deal with variables with NL, TD and TEL effects 
        # (Iteration 1: PH assumption imposed and no TEL effect assumed)
      for (k in 1:nbNLTDTEL){
        covp <- paste("QWR[,", (dim(data)[2] + (pNLTDTEL[k]-1)*(m+p) +1):
                        (dim(data)[2] + pNLTDTEL[k]*(m+p)), "]", sep="")
        tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
      }

      # Deal with variables with NL and TD effects only (Iteration 1: PH assumption imposed)
      if (nbNLTD != 0){ 
        for (k in 1:nbNLTD){
          covp <- paste("QWR[,", (dim(data)[2] + (pNLTD[k]-1)*(m+p) +1):
                          (dim(data)[2] + (pNLTD[k])*(m+p)), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
      }

      # Deal with variables with NL and TEL effects only (Iteration 1: no TEL effect assumed)
      if (nbNLTEL != 0){ 
        for (k in 1:nbNLTEL){
          covp <- paste("QWR[,", (dim(data)[2] + (pNLTEL[k]-1)*(m+p) +1):
                          (dim(data)[2] + (pNLTEL[k])*(m + p)), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
      }

      # Deal with variables with TD and TEL effects only (PH assumption imposed and no TEL effect assumed)
      if (nbTDTEL != 0){ 
        for (k in 1:nbTDTEL){
          covp <- paste("QWR[,", variablesNEW[NL==0 & TD==1 & TEL==1][k], "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
      }      
      
      # Model estimation
      modX <- coxph(eval(parse(text = substr(tt2, 15, nchar(tt2) - 1))), method="efron")
      vrais <- c(vrais, modX$loglik[2]) # save likelihood
      
      ## Iteration 1, step 2 (estimate TD effects) ##

      modNL <- modX
      tt2 <- tt 
      VV <- matrix(ncol=((nbNLTDTEL+nbNLTD+nbTDTEL)*(m+p+1)+nbNLTEL), nrow=dim(QWR)[1]) 

      # Deal with variables with NL, TD and TEL effects 
        # (conditional on NL effect estimate from Step 1 and no TEL effect assumed)
      for (k in 1:nbNLTDTEL){
        kal <- QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))] %*% # Splines NL  
                modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]                                      # NL coef of Iter1 Step1 
        VV[,((m+p+1)*(k-1)+1):(k*(m+p+1))] <- as.vector(kal) * 
          QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+1+m+p)] # Splines for time 
        covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep="")
        tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
      }

      # Deal with variables with NL and TD effects only (conditional on NL effect estimate from Step 1)
      if (nbNLTD != 0){
        for (k in 1:nbNLTD){
          VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *                 # Splines for time
                as.vector(QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL  
                modNL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+k)*(m+p))])                    # NL coef of Iter1 Step1 
          covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        } 
      }

      # Deal with variables with TD and TEL effects only (no TEL effect assumed)
      if (nbTDTEL != 0){
        for (k in 1:nbTDTEL){
          VV[,((1+m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTD+k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] * # Splines for time  
                QWR[,variablesNEW[NL==0 & TD==1 & TEL==1][k]]                         # Variable 
          covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTD+k-1)+1):((k+nbNLTD+nbNLTDTEL)*(m+p+1)), "]", 
                         sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
      }

      # Deal with variables with NL and TEL effects only 
       # (conditional on NL effect estimate from Step 1 and no TEL effect assumed)
      if (nbNLTEL != 0){
        for (k in 1:nbNLTEL){
          VV[,((1+m+p)*(nbNLTDTEL+nbNLTD+nbTDTEL)+k)] <- 
              as.vector(QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL   
              modNL$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p))])        # NL coef of Iter1 Step1
          covp <- paste( "VV[,", ((1+m+p)*(nbNLTDTEL+nbNLTD+nbTDTEL)+k), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
        } 
      }

      modX <- coxph(eval(parse(text=substr(tt2, 15, nchar(tt2)-1))), method="efron")
      vrais <- c(vrais,modX$loglik[2]) # estimate the TD effects conditional on the NL effects

      ## Iteration 1, step 3 (estimate TEL effects) ##

      modTD <- modX
      tt2 <- tt 
      VV <- matrix(ncol=((nbNLTDTEL+nbNLTEL+nbTDTEL)*(m+p+1)+nbNLTD), nrow=dim(QWR)[1]) 

      # Deal with variables with NL, TD and TEL effects 
       # (conditional on NL and TD effects estimates from Step 1 and Step 2 respectively)
      for (k in 1:nbNLTDTEL) {
        kal <- QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))] %*% # Splines NL 
               modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]                                       # NL coef Iter1 Step1 
        kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))] %*%         # Splines for time  
                  modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]                                # TD coef Iter 1 Step 2
        VV[,((m+p+1)*(k-1)+1):(k*(m+p+1))] <- as.vector(kal) * as.vector(kalTD) *                                                            
                  QWR[,(dim(data)[2]+(m+p)*nbNL+1+m+p+(posTELNLTD[k]-1)*(m+p+1)+1):
                        (dim(data)[2]+(m+p)*nbNL+1+m+p+posTELNLTD[k]*(m+p+1))]                   # Splines TEL 
        covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep="")  
        tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
      }

      # Deal with variables with NL and TEL effects only (conditional on NL effect estimate from Step 1)
      if (nbNLTEL != 0) {
        for (k in 1:nbNLTEL) {
          VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):                         # Splines TEL
                      (dim(data)[2]+(m+p)*nbNL+(m+p+1)+posTELNL[k]*(m+p+1))] *                   
                as.vector(QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL 
                modNL$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p))])        # NL coef Iter1 Step1 
          covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse= "+") )
          }
      }

      # Deal with variables with TD and TEL effects only (conditional on TD effect estimate from Step 2)
      if (nbTDTEL != 0){
        for (k in 1:nbTDTEL){
          VV[,((1+m+p)*(nbNLTDTEL+nbNLTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTEL+k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELTD[k]-1)*(m+p+1)+1):                       # Spline TEL
                      (dim(data)[2]+(m+p)*nbNL+(m+p+1)+posTELTD[k]*(m+p+1))] *  
                as.vector(QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))] %*%     # Spline time
                modTD$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p+1))] * # TD coef Itel1 Step 2 
                QWR[, variablesNEW[NL==0 & TD==1 & TEL==1][k]])                                       # variable
          covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTEL+k-1)+1):((k+nbNLTEL+nbNLTDTEL)*(m+p+1)), "]", 
                         sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
      }

      # Deal with variables with NL and TD effects only 
        # (conditional on NL and TD effects estimate from Steps 1 and 2)
      if (nbNLTD != 0) {      
        for (k in 1:nbNLTD) {
          kal <- QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL 
                  modNL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+k)*(m+p))]          # NL coef Iter1 Step1 
          kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))] %*%   # Splines for time 
                    modTD$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))]    # TD coef Iter 1 Step 2
          VV[,((1+m+p)*(nbNLTDTEL+nbNLTEL+nbTDTEL))+k] <- as.vector(kal) * as.vector(kalTD)  
          covp <- paste( "VV[,", ((1+m+p)*(nbNLTDTEL+nbNLTEL+nbTDTEL)+k), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
      }
      
      modX <- coxph(eval(parse(text=substr(tt2,15, nchar(tt2)-1))), method="efron")
      vrais <- c(vrais,modX$loglik[2])

      diff <- 1

      ## Subsequent iterations until convergence ##
      
      while(diff > 0.0001){

        ## Step 1 (estimate NL effects) ##
        
        modTEL <- modX
        tt2 <- tt
        VV <- matrix(ncol=((nbNLTDTEL+nbNLTD+nbNLTEL)*(m+p)+nbTDTEL), nrow=dim(QWR)[1])
        
        for (k in 1:nbNLTDTEL){
          kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*%          # Splines time
                    modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]                               # TD coef previous Step 2 
          kalTEL <- QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNLTD[k]-1)*(m+p+1)+1):             # Splines TEL
                          (dim(data)[2]+nbNL*(m+p)+(m+p+1)+(posTELNLTD[k])*(m+p+1))] %*%  
                    modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]                              # TEL coef previous Step 3
          VV[,((m+p)*(k-1)+1):(k*(m+p))] <- as.vector(kalTD) * as.vector(kalTEL) * 
                    QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))] # Splines NL                       
          covp <- paste( "VV[,", ((m+p)*(k-1)+1):(k*(m+p)), "]", sep = "")
          tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
        }
        
        if(nbNLTD != 0){
          for (k in 1:nbNLTD){
            VV[,((m+p)*(nbNLTDTEL+k-1)+1):((m+p)*(nbNLTDTEL+k))] <- 
                  QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] *         # Splines NL 
                  as.vector(QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] %*% # Splines time 
                  modTD$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))])            # TD coef previous Step2 
            covp <- paste( "VV[,", ((m+p)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p)), "]", sep = "")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            
          }
        }
        
        if(nbNLTEL != 0){
          for (k in 1:nbNLTEL){
            VV[,((m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((m+p)*(nbNLTDTEL+nbNLTD+k))] <- 
                QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] * # Splines NL 
                as.vector(QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):   # Spline TEL 
                                (dim(data)[2]+nbNL*(m+p)+(m+p+1)+posTELNL[k]*(m+p+1))] %*%  
                modTEL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))])     # TEL coef previous Step 3 
            covp <- paste( "VV[,", ((m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((k+nbNLTD+nbNLTDTEL)*(m+p)), "]", 
                           sep = "")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+") )
          }
        }
 
        if(nbTDTEL != 0){
          for (k in 1:nbTDTEL){
            kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*%    # Splines time 
                      modTD$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p+1)+1):                       # TD coef previous Step 2
                                   (Nn+(nbNLTDTEL+nbNLTD+k)*(m+p+1))]  
            kalTEL <- QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELTD[k]-1)*(m+p+1)+1):         # Splines TEL
                            (dim(data)[2]+nbNL*(m+p)+(m+p+1)+(posTELTD[k])*(m+p+1))] %*%  
                      modTEL$coef[(Nn+(nbNLTDTEL+nbNLTEL+k-1)*(m+p+1)+1):                     # TEL coef previous Step 3
                                    (Nn+(nbNLTDTEL+nbNLTEL+k)*(m+p+1))] 
            VV[,(nbNLTDTEL+nbNLTD+nbNLTEL)*(m+p)+k] <- as.vector(kalTD) * as.vector(kalTEL) * 
                      QWR[, variablesNEW[NL==0 & TD==1 & TEL==1][k]]                          # Variable
            covp <- paste( "VV[,", (nbNLTDTEL+nbNLTD+nbNLTEL)*(m+p)+k, "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }

        modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
        vrais <- c(vrais, modX$loglik[2]) 

        ## Step 2 (estimate TD effects) ##

        modNL <- modX
        tt2 <- tt
        VV <- matrix(ncol=((nbNLTDTEL+nbNLTD+nbTDTEL)*(m+p+1)+nbNLTEL),nrow=dim(QWR)[1])
        
        for (k in 1:nbNLTDTEL){
          kal <- QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))] %*% # Splines NL  
                  modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]                                      # NL coef of previous Step1 
          kalTEL <- QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNLTD[k]-1)*(m+p+1)+1):              # Splines TEL
                          (dim(data)[2]+nbNL*(m+p)+(m+p+1)+posTELNLTD[k]*(m+p+1))] %*% 
                    modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]                               # TEL coef of previous Step3 
          VV[,((m+p+1)*(k-1)+1):(k*(m+p+1))] <- as.vector(kal)*as.vector(kalTEL) * 
                    QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+1+m+p)]              # Splines time  
          covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }

        if(nbNLTD != 0){
          for (k in 1:nbNLTD){
            VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))] <- 
                  QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *                 # Splines for time 
                  as.vector(QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL 
                  modNL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+k)*(m+p))])                    # NL coef of previous Step1 
            covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
          }
        }
        
        if(nbTDTEL != 0){
          for (k in 1:nbTDTEL){
            VV[,((1+m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTD+k))] <-
                    QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *     # Splines time 
                    as.vector(QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELTD[k]-1)*(m+p+1)+1): # Splines TEL 
                                (dim(data)[2]+nbNL*(m+p)+(posTELTD[k]+1)*(m+p+1))] %*%        
                    modTEL$coef[(Nn+(nbNLTDTEL+nbNLTEL+k-1)*(m+p+1)+1):                       # TEL coef of previous Step3 
                                (Nn+(nbNLTDTEL+nbNLTEL+k)*(m+p+1))] *       
                  QWR[,variablesNEW[NL==0 & TD==1 & TEL==1][k]])                              # Variable
            covp <- paste("VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTD+k-1)+1):((k+nbNLTD+nbNLTDTEL)*(m+p+1)), "]", 
                          sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
          }
        }

        if(nbNLTEL != 0){
          for (k in 1:nbNLTEL){        
            kal <- QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*%   # Splines NL  
                   modNL$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p))] # NL coef of previous Step1 
            kalTEL <- QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):              # Splines TEL
                            (dim(data)[2]+nbNL*(m+p)+(m+p+1)+posTELNL[k]*(m+p+1))] %*%          
                      modTEL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))]       # TEL coef of previous Step3 
            VV[, (nbNLTDTEL+nbNLTD+nbTDTEL)*(m+p+1)+k] <- as.vector(kal) * as.vector(kalTEL)  
            covp <- paste( "VV[,", (nbNLTDTEL+nbNLTD+nbTDTEL)*(m+p+1)+k, "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
          }
        }
                          
        modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
        vrais <- c(vrais, modX$loglik[2])

        ## Step 3 (estimate TEL effects) ##

        modTD <- modX
        tt2 <- tt   
        VV <- matrix(ncol=((nbNLTDTEL+nbNLTEL+nbTDTEL)*(m+p+1)+nbNLTD), nrow=dim(QWR)[1]) 
        
        for (k in 1:nbNLTDTEL){
          kal <- QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))] %*% # Splines NL
                  modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]                                      # NL coef previous Step1  
          kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*%           # Splines for time
                    modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]                                # TD coef previous Step 2
          VV[,((m+p+1)*(k-1)+1):(k*(m+p+1))]<-as.vector(kal) * as.vector(kalTD) * 
                QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNLTD[k]-1)*(m+p+1)+1):
                      (dim(data)[2]+(m+p)*nbNL+(posTELNLTD[k]+1)*(m+p+1))]                         # Splines TEL
          covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep="")
          tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
        }
        
        if (nbNLTEL != 0){
          for (k in 1:nbNLTEL){
            VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))] <-
              QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):                         # Splines TEL
                    (dim(data)[2]+(m+p)*nbNL+(posTELNL[k]+1)*(m+p+1))] *   
              as.vector(QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL 
              modNL$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p))])        # NL coef previous Step1 
            covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            }
        }
        
        if (nbTDTEL != 0){
          for (k in 1:nbTDTEL){
            kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*%                  # Splines TEL 
                      modTD$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p+1))] * # Splines time 
                      QWR[,variablesNEW[NL==0 & TD==1 & TEL==1][k]]                                         # TD previous Step2 
            VV[,((1+m+p)*(nbNLTDTEL+nbNLTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTEL+k))] <- as.vector(kalTD) *                          
                      QWR[,(dim(data)[2]+(m+p)*nbNL+(m+p+1)+(posTELTD[k]-1)*(m+p+1)+1):                     # Variable
                            (dim(data)[2]+(m+p)*nbNL+(posTELTD[k]+1)*(m+p+1))]
            covp <- paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTEL+k-1)+1):((k+nbNLTEL+nbNLTDTEL)*(m+p+1)), "]", 
                           sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            }
        }

        if (nbNLTD != 0){      
          for (k in 1:nbNLTD){
            kal <- QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL 
              modNL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+k)*(m+p))]                # NL coef previous Step1 
            kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))] %*%   # Splines for time 
              modTD$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))]            # TD coef previous Step 2
            VV[,((1+m+p)*(nbNLTDTEL+nbNLTEL+nbTDTEL))+k] <- as.vector(kal) * as.vector(kalTD)  
            covp <- paste( "VV[,", ((1+m+p)*(nbNLTDTEL+nbNLTEL+nbTDTEL)+k), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }
        
        modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
        vrais <- c(vrais, modX$loglik[2])

        diff <- vrais[length(vrais)] - vrais[length(vrais) - 3]
      }

      
    ### ELSE IF (at least one variable has 2 effects, but none have 3 effects) ###

    } else if (nbNLTD+nbNLTEL+nbTDTEL != 0){ 

      ## Iteration 1, step 1 (estimate NL effects) ##
      
      if (nbNLTD+nbNLTEL != 0){ #if there is a NL effect 
        
        if (nbNLTD != 0){
          for (k in 1:nbNLTD){
            covp <- paste("QWR[,", (dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+(pNLTD[k])*(m+p)), "]",# Splines for variable
                          sep="")  
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }
        
        if (nbNLTEL != 0){
          for (k in 1:nbNLTEL){
            covp <- paste("QWR[,", (dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+(pNLTEL[k])*(m+p)), "]", # Splines for variable
                          sep="") 
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }

        # Deal with variables with TD and TEL effects only (PH assumption imposed and no TEL effect assumed)
        if (nbTDTEL != 0){ 
          for (k in 1:nbTDTEL){
            covp <- paste("QWR[,", variablesNEW[NL==0 & TD==1 & TEL==1][k], "]", sep="") # Variable
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }
        
        modX <- coxph(eval(parse(text = substr(tt2, 15, nchar(tt2)-1))), method="efron")
        vrais <- c(vrais, modX$loglik[2])

        modNL <- modX
        tt2 <- tt
      } 

      ## Iteration 1, step 2 (estimate TD effects) ## 
      
      if (nbNLTD+nbTDTEL != 0){ 
        
        VV <- matrix(ncol=((nbNLTD+nbTDTEL)*(m+p+1)+nbNLTEL), nrow=dim(QWR)[1]) 
        
        if (nbNLTD != 0){
          for (k in 1:nbNLTD){
            VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *                 # Splines for time 
                as.vector(QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL 
                modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])                                        # NL coef of Iter1 Step1 
            covp <- paste("VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            }
        }
        
        if (nbTDTEL != 0){ 
          for (k in 1:nbTDTEL){
            VV[,((1+m+p)*(nbNLTD+k-1)+1):((1+m+p)*(nbNLTD+k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] * # Splines for time 
                QWR[,variablesNEW[NL==0 & TD==1 & TEL==1][k]]                         # Variable
            
            covp <- paste("VV[,", ((m+p+1)*(nbNLTD+k-1)+1):((k+nbNLTD)*(m+p+1)), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
           }
        }
        
        # Deal with variables with NL and TEL effects only (conditional on NL effect estimate from Step 1 and no TEL effect assumed)
        if (nbNLTEL != 0){
          for (k in 1:nbNLTEL){
            VV[,((1+m+p)*(nbNLTD+nbTDTEL)+k)] <- 
              as.vector(QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL   
              modNL$coef[(Nn+(nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTD+k)*(m+p))])                            # NL coef of Iter1 Step1
            covp <- paste("VV[,", ((1+m+p)*(nbNLTD+nbTDTEL)+k), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
          } 
        }      
      
        modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
        vrais <- c(vrais, modX$loglik[2]) 

        modTD <- modX
        tt2 <- tt
      }

      ## Iteration 1, step 3 (estimate TEL effects) ##

      if (nbNLTEL+nbTDTEL != 0){ 
        
        VV <- matrix(ncol=((nbNLTEL+nbTDTEL)*(m+p+1)+nbNLTD), nrow=dim(QWR)[1]) 

        if (nbNLTEL != 0){
          for (k in 1:nbNLTEL){
            VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))] <- 
              QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):                         # Splines TEL
                    (dim(data)[2]+(m+p)*nbNL+(posTELNL[k]+1)*(m+p+1))] * 
              as.vector(QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL 
              modNL$coef[(Nn+(nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTD+k)*(m+p))])                            # NL coef Iter1 Step1 
            covp <- paste("VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse="+") )
          }
        }
                
        if (nbTDTEL != 0){
          for (k in 1:nbTDTEL){
            kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*% # Splines time
                      modTD$coef[(Nn+(nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTD+k)*(m+p+1))] *    # TD coef Itel1 Step2 
                      QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]                    # Variable
            VV[,((1+m+p)*(nbNLTEL+k-1)+1):((1+m+p)*(nbNLTEL+k))] <- as.vector(kalTD) *
                      QWR[,(dim(data)[2]+(m+p)*nbNL+(m+p+1)+(posTELTD[k]-1)*(m+p+1)+1):    # Spline TEL
                            (dim(data)[2]+(m+p)*nbNL+(posTELTD[k]+1)*(m+p+1))]  
            covp <- paste("VV[,", ((m+p+1)*(nbNLTEL+k-1)+1):((k+nbNLTEL)*(m+p+1)), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
          }
        }

        if (nbNLTD != 0){      
          for (k in 1:nbNLTD){
            kal <- QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL 
              modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))]                                    # NL coef Iter1 Step1 
            kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))] %*%   # Splines for time 
              modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))]                                # TD coef Iter1 Step2
            VV[,((1+m+p)*(nbNLTEL+nbTDTEL))+k] <- as.vector(kal) * as.vector(kalTD)  
            covp <- paste("VV[,", ((1+m+p)*(nbNLTEL+nbTDTEL)+k), "]", sep="")
            tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
          }
        }        
                
        modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
        vrais <- c(vrais, modX$loglik[2])

        modTEL <- modX
        tt2 <- tt
      }
      
      diff <- 1

      ## Subsequent iterations until convergence ##
      
      while(diff > 0.0001){ 

        ## Step 1 (estimate NL effects) ##
        
        if (nbNLTD+nbNLTEL != 0){
          
          VV <- matrix(ncol=((nbNLTD+nbNLTEL)*(m+p)+nbTDTEL), nrow=dim(QWR)[1])
          
          if(nbNLTD != 0){
            for (k in 1:nbNLTD){
              VV[,((m+p)*(k-1)+1):((m+p)*(k))] <- 
                    QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] *         # Splines NL 
                    as.vector(QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] %*% # Splines time 
                    modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])                                # TD coef Iter1 Step2               
              covp <- paste("VV[,", ((m+p)*(k-1)+1):((k)*(m+p)), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
            }
          }
          
          if(nbNLTEL != 0){
            for (k in 1:nbNLTEL){
              VV[,((m+p)*(nbNLTD+k-1)+1):((m+p)*(nbNLTD+k))] <- 
                    QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] * # Splines NL 
                    as.vector(QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):   # Spline TEL
                                    (dim(data)[2]+nbNL*(m+p)+(posTELNL[k]+1)*(m+p+1))] %*%  
                    modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])                         # TEL coef Itel1 Step3 
              covp <- paste("VV[,", ((m+p)*(nbNLTD+k-1)+1):((k+nbNLTD)*(m+p)), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            }
          }

          if(nbTDTEL != 0){
            for (k in 1:nbTDTEL){
              kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*% # Splines time 
                  modTD$coef[(Nn+(nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTD+k)*(m+p+1))]            # TD coef previous Step 2 
              kalTEL <- QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELTD[k]-1)*(m+p+1)+1):      # Splines TEL
                              (dim(data)[2]+nbNL*(m+p)+(m+p+1)+(posTELTD[k])*(m+p+1))] %*%  
                  modTEL$coef[(Nn+(nbNLTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTEL+k)*(m+p+1))]         # TEL coef previous Step 3
              VV[,(nbNLTD+nbNLTEL)*(m+p)+k] <- as.vector(kalTD) * as.vector(kalTEL) * 
                  QWR[, variablesNEW[NL==0 & TD==1 & TEL==1][k]]                             # Variable
              covp <- paste("VV[,", (nbNLTD+nbNLTEL)*(m+p)+k, "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            }
          }
            
          modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
          vrais <- c(vrais, modX$loglik[2])

          modNL <- modX
          tt2 <- tt
        }

          ## Step 2 (estimate TD effects) ##  
        
          if (nbNLTD+nbTDTEL != 0){
            
          VV <- matrix(ncol=((nbNLTD+nbTDTEL)*(m+p+1)+nbNLTEL), nrow=dim(QWR)[1])
          
          if(nbNLTD != 0){
            for (k in 1:nbNLTD){
              VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))] <- 
                  QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *
                  as.vector(QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*%
                  modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
              covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
            }
          }
          
          if(nbTDTEL != 0){
            for (k in 1:nbTDTEL){
              VV[,((1+m+p)*(nbNLTD+k-1)+1):((1+m+p)*(nbNLTD+k))] <- 
                  QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *     # Splines time
                  as.vector(QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELTD[k]-1)*(m+p+1)+1): # Splines TEL
                                  (dim(data)[2]+nbNL*(m+p)+(posTELTD[k]+1)*(m+p+1))] %*%    
                  modTEL$coef[(Nn+(nbNLTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTEL+k)*(m+p+1))] *      # TEL coef of Iter1 Step3
                  QWR[,variablesNEW[NL==0 & TD==1 & TEL==1][k]])                            # Variable 
              covp <- paste( "VV[,", ((m+p+1)*(nbNLTD+k-1)+1):((k+nbNLTD)*(m+p+1)), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
            }
          }
          
          if(nbNLTEL != 0){
            for (k in 1:nbNLTEL){        
              kal <- QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL  
                     modNL$coef[(Nn+(nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTD+k)*(m+p))]                   # NL coef of previous Step1 
              kalTEL <- QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):            # Splines TEL
                              (dim(data)[2]+nbNL*(m+p)+(m+p+1)+posTELNL[k]*(m+p+1))] %*%  
                modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))]                                 # TEL coef of previous Step3 
              VV[, (nbNLTD+nbTDTEL)*(m+p+1)+k] <- as.vector(kal) * as.vector(kalTEL)  
              covp <- paste("VV[,", (nbNLTD+nbTDTEL)*(m+p+1)+k, "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
            }
          }          
          
          modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
          vrais <- c(vrais, modX$loglik[2])
          
          modTD <- modX
          tt2 <- tt 
          
        }
        
        ## Step 3 (estimate TEL effects) ##
        
        if (nbNLTEL+nbTDTEL != 0){ 
          
          VV <- matrix(ncol=((nbNLTEL+nbTDTEL)*(m+p+1)+nbNLTD), nrow=dim(QWR)[1]) 
          
          if (nbNLTEL != 0){
            for (k in 1:nbNLTEL){
              VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))] <- 
                QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):                         # Splines TEL 
                     (dim(data)[2]+(m+p)*nbNL+(posTELNL[k]+1)*(m+p+1))] *                               # Splines TEL 
                as.vector(QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))] %*% # Splines NL
                modNL$coef[(Nn+(nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTD+k)*(m+p))])                            # NL coef previous Step1
              covp <- paste("VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse="+") )
            }
          }
          
          if (nbTDTEL != 0){
            for (k in 1:nbTDTEL){
              kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)] %*% # Splines time
                        modTD$coef[(Nn+(nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTD+k)*(m+p+1))] *    # TD coef previous Step2 
                        QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]                    # Variable 
              VV[,((1+m+p)*(nbNLTEL+k-1)+1):((1+m+p)*(nbNLTEL+k))] <- as.vector(kalTD ) *
                        QWR[,(dim(data)[2]+(m+p)*nbNL+(m+p+1)+(posTELTD[k]-1)*(m+p+1)+1):    # Splines TEL
                              (dim(data)[2]+(m+p)*nbNL+(posTELTD[k]+1)*(m+p+1))]          
              
              covp <- paste("VV[,", ((m+p+1)*(nbNLTEL+k-1)+1):((k+nbNLTEL)*(m+p+1)), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse= "+") )
            }
          }

          if (nbNLTD != 0){      
            for (k in 1:nbNLTD){
              kal <- QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))] %*% # Splines NL 
                modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))]                                    # NL coef previous Step1 
              kalTD <- QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))] %*%   # Splines for time 
                modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))]                                # TD coef previous Step 2
              VV[,((1+m+p)*(nbNLTEL+nbTDTEL))+k] <- as.vector(kal) * as.vector(kalTD)  
              covp <- paste( "VV[,", ((1+m+p)*(nbNLTEL+nbTDTEL)+k), "]", sep="")
              tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
            }
          }
          
          modX <- coxph(eval(parse(text=substr(tt2,15,nchar(tt2)-1))), method="efron")
          vrais <- c(vrais, modX$loglik[2]) 
          
          modTEL <- modX 
          tt2 <- tt 
        }
        
        if ((i4!=0 & i6!=0) | (i4!=0 & i7!=0) | (i6!=0 & i7!=0)){
          diff <- vrais[length(vrais)] - vrais[length(vrais) - 3]
        } else { 
          diff <- vrais[length(vrais)] - vrais[length(vrais) - 2]
        }
      }
      
        
    } else {
      
    ### ELSE IF (All variables have 0 or 1 effect)
      
      modX <- coxph(eval(parse(text=substr(tt2, 15, nchar(tt2) - 1))), method="efron")
      vrais <- c(modX$loglik[2])
    } 

    rm(QWR)
    gc()

    ### Store parameters estimated ###
    
    MAT <- matrix(ncol=V, nrow= 1 + m+p + m+p+1 + m+p+1) # 1 variable per column
                                                         # 1 line for coefficients for variables without effects
                                                         # m+p lines for spline coefficients for NL effects
                                                         # m+p+1 lines for spline coefficients for TD effects
                                                         # m+p+1 lines for spline coefficients for TEL effects
    if (i1 != 0) { # for variables without effects
      for (j in 1:i1) { 
        MAT[1, j] <- modX$coef[j]
      }
    }
    if (i2 != 0) { # for variables with only NL effect, m+p coefficients per variable 
      for (j in 1:i2) { 
        MAT[2:(m+p+1), i1+j] <- modX$coef[(i1 + (j-1)*(m+p) + 1):(i1 + (j-1)*(m+p) + m+p)]
      }
    }
    if (i3 != 0) { # for variables with only TD effect, m+p+1 coefficients per variable 
      for (j in 1:i3) { 
        MAT[(m+p+2):(2*m+2*p+2), i1+i2+j] <- modX$coef[(i1 + i2*(m+p) + (j-1)*(m+p+1) + 1):
                                                       (i1 + i2*(m+p) + (j-1)*(m+p+1) + m+p+1)]
      }
    }
    if (i5 != 0) { # for variables with only TEL effect, m+p+1 coefficients per variable 
      for (j in 1:i5) { 
        MAT[(2*m + 2*p + 3):(3*m + 3*p + 3), i1+i2+i3+j] <- 
            modX$coef[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + 1):
                      (i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + m+p+1)]
      }
    }
    if (i8!=0) { # for variables with NL+TD+TEL effects
      for (j in 1:i8) {
        MAT[2:(m+p+1), i1+i2+i3+i5+j] <- modNL$coef[(i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p) + 1):
                                                    (i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p) + m+p)]
        MAT[(m+p+2):(2*(m+p+1)), i1+i2+i3+i5+j] <- 
            modTD$coef[(i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1) * (m+p+1) + 1):
                       (i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1) * (m+p+1) + m+p+1)]
        MAT[(2*m+2*p+3):(3*(m+p+1)), i1+i2+i3+i5+j] <- 
            modX$coef[(i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p+1) + 1):
                      (i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p+1) + m+p+1)]
      }
    }
    if (i4!=0) { # for variables with NL+TD effects
      for (j in 1:i4){
        MAT[2:(m+p+1), i1+i2+i3+i5+i8+j] <- modNL$coef[(i1 + (i2+i8)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + 1):
                                                       (i1 + (i2+i8)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + m+p)]
        MAT[(m+p+2):(2*(m+p+1)), i1+i2+i3+i5+i8+j] <- 
            modTD$coef[(i1 + i2*(m+p) + (i3+i5+i8) * (m+p+1) + (j-1)*(m+p+1) + 1):
                       (i1 + i2*(m+p) + (i3+i5+i8) * (m+p+1) + (j-1)*(m+p+1) + m+p+1)] *
            modX$coef[i1 + i2*(m+p) + (i3+i5+i8+i6+i7)*(m+p+1) + j] # Coef from step3 for fixed values of NL.TD effect
      }                     
    }
    if (i6!=0) { # for variables with NL+TEL effects
      for(j in 1:i6){
        MAT[2:(m+p+1), i1+i2+i3+i5+i8+i4+j] <- 
          modNL$coef[(i1 + (i2+i8+i4)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + 1):
                     (i1 + (i2+i8+i4)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + m+p)]
        MAT[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+i4+ j] <- 
          modX$coef[(i1 + i2 * (m + p) + (i3+i5+i8) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):
                    (i1 + i2 * (m + p) + (i3+i5+i8) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)]
      }
    }
    if (i7!=0) { # for variables with TD+TEL effects
      for(j in 1:i7){
        MAT[(m+p+2):(2*(m+p+1)), i1+i2+i3+i5+i8+i6+i4+j] <- 
              modTD$coef[(i1 + (i2+i8)*(m+p) + (i3+i5+i4)*(m+p+1) + (j-1)*(m+p+1) + 1):
                         (i1 + (i2+i8)*(m+p) + (i3+i5+i4)*(m+p+1) + (j-1)*(m+p+1) + m+p+1)]
        MAT[(2*m+2*p+3):(3*(m+p+1)), i1+i2+i3+i5+i8+i6+i4+j] <- 
              modX$coef[(i1 + i2*(m+p) + (i3+i5+i8+i6)*(m+p+1) + (j-1)*(m+p+1) + 1):
                        (i1 + i2*(m+p) + (i3+i5+i8+i6)*(m+p+1) + (j-1)*(m+p+1) + m+p+1)]
      }
    }
   
    ### Store SE estimated ###
    
    MATse <- matrix(ncol=V, nrow=3*(1+m+p)) 
    if (i1 != 0){ # for variables without effects
      for (j in 1:i1){
        MATse[1, j] <- sqrt(diag(modX$var)[j])
      }
    }
    if (i2 != 0){
      for (j in 1:i2){ # for variables with only NL effect, m+p coefficients per variable
        MATse[2:(m+p+1), i1+j] <- sqrt(diag(modX$var)[(i1 + (j-1)*(m+p) + 1):(i1 + (j-1)*(m+p) + m+p)])
      }
    }
    if (i3 != 0){ # for variables with only TD effect, m+p+1 coefficients per variable
      for (j in 1:i3){
        MATse[(m+p+2):(2*m+2*p+2), i1+i2+j] <- sqrt(diag(modX$var)[(i1 + i2*(m+p) + (j-1)*(m+p+1) + 1):
                                                                   (i1 + i2*(m+p) + (j-1)*(m+p+1) + m+p+1)])
      }
    }
    if (i5 != 0){ # for variables with only TEL effect, m+p+1 coefficients per variable
      for (j in 1:i5){ # row m+p+2 to 2*(m+p+1) stores TD effects
        MATse[(2*m+2*p+3):(3*m+3*p+3), i1+i2+i3+j] <- 
          sqrt(diag(modX$var)[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + 1):
                              (i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
      }
    }
    if(i8 != 0){ # for variables with NL+TD+TEL effects
      for (j in 1:i8){
        MATse[2:(m+p+1), i1+i2+i3+i5+j] <- 
          sqrt(diag(modNL$var)[(i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p) + 1):
                               (i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p) + m+p)])
        MATse[(m+p+2):(2*(m+p+1)), i1+i2+i3+i5+j] <- 
          sqrt(diag(modTD$var)[(i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p+1) + 1):
                               (i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
        MATse[(2*m+2*p+3):(3*(m+p+1)), i1+i2+i3+i5+j] <- 
          sqrt(diag(modX$var)[(i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p+1) + 1):
                              (i1 + i2*(m+p) + i3*(m+p+1) + i5*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
      }
    }
    if(i4 != 0){ # for variables with NL+TD effects
      for (j in 1:i4){
        MATse[2:(m+p+1), i1+i2+i3+i5+i8+j] <- 
          sqrt(diag(modNL$var)[(i1 + (i2+i8)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + 1):
                               (i1 + (i2+i8)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + m+p)])
        MATse[(m+p+2):(2*(m+p+1)), i1+i2+i3+i5+i8+j] <- 
          sqrt(diag(modTD$var)[(i1 + i2*(m+p) + (i3+i5+i8)*(m+p+1) + (j-1)*(m+p+1) + 1):
                               (i1 + i2*(m+p) + (i3+i5+i8)*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
      }                     
    }
    if(i6 != 0){ # for variables with NL+TEL effects
      for(j in 1:i6){
        MATse[2:(m+p+1), i1+i2+i3+i5+i8+i4+j] <-
          sqrt(diag(modNL$var)[(i1 + (i2+i8+i4)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + 1):
                               (i1 + (i2+i8+i4)*(m+p) + (i3+i5)*(m+p+1) + (j-1)*(m+p) + m+p)])
        MATse[(2*m+2*p+3):(3*(m+p+1)), i1+i2+i3+i5+i8+i4+j] <-
          sqrt(diag(modX$var)[(i1 + i2*(m+p) + (i3+i5+i8)*(m+p+1) + (j-1)*(m+p+1) + 1):
                              (i1 + i2*(m+p) + (i3+i5+i8)*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
      }
    }
    if(i7 != 0){ # for variables with TD+TEL effects
      for(j in 1:i7){
        MATse[(m+p+2):(2*(m+p+1)), i1+i2+i3+i5+i8+i4+i6+j] <-
          sqrt(diag(modTD$var)[(i1 + (i2+i8)*(m+p) + (i3+i5+i4)*(m+p+1) + (j-1)*(m+p+1) + 1):
                               (i1 + (i2+i8)*(m+p) + (i3+i5+i4)*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
        MATse[(2*m+2*p+3):(3*(m+p+1)), i1+i2+i3+i5+i8+i6+i4+j] <-
          sqrt(diag(modX$var)[(i1 + i2*(m+p) + (i3+i5+i8+i6)*(m+p+1) + (j-1)*(m+p+1) + 1):
                              (i1 + i2*(m+p) + (i3+i5+i8+i6)*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
      }
    }

    ### Prepare output ###
    
    var_order <- c(variablesNEW[NL==0 & TD==0 &TEL==0], variablesNEW[NL==1 & TD==0 & TEL==0], 
                   variablesNEW[NL==0 & TD==1 &TEL==0], variablesNEW[NL==0 & TD==0 & TEL==1],
                   variablesNEW[NL==1 & TD==1 &TEL==1], variablesNEW[NL==1 & TD==1 & TEL==0],
                   variablesNEW[NL==1 & TD==0 &TEL==1], variablesNEW[NL==0 & TD==1 & TEL==1])

    coefficients <- MAT[1, match(variablesNEW, var_order)] # Coefficient for variables without effects 
    names(coefficients) <- variables
    se_coef <- MATse[1, match(variablesNEW, var_order)]    # SE for variables without effects
    
    coefficients_splines_NL <- as.matrix(MAT[2:(m+p+1), match(variablesNEW,var_order)])
    coefficients_splines_NL <- rbind(rep(0,V), coefficients_splines_NL)
    coefficients_splines_NL[1, (NL==0)] <- NA
    colnames(coefficients_splines_NL) <- variables
   
    coefficients_splines_TD <- as.matrix(MAT[(m+p+2):(2*(m+p+1)), match(variablesNEW,var_order)])
    colnames(coefficients_splines_TD) <- variables
    
    coefficients_splines_TEL <- as.matrix(MAT[(2*(m+p)+3):(3*(m+p+1)), match(variablesNEW,var_order)])
    colnames(coefficients_splines_TEL) <- variables

    knots_covariates <- knotsNEW[1:V,]
    if (V>1){rownames(knots_covariates) <- variables}
    if (V>1){
      knots_covariates[(NL==0),] <- rep(NA, p+1+m+p+1)
    } else { 
      if (NL==0) {knots_covariates <- rep(NA, p+1+m+p+1)} 
    }
    
    knots_time <- knotsNEW[V+1,]
    
    if(nTEL>1){
      knots_TEL <- knotsNEW[(V+2):(V+1+nTEL),]
      rownames(knots_TEL) <- variables[TEL==1]
    } else {
      knots_TEL <- knotsNEW[(V+2),]
      names(knots_TEL) <- variables[TEL==1]
    }
     
    nEvents <- sum(data[,TypeNEW[3]]==1)
    rm(data)
    gc()

    list(Partial_Log_Likelihood = vrais[length(vrais)],
         log_likelihood_history = vrais,
         Number_of_parameters = nbpara + nbonlyNL*(m+p) + (nbonlyTD+nbonlyTEL+2*nbTDTEL)*(1+m+p) + 
                                (nbNLTD+nbNLTEL)*(m+p+m+p+1) + nbNLTDTEL*(3*(m+p+1)-1), 
         Number_events = nEvents, 
         Number_knots = m, 
         Degree_of_splines = p, 
         knots_covariates = knots_covariates, 
         knots_time = knots_time, 
         knots_TEL = knots_TEL,
         coefficients = coefficients, 
         Standard_Error=se_coef,
         coefficients_splines_NL = coefficients_splines_NL,
         coefficients_splines_TD = coefficients_splines_TD, 
         coefficients_splines_TEL = coefficients_splines_TEL,
         variables = variables)
  }
}



last_prog <- function(data, Type, variables, TD, NL, m, p, knots=-999){
  
  if (length(Type) == 2){
    Type2 <- Type
    data$StartV0 <- rep(0, dim(data)[1])
    Type <- c("StartV0", Type2[1], Type2[2])
  } 
  
  i1 <- sum((NL+TD)==0)          # number of non-nl non-td variables
  i2 <- sum(((NL==1) & (TD==0))) # number of NL non-td variables
  i3 <- sum(((NL==0) & (TD==1))) # number of non-nl TD variables 
  i4 <- sum((NL+TD)==2)          # number of NL and TD variables
  nonpara <- TD+NL 
    
  V <- length(variables)         # number of covariates
  
  variablesNEW <- match(variables, names(data)) # positions in data of items in argument variables 
  TypeNEW <- match(Type, names(data))           # positions in data of items in argument Type
  
  ## Store interior and exterior knots in (2*(p+1)+m) columns, for variables (rows 1:V) and time (row V+1)

  listeprobaquantile <- seq(1, m)/(m+1) 
  knotsNEW <- matrix(nrow=V+1, ncol= p+1 + m + p+1) 

  # If no user-defined knots
  if (is.matrix(knots)==FALSE){ 
    if (is.numeric(knots)==TRUE & knots==-999){ # use default knots 
      
      # Knots for variables:
        # p+1 exterior knots at min(variables[i])
        # interior knots at quantiles
        # p+1 exterior knots equally spaced between max(variable) and max(variable)+p
      for (i in 1:V){
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p+1), 
                           quantile(data[, variables[i]], probs=listeprobaquantile),  
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1)) 
      }
      
      # Knots for time:
        # p+1 exterior knots at 0
        # interior knots at quantiles of event times
        # p+1 exterior knots equally spaced between max(time) and max(time)+p
      knotsNEW[V+1, ] <- c(rep(0, p+1), 
                             quantile(data[data[, TypeNEW[3]]==1, TypeNEW[2]], probs=listeprobaquantile), 
                             seq(max(data[data[,TypeNEW[3]]==1, TypeNEW[2]]), 
                                 max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1)) 
    }                                                                          
  }

  # If (some) user-defined interior knots are specified 
  if (is.matrix(knots) == TRUE){ 
    if (dim(knots)[1]!=(length(variables)+1) | dim(knots)[2]!=m) 
      stop("Error Message: variable knots should be a matrix of dimension (length(variables)+1)*m")

    # Knots for variables    
    for (i in 1:V){
      if (is.na(knots[i,1]) == TRUE){ # if user-defined knots for that variable is NA, use default
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p+1), 
                           quantile(data[, variables[i]], probs=listeprobaquantile), 
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
      } else {
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),   
                           knots[i,], 
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
        ## else use user defined interior knots 
      }
    }

    # Knots for time    
    if (is.na(knots[(V+1),1]) == TRUE){ # if user-defined knots for time is NA, use default
      knotsNEW[(V+1), ] <- c(rep(0, p+1), 
                             quantile(data[data[, TypeNEW[3]]==1, TypeNEW[2]], probs=listeprobaquantile), 
                             seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), 
                                 max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1))
     } else {
      knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                             knots[(V+1),], 
                             seq(max(data[data[,TypeNEW[3]]==1, TypeNEW[2]]), 
                                 max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1))
    }
  }

  ## Prepare data (QWR) for estimation
  
  data <- as.matrix(data)
  listeT <- c(0, sort(unique(data[data[, TypeNEW[3]]==1, TypeNEW[2]]))) # list of unique event times
  ncol <- dim(data)[2]  
 
  X <- split(data, data[, 1]) 
  matX <- sapply(X, DvlpMatrix, listeT=listeT, ncol=ncol, TypeNEW=TypeNEW) 
  QWR <- do.call(rbind, matX) # delete rows not in any risk sets to improve efficiency 

  # Add splines for NL effects in QWR   
  nbNL <- sum(NL)
  if (nbNL != 0){  
    for (i in 1:nbNL){ 
      QWR <- cbind(QWR, 
              splineDesign(knotsNEW[seq(1,V, 1)[NL==1][i],], x=QWR[, variablesNEW[NL==1][i]], ord=p+1)[,-1])
    }
  }

  # Add splines for TD effects in QWR (splines for time)  
  QWR <- cbind(QWR, splineDesign(knotsNEW[V+1,], x=QWR[, TypeNEW[2]], ord=p+1))

  ## Construct equation (tt) of Cox model (part to be estimated at all steps of all iterations)  
  tt <- paste("modX<-coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,", TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", 
              sep="")
 
  Nn <- 0 # number of parameters added to model up to now

  # Add parameters for variables without NL, TD and TEL effects 
  nbpara <- sum((NL==0 & TD==0)) 
  if (nbpara != 0) { 
    for (k in 1:nbpara) {
      tt <- paste(tt, "QWR[,", variablesNEW[NL==0 & TD==0][k], "]+", sep="")
      Nn <- Nn + 1
    }
  }

  # Add parameters for variables with NL effect only   
  nbonlyNL <- sum((NL==1 & TD==0)) 
  if (nbonlyNL != 0){
    onlyNL <- match(variablesNEW[NL==1 & TD==0], variablesNEW[NL==1]) 
    for (k in 1:nbonlyNL){
      covp <- paste("QWR[,", (dim(data)[2] + (onlyNL[k]-1)*(m+p)+1):(dim(data)[2] + onlyNL[k]*(m+p)), "]",
                     sep="")
      tt <- paste(tt, paste(c(covp,""), collapse="+"))
      Nn <- Nn+m+p
    }
  }

  # Add parameters for variables with TD effect only   
  nbonlyTD <- sum((NL==0 & TD==1)) # no NL only TD
  if (nbonlyTD != 0){
    for (k in 1:nbonlyTD){
      flag<-dim(QWR)[2]+1
      QWR <- cbind(QWR, QWR[, variablesNEW[NL==0 & TD==1][k]] * 
                     QWR[, (dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*(nbNL+1)+1)]) 
      covp <- paste("QWR[,", flag: dim(QWR)[2], "]", sep="")
      tt <- paste(tt, paste(c(covp,""), collapse="+"))  
      Nn <- Nn + m+p+1
    }#
  }
  
  tt2 <- tt

  nbNLTD <- sum((NL==1 & TD==1)) 
  NOTonlyNL <- match(variablesNEW[NL==1 & TD==1], variablesNEW[NL==1])

  ### If some variables have up to 2 effects to be estimated ###

  if (nbNLTD != 0){ 

    ## Iteration 1, step 1 (estimate NL effects) ##
    
    for (k in 1:nbNLTD){
      covp <- paste("QWR[,",(dim(data)[2] + (NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2] + NOTonlyNL[k]*(m+p)), "]",
                    sep="")
      tt2 <- paste(tt2, paste(c(covp,""), collapse= "+")) 
    }
    
    vrais <- c()
    modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 1))), method="efron")
    vrais <- c(vrais, modX$loglik[2])

    ## Iteration 1, step 2 (estimate TD effects) ##
    
    mod <- modX 
    tt2 <- tt
    VV <- matrix(ncol=nbNLTD*(m+p+1), nrow=dim(QWR)[1])
    
    for (k in 1:nbNLTD){
      VV[,((1+m+p)*(k-1)+1):((1+m+p)*k)] <- QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] * 
                as.vector(QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))] %*%
                mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))])
      covp <- paste("VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
      tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
    }

    modX <- coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
    vrais <- c(vrais,modX$loglik[2])
     
    diff <- 1

    ## Subsequent iterations until convergence ##
    
    while (diff > 0.00001){ 
 
      ## Step 1 (estimate NL effects) ##
      
      mod <- modX
      tt2 <- tt
      VV <- matrix(ncol=nbNLTD*(m+p),nrow=dim(QWR)[1])
      
      for (k in 1:nbNLTD){
        VV[,((m+p)*(k-1)+1):((m+p)*(k))] <- 
              QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))] *
              as.vector(QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] %*%
              mod$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])
        covp <- paste("VV[,", ((m+p)*(k-1)+1):((k)*(m+p)), "]", sep="")
        tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
      }
      
      modX <- coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))), method="efron")
      vrais <- c(vrais, modX$loglik[2])

      ## Step 2 (estimate TD effects) ##
      
      mod <- modX
      tt2 <- tt
      VV <- matrix(ncol=nbNLTD*(m+p+1), nrow=dim(QWR)[1])

      for (k in 1:nbNLTD){
        VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))] <- 
              QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] *
              as.vector(QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))] %*%
              mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
        covp <- paste("VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
        tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
      }
     
      modX <- coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))), method="efron")
      vrais <- c(vrais, modX$loglik[2])
      diff <- abs(vrais[length(vrais)] - vrais[length(vrais) - 2])
    }
  
  } else {

    ### ELSE IF (All variables have 0 or 1 effect)
    
    modX <- coxph(eval(parse(text=substr(tt2, 13, nchar(tt2)-1))), method="efron")
    vrais <- c(modX$loglik[2])
  }
  rm(QWR, X, matX)
  gc() 

  ### Store parameters estimated ###

  MAT <- matrix(ncol=V, nrow=1+m+p+m+p+1) 
  if (i1 != 0) { #if no NL and TD
    for (j in 1:i1){
      MAT[1, j] <- modX$coef[j]
    }
  }
  if (i2 != 0){
    for (j in 1:i2){
      MAT[2:(m+p+1), i1+j] <- modX$coef[(i1 + (j-1)*(m+p) + 1):(i1 + (j-1)*(m+p) + m+p)]
    }
  }
  if (i3 != 0){
    for (j in 1:i3){
      MAT[(m+p+2):(2*m+2*p+2), i1+i2+j] <- modX$coef[(i1 + i2*(m+p) + (j-1)*(m+p+1) + 1):
                                                     (i1 + i2*(m+p) + (j-1)*(m+p+1) + m+p+1)]
    }
  }
  if (i4 != 0){
    for (j in 1:i4){
      MAT[2:(m+p+1), i1+i2+i3+j] <- mod$coef[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p) + 1):
                                             (i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p) + m+p)]
      MAT[(m+p+2):(2*(m+p+1)), i1+i2+i3+j] <- modX$coef[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + 1):
                                                        (i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + m+p+1)]
    }
  }

  ### Store SE estimated ###
  
  MATse <- matrix(ncol=V, nrow=1+m+p+m+p+1) 
  if (i1 != 0){
    for (j in 1:i1){
      MATse[1, j] <- sqrt(diag(modX$var)[j])
    }
  }
  if (i2 != 0){
    for (j in 1:i2){
      MATse[2:(m+p+1), i1+j] <- sqrt(diag(modX$var)[(i1 + (j-1)*(m+p) + 1):(i1 + (j-1)*(m+p) + m+p)])
    }
  }
  if (i3 != 0){
    for (j in 1:i3){
      MATse[(m+p+2):(2*m+2*p+2), i1+i2+j] <- sqrt(diag(modX$var)[(i1 + i2*(m+p) + (j-1)*(m+p+1) + 1):
                                                                 (i1 + i2*(m+p) + (j-1)*(m+p+1) + m+p+1)])
    }
  }
  if (i4 != 0){
    for (j in 1:i4){
      MATse[2:(m+p+1), i1+i2+i3+j] <- sqrt(diag(mod$var)[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p) + 1):
                                                         (i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p) + m+p)])
      MATse[(m+p+2):(2*(m+p+1)), i1+i2+i3+j] <- sqrt(diag(modX$var)[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + 1):
                                                                    (i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + m+p+1)])
    }
  }

  ### Prepare output ###
  
  var_order <- c(variablesNEW[NL==0 & TD==0], variablesNEW[NL==1 & TD==0], variablesNEW[NL==0 & TD==1], 
                 variablesNEW[NL==1 & TD==1])
  
  coefficients <- MAT[1, match(variablesNEW, var_order)] 
  names(coefficients) <- variables
  
  se_coef <- MATse[1, match(variablesNEW,var_order)] 
  
  coefficients_splines_NL <- as.matrix(MAT[2:(m+p+1), match(variablesNEW,var_order)])
  coefficients_splines_NL <- rbind(rep(0,V), coefficients_splines_NL)
  coefficients_splines_NL[1,(NL==0)] <- NA
  colnames(coefficients_splines_NL) <- variables
  
  coefficients_splines_TD <- as.matrix(MAT[(m+p+2):(2*(m+p+1)), match(variablesNEW,var_order)])
  colnames(coefficients_splines_TD) <- variables

  knots_covariates <- knotsNEW[1:V,]
  if (V>1){rownames(knots_covariates) <- variables}
  if (V>1){
    knots_covariates[(NL==0),] <- rep(NA, p+1+m+p+1)
  } else { 
    if (NL==0) {knots_covariates <- rep(NA, p+1+m+p+1)}
  }
  knots_time <- knotsNEW[V+1,]
  
  nEvents <- sum(data[, TypeNEW[3]]==1) 
  
  rm(data, modX)
  gc()
  
  list(Partial_Log_Likelihood = vrais[length(vrais)], 
       Number_of_parameters = nbpara + nbonlyNL*(m+p) + nbonlyTD*(1+m+p) + nbNLTD * (m+p+m+p+1), 
       Number_events = nEvents, 
       Number_knots = m, 
       Degree_of_splines = p, 
       knots_covariates = knots_covariates, 
       knots_time = knots_time, 
       coefficients = coefficients, 
       Standard_Error = se_coef, 
       coefficients_splines_NL = coefficients_splines_NL, 
       coefficients_splines_TD = coefficients_splines_TD, 
       variables = variables)
}



TDestim.sim <- function(x, m, p, coefs, knots){
  tdfct <- 0
  for (k in 1:(m+p+1)){
    tdfct <- tdfct + coefs[k] * spli(x, k, p, knots)
  }
  return(tdfct)
}



NLestim.sim <- function(x, m, p, coefs, knots){
  nlfct <- 0
  for (k in 1:(m+p+1)){
    nlfct <- nlfct + coefs[k] * spli(x, k, p, knots)
  }
  return(nlfct)
}



TELestim.sim <- function(x, m, p, coefs, knots){
  telfct <- 0
  for (k in 1:(m+p+1)){
    telfct <- telfct + coefs[k] * spli(x, k, p, knots)
  }
  return(telfct)
}



spli <- function(x, j, p, knots){
  if (p == 0){
    b <- ifelse(x >= knots[j] & x < knots[j+1], 1, 0)
    return(b)
  } else {
    a1 <- ifelse(rep(knots[j] != knots[j+p], length(x)), (x - knots[j])/(knots[j+p] - knots[j]), 0)
    a2 <- ifelse(rep(knots[j+p+1] != knots[j+1], length(x)), 
                 (knots[j+p+1] - x)/(knots[j+p+1] - knots[j+1]), 0)
    return(a1 * spli(x, j, p-1, knots) + a2 * spli(x, j+1, p-1, knots))
  }
}


DvlpMatrix <- function(data, listeT, ncol, TypeNEW){
  data <- matrix(data, ncol=ncol) # data is a list of full data on each individual, first transform this list in a matrix
  
  if (max(data[, TypeNEW[2]]) < min(listeT[listeT != 0])){  
    XX <- data  # if max(stop time) is less then min(event time), then data remain the same

  } else { # if max(stop time of certain individual) > min(event time)
    aindex <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) - 1))
    for (i in 1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1)){ # for the length of unique event times 
      # for each row of an individual, if start time < ith order of event time<=stop time of row j, stores
      for (j in 1:(dim(data)[1])){ 
        # the largest row satisfy this condition in aindex[i]
        if (as.numeric(data[j, TypeNEW[1]]) < as.numeric(listeT[1+i]) & as.numeric(data[j, TypeNEW[2]]) >= as.numeric(listeT[1+i]))
          aindex[i] <- j
      }
    }   

    # XX: first column repeats ID, column for start day assigns c(0, unique event time-last one), 
      #   column for stop day assigns unique event time 
    XX <- matrix(nrow = sum(listeT <= max(data[, TypeNEW[2]])) - 1, ncol=ncol)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), 1] <- 
        rep(data[1, 1], sum(listeT <= max(data[, TypeNEW[2]])) - 1)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[1]] <- 
        listeT[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1)]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[2]] <- 
        listeT[2:(sum(listeT <= max(data[, TypeNEW[2]])))]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[3]] <- 
        c(rep(0,(sum(listeT <= max(data[, TypeNEW[2]])) - 2)), data[dim(data)[1],TypeNEW[3]])
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])] <- 
        as.matrix(data[aindex, -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])])
  }
  X <- XX
  list(X)
}

