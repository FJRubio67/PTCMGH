## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())

#library(devtools)
#install_github("FJRubio67/HazReg")
library(HazReg)
#library(devtools)
#install_github("FJRubio67/PTCMGH")
library(PTCMGH)
library(survival)
library(knitr)

# Data simulation
n = 10000
seed = 123
set.seed(seed)
# Design matrix
des0 <- cbind(1, rnorm(n), rnorm(n))

# Simulation: PTCM with GH structure and LogNormal baseline hazard, and log-link
sim = simPTCMGH(n = n,
          seed = seed,
          hstr = "GH",
          dist = "LogNormal",
          des_theta = des0,
          des_t = des0[, -1],
          des_h = des0[, -1],
          des = NULL,
          par_base = c(-0.25, 0.1),
          alpha = c(0.5, 0.5, 0.5),
          beta_t = c(0.5, 0.5),
          beta_h = c(1.0, 1.0),
          beta = NULL)

# Inducing censoring
cens = 6
status <- ifelse(sim < cens, 1, 0)
mean(status)

times <- ifelse(sim< cens, sim, cens)

# Kaplan-Meier estimator for the survival times
km <- survfit(Surv(times, status) ~ 1)
plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Model fitting
OPT <- PTCMMLE(init = c(0,0,0),
                    times = times,
                    status = status,
                    hstr = "baseline",
                    dist = "LogNormal",
                    des_theta = NULL,
                    des_t = NULL,
                    des_h = NULL,
                    des = NULL,
                    method = "nlminb", 
                    maxit = 10000) 

# MLE
MLE <- c(OPT$OPT$par[1],exp(OPT$OPT$par[2]),exp(OPT$OPT$par[3]))

# Fitted survival function 
spt <- Vectorize(function(t) exp(- MLE[3]*(1-exp(-chlnorm(t,MLE[1],MLE[2])))) )

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt,0,6, col = "red", add= T) 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Model fitting
OPT_PH <- PTCMMLE(init = c(0,0,0,0,0,0,0),
               times = times,
               status = status,
               hstr = "PH",
               dist = "LogNormal",
               des_theta = des0,
               des_t = NULL,
               des_h = NULL,
               des = des0[,-1],
               method = "nlminb", 
               maxit = 10000) 

# MLE
MLE_PH <- c(OPT_PH$OPT$par[1],exp(OPT_PH$OPT$par[2]),OPT_PH$OPT$par[-c(1:2)])

# Population survival function 
spt_ph <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_PH[3:5]))
  F_i <- 1 - exp(-chlnorm(t,MLE_PH[1],MLE_PH[2])*exp(des0[,-1]%*%MLE_PH[6:7]))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
  }) 

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_ph,0,6, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Model fitting
OPT_AFT <- PTCMMLE(init = c(0,0,0,0,0,0,0),
                  times = times,
                  status = status,
                  hstr = "AFT",
                  dist = "LogNormal",
                  des_theta = des0,
                  des_t = NULL,
                  des_h = NULL,
                  des = des0[,-1],
                  method = "nlminb", 
                  maxit = 10000) 

# MLE
MLE_AFT <- c(OPT_AFT$OPT$par[1],exp(OPT_AFT$OPT$par[2]),OPT_AFT$OPT$par[-c(1:2)])

# Population survival function
spt_aft <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_AFT[3:5]))
  F_i <- 1 - exp(-chlnorm(as.vector(t*exp(des0[,-1]%*%MLE_AFT[6:7])),MLE_AFT[1],MLE_AFT[2]))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
}) 


plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_aft,0,6, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Model fitting
OPT_AH <- PTCMMLE(init = c(0,0,0,0,0,0,0),
                   times = times,
                   status = status,
                   hstr = "AH",
                   dist = "LogNormal",
                   des_theta = des0,
                   des_t = des0[,-1],
                   des_h = NULL,
                   des = NULL,
                   method = "nlminb", 
                   maxit = 10000) 

MLE
MLE_AH <- c(OPT_AH$OPT$par[1],exp(OPT_AH$OPT$par[2]),OPT_AH$OPT$par[-c(1:2)])

# Population survival function
spt_ah <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_AH[3:5]))
  F_i <- 1 - exp(-chlnorm(as.vector(t*exp(des0[,-1]%*%MLE_AH[6:7])),MLE_AH[1],MLE_AH[2])*as.vector(exp(-des0[,-1]%*%MLE_AH[6:7])))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
}) 

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_ah,0,6, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Model fitting
OPT_GH <- PTCMMLE(init = c(0,0,0,0,0,0,0,0,0),
                   times = times,
                   status = status,
                   hstr = "GH",
                   dist = "LogNormal",
                   des_theta = des0,
                   des_t = des0[,-1],
                   des_h = des0[,-1],
                   des = NULL,
                   method = "nlminb", 
                   maxit = 10000) 

# MLE
MLE_GH <- c(OPT_GH$OPT$par[1],exp(OPT_GH$OPT$par[2]),OPT_GH$OPT$par[-c(1:2)])

# Comparison: MLE vs true parameter values
cbind(MLE_GH, c(-0.25, 0.1,0.5, 0.5, 0.5,0.5, 0.5, 1.0, 1.0))


# Population survival function
spt_gh <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_GH[3:5]))
  F_i <- 1 - exp(-chlnorm(as.vector(t*as.vector(exp(des0[,-1]%*%MLE_GH[6:7]))),MLE_GH[1],MLE_GH[2])*as.vector(exp(des0[,-1]%*%MLE_GH[8:9] - des0[,-1]%*%MLE_GH[6:7])))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
}) 

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_gh,0,6, col = "red", add= T)  

# Confidence intervals for the parameters (positive parameters transformed to log-scale)
CI_GH <- Conf_Int(OPT_GH$log_lik, OPT_GH$OPT$par, level = 0.95)
kable(CI_GH, digits = 3)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Akaike information criterion
AIC <- 2*OPT$OPT$objective + 2*length(OPT$OPT$par)
AIC_PH <- 2*OPT_PH$OPT$objective + 2*length(OPT_PH$OPT$par)
AIC_AFT <- 2*OPT_AFT$OPT$objective + 2*length(OPT_AFT$OPT$par)
AIC_AH <- 2*OPT_AH$OPT$objective + 2*length(OPT_AH$OPT$par)
AIC_GH <- 2*OPT_GH$OPT$objective + 2*length(OPT_GH$OPT$par)

# AICs: baseline, PH, AFT, AH, GH
AICs <- c(AIC, AIC_PH, AIC_AFT, AIC_AH, AIC_GH)

AICs

which.min(AICs)


