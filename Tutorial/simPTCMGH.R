## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())

#library(devtools)
#install_github("FJRubio67/HazReg")
library(HazReg)
library(PTCMGH)

n = 10000
seed = 123
set.seed(seed)
des0 <- cbind(1, rnorm(n), rnorm(n))


sim = simPTCMGH(n = n,
          seed = seed,
          hstr = "AFT",
          dist = "LogNormal",
          des_theta = des0,
          des_t = NULL,
          des_h = NULL,
          des = des0[, -1],
          par_base = c(-0.25, 0.1),
          alpha = c(0.5, 0.5, 0.5),
          beta_t = NULL,
          beta_h = NULL,
          beta = c(0.5, 0.5))


cens = 4
status <- ifelse(sim < cens, 1, 0)
mean(status)

times <- ifelse(sim< cens, sim, cens)

library(survival)
# Kaplan-Meier estimator for the survival times
km <- survfit(Surv(times, status) ~ 1)

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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


MLE <- c(OPT$OPT$par[1],exp(OPT$OPT$par[2]),exp(OPT$OPT$par[3]))

spt <- Vectorize(function(t) exp(- MLE[3]*(1-exp(-chlnorm(t,MLE[1],MLE[2])))) )

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt,0,4, col = "red", add= T) 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


MLE_PH <- c(OPT_PH$OPT$par[1],exp(OPT_PH$OPT$par[2]),OPT_PH$OPT$par[-c(1:2)])

spt_ph <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_PH[3:5]))
  F_i <- 1 - exp(-chlnorm(t,MLE_PH[1],MLE_PH[2])*exp(des0[,-1]%*%MLE_PH[6:7]))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
  }) 

plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_ph,0,4, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


MLE_AFT <- c(OPT_AFT$OPT$par[1],exp(OPT_AFT$OPT$par[2]),OPT_AFT$OPT$par[-c(1:2)])

cbind(MLE_AFT, c(-0.25, 0.1,0.5, 0.5, 0.5,0.5, 0.5))


spt_aft <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_AFT[3:5]))
  F_i <- 1 - exp(-chlnorm(as.vector(t*exp(des0[,-1]%*%MLE_AFT[6:7])),MLE_AFT[1],MLE_AFT[2]))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
}) 




plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_aft,0,4, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


MLE_AH <- c(OPT_AH$OPT$par[1],exp(OPT_AH$OPT$par[2]),OPT_AH$OPT$par[-c(1:2)])

cbind(MLE_AH, c(-0.25, 0.1,0.5, 0.5, 0.5,0.5, 0.5))


spt_ah <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_AH[3:5]))
  F_i <- 1 - exp(-chlnorm(as.vector(t*exp(des0[,-1]%*%MLE_AH[6:7])),MLE_AH[1],MLE_AH[2])*as.vector(exp(-des0[,-1]%*%MLE_AH[6:7])))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
}) 




plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_ah,0,4, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


MLE_GH <- c(OPT_GH$OPT$par[1],exp(OPT_GH$OPT$par[2]),OPT_GH$OPT$par[-c(1:2)])

cbind(MLE_GH, c(-0.25, 0.1,0.5, 0.5, 0.5,0.5, 0.5, 0.5, 0.5))


spt_gh <- Vectorize(function(t){
  
  theta_i <- as.vector(exp(des0 %*% MLE_GH[3:5]))
  F_i <- 1 - exp(-chlnorm(as.vector(t*as.vector(exp(des0[,-1]%*%MLE_GH[6:7]))),MLE_GH[1],MLE_GH[2])*as.vector(exp(des0[,-1]%*%MLE_GH[8:9] - des0[,-1]%*%MLE_GH[6:7])))
  survs <- exp(-theta_i*F_i)
  return(mean(survs))
}) 




plot(km$time, km$surv, type = "l", col = "black", lwd = 2, lty = 1, 
     ylim = c(0,1), xlab = "Time", ylab = "Survival")
curve(spt_gh,0,4, col = "red", add= T)  


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
AIC <- 2*OPT$OPT$objective + 2*length(OPT$OPT$par)
AIC_PH <- 2*OPT_PH$OPT$objective + 2*length(OPT_PH$OPT$par)
AIC_AFT <- 2*OPT_AFT$OPT$objective + 2*length(OPT_AFT$OPT$par)
AIC_AH <- 2*OPT_AH$OPT$objective + 2*length(OPT_AH$OPT$par)
AIC_GH <- 2*OPT_GH$OPT$objective + 2*length(OPT_GH$OPT$par)

c(AIC, AIC_PH, AIC_AFT, AIC_AH, AIC_GH)

