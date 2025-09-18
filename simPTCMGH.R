rm(list=ls())


library(HazReg)

source("routines.R")

#hstr = "AFT"
#baseline = "LN"
#des_theta = cbind(1, rnorm(n), rnorm(n))
#des_t = NULL
#des_h = NULL
#des = des_theta[, -1]
#par_base = c(0, 0.1)
#alpha = c(0.5, 0.5, 0.5)
#beta_t = NULL
#beta_h = NULL
#beta = c(0.5, 0.5)

n = 1000
seed = 123
set.seed(seed)
des0 <- cbind(1, rnorm(n), rnorm(n))
sim = simPTCMGH(n = n,
          seed = seed,
          hstr = "AFT",
          baseline = "LN",
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
