simPTCMGH <- function(n,
                      seed,
                      hstr,
                      baseline,
                      des_theta = NULL,
                      des_t = NULL,
                      des_h = NULL,
                      des = NULL,
                      par_base = NULL,
                      alpha = NULL,
                      beta_t = NULL,
                      beta_h = NULL,
                      beta = NULL) { 
  if (!is.null(des)) 
    des <- as.matrix(des)
  if (!is.null(des_h)) 
    des_h <- as.matrix(des_h)
  if (!is.null(des_t)) 
    des_t <- as.matrix(des_t)
  # Baseline hazards
  if (baseline == "LN") 
    quantf <- function(p) qlnorm(p, par_base[1], par_base[2])
  if (baseline == "LL") 
    quantf <- function(p) qllogis(p, par_base[1], par_base[2])
  if (baseline == "G") 
    quantf <- function(p) qgamma(p, par_base[1], par_base[2])
  if (baseline == "W") 
    quantf <- function(p) qweibull(p, par_base[1], par_base[2])
  if (baseline == "PGW") 
    quantf <- function(p) qpgw(p, par_base[1], par_base[2], par_base[3])
  if (baseline == "EW") 
    quantf <- function(p) qew(p, par_base[1], par_base[2], par_base[3])
  if (baseline == "GG") 
    quantf <- function(p) qggamma(p, par_base[1], par_base[2], par_base[3])
  
  
  # Simulated uniform random variates
  set.seed(seed)
  u0 = runif(n)
  # Transformed random variates  
  u  = -as.vector(exp(-des_theta%*%alpha))*log(u0)
  
  times <- rep(Inf, n)
  idx <- which(u < 1)
  
  if (length(idx) > 0) {
    # General Hazards model
    if (hstr == "GH") {
      exp.xbeta_t <- exp(des_t %*% beta_t)
      exp.dif <- exp(des_t %*% beta_t - des_h %*% beta_h)
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif[idx]))
      times[idx] <- as.vector(quantf(p0)/exp.xbeta_t[idx])
    }
    # Proportional Hazards model
    if (hstr == "PH") {
      exp.xbeta_t <- 1
      exp.dif <- exp(-des %*% beta)
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif[idx]))
      times[idx] <- as.vector(quantf(p0)/exp.xbeta_t[idx])
    }
    # Accelerated Failure Time model
    if (hstr == "AFT") {
      exp.xbeta_t <- exp(des %*% beta)
      exp.dif <- 1
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif))
      times[idx] <- as.vector(quantf(p0)/exp.xbeta_t[idx])
    }
    # Accelerated Hazards model
    if (hstr == "AH") {
      exp.xbeta_t <- exp(des %*% beta)
      exp.dif <- exp(des %*% beta)
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif[idx]))
      times[idx] <- as.vector(quantf(p0)/exp.xbeta_t)
    }
  }
  
  # Output
  return(as.vector(times))
}




library(HazReg)


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
set.seed(seed)
des0 <- cbind(1, rnorm(n), rnorm(n))
seed = 123
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
