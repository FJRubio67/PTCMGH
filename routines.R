# Function to simulate from a promotion-time cure model with GH structure

simPTCMGH <- function(n,
                      seed = 123,
                      hstr = NULL,
                      dist = NULL,
                      des_theta = NULL,
                      des_t = NULL,
                      des_h = NULL,
                      des = NULL,
                      par_base = NULL,
                      alpha0 = NULL,
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
  if (dist == "LN") 
    quantf <- function(p) qlnorm(p, par_base[1], par_base[2])
  if (dist == "LL") 
    quantf <- function(p) qllogis(p, par_base[1], par_base[2])
  if (dist == "G") 
    quantf <- function(p) qgamma(p, par_base[1], par_base[2])
  if (dist == "W") 
    quantf <- function(p) qweibull(p, par_base[1], par_base[2])
  if (dist == "PGW") 
    quantf <- function(p) qpgw(p, par_base[1], par_base[2], par_base[3])
  if (dist == "EW") 
    quantf <- function(p) qew(p, par_base[1], par_base[2], par_base[3])
  if (dist == "GG") 
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
      exp_xbeta_t <- exp(des_t %*% beta_t)
      exp.dif <- exp(des_t %*% beta_t - des_h %*% beta_h)
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif[idx]))
      times[idx] <- as.vector(quantf(p0)/exp_xbeta_t[idx])
    }
    # Proportional Hazards model
    if (hstr == "PH") {
      exp_xbeta_t <- 1
      exp.dif <- exp(-des %*% beta)
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif[idx]))
      times[idx] <- as.vector(quantf(p0)/exp_xbeta_t[idx])
    }
    # Accelerated Failure Time model
    if (hstr == "AFT") {
      exp_xbeta_t <- exp(des %*% beta)
      exp.dif <- 1
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif))
      times[idx] <- as.vector(quantf(p0)/exp_xbeta_t[idx])
    }
    # Accelerated Hazards model
    if (hstr == "AH") {
      exp_xbeta_t <- exp(des %*% beta)
      exp.dif <- exp(des %*% beta)
      p0 <- as.vector(1 - exp(log(1 - u[idx]) * exp.dif[idx]))
      times[idx] <- as.vector(quantf(p0)/exp_xbeta_t)
    }
  }
  
  # Output
  return(as.vector(times))
}


#################################################################################
#################################################################################
# IN PROGRESS
# Function to find the MLE
#################################################################################
#################################################################################

PTCMMLE <- function(init,
                    times,
                    status,
                    hstr = NULL,
                    dist =NULL,
                    des_theta = NULL,
                    des_t = NULL,
                    des_h = NULL,
                    des = NULL,
                    method = "Nelder-Mead", 
                    maxit = 100) 
{
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  times_obs <- times[status]
  n_obs <- sum(status)
  n <- length(times)
  
  if (!is.null(des_h)) 
    des_h <- as.matrix(des_h)

  if (!is.null(des)) {
    des <- as.matrix(des)
    des_obs <- des[status, ]
  }
  
  if (!is.null(des_theta)) {
    des_theta <- as.matrix(des_theta)
    des_theta_obs <- des_theta[status, ]
  } 
  
  if (!is.null(des_t)) {
    des_t <- as.matrix(des_t)
    des_t_obs <- des_t[status, ]
  }
  
  #----------------------------------------------------------------------------
  # Baseline
  #----------------------------------------------------------------------------
  if (hstr == "baseline") {
    if (dist == "PGW") {
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        theta = exp(par[4])
        alpha0 = par[4]
        lhaz0_obs <- hpgw(times_obs, ae0, be0, ce0, log = TRUE)
        chaz0 <- chpgw(times, ae0, be0, ce0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)
      }
    }
    if (dist == "EW") {
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        theta = exp(par[4])
        alpha0 = par[4]
        lhaz0_obs <- hew(times_obs, ae0, be0, ce0, log = TRUE)
        chaz0 <- chew(times, ae0, be0, ce0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)

      }
    }
    if (dist == "GenGamma") {
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        theta = exp(par[4])
        alpha0 = par[4]
        lhaz0_obs <- hggamma(times_obs, ae0, be0, ce0, log = TRUE)
        chaz0 <- chggamma(times, ae0, be0, ce0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)
      }
    }
    if (dist == "LogNormal") {
      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        theta = exp(par[3])
        alpha0 = par[3]
        lhaz0_obs <- hlnorm(times_obs, ae0, be0, log = TRUE)
        chaz0 <- chlnorm(times, ae0, be0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)
      }
    }
    if (dist == "LogLogistic") {
      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        theta = exp(par[3])
        alpha0 = par[3]
        lhaz0_obs <- hllogis(times_obs, ae0, be0, log = TRUE)
        chaz0 <- chllogis(times, ae0, be0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)
      }
    }
    if (dist == "Gamma") {
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        theta = exp(par[3])
        alpha0 = par[3]
        lhaz0_obs <- hgamma(times_obs, ae0, be0, log = TRUE)
        chaz0 <- chgamma(times, ae0, be0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)
      }
    }
    if (dist == "Weibull") {
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        theta = exp(par[3])
        alpha0 = par[3]
        lhaz0_obs <- hweibull(times_obs, ae0, be0, log = TRUE)
        chaz0 <- chweibull(times, ae0, be0)
        val <- n_obs*alpha0 - n*theta + sum(lhaz0_obs) - sum(chaz0[status]) + 
          theta*sum(exp(-chaz0))
        return(-val)
      }
    }
  }
  
  #----------------------------------------------------------------------------
  # Proportional Hazards
  #----------------------------------------------------------------------------
  if (hstr == "PH") {    
    
    p <- ncol(des)
    p_theta <- ncol(des_theta)
    
    if (dist == "PGW") {
      
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        beta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hpgw(times_obs, ae0, be0, ce0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chpgw(times, ae0, be0, ce0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
          
        return(-val)
      }
    }
    if (dist == "EW") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        beta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hew(times_obs, ae0, be0, ce0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chew(times, ae0, be0, ce0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "GenGamma") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        beta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hggamma(times_obs, ae0, be0, ce0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chggamma(times, ae0, be0, ce0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
       }
    }
    if (dist == "LogNormal") {

      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hlnorm(times_obs, ae0, be0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chlnorm(times, ae0, be0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
        
      }
    }
    if (dist == "LogLogistic") {

      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hllogis(times_obs, ae0, be0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chllogis(times, ae0, be0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "Gamma") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hgamma(times_obs, ae0, be0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chgamma(times, ae0, be0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "Weibull") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hweibull(times_obs, ae0, be0, log = TRUE) + 
          x_beta_obs
        chaz0 <- chweibull(times, ae0, be0)*exp_x_beta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
  }
  
  #----------------------------------------------------------------------------
  # Accelerated Failure Time (IN PROGRESS)
  #----------------------------------------------------------------------------
  if (hstr == "AFT") {
    
    p <- ncol(des)
    p_theta <- ncol(des_theta)
    
    if (dist == "PGW") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        beta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hpgw(times_obs * exp_x_beta_obs, ae0, 
                         be0, ce0, log = TRUE) + x_beta_obs
        chaz0 <- chpgw(times * exp_x_beta, ae0, be0, ce0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "EW") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        beta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hew(times_obs * exp_x_beta_obs, ae0, 
                          be0, ce0, log = TRUE) + x_beta_obs
        chaz0 <- chew(times * exp_x_beta, ae0, be0, ce0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "GenGamma") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        beta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hggamma(times_obs * exp_x_beta_obs, ae0, 
                         be0, ce0, log = TRUE) + x_beta_obs
        chaz0 <- chggamma(times * exp_x_beta, ae0, be0, ce0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "LogNormal") {

      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hlnorm(times_obs * exp_x_beta_obs, ae0, 
                             be0, log = TRUE) + x_beta_obs
        chaz0 <- chlnorm(times * exp_x_beta, ae0, be0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "LogLogistic") {

      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hllogis(times_obs * exp_x_beta_obs, ae0, 
                            be0, log = TRUE) + x_beta_obs
        chaz0 <- chllogis(times * exp_x_beta, ae0, be0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "Gamma") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hgamma(times_obs * exp_x_beta_obs, ae0, 
                             be0, log = TRUE) + x_beta_obs
        chaz0 <- chgamma(times * exp_x_beta, ae0, be0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "Weibull") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        beta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_beta <- des %*% beta
        x_beta_obs <- x_beta[status]
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta_obs <- exp_x_beta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hweibull(times_obs * exp_x_beta_obs, ae0, 
                            be0, log = TRUE) + x_beta_obs
        chaz0 <- chweibull(times * exp_x_beta, ae0, be0)
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
  }
  
  #----------------------------------------------------------------------------
  # Accelerated Hazards
  #----------------------------------------------------------------------------
  if (hstr == "AH") {
    
    p <- ncol(des_t)
    p_theta <- ncol(des_theta)
    
    if (dist == "PGW") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        eta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hpgw(times_obs * exp_x_eta_obs, ae0, 
                          be0, ce0, log = TRUE)
        chaz0 <- chpgw(times * exp_x_eta, 
                       ae0, be0, ce0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "EW") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        eta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hew(times_obs * exp_x_eta_obs, ae0, 
                          be0, ce0, log = TRUE)
        chaz0 <- chew(times * exp_x_eta, 
                       ae0, be0, ce0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
     }
    }
    if (dist == "GenGamma") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha0 = par[4:(3 + p_theta)]
        eta <- par[(4 + p_theta):(3 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hggamma(times_obs * exp_x_eta_obs, ae0, 
                          be0, ce0, log = TRUE)
        chaz0 <- chggamma(times * exp_x_eta, 
                       ae0, be0, ce0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "LogNormal") {

      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        eta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hlnorm(times_obs * exp_x_eta_obs, ae0, 
                          be0, log = TRUE)
        chaz0 <- chlnorm(times * exp_x_eta, 
                       ae0, be0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "LogLogistic") {

      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        eta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hllogis(times_obs * exp_x_eta_obs, ae0, 
                          be0, log = TRUE)
        chaz0 <- chllogis(times * exp_x_eta, 
                       ae0, be0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))

        return(-val)
      }
    }
    if (dist == "Gamma") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        eta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hgamma(times_obs * exp_x_eta_obs, ae0, 
                          be0, log = TRUE)
        chaz0 <- chgamma(times * exp_x_eta, 
                       ae0, be0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
    if (dist == "Weibull") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha0 = par[3:(2 + p_theta)]
        eta <- par[(3 + p_theta):(2 + p_theta + p)]
        x_eta <- des_t %*% eta
        x_eta_obs <- x_eta[status]
        exp_x_eta <- as.vector(exp(x_eta))
        exp_x_eta_obs <- exp_x_eta[status]
        w_theta <- des_theta %*% alpha0
        theta <- as.vector(exp(w_theta))
        theta_obs <- theta[status]
        lhaz0_obs <- hweibull(times_obs * exp_x_eta_obs, ae0, 
                          be0, log = TRUE)
        chaz0 <- chweibull(times * exp_x_eta, 
                       ae0, be0)/exp_x_eta
        val <- sum(w_theta[status]) + sum(lhaz0_obs) - sum(chaz0[status]) - 
          sum(theta) + sum(exp( w_theta - chaz0 ))
        
        return(-val)
      }
    }
  }
  
  #----------------------------------------------------------------------------
  # General Hazards
  #----------------------------------------------------------------------------
  if (hstr == "GH") {

    p0 <- dim(des_t)[2]
    p1 <- dim(des)[2]
    p_theta <- ncol(des_theta)
    
    if (dist == "PGW") {

      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        
        
        
        alpha <- par[4:(3 + p0)]
        beta <- par[(4 + p0):(3 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x_beta <- des %*% beta
        x.alpha_obs <- x.alpha[status]
        x_beta_obs <- x_beta[status]
        exp_x.alpha <- as.vector(exp(x.alpha))
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta.dif <- as.vector(exp(x_beta - x.alpha))
        exp_x.alpha_obs <- as.vector(exp_x.alpha[status])
        exp_x_beta_obs <- as.vector(exp_x_beta[status])
        lhaz0 <- hpgw(times_obs * exp_x.alpha_obs, ae0, 
                      be0, ce0, log = TRUE) + x_beta_obs
        val <- -sum(lhaz0) + sum(chpgw(times * exp_x.alpha, 
                                       ae0, be0, ce0) * exp_x_beta.dif)
        return(sum(val))
      }
    }
    if (dist == "EW") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p0)]
        beta <- par[(4 + p0):(3 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x_beta <- des %*% beta
        x.alpha_obs <- x.alpha[status]
        x_beta_obs <- x_beta[status]
        exp_x.alpha <- as.vector(exp(x.alpha))
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta.dif <- as.vector(exp(x_beta - x.alpha))
        exp_x.alpha_obs <- as.vector(exp_x.alpha[status])
        exp_x_beta_obs <- as.vector(exp_x_beta[status])
        lhaz0 <- hew(times_obs * exp_x.alpha_obs, ae0, 
                     be0, ce0, log = TRUE) + x_beta_obs
        val <- -sum(lhaz0) + sum(chew(times * exp_x.alpha, 
                                      ae0, be0, ce0) * exp_x_beta.dif)
        return(sum(val))
      }
    }
    if (dist == "GenGamma") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p0)]
        beta <- par[(4 + p0):(3 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x_beta <- des %*% beta
        x.alpha_obs <- x.alpha[status]
        x_beta_obs <- x_beta[status]
        exp_x.alpha <- as.vector(exp(x.alpha))
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta.dif <- as.vector(exp(x_beta - x.alpha))
        exp_x.alpha_obs <- as.vector(exp_x.alpha[status])
        exp_x_beta_obs <- as.vector(exp_x_beta[status])
        lhaz0 <- hggamma(times_obs * exp_x.alpha_obs, 
                         ae0, be0, ce0, log = TRUE) + x_beta_obs
        val <- -sum(lhaz0) + sum(chggamma(times * exp_x.alpha, 
                                          ae0, be0, ce0) * exp_x_beta.dif)
        return(sum(val))
      }
    }
    if (dist == "LogNormal") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p0)]
        beta <- par[(3 + p0):(2 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x_beta <- des %*% beta
        x.alpha_obs <- x.alpha[status]
        x_beta_obs <- x_beta[status]
        exp_x.alpha <- as.vector(exp(x.alpha))
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta.dif <- as.vector(exp(x_beta - x.alpha))
        exp_x.alpha_obs <- as.vector(exp_x.alpha[status])
        exp_x_beta_obs <- as.vector(exp_x_beta[status])
        lhaz0 <- hlnorm(times_obs * exp_x.alpha_obs, 
                        ae0, be0, log = TRUE) + x_beta_obs
        val <- -sum(lhaz0) + sum(chlnorm(times * exp_x.alpha, 
                                         ae0, be0) * exp_x_beta.dif)
        return(sum(val))
      }
    }
    if (dist == "LogLogistic") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log_lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p0)]
        beta <- par[(3 + p0):(2 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x_beta <- des %*% beta
        x.alpha_obs <- x.alpha[status]
        x_beta_obs <- x_beta[status]
        exp_x.alpha <- as.vector(exp(x.alpha))
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta.dif <- as.vector(exp(x_beta - x.alpha))
        exp_x.alpha_obs <- as.vector(exp_x.alpha[status])
        exp_x_beta_obs <- as.vector(exp_x_beta[status])
        lhaz0 <- hllogis(times_obs * exp_x.alpha_obs, 
                         ae0, be0, log = TRUE) + x_beta_obs
        val <- -sum(lhaz0) + sum(chllogis(times * exp_x.alpha, 
                                          ae0, be0) * exp_x_beta.dif)
        return(sum(val))
      }
    }
    if (dist == "Gamma") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log_lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p0)]
        beta <- par[(3 + p0):(2 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x_beta <- des %*% beta
        x.alpha_obs <- x.alpha[status]
        x_beta_obs <- x_beta[status]
        exp_x.alpha <- as.vector(exp(x.alpha))
        exp_x_beta <- as.vector(exp(x_beta))
        exp_x_beta.dif <- as.vector(exp(x_beta - x.alpha))
        exp_x.alpha_obs <- as.vector(exp_x.alpha[status])
        exp_x_beta_obs <- as.vector(exp_x_beta[status])
        lhaz0 <- hgamma(times_obs * exp_x.alpha_obs, 
                        ae0, be0, log = TRUE) + x_beta_obs
        val <- -sum(lhaz0) + sum(chgamma(times * exp_x.alpha, 
                                         ae0, be0) * exp_x_beta.dif)
        return(sum(val))
      }
    }
  }
  if (method != "nlminb") {
    OPT <- optim(init, log_lik, control = list(maxit = maxit), 
                 method = method)
  }
  if (method == "nlminb") {
    OPT <- nlminb(init, log_lik, control = list(iter.max = maxit))
  }
  OUT <- list(log_lik = log_lik, OPT = OPT)
  return(OUT)
}

