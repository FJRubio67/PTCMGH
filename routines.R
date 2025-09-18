# Function to simulate from a promotion-time cure model with GH structure

simPTCMGH <- function(n,
                      seed = 123,
                      hstr = NULL,
                      baseline = NULL,
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



# IN PROGRESS
# Function to find the MLE

PTCMMLE <- function(init,
                    times,
                    status,
                    hstr = NULL,
                    baseline =NULL,
                    des_theta = NULL,
                    des_t = NULL,
                    des_h = NULL,
                    des = NULL,
                    method = "Nelder-Mead", 
                    maxit = 100) 
{
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  times.obs <- times[status]
  if (!is.null(des)) 
    des <- as.matrix(des)
  if (!is.null(des_h)) 
    des_h <- as.matrix(des_h)
  if (!is.null(des_t)) 
    des_t <- as.matrix(des_t)
  
  
  if (!is.null(des)) {
    des <- as.matrix(des)
    des.obs <- des[status, ]
  }
  
  
  if (!is.null(des_t)) {
    des_t <- as.matrix(des_t)
    des_t.obs <- des_t[status, ]
  }
  
  
  
  if (hstr == "baseline") {
    if (dist == "PGW") {
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE)
        val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, 
                                       ce0))
        return(val)
      }
    }
    if (dist == "EW") {
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE)
        val <- -sum(lhaz0) + sum(chew(times, ae0, be0, 
                                      ce0))
        return(val)
      }
    }
    if (dist == "GenGamma") {
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE)
        val <- -sum(lhaz0) + sum(chggamma(times, ae0, 
                                          be0, ce0))
        return(val)
      }
    }
    if (dist == "LogNormal") {
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chlnorm(times, ae0, 
                                         be0))
        return(val)
      }
    }
    if (dist == "LogLogistic") {
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chllogis(times, ae0, 
                                          be0))
        return(val)
      }
    }
    if (dist == "Gamma") {
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chgamma(times, ae0, 
                                         be0))
        return(val)
      }
    }
    if (dist == "Weibull") {
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chweibull(times, ae0, 
                                           be0))
        return(val)
      }
    }
  }
  if (hstr == "PH") {
    if (dist == "PGW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hpgw(times.obs, ae0, be0, ce0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chpgw(times, ae0, be0, 
                                       ce0) * exp.x.beta)
        return(val)
      }
    }
    if (dist == "EW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hew(times.obs, ae0, be0, ce0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chew(times, ae0, be0, 
                                      ce0) * exp.x.beta)
        return(val)
      }
    }
    if (dist == "GenGamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hggamma(times.obs, ae0, be0, ce0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chggamma(times, ae0, 
                                          be0, ce0) * exp.x.beta)
        return(val)
      }
    }
    if (dist == "LogNormal") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hlnorm(times.obs, ae0, be0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chlnorm(times, ae0, 
                                         be0) * exp.x.beta)
        return(val)
      }
    }
    if (dist == "LogLogistic") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hllogis(times.obs, ae0, be0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chllogis(times, ae0, 
                                          be0) * exp.x.beta)
        return(val)
      }
    }
    if (dist == "Gamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hgamma(times.obs, ae0, be0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chgamma(times, ae0, 
                                         be0) * exp.x.beta)
        return(val)
      }
    }
    if (dist == "Weibull") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        lhaz0 <- hweibull(times.obs, ae0, be0, log = TRUE) + 
          x.beta.obs
        val <- -sum(lhaz0) + sum(chweibull(times, ae0, 
                                           be0) * exp.x.beta)
        return(val)
      }
    }
  }
  if (hstr == "AFT") {
    if (dist == "PGW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hpgw(times.obs * exp.x.beta.obs, ae0, 
                      be0, ce0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chpgw(times * exp.x.beta, 
                                       ae0, be0, ce0))
        return(val)
      }
    }
    if (dist == "EW") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hew(times.obs * exp.x.beta.obs, ae0, 
                     be0, ce0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chew(times * exp.x.beta, 
                                      ae0, be0, ce0))
        return(val)
      }
    }
    if (dist == "GenGamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        beta <- par[4:(3 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hggamma(times.obs * exp.x.beta.obs, 
                         ae0, be0, ce0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chggamma(times * exp.x.beta, 
                                          ae0, be0, ce0))
        return(val)
      }
    }
    if (dist == "LogNormal") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hlnorm(times.obs * exp.x.beta.obs, ae0, 
                        be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chlnorm(times * exp.x.beta, 
                                         ae0, be0))
        return(val)
      }
    }
    if (dist == "LogLogistic") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hllogis(times.obs * exp.x.beta.obs, 
                         ae0, be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chllogis(times * exp.x.beta, 
                                          ae0, be0))
        return(val)
      }
    }
    if (dist == "Gamma") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hgamma(times.obs * exp.x.beta.obs, ae0, 
                        be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chgamma(times * exp.x.beta, 
                                         ae0, be0))
        return(val)
      }
    }
    if (dist == "Weibull") {
      p <- ncol(des)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        beta <- par[3:(2 + p)]
        x.beta <- des %*% beta
        x.beta.obs <- x.beta[status]
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.obs <- exp.x.beta[status]
        lhaz0 <- hweibull(times.obs * exp.x.beta.obs, 
                          ae0, be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chweibull(times * exp.x.beta, 
                                           ae0, be0))
        return(val)
      }
    }
  }
  if (hstr == "AH") {
    if (dist == "PGW") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hpgw(times.obs * exp.x.alpha.obs, ae0, 
                      be0, ce0, log = TRUE)
        val <- -sum(lhaz0) + sum(chpgw(times * exp.x.alpha, 
                                       ae0, be0, ce0)/exp.x.alpha)
        return(val)
      }
    }
    if (dist == "EW") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hew(times.obs * exp.x.alpha.obs, ae0, 
                     be0, ce0, log = TRUE)
        val <- -sum(lhaz0) + sum(chew(times * exp.x.alpha, 
                                      ae0, be0, ce0)/exp.x.alpha)
        return(val)
      }
    }
    if (dist == "GenGamma") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hggamma(times.obs * exp.x.alpha.obs, 
                         ae0, be0, ce0, log = TRUE)
        val <- -sum(lhaz0) + sum(chggamma(times * exp.x.alpha, 
                                          ae0, be0, ce0)/exp.x.alpha)
        return(val)
      }
    }
    if (dist == "LogNormal") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hlnorm(times.obs * exp.x.alpha.obs, 
                        ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chlnorm(times * exp.x.alpha, 
                                         ae0, be0)/exp.x.alpha)
        return(val)
      }
    }
    if (dist == "LogLogistic") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hllogis(times.obs * exp.x.alpha.obs, 
                         ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chllogis(times * exp.x.alpha, 
                                          ae0, be0)/exp.x.alpha)
        return(val)
      }
    }
    if (dist == "Gamma") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hgamma(times.obs * exp.x.alpha.obs, 
                        ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chgamma(times * exp.x.alpha, 
                                         ae0, be0)/exp.x.alpha)
        return(val)
      }
    }
    if (dist == "Weibull") {
      p <- ncol(des_t)
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p)]
        exp.x.alpha <- as.vector(exp(des_t %*% alpha))
        exp.x.alpha.obs <- exp.x.alpha[status]
        lhaz0 <- hweibull(times.obs * exp.x.alpha.obs, 
                          ae0, be0, log = TRUE)
        val <- -sum(lhaz0) + sum(chweibull(times * exp.x.alpha, 
                                           ae0, be0)/exp.x.alpha)
        return(val)
      }
    }
  }
  if (hstr == "GH") {
    if (dist == "PGW") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p0)]
        beta <- par[(4 + p0):(3 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x.beta <- des %*% beta
        x.alpha.obs <- x.alpha[status]
        x.beta.obs <- x.beta[status]
        exp.x.alpha <- as.vector(exp(x.alpha))
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.dif <- as.vector(exp(x.beta - x.alpha))
        exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
        exp.x.beta.obs <- as.vector(exp.x.beta[status])
        lhaz0 <- hpgw(times.obs * exp.x.alpha.obs, ae0, 
                      be0, ce0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chpgw(times * exp.x.alpha, 
                                       ae0, be0, ce0) * exp.x.beta.dif)
        return(sum(val))
      }
    }
    if (dist == "EW") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p0)]
        beta <- par[(4 + p0):(3 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x.beta <- des %*% beta
        x.alpha.obs <- x.alpha[status]
        x.beta.obs <- x.beta[status]
        exp.x.alpha <- as.vector(exp(x.alpha))
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.dif <- as.vector(exp(x.beta - x.alpha))
        exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
        exp.x.beta.obs <- as.vector(exp.x.beta[status])
        lhaz0 <- hew(times.obs * exp.x.alpha.obs, ae0, 
                     be0, ce0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chew(times * exp.x.alpha, 
                                      ae0, be0, ce0) * exp.x.beta.dif)
        return(sum(val))
      }
    }
    if (dist == "GenGamma") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        ce0 <- exp(par[3])
        alpha <- par[4:(3 + p0)]
        beta <- par[(4 + p0):(3 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x.beta <- des %*% beta
        x.alpha.obs <- x.alpha[status]
        x.beta.obs <- x.beta[status]
        exp.x.alpha <- as.vector(exp(x.alpha))
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.dif <- as.vector(exp(x.beta - x.alpha))
        exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
        exp.x.beta.obs <- as.vector(exp.x.beta[status])
        lhaz0 <- hggamma(times.obs * exp.x.alpha.obs, 
                         ae0, be0, ce0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chggamma(times * exp.x.alpha, 
                                          ae0, be0, ce0) * exp.x.beta.dif)
        return(sum(val))
      }
    }
    if (dist == "LogNormal") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p0)]
        beta <- par[(3 + p0):(2 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x.beta <- des %*% beta
        x.alpha.obs <- x.alpha[status]
        x.beta.obs <- x.beta[status]
        exp.x.alpha <- as.vector(exp(x.alpha))
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.dif <- as.vector(exp(x.beta - x.alpha))
        exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
        exp.x.beta.obs <- as.vector(exp.x.beta[status])
        lhaz0 <- hlnorm(times.obs * exp.x.alpha.obs, 
                        ae0, be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chlnorm(times * exp.x.alpha, 
                                         ae0, be0) * exp.x.beta.dif)
        return(sum(val))
      }
    }
    if (dist == "LogLogistic") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log.lik <- function(par) {
        ae0 <- par[1]
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p0)]
        beta <- par[(3 + p0):(2 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x.beta <- des %*% beta
        x.alpha.obs <- x.alpha[status]
        x.beta.obs <- x.beta[status]
        exp.x.alpha <- as.vector(exp(x.alpha))
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.dif <- as.vector(exp(x.beta - x.alpha))
        exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
        exp.x.beta.obs <- as.vector(exp.x.beta[status])
        lhaz0 <- hllogis(times.obs * exp.x.alpha.obs, 
                         ae0, be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chllogis(times * exp.x.alpha, 
                                          ae0, be0) * exp.x.beta.dif)
        return(sum(val))
      }
    }
    if (dist == "Gamma") {
      p0 <- dim(des_t)[2]
      p1 <- dim(des)[2]
      log.lik <- function(par) {
        ae0 <- exp(par[1])
        be0 <- exp(par[2])
        alpha <- par[3:(2 + p0)]
        beta <- par[(3 + p0):(2 + p0 + p1)]
        x.alpha <- des_t %*% alpha
        x.beta <- des %*% beta
        x.alpha.obs <- x.alpha[status]
        x.beta.obs <- x.beta[status]
        exp.x.alpha <- as.vector(exp(x.alpha))
        exp.x.beta <- as.vector(exp(x.beta))
        exp.x.beta.dif <- as.vector(exp(x.beta - x.alpha))
        exp.x.alpha.obs <- as.vector(exp.x.alpha[status])
        exp.x.beta.obs <- as.vector(exp.x.beta[status])
        lhaz0 <- hgamma(times.obs * exp.x.alpha.obs, 
                        ae0, be0, log = TRUE) + x.beta.obs
        val <- -sum(lhaz0) + sum(chgamma(times * exp.x.alpha, 
                                         ae0, be0) * exp.x.beta.dif)
        return(sum(val))
      }
    }
  }
  if (method != "nlminb") {
    OPT <- optim(init, log.lik, control = list(maxit = maxit), 
                 method = method)
  }
  if (method == "nlminb") {
    OPT <- nlminb(init, log.lik, control = list(iter.max = maxit))
  }
  OUT <- list(log_lik = log.lik, OPT = OPT)
  return(OUT)
}

