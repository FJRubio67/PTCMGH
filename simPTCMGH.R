library(HazReg)

simPTCMGH <- function(n, seed, hstr, baseline, des_theta = NULL, des_t = NULL, des_h = NULL, des = NULL, 
                      theta = NULL, beta_t = NULL, beta_h = NULL, beta = NULL){
  
  set.seed(seed)
  u <- runif(n)
  
  if(!is.null(des_theta)) des_theta = as.matrix(des_theta)
  if(!is.null(des_t)) des_t = as.matrix(des_t)
  if(!is.null(des_h)) des_h = as.matrix(des_h)
  
  theta = as.vector(theta)
  alpha = as.vector(alpha)
  beta_t = as.vector(beta_t)
  beta_h = as.vector(beta_h)
  
  if (baseline == "LN") 
    quantf <- function(p) qlnorm(p, theta[1], theta[2])
  if (baseline == "LL") 
    quantf <- function(p) qllogis(p, theta[1], theta[2])
  if (baseline == "G") 
    quantf <- function(p) qgamma(p, theta[1], theta[2])
  if (baseline == "W") 
    quantf <- function(p) qweibull(p, theta[1], theta[2])
  if (baseline == "PGW") 
    quantf <- function(p) qpgw(p, theta[1], theta[2], theta[3])
  if (baseline == "EW") 
    quantf <- function(p) qew(p, theta[1], theta[2], theta[3])
  if (baseline == "GG") 
    quantf <- function(p) qggamma(p, theta[1], theta[2], theta[3])
  
  
  if(hstr == "PH"){
    theta_i <- as.vector(exp(des_theta%*%theta))
    des_beta_h <- as.vector(des_h%*%beta_h)
    F_i <- 1 - exp(0)
  }
  
  
  
}