library(HazReg)

simPTCMGH <- function(n, seed, hstr, baseline, des_theta = NULL, des_t = NULL, des_h = NULL, des = NULL, 
                      par_base = NULL, alpha = NULL, beta_t = NULL, beta_h = NULL, beta = NULL){
  
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
  
  
  if(hstr == "PH"){
    mexp_alpha <- as.vector(exp(-des_theta%*%alpha))
    mexp_betah <- as.vector(exp(-des_h%*%beta_h))
    test_surv <- 1 + log(u)*mexp_alpha
    
    sim <- vector()
    for(i in 1:n){
      if(test_surv[i] > 0){
        arg <- 1 - exp( mexp_betah[i]*log( 1 + log(u[i])*mexp_alpha[i] ) )
        sim[i] <- quantf(arg)
      }
      else sim[i]  <- Inf
    }
    
  }
  
  
  
}