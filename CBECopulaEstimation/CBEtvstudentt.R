library(stats)

# Compute the negative copula log-likelihood of a bivariate t distribution
# with time-varying correlation
bivt_tvp1_CL <- function(theta, Zdata, rhobar) {
  
  T <- nrow(Zdata)
  nu <- theta[4]
  x <- qt(Zdata[, 1], nu)
  y <- qt(Zdata[, 2], nu)
  
  kappa <- rep(-999.99, T)
  kappa[1] <- rhobar  # this is the MLE of kappa in the time-invariant version of this model
  for (jj in 2:T) {
    if (jj <= 10) {
      psi <- theta[1] + theta[2] * mean(x[1:(jj - 1)] * y[1:(jj - 1)]) + theta[3] * kappa[jj - 1]
    } else {
      psi <- theta[1] + theta[2] * mean(x[(jj - 10):(jj - 1)] * y[(jj - 10):(jj - 1)]) + theta[3] * kappa[jj - 1]
    }
    kappa[jj] <- 1.998 / (1 + exp(-psi)) - 0.999  # a modified logistic transformation
  }
  rhohat <- kappa  # time-path of conditional copula parameter
  
  K = lgamma((nu+2)/2) + lgamma(nu/2) - 2*lgamma((nu+1)/2) 
  CL = K - 0.5*log((1-kappa^2))
  CL = CL - (nu+2)/2*log(1+(x^2 + y^2 - 2*kappa*x*y)/(nu*(1-kappa^2)));
  CL = CL + (nu+1)/2*log(1+x^2/nu) + (nu+1)/2*log(1+y^2/nu);
  CL = sum(CL);
  CL = -CL;
  
  return(list(CL,rhohat,nu))
  
}

# Create a function equal to the one before (bivt_tvp1_CL) to be used for optimization
bivt_tvp1_CL_only <- function(theta, Zdata, rhobar) { 
  CL <- bivt_tvp1_CL(theta, Zdata, rhobar)[1]
  return(CL)
}