library(stats)

# Compute the negative copula log-likelihood of a bivariate Normal distribution
# with time-varying correlation
bivnorm_tvp1_CL <- function(theta, Zdata, rhobar) {
  
  ### The following function is adapted in R from Andrew J. Patton's Matlab code
  ### URL: https://public.econ.duke.edu/~ap172/code.html
  
  ### Written for the following papers:
  ### Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
  ### Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
  ### Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
  
  T <- nrow(Zdata)
  x <- qnorm(Zdata[, 1], 0, 1)
  y <- qnorm(Zdata[, 2], 0, 1)
  
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
  rhohat <- kappa  # time-path of copula parameter
  
  CL <- -1 * (2 * (1 - kappa^2))^(-1) * (x^2 + y^2 - 2 * kappa * x * y)
  CL <- CL + 0.5 * (x^2 + y^2)
  CL <- CL - 0.5 * log(1 - kappa^2)
  CL <- sum(CL)
  CL <- -CL
  
  return(list(CL,rhohat))
}

# Create a function equal to the one before (bivnorm_tvp1_CL) to be used for optimization
bivnorm_tvp1_CL_only <- function(theta, Zdata, rhobar) { 
  CL <- bivnorm_tvp1_CL(theta, Zdata, rhobar)[1]
  return(CL)
}
  