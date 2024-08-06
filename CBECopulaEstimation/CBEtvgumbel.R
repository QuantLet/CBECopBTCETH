# Compute the negative copula log-likelihood of Gumbel copula with time-varying parameter 
gumbel_tvp1_CL <- function(theta, data, kappabar) {
  
  ### The following function is adapted in R from Andrew J. Patton's Matlab code
  ### URL: https://public.econ.duke.edu/~ap172/code.html
  
  ### Written for the following papers:
  ### Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence, International Economic Review, 47(2), 527-556. 
  ### Patton, A.J., 2006, Estimation of Multivariate Models for Time Series of Possibly Different Lengths, Journal of Applied Econometrics, 21(2), 147-173.  
  ### Patton, A.J., 2004, On the Out-of-Sample Importance of Skewness and Asymmetric Dependence for Asset Allocation, Journal of Financial Econometrics, 2(1), 130-168. 
  
  T <- nrow(data)
  u <- data[, 1]
  v <- data[, 2]
  
  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]

  kappa <- rep(-999.99, T)
  kappa[1] <- kappabar
  for (jj in 2:T) {
    if (jj <= 10) {
      psi1 <- theta1 + theta2 * kappa[jj - 1] + theta3 * mean(abs(u[1:(jj - 1)] - v[1:(jj - 1)]))
    } else {
      psi1 <- theta1 + theta2 * kappa[jj - 1] + theta3 * mean(abs(u[(jj - 10):(jj - 1)] - v[(jj - 10):(jj - 1)]))
    }
    kappa[jj] <- 1.0001 + psi1^2
  }    
      
  ut <- -log(data[, 1])
  vt <- -log(data[, 2])
  
  CL <- -(ut^kappa + vt^kappa)^(1 / kappa) - log(u) - log(v)
  CL <- CL + (kappa - 1) * (log(ut) + log(vt)) - (2 - 1 / kappa) * (log(ut^kappa + vt^kappa))
  CL <- CL + log((ut^kappa + vt^kappa)^(1 / kappa) + kappa - 1)
  CL <- sum(CL)
  CL <- -CL
  
  if (!is.double(theta)) {
    CL <- 1e6
  } else if (!is.double(CL)) {
    CL <- 1e7
  } else if (is.nan(CL)) {
    CL <- 1e8
  } else if (is.infinite(CL)) {
    CL <- 1e9
  }
  
  return(list(CL = CL, kappa = kappa))
}

# Create a function equal to the one before (gumbel_tvp1_CL) to be used for optimization
gumbel_tvp1_CL_only <- function(theta, data, kappabar) { 
  CL <- gumbel_tvp1_CL(theta, data, kappabar)[1]
  return(CL)
}

