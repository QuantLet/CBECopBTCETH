# Compute the value of the symmetrized Joe-Clayton copula pdf at a specified point
sym_jc_pdf <- function(u, v, tauU, tauL) { # same function found in sjcLL.R
  
  ### The following function is adapted in R from Andrew J. Patton's Matlab code
  ### URL: https://public.econ.duke.edu/~ap172/code.html
  
  ### Published in:
  ### Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence,
  ### International Economic Review, 47(2), 527-556.
  
  T <- max(c(length(u), length(v), length(tauU), length(tauL)))
  
  # Stretching the input vectors to match
  if (length(u) < T) {
    u <- rep(u, T)
  }
  if (length(v) < T) {
    v <- rep(v, T)
  }
  if (length(tauU) < T) {
    tauU <- rep(tauU[1], T)
  }
  if (length(tauL) < T) {
    tauL <- rep(tauL[1], T)
  }
  
  k1 <- 1 / log2(2 - tauU)
  k2 <- -1 / log2(tauL)
  CL1 <- ((1 - (1 - u)^k1)^(k2 - 1) * (1 - u)^(k1 - 1) * (-1 + k1 * (k2 * (-1 + (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(k2^(-1))) + (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(k2^(-1)))) * (1 - (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(-k2^(-1)))^(k1^(-1)) * (1 - (1 - v)^k1)^(k2 - 1) * (1 - v)^(k1 - 1))
  CL2 <- (((-1 + (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(k2^(-1)))^2) * ((1 - (1 - u)^k1)^k2 + (1 - (1 - v)^k1)^k2 - (1 - (1 - u)^k1)^k2 * (1 - (1 - v)^k1)^k2)^2)
  CL1 <- CL1 / CL2
  
  k1 <- 1 / log2(2 - tauL)
  k2 <- -1 / log2(tauU)
  u <- 1 - u
  v <- 1 - v
  CL3 <- ((1 - (1 - u)^k1)^(k2 - 1) * (1 - u)^(k1 - 1) * (-1 + k1 * (k2 * (-1 + (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(k2^(-1))) + (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(k2^(-1)))) * (1 - (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(-k2^(-1)))^(k1^(-1)) * (1 - (1 - v)^k1)^(k2 - 1) * (1 - v)^(k1 - 1))
  CL4 <- (((-1 + (-1 + (1 - (1 - u)^k1)^(-k2) + (1 - (1 - v)^k1)^(-k2))^(k2^(-1)))^2) * ((1 - (1 - u)^k1)^k2 + (1 - (1 - v)^k1)^k2 - (1 - (1 - u)^k1)^k2 * (1 - (1 - v)^k1)^k2)^2)
  CL3 <- CL3 / CL4
  CL <- 0.5 * (CL1 + CL3)
  
  return(CL)
}

# Compute the negative copula log-likelihood of the symmetrized Joe-Clayton copula 
# with time-varying tail dependence
sym_jc_tvp_CL <- function(theta, Z, thetabar) {
  
  ### The following function is adapted in R from Andrew J. Patton's Matlab code
  ### URL: https://public.econ.duke.edu/~ap172/code.html
  
  ### Published in:
  ### Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence,
  ### International Economic Review, 47(2), 527-556.
  
  w1 <- theta[1]
  a1 <- theta[2]
  b1 <- theta[3]
  w2 <- theta[4]
  a2 <- theta[5]
  b2 <- theta[6]
  
  u <- Z[, 1]
  v <- Z[, 2]
  T <- nrow(Z)
  
  TAU1 <- rep(-999.99, T)
  TAU2 <- rep(-999.99, T)
  TAU1[1] <- thetabar[1] # tau1 and tau2, based on time-invariant version of this model
  TAU2[1] <- thetabar[2]
  psi <- matrix(0, nrow = T, ncol = 2)
  
  for (jj in 2:T) {
    if (jj <= 10) {
      psi1 <- w1 + b1 * TAU1[jj - 1] + a1 * mean(abs(u[1:(jj - 1)] - v[1:(jj - 1)]))
      psi2 <- w2 + b2 * TAU2[jj - 1] + a2 * mean(abs(u[1:(jj - 1)] - v[1:(jj - 1)]))
    } else {
      psi1 <- w1 + b1 * TAU1[jj - 1] + a1 * mean(abs(u[(jj - 10):(jj - 1)] - v[(jj - 10):(jj - 1)]))
      psi2 <- w2 + b2 * TAU2[jj - 1] + a2 * mean(abs(u[(jj - 10):(jj - 1)] - v[(jj - 10):(jj - 1)]))
    }
    psi[jj, ] <- c(psi1, psi2)
    TAU1[jj] <- 0.998 / (1 + exp(-psi1)) + 0.001  # tail dependence parameters
    TAU2[jj] <- 0.998 / (1 + exp(-psi2)) + 0.001
  }
  
  CL <- sym_jc_pdf(u, v, TAU1, TAU2)
  CL <- log(CL)
  CL <- sum(CL)
  CL <- -CL
  
  if (is.nan(CL)) {
    CL <- 1e6
  }
  if (!is.double(CL)) {
    CL <- 1e7
  }
  if (!is.double(theta)) {
    CL <- 1e8
  }
  
  return(list(CL = CL, TAU1 = TAU1, TAU2 = TAU2))
}

# Create a function equal to the one before (sym_jc_tvp_CL) to be used for optimization
sym_jc_tvp_CL_only <- function(theta, Z, thetabar) { 
  CL <- sym_jc_tvp_CL(theta, Z, thetabar)[1]
  return(CL)
}
