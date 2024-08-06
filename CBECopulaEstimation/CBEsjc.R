# Compute the value of the symmetrized Joe-Clayton copula pdf at a specified point
sym_jc_pdf <- function(u, v, tauU, tauL) {
  
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

# Calculate the negative copula log-likelihood of a member of the 
# symmetrized Joe-Clayton family. From Joe(1997), p. 153. 
bb7CL <- function(theta, data) {
  
  ### The following function is adapted in R from Andrew J. Patton's Matlab code
  ### URL: https://public.econ.duke.edu/~ap172/code.html

  ### Published in:
  ### Patton, A.J., 2006, Modelling Asymmetric Exchange Rate Dependence,
  ### International Economic Review, 47(2), 527-556.
  
  CL <- sym_jc_pdf(data[,1], data[,2], theta[1], theta[2])
  CL <- log(CL)
  
  CL <- sum(CL)
  CL <- -CL
  
  if (is.nan(CL)) {
    CL <- 1e6
  }
  if (is.complex(CL)) {
    CL <- 1e7
  }
  if (is.complex(theta)) {
    CL <- 1e8
  }
  
  return(CL)
}

