library(quantmod)
library(qrmtools)
library(rugarch)
library(copula)
library(ggplot2)
library(gridExtra)

# Obtain daily prices of Bitcoin (BTC) and Ethereum (ETH) from Yahoo Finance 
df_BTC <- getSymbols('BTC-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
df_ETH <- getSymbols('ETH-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
BTC_ETH <- cbind(df_BTC$`BTC-EUR.Adjusted`,df_ETH$`ETH-EUR.Adjusted`)
colnames(BTC_ETH) <- c("BTC-EUR","ETH-EUR") 

# Calculate daily (log-)returns
R_BTCETH <- returns(BTC_ETH)

# Marginal specification

# For BTC:
# GARCH Model	: EGARCH(3,4) 
# Mean Model	: ARFIMA(0,0,0) with mean
# Distribution	: Skew t distribution for standardized residuals

ctrl = list(RHO = 1,DELTA = 1e-8,MAJIT = 100,MINIT = 650,TOL = 1e-6)
spec_BTC = ugarchspec(variance.model = list(model="eGARCH",garchOrder = c(3, 4)),mean.model = list(armaOrder = c(0,0),include.mean=TRUE),distribution.model = "sstd") 
garch.fit_BTC = ugarchfit(data = R_BTCETH$`BTC-EUR`, spec = spec_BTC, solver = "solnp", solver.control = ctrl)
garch.fit_BTC

# For ETH:
# GARCH Model	: EGARCH(4,3) 
# Mean Model	: ARFIMA(0,0,0) with zero-mean
# Distribution	: Skew t distribution for standardized residuals

spec_ETH = ugarchspec(variance.model = list(model="eGARCH",garchOrder = c(4, 3)),mean.model = list(armaOrder = c(0,0),include.mean=FALSE),distribution.model = "sstd") 
garch.fit_ETH = ugarchfit(data = R_BTCETH$`ETH-EUR`, spec = spec_ETH, solver = "solnp", solver.control = ctrl)
garch.fit_ETH

# Put together standardized residuals from ARMA-EGARCH models of BTC and ETH
residuals_BTC <- residuals(garch.fit_BTC,standardize=TRUE)
residuals_ETH <- residuals(garch.fit_ETH,standardize=TRUE)
residuals_BTCETH <- cbind(residuals_BTC,residuals_ETH)
colnames(residuals_BTCETH) <- c("residuals_BTC","residuals_ETH")
df_residuals_BTCETH <- fortify.zoo(residuals_BTCETH) # convert to data frame for plotting

# Transform standardized residuals into pseudo-observations
U <- pobs(df_residuals_BTCETH[,-1]) # eliminate first column (Date)
colnames(U)<- c("U_BTC","U_ETH")

# Time-invariant SJC Copula

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
# symmetrised Joe-Clayton family. From Joe(1997), p. 153. 
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

# Maximize the copula-likelihood function
theta0 <- c(0.25,0.25) # initialize parameters
result_sjc <- optim( # maximum likelihood estimation
  par = theta0,
  fn = bb7CL,
  hessian = TRUE,
  control=list(trace=TRUE, maxit=1000),
  data = U
)

# Obtain the time-invariant tail dependence coefficients
theta_sjc <- result_sjc$par # (tauU,tauL)

# Time-varying SJC Copula

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

# Maximize the copula-likelihood function
kappa_tvsjc <- theta_sjc # parameters of the SJC copula without time-variation
theta0 <- c(log(kappa_tvsjc[1]/(1-kappa_tvsjc[1])), 0, 0, # initialize the parameters of both evolution equations
            log(kappa_tvsjc[2]/(1-kappa_tvsjc[2])), 0, 0)
result_tvsjc <- optim( # maximum likelihood estimation
  par = theta0,
  fn = sym_jc_tvp_CL_only, 
  hessian = TRUE,
  control=list(trace=TRUE, maxit=1000),
  Z = U,
  thetabar = kappa_tvsjc
)

# Extract the parameters of both evolution equations that maximize the copula-likelihood function  
theta_tvsjc <- result_tvsjc$par

# Obtain the time-varying tail dependence coefficients (tauU and tauL)
tauU_tvsjc <- sym_jc_tvp_CL(theta_tvsjc, U, kappa_tvsjc)[[2]]
tauL_tvsjc <- sym_jc_tvp_CL(theta_tvsjc, U, kappa_tvsjc)[[3]]

# Plot time-path of the tail dependence coefficients
T_tauLU_tvsjc <- length(tauL_tvsjc)
data_tvsjc <- data.frame(
  time = df_residuals_BTCETH[,1],
  tauL_tvsjc = tauL_tvsjc,
  tauU_tvsjc = tauU_tvsjc,
  constant_lower = theta_sjc[2] * rep(1, T_tauLU_tvsjc),
  constant_upper = theta_sjc[1] * rep(1, T_tauLU_tvsjc)
)
p1 <- ggplot(data_tvsjc, aes(x = time)) +
  geom_line(aes(y = tauL_tvsjc, color = "Time-varying")) +
  geom_line(aes(y = constant_lower, color = "Constant"), linetype = "dashed") +
  labs(title = "SJC copula - lower tail", 
       x = "time", y = expression(lambda^{"L"}), color = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  scale_color_manual(values = c("red","black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.background = element_blank(),
        legend.position = c(0.9, 0.95))
p2 <- ggplot(data_tvsjc, aes(x = time)) +
  geom_line(aes(y = tauU_tvsjc, color = "Time-varying")) +
  geom_line(aes(y = constant_upper, color = "Constant"), linetype = "dashed") +
  labs(title = "SJC copula - upper tail", 
       x = "time", y = expression(lambda^{"U"}), color = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  scale_color_manual(values = c("red","black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.background = element_blank(),
        legend.position = c(0.9, 0.95))
grid.arrange(p1, p2, nrow = 2)
