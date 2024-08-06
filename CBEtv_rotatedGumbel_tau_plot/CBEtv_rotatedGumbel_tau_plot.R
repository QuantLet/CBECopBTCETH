library(quantmod)
library(qrmtools)
library(rugarch)
library(copula)
library(ggplot2)

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

# Time-invariant rotated Gumbel copula (180 degrees) or Gumbel survival copula

fitcop_gumbel_s <- fitCopula(rotCopula(gumbelCopula(dim = 2)), data = U, method = "mpl")

# Time-varying rotated Gumbel copula (180 degrees) or Gumbel survival copula

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

# Maximize the copula-likelihood function
alpha1_tvgumbel_s <- coef(fitcop_gumbel_s) # parameter of the rotated Gumbel copula without time-variation
theta0 <- c(sqrt(alpha1_tvgumbel_s-1), 0, 0) # initialize the parameters of the evolution equation
result_tvgumbel_s <- optim( # maximum likelihood estimation
  par = theta0,
  fn = gumbel_tvp1_CL_only,
  hessian = TRUE,
  control=list(trace=TRUE, maxit=1000),
  data = 1-U,  # rotated copula
  kappabar = alpha1_tvgumbel_s
)

# Extract the parameters of the evolution equation that maximize the copula-likelihood function  
theta_tvgumbel_s <- result_tvgumbel_s$par

# Obtain the time-varying copula parameter
alphat_tvgumbel_s <- gumbel_tvp1_CL(theta_tvgumbel_s, 1-U, alpha1_tvgumbel_s)[[2]]

# Plot time-path of Kendall's tau derived from the copula parameter
kendallt_tvgumbel_s <- 1-1/alphat_tvgumbel_s
kendall1_tvgumbel_s <- 1-1/alpha1_tvgumbel_s
T_kendallt_tvgumbel_s <- length(kendallt_tvgumbel_s)
data_kendall_tvgumbel_s <- data.frame(time = df_residuals_BTCETH[,1], 
                                      kendallt_tvgumbel_s = kendallt_tvgumbel_s, 
                                      kendall1_tvgumbel_s = kendall1_tvgumbel_s * rep(1,T_kendallt_tvgumbel_s))
ggplot(data_kendall_tvgumbel_s, aes(x = time)) +
  geom_line(aes(y = kendallt_tvgumbel_s, color = "Time-varying")) +
  geom_line(aes(y = kendall1_tvgumbel_s, color = "Constant"), linetype = "dashed") +
  ggtitle("Time-path of Kendall's tau derived from the rotated Gumbel copula parameter") +
  labs(x = "time", y = expression(tau), color = "") +
  theme_classic() +
  scale_color_manual(values = c("red","black")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.background = element_blank(),
        legend.position = c(0.9,0.1))
