library(quantmod)
library(qrmtools)
library(rugarch)
library(copula)
library(ggpubr)
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

# Fit a rotated Gumbel copula (180 degrees) or Gumbel survival copula
fitcop_gumbel_s <- fitCopula(rotCopula(gumbelCopula(dim = 2)), data = U, method = "mpl")
alpha_gumbel_s <- coef(fitcop_gumbel_s)

# Wireframe and contour plot of the copula density for theta = 2.69
p1 <- wireframe2(rotCopula(gumbelCopula(param = alpha_gumbel_s, dim = 2)),
                  FUN = dCopula, 
                  delta = 0.025,
                  xlab = expression(u["1"]),
                  ylab = expression(u["2"]), 
                  zlab = "")
p2 <- contourplot2(rotCopula(gumbelCopula(param = alpha_gumbel_s, dim = 2)), 
                  FUN = dCopula, 
                  n.grid = 42, 
                  cuts = 33, 
                  lwd = 1/2,
                  xlab = expression(u["1"]),
                  ylab = expression(u["2"]))
title1 = text_grob(bquote(bold("Rotated Gumbel copula density" ~ theta == .(round(alpha_gumbel_s,2)))))
grid.arrange(p1, p2, nrow = 1, top = title1)
