library(quantmod)
library(qrmtools)
library(moments)
library(tseries)
library(fDMA)

# Obtain daily prices of Bitcoin (BTC) and Ethereum (ETH) from Yahoo Finance 
df_BTC <- getSymbols('BTC-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
df_ETH <- getSymbols('ETH-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
BTC_ETH <- cbind(df_BTC$`BTC-EUR.Adjusted`,df_ETH$`ETH-EUR.Adjusted`)
colnames(BTC_ETH) <- c("BTC-EUR","ETH-EUR") 

# Calculate daily (log-)returns
R_BTCETH <- returns(BTC_ETH)

# Summary Statistics

# Calculate Min, Median, Mean, Max value
summary(R_BTCETH)[,c(2,3)]

# Calculate Standard deviation (Std. dev.)
apply(R_BTCETH, 2, sd)

# Calculate skewness
skewness(R_BTCETH$`BTC-EUR`) # BTC 
skewness(R_BTCETH$`ETH-EUR`) # ETH 

# Calculate kurtosis
kurtosis(R_BTCETH$`BTC-EUR`) # BTC
kurtosis(R_BTCETH$`ETH-EUR`) # ETH 

# Jarque–Bera test statistic of normality
jarque.bera.test(R_BTCETH$`BTC-EUR`) # BTC
jarque.bera.test(R_BTCETH$`ETH-EUR`) # ETH

# Ljung-Box test statistic for the presence of no autocorrelation in returns
Box.test(R_BTCETH$`BTC-EUR`, lag = 20, type = "Ljung-Box") # BTC
Box.test(R_BTCETH$`ETH-EUR`, lag = 20, type = "Ljung-Box") # ETH

# Ljung-Box test statistic for the presence of no autocorrelation in squared returns
Box.test(R_BTCETH$`BTC-EUR`^2, lag = 20, type = "Ljung-Box") # BTC
Box.test(R_BTCETH$`ETH-EUR`^2, lag = 20, type = "Ljung-Box") # ETH

# ARCH-LM test statistic for the presence of heteroscedasticity 
archtest(as.vector(R_BTCETH$`BTC-EUR`),lag=20) # BTC
archtest(as.vector(R_BTCETH$`ETH-EUR`),lag=20) # ETH

# Augmented Dickey-Fuller test for stationarity
tseries::adf.test(R_BTCETH$`BTC-EUR`) # BTC 
tseries::adf.test(R_BTCETH$`ETH-EUR`) # ETH
