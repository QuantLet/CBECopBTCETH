library(quantmod)
library(qrmtools)

# Obtain daily prices of Bitcoin (BTC) and Ethereum (ETH) from Yahoo Finance 
df_BTC <- getSymbols('BTC-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
df_ETH <- getSymbols('ETH-EUR',src='yahoo',auto.assign=FALSE, from = "2017-11-11",
                     to = "2024-05-10")
BTC_ETH <- cbind(df_BTC$`BTC-EUR.Adjusted`,df_ETH$`ETH-EUR.Adjusted`)
colnames(BTC_ETH) <- c("BTC-EUR","ETH-EUR") 

# Calculate daily (log-)returns
R_BTCETH <- returns(BTC_ETH)

# Plot daily return series
plot.zoo(R_BTCETH, main = "Returns of cryptocurrencies data", xlab = "", col = c("blue3","red3"))

# Plot histogram of returns
hist(R_BTCETH$`BTC-EUR`, main = "Distribution of BTC returns", xlab = "", col = "blue3") # BTC
hist(R_BTCETH$`ETH-EUR`, main = "Distribution of ETH returns", xlab = "", col = "red3") # ETH
