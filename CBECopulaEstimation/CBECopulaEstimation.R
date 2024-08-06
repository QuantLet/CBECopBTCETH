library(copula)
library(stringr)

# Load uniform marginals (pseudo-observations)
U <- as.matrix(read.table("PseudoObservations.txt", header = TRUE, sep = " ", dec = "."))

# Copula-parameter estimation

# Create table to summarize results
summarycop <- data.frame()

# Normal Copula
fitcop_normal <- fitCopula(normalCopula(dim = 2), data = U, method = "mpl", 
                           estimate.variance = TRUE)
summary(fitcop_normal)
rho_normal <- coef(fitcop_normal) # parameters
LL_normal <- logLik(fitcop_normal) # log-likelihood
k <- length(coef(fitcop_normal))  # number of parameters
T <- nrow(U) # number of observations
AIC_normal <- 2*k - 2*LL_normal # AIC
BIC_normal <- k*log(T) - 2*LL_normal # BIC
summary_normal <- c("Normal",
                    round(LL_normal,2),
                    round(AIC_normal,2),
                    round(BIC_normal,2),
                    str_flatten_comma(lapply(coef(fitcop_normal),round,2)))
summarycop <- rbind(summarycop,summary_normal)
colnames(summarycop) <- c("Copula",
                          "Log-likelihood",
                          "AIC",
                          "BIC",
                          "Parameters")

# Estimate z value and p-value for Normal copula parameter
stat_normal <- rho_normal[[1]]/sqrt(fitcop_normal@var.est[1,1]) # z value
pvalue_normal <- 2*(1-pnorm(stat_normal)) # p-value

# Student t Copula
fitcop_t <- fitCopula(tCopula(dim = 2,dispstr = "un"), data = U, method = "ml",
                      estimate.variance = TRUE)
summary(fitcop_t)
rho_t <- coef(fitcop_t)[1] # parameters
df_t <- coef(fitcop_t)[2]
LL_t <- logLik(fitcop_t) # log-likelihood
k <- length(coef(fitcop_t))  # number of parameters
T <- nrow(U) # number of observations
AIC_t <- 2*k - 2*LL_t # AIC
BIC_t <- k*log(T) - 2*LL_t # BIC
summary_t <- c("Student t",
               round(LL_t,2),
               round(AIC_t,2),
               round(BIC_t,2),
               str_flatten_comma(lapply(coef(fitcop_t),round,2)))
summarycop <- rbind(summarycop,summary_t)

# Estimate z value and p-value for t copula parameters
stat_t <- c(rho_t[[1]]/sqrt(fitcop_t@var.est[1,1]),df_t[[1]]/sqrt(fitcop_t@var.est[2,2])) # z value
pvalue_t <- c(2*(1-pnorm(stat_t[1])),2*(1-pnorm(stat_t[2]))) # p-value

# Skew t Copula
source("CBEskewt.R")
dim <- ncol(U);
start <- list(rho=numeric(dim*(dim-1)/2),delta=numeric(dim),nu=6);
system.time(stcopmle<-stcop.mle(U, start=start, # maximum likelihood estimation
                                control=list(reltol=1e-4,trace=TRUE)))
showResult(stcopmle)
rho_skewt <- showResult(stcopmle)[[1]] # parameters
delta_skewt <- showResult(stcopmle)[[2]]  
df_skewt <- showResult(stcopmle)[[3]]
LL_skewt <- showResult(stcopmle)[[4]] # log-likelihood
k <- 4  # number of parameters
T <- nrow(U) # number of observations
AIC_skewt <- 2*k - 2*LL_skewt # AIC
BIC_skewt <- k*log(T) - 2*LL_skewt # BIC
summary_skewt <- c("Skew t",
                   round(LL_skewt,2),
                   round(AIC_skewt,2),
                   round(BIC_skewt,2),
                   str_flatten_comma(lapply(c(rho_skewt,delta_skewt,df_skewt),round,2)))
summarycop <- rbind(summarycop,summary_skewt)

# Estimate z value and p-value for skew t copula parameters
parameters_skewt<-c(showResult(stcopmle)$rho,
                    showResult(stcopmle)$delta[1],
                    showResult(stcopmle)$delta[2],
                    showResult(stcopmle)$nu)
fisher_info_skewt<- solve(stcopmle$hessian)
stderror_skewt<-sqrt(diag(fisher_info_skewt))
stat_skewt <- parameters_skewt/stderror_skewt # z value
pvalue_skewt <- c(2*(1-pnorm(stat_skewt[1])), # p-value
              2*(pnorm(stat_skewt[2])),
              2*(pnorm(stat_skewt[3])),
              2*(1-pnorm(stat_skewt[4])))

# Gumbel Copula
fitcop_gumbel <- fitCopula(gumbelCopula(dim = 2), data = U, method = "mpl")
summary(fitcop_gumbel)
alpha_gumbel <- coef(fitcop_gumbel) # parameters
LL_gumbel <- logLik(fitcop_gumbel) # log-likelihood
k <- length(coef(fitcop_gumbel))  # number of parameters
T <- nrow(U) # number of observations
AIC_gumbel <- 2*k - 2*LL_gumbel # AIC
BIC_gumbel <- k*log(T) - 2*LL_gumbel # BIC
summary_gumbel <- c("Gumbel",
                    round(LL_gumbel,2),
                    round(AIC_gumbel,2),
                    round(BIC_gumbel,2),
                    str_flatten_comma(lapply(coef(fitcop_gumbel),round,2)))
summarycop <- rbind(summarycop,summary_gumbel)

# Estimate z value and p-value for Gumbel copula parameters
stat_gumbel <- alpha_gumbel[[1]]/sqrt(fitcop_gumbel@var.est[1,1]) # z value
pvalue_gumbel <- 2*(1-pnorm(stat_gumbel))

# Clayton Copula
fitcop_clayton <- fitCopula(claytonCopula(dim = 2), data = U, method = "mpl",
                            estimate.variance = TRUE)
summary(fitcop_clayton)
alpha_clayton <- coef(fitcop_clayton) # parameters
LL_clayton <- logLik(fitcop_clayton) # log-likelihood
k <- length(coef(fitcop_clayton))  # number of parameters
T <- nrow(U) # number of observations
AIC_clayton <- 2*k - 2*LL_clayton # AIC
BIC_clayton <- k*log(T) - 2*LL_clayton # BIC
summary_clayton <- c("Clayton",
                     round(LL_clayton,2),
                     round(AIC_clayton,2),
                     round(BIC_clayton,2),
                     str_flatten_comma(lapply(coef(fitcop_clayton),round,2)))
summarycop <- rbind(summarycop,summary_clayton)

# Estimate z value and p-value for Clayton copula parameters
stat_clayton <- alpha_clayton[[1]]/sqrt(fitcop_clayton@var.est[1,1]) # z value
pvalue_clayton <- 2*(1-pnorm(stat_clayton)) # p-value

# Frank Copula
fitcop_frank <- fitCopula(frankCopula(dim = 2), data = U, method = "mpl",
                          estimate.variance = TRUE)
summary(fitcop_frank)
alpha_frank <- coef(fitcop_frank) # parameters
LL_frank <- logLik(fitcop_frank) # log-likelihood
k <- length(coef(fitcop_frank))  # number of parameters
T <- nrow(U) # number of observations
AIC_frank <- 2*k - 2*LL_frank # AIC
BIC_frank <- k*log(T) - 2*LL_frank # BIC
summary_frank <- c("Frank",
                   round(LL_frank,2),
                   round(AIC_frank,2),
                   round(BIC_frank,2),
                   str_flatten_comma(lapply(coef(fitcop_frank),round,2)))
summarycop <- rbind(summarycop,summary_frank)

# Estimate z value and p-value for Frank copula parameters
stat_frank <- alpha_frank[[1]]/sqrt(fitcop_frank@var.est[1,1]) # z value
pvalue_frank <- 2*(1-pnorm(stat_frank)) # p-value

# Joe Copula
fitcop_joe <- fitCopula(joeCopula(dim = 2), data = U, method = "mpl",
                        estimate.variance = TRUE)
summary(fitcop_joe)
alpha_joe <- coef(fitcop_joe) # parameters
LL_joe <- logLik(fitcop_joe) # log-likelihood
k <- length(coef(fitcop_joe)) # number of parameters
T <- nrow(U) # number of observations
AIC_joe <- 2*k - 2*LL_joe # AIC
BIC_joe <- k*log(T) - 2*LL_joe # BIC
summary_joe <- c("Joe",
                 round(LL_joe,2),
                 round(AIC_joe,2),
                 round(BIC_joe,2),
                 str_flatten_comma(lapply(coef(fitcop_joe),round,2)))
summarycop <- rbind(summarycop,summary_joe)

# Estimate z value and p-value for Joe copula parameters
stat_joe <- alpha_joe[[1]]/sqrt(fitcop_joe@var.est[1,1]) # z value
pvalue_joe <- 2*(1-pnorm(stat_joe))

# Symmetrized Joe-Clayton Copula
source("CBEsjc.R")
theta0 <- c(0.25,0.25) # initialize parameters
result_sjc <- optim( # maximum likelihood estimation
  par = theta0,
  fn = bb7CL,
  hessian = TRUE,
  control=list(trace=TRUE, maxit=1000),
  data = U
)
theta_sjc <- result_sjc$par # parameters (tauU,tauL)
LL_sjc <- -result_sjc$value # log-likelihood
k <- length(theta_sjc)  # number of parameters
T <- nrow(U) # number of observations
AIC_sjc <- 2*k - 2*LL_sjc # AIC
BIC_sjc <- k*log(T) - 2*LL_sjc # BIC
summary_sjc <- c("Symmetrized Joe-Clayton (SJC)",
                 round(LL_sjc,2),
                 round(AIC_sjc,2),
                 round(BIC_sjc,2),
                 str_flatten_comma(lapply(theta_sjc,round,2)))
summarycop <- rbind(summarycop,summary_sjc)

# Estimate z value and p-value for Symmetrized Joe-Clayton copula  parameters
fisher_info_sjc <- solve(result_sjc$hessian)
stderror_sjc <- sqrt(diag(fisher_info_sjc))
stat_sjc <- result_sjc$par/stderror_sjc # z value
pvalue_sjc <- 2*(1-pnorm(stat_sjc)) # p-value

# Rotated Gumbel copula (180 degrees) or Gumbel survival copula
fitcop_gumbel_s <- fitCopula(rotCopula(gumbelCopula(dim = 2)), data = U, method = "mpl",
                             estimate.variance = TRUE)
summary(fitcop_gumbel_s)
alpha_gumbel_s <- coef(fitcop_gumbel_s) # parameters
LL_gumbel_s <- logLik(fitcop_gumbel_s) # log-likelihood
k <- length(coef(fitcop_gumbel_s))  # number of parameters
T <- nrow(U) # number of observations
AIC_gumbel_s <- 2*k - 2*LL_gumbel_s # AIC
BIC_gumbel_s <- k*log(T) - 2*LL_gumbel_s # BIC
summary_gumbel_s <- c("Rotated Gumbel (180 degrees)",
                      round(LL_gumbel_s,2),
                      round(AIC_gumbel_s,2),
                      round(BIC_gumbel_s,2),
                      str_flatten_comma(lapply(coef(fitcop_gumbel_s),round,2)))
summarycop <- rbind(summarycop,summary_gumbel_s)

# Estimate z value and p-value for rotated Gumbel (180 degrees) copula  parameters
stat_gumbel_s <- alpha_gumbel_s[[1]]/sqrt(fitcop_gumbel_s@var.est[1,1]) # z value
pvalue_gumbel_s <- 2*(1-pnorm(stat_gumbel_s))

# Results for time-invariant (static) copula models
summarycop

# Time-varying Normal copula
source("CBEtvnormal.R")
rho1_tvnormal <- coef(fitcop_normal) # parameter of the Normal copula without time-variation
theta0 <- c(log((1 + rho1_tvnormal) / (1 - rho1_tvnormal)), 0, 0) # initialize parameters of the evolution equation
result_tvnormal <- optim( # maximum likelihood estimation
  par = theta0,
  fn = bivnorm_tvp1_CL_only,
  hessian = TRUE,
  control=list(trace=TRUE, maxit=10000),
  Zdata = U,
  rhobar = rho1_tvnormal
)
theta_tvnormal <- result_tvnormal$par # parameters of the evolution equation
LL_tvnormal <- -result_tvnormal$value # log-likelihood
rhot_tvnormal <- bivnorm_tvp1_CL(theta_tvnormal, U, rho1_tvnormal)[[2]] # time-path of copula parameter
k <- length(theta_tvnormal)  # number of parameters
T <- nrow(U) # number of observations
AIC_tvnormal <- 2*k - 2*LL_tvnormal # AIC
BIC_tvnormal <- k*log(T) - 2*LL_tvnormal # BIC
summary_tvnormal <- c("Time-varying Normal",
                      round(LL_tvnormal,2),
                      round(AIC_tvnormal,2),
                      round(BIC_tvnormal,2),
                      str_flatten_comma(lapply(theta_tvnormal,round,2)))
summarycop <- rbind(summarycop,summary_tvnormal)

# Estimate z value and p-value for the parameters of the evolution equation corresponding to the time-varying Normal copula
fisher_info_tvnormal <- solve(result_tvnormal$hessian)
stderror_tvnormal <- sqrt(diag(fisher_info_tvnormal))
stat_tvnormal <- result_tvnormal$par/stderror_tvnormal # z value
pvalue_tvnormal <- c(2*pnorm(stat_tvnormal[1]), # p-value
                     2*(1-pnorm(stat_tvnormal[2])),
                     2*(1-pnorm(stat_tvnormal[3])))

# Time-varying Student t Copula
source("CBEtvstudentt.R")
rho1_tvt <- coef(fitcop_t)[1] # parameters of the Student t copula without time-variation
nubar <- coef(fitcop_t)[2]
theta0 <- c(log((1 + rho1_tvt) / (1 - rho1_tvt)), 0, 0,nubar) # initialize parameters of the evolution equation and degrees of freedom
result_tvt <- optim( # maximum likelihood estimation
  par = theta0,
  fn = bivt_tvp1_CL_only,
  hessian = TRUE,
  control=list(trace=TRUE, maxit=1000),
  Zdata = U,
  rhobar = rho1_tvt
)
theta_tvt <- result_tvt$par # parameters of the evolution equation and degrees of freedom
LL_tvt <- -result_tvt$value # log-likelihood
rhot_tvt <- bivt_tvp1_CL(theta_tvt, U, rho1_tvt)[[2]] # time-path of copula parameter
k <- length(theta_tvt) # number of parameters
T <- nrow(U) # number of observations
AIC_tvt <- 2*k - 2*LL_tvt # AIC
BIC_tvt <- k*log(T) - 2*LL_tvt # BIC
summary_tvt <- c("Time-varying Student t",
                 round(LL_tvt,2),
                 round(AIC_tvt,2),
                 round(BIC_tvt,2),
                 str_flatten_comma(lapply(theta_tvt,round,2)))
summarycop <- rbind(summarycop,summary_tvt)

# Estimate z value and p-value for the parameters of the evolution equation corresponding to the time-varying Student t copula
fisher_info_tvt <- solve(result_tvt$hessian)
stderror_tvt <- sqrt(diag(fisher_info_tvt))
stat_tvt <- result_tvt$par/stderror_tvt # z value
pvalue_tvt <- c(2*pnorm(stat_tvt[1]), # p-value
            2*(1-pnorm(stat_tvt[2])),
            2*(1-pnorm(stat_tvt[3])),
            2*(1-pnorm(stat_tvt[4])))

# Time-varying rotated Gumbel (180 degrees) copula 
source("CBEtvgumbel.R")
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
theta_tvgumbel_s <- result_tvgumbel_s$par # parameters of the evolution equation
LL_tvgumbel_s <- -result_tvgumbel_s$value # log-likelihood
alphat_tvgumbel_s <- gumbel_tvp1_CL(theta_tvgumbel_s, 1-U, alpha1_tvgumbel_s)[[2]] # time-path of copula parameter
k <- length(theta_tvgumbel_s)  # number of parameters
T <- nrow(U) # number of observations
AIC_tvgumbel_s <- 2*k - 2*LL_tvgumbel_s # AIC
BIC_tvgumbel_s <- k*log(T) - 2*LL_tvgumbel_s # BIC
summary_tvgumbel_s <- c("Time-varying rotated Gumbel (180 degrees)",
                        round(LL_tvgumbel_s,2),
                        round(AIC_tvgumbel_s,2),
                        round(BIC_tvgumbel_s,2),
                        str_flatten_comma(lapply(theta_tvgumbel_s,round,2)))
summarycop <- rbind(summarycop,summary_tvgumbel_s)

# Estimate z value and p-value for the parameters of the evolution equation corresponding to time-varying rotated Gumbel copula
fisher_info_tvgumbel_s <- solve(result_tvgumbel_s$hessian)
stderror_tvgumbel_s <-sqrt(diag(fisher_info_tvgumbel_s))
stat_tvgumbel_s <- result_tvgumbel_s$par/stderror_tvgumbel_s # z value
pvalue_tvgumbel_s <- c(2*(1-pnorm(stat_tvgumbel_s[1])), # p-value
                       2*pnorm(stat_tvgumbel_s[2]),
                       2*pnorm(stat_tvgumbel_s[3]))

# Time-varying Symmetrized Joe-Clayton (SJC) copula
source("CBEtvsjc.R")
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
theta_tvsjc <- result_tvsjc$par # parameters of both evolution equations
LL_tvsjc <- -result_tvsjc$value # log-likelihood
tauU_tvsjc <- sym_jc_tvp_CL(theta_tvsjc, U, kappa_tvsjc)[[2]] # time-path of copula parameter tauU
tauL_tvsjc <- sym_jc_tvp_CL(theta_tvsjc, U, kappa_tvsjc)[[3]] # time-path of copula parameter tauL
k <- length(theta_tvsjc)  # number of parameters
T <- nrow(U) # number of observations
AIC_tvsjc <- 2*k - 2*LL_tvsjc # AIC
BIC_tvsjc <- k*log(T) - 2*LL_tvsjc # BIC
summary_tvsjc <- c("Time-varying Symmetrized Joe-Clayton (SJC)",
                   round(LL_tvsjc,2),
                   round(AIC_tvsjc,2),
                   round(BIC_tvsjc,2),
                   str_flatten_comma(lapply(theta_tvsjc,round,2)))
summarycop <- rbind(summarycop,summary_tvsjc)

# Estimate z value and p-value for the parameters of the evolution equation corresponding to time-varying SJC copula
fisher_info_tvsjc <- solve(result_tvsjc$hessian)
stderror_tvsjc <- sqrt(diag(fisher_info_tvsjc))
stat_tvsjc <- result_tvsjc$par/stderror_tvsjc # z value
pvalue_tvsjc <- c(2*(1-pnorm(stat_tvsjc[1])),
                  2*(pnorm(stat_tvsjc[2])),
                  2*(pnorm(stat_tvsjc[3])),
                  2*(1-pnorm(stat_tvsjc[4])),
                  2*(pnorm(stat_tvsjc[5])),
                  2*(pnorm(stat_tvsjc[6])))

# Results for time-varying copula models
summarycop[10:13,]

# Results for static and time-varying copula models ordered by AIC values
summarycop$AIC <- as.numeric(as.character(summarycop$AIC)) # convert column to numeric
summarycop$BIC <- as.numeric(as.character(summarycop$BIC))
summarycop$`Log-likelihood` <- as.numeric(as.character(summarycop$`Log-likelihood`))
summarycop <- summarycop[order(summarycop$AIC, decreasing = FALSE),] # put in ascending order according to AIC
rownames(summarycop) <- 1:nrow(summarycop) # reset index
View(summarycop)
