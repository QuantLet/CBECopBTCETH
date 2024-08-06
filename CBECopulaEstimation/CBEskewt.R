### The following code for Skew t copula estimation is adapted from
### the paper "Maximum likelihood estimation of skew t-copula" by Yoshiba, T. (2014).
### URL: https://www.ism.ac.jp/editsec/resmemo/resmemo-file/resm1183.pdf

# skew-t coplua estimation (MLE) using "sn" ver.1.0-0
library(sn)
library(signal)
## redefine qst on "sn" ver.1.0-0 here
## negative log-likelihood for multivariate skew-t copula
## udat[1:n,1:dim] : pseudo sample (N observations for [0,1]^dim)
stcopnll <- function(para, udat=NULL){
  mpoints <- 150;
  dp <- stOrgPara(para);
  delta <- dp$delta;
  zeta <- delta/sqrt(1-delta*delta);
  dim <- length(delta);
  Omega <- diag(dim);
  Omega[upper.tri(Omega)] <- Omega[lower.tri(Omega)] <- dp$rho;
  iOmega <- solve(Omega);
  alpha <- iOmega %*% delta /sqrt(1-(t(delta) %*% iOmega %*% delta)[1,1]);
  nu <- dp$nu;
  ix <- ipqst(udat,zeta,nu,mpoints,rel.tol=1e-6);
  ## Activate the following line instead of monotone interpolating quantile
  ## function ipqst() to use accurate quantile function aqst()
  ## ix <- aqst(udat,zeta,nu,mpoints);
  lm <- matrix(0,nrow=nrow(udat),ncol=dim);
  for(j in 1:dim){ lm[,j] <- dst(ix[,j], alpha=zeta[j], nu=nu, log=TRUE); }
  lc <- dmst(ix,Omega=Omega,alpha=alpha,nu=nu,log=TRUE);
  -sum(lc)+sum(lm)
}
stcop.mle <- function (udat, start = NULL, gr = NULL, ...,
                       method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                       lower = -Inf, upper = Inf,
                       control = list(), hessian = TRUE)
{
  iniPar <- stIntPara(start$rho,start$delta,start$nu);
  method <- match.arg(method);
  fit <- optim(iniPar, stcopnll, method=method, control=control, udat=udat, hessian=hessian);
  list(call = match.call(), dp = stOrgPara(fit$par), logL = -fit$value,
       details=fit, nobs = nrow(udat), method = method, hessianX= fit$hessian);
}
## show estimated parameters and the log-likelihood ##
showResult <- function(fit){
  dp <- fit$dp;
  list(rho=dp$rho,delta=dp$delta,nu=dp$nu,logL=fit$logL,hessian=fit$hessianX);
}




## random number generator of skew t-copula
rstcop <- function(n,rho,delta,nu,...){
  dim <- length(delta);
  zeta <- delta/sqrt(1-delta*delta);
  Omega <- diag(dim);
  Omega[upper.tri(Omega)] <- Omega[lower.tri(Omega)] <- rho;
  iOmega <- solve(Omega);
  alpha <- iOmega %*% delta /sqrt(1-(t(delta) %*% iOmega %*% delta)[1,1]);
  x <- rmst(n=n,Omega=Omega,alpha=alpha,nu=nu);
  u <- matrix(0,nrow=n,ncol=dim);
  for(j in 1:dim){ u[,j] <- pst(x[,j], alpha=zeta[j], nu=nu,...); }
  list(x=x,u=u);
}




## transforming original parameters to internal parameters ##
stIntPara <- function(rho,delta,nu){
  ndim <- length(delta)+1;
  R <- diag(ndim);
  for(i in 2:ndim){
    R[i,1] <- R[1,i] <- delta[i-1];
    if(i>=3){ for(j in 2:(i-1)){ R[i,j] <- R[j,i] <-
      rho[i-ndim+(j-1)*(ndim-2-(j-2)/2)]; } }
  }
  LTR <- t(chol(R));
  Mtheta <- matrix(0,nrow=ndim,ncol=ndim);
  for(i in 2:ndim){
    Mtheta[i,1] <- acos(LTR[i,1]);
    cumsin <- sin(Mtheta[i,1]);
    if(i>=3){ for(j in 2:(i-1)){
      Mtheta[i,j] <- acos(LTR[i,j]/cumsin);
      cumsin <- cumsin*sin(Mtheta[i,j]); }
    }
  }
  c(Mtheta[lower.tri(Mtheta)],log(nu-2.0));
}
## transforming internal parameters to original parameters ##
stOrgPara <- function(para){
  ntheta <- length(para)-1;
  theta <- para[1:ntheta];
  ndim <- (1+sqrt(1+8*ntheta))/2;
  LTR <- diag(ndim);
  for(i in 2:ndim){
    LTR[i,1] <- cos(theta[i-1]);
    cumsin <- sin(theta[i-1]);
    if(i>=3){ for(j in 2:(i-1)){
      k <- i+ndim*(j-1)-j*(j+1)/2;
      LTR[i,j] <- cumsin*cos(theta[k]);
      cumsin <- cumsin*sin(theta[k]); }
    }
    LTR[i,i] <- cumsin;
  }
  R <- LTR %*% t(LTR);
  Omega <- R[-1,-1];
  delta <- R[1,-1];
  nu <- exp(para[ntheta+1])+2.0;
  list(rho = Omega[lower.tri(Omega)], delta = delta, nu = nu);
}



## redefine qst on "sn" ver.1.0-0
qst <- function (p, xi = 0, omega = 1, alpha = 0, nu = Inf, tol = 1e-08, maxit
                 = 30, ...)
{
  if (length(alpha) > 1)
    stop("'alpha' must be a single value")
  if (length(nu) > 1)
    stop("'nu' must be a single value")
  if (nu <= 0)
    stop("nu must be non-negative")
  if (nu == Inf)
    return(qsn(p, xi, omega, alpha))
  if (nu == 1)
    return(qsc(p, xi, omega, alpha))
  if (alpha == Inf)
    return(xi + omega * sqrt(qf(p, 1, nu)))
  if (alpha == -Inf)
    return(xi - omega * sqrt(qf(1 - p, 1, nu)))
  na <- is.na(p) | (p < 0) | (p > 1)
  abs.alpha <- abs(alpha)
  if (alpha < 0)
    p <- (1 - p)
  zero <- (p == 0)
  one <- (p == 1)
  x <- xa <- xb <- xc <- fa <- fb <- fc <- rep(NA, length(p))
  nc <- rep(TRUE, length(p))
  nc[(na | zero | one)] <- FALSE
  fc[!nc] <- 0
  xa[nc] <- qt(p[nc], nu)
  xb[nc] <- sqrt(qf(p[nc], 1, nu))
  fa[nc] <- pst(xa[nc], 0, 1, abs.alpha, nu, ...) - p[nc]
  fb[nc] <- pst(xb[nc], 0, 1, abs.alpha, nu, ...) - p[nc]
  regula.falsi <- FALSE
  it <- 0
  while (sum(nc) > 0 & it < maxit) {
    xc[nc] <- if (regula.falsi)
      xb[nc] - fb[nc] * (xb[nc] - xa[nc])/(fb[nc] - fa[nc])
    else (xb[nc] + xa[nc])/2
    fc[nc] <- pst(xc[nc], 0, 1, abs.alpha, nu, ...) - p[nc]
    pos <- (fc[nc] > 0)
    xa[nc][!pos] <- xc[nc][!pos]
    fa[nc][!pos] <- fc[nc][!pos] 
    xb[nc][pos] <- xc[nc][pos]
    fb[nc][pos] <- fc[nc][pos]
    x[nc] <- xc[nc]
    nc[(abs(fc) < tol)] <- FALSE
    regula.falsi <- !regula.falsi
    it <- it + 1
  }
  x <- replace(x, zero, -Inf)
  x <- replace(x, one, Inf)
  Sign <- function(x) sign(x)+ as.numeric(x==0)
  q <- as.numeric(xi + omega * Sign(alpha)* x)
  names(q) <- names(p)
  return(q)
}



## accurate quantiles ##
aqst <- function(udat,zeta,nu, ...){
  dim <- ncol(udat);
  ax <- matrix(0,nrow=nrow(udat),ncol=dim);
  for(j in 1:dim){
    ax[,j] <- qst(udat[,j], alpha=zeta[j], nu=nu, ...);
  }
  ax
}
## empirical quantiles with random sampling ##
rsqst <- function(udat,zeta,nu,simNum){
  dim <- ncol(udat);
  sx <- matrix(0,nrow= nrow(udat),ncol=dim);
  sy <- matrix(0,nrow=simNum,ncol=dim);
  for(j in 1:dim){
    sy[,j] <- sort(rst(simNum, alpha=zeta[j], nu=nu));
    sx[,j] <- sy[udat[,j]*(simNum-1)+1,j];
  }
  sx
}
## interpolating quantiles ##
ipqst <- function(udat,zeta,nu,mpoints, ...){
  dim <- ncol(udat);
  ix <- matrix(0,nrow=nrow(udat),ncol=dim);
  for(j in 1:dim){
    minx <- qst(min(udat[,j]), alpha=zeta[j], nu=nu, ...);
    maxx <- qst(max(udat[,j]), alpha=zeta[j], nu=nu, ...); 
    xx <- seq(minx,maxx,length.out=mpoints);
    px <- sort(pst(xx, alpha=zeta[j], nu=nu, ...));
    ix[,j] <- pchip(px, xx, udat[,j]);
  }
  ix
}

