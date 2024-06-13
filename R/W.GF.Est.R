# weighted G-Formula RMST estimator
#@@df: data.frame
#@@p: weights (e.g., calibration weights, inverse propensity score weights)
#@@varX: names of calibrated covariates
#@@gf: vector function of calibrated covariates (e.g, moment function)
#@@tau: truncation time t* for RMST
W.GF.Est <- function(df, p, varX, gf, tau){
  EW.mu <- EW.sd <- NULL
  df$A <- as.numeric(df$A) - 1
  fit <- fit.rmst.reg(df, varX, gf, tau)
  gamma <- coef(fit)
  vgamma <- vcov(fit)
  dfx <- t(apply(df[,varX], 1, gf))
  mdf0 <- cbind(1,0,dfx, 0*dfx)
  mdf1 <- cbind(1,1,dfx, 1*dfx)
  mu0 <- sum(p*(as.matrix(mdf0)%*%matrix(gamma,ncol = 1)))/sum(p)
  mu1 <- sum(p*(as.matrix(mdf1)%*%matrix(gamma,ncol = 1)))/sum(p)
  delta <- mu0 - mu1
  J0 <- t(as.matrix(mdf0)) %*% p/sum(p)
  J1 <- t(as.matrix(mdf1)) %*% p/sum(p)
  delta.J <- t(as.matrix(mdf0)-as.matrix(mdf1)) %*% p/sum(p)
  sd0 <- as.numeric(sqrt(t(J0) %*% vgamma %*% J0))
  sd1 <- as.numeric(sqrt(t(J1) %*% vgamma %*% J1))
  delta.sd <- as.numeric(sqrt(t(delta.J) %*% vgamma %*% delta.J))
  EW.mu <- c(mu0,mu1,delta)
  EW.sd <- c(sd0,sd1,delta.sd)

  return(list(mu = EW.mu, sd = EW.sd))
}
