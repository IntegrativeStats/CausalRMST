# weighted Kaplan Meier RMST estimator
#@@df: data.frame
#@@p: weights (e.g., calibration weights, inverse propensity score weights)
#@@tau: truncation time t* for RMST
W.KM.Est <- function(df,p,tau){
  EW.mu <- EW.sd <- NULL
  df$A <- as.numeric(df$A) - 1
  mu1 <- my_akm_rmst(df$Y[df$A==1], df$status[df$A==1],p[df$A==1],tau)
  mu0 <- my_akm_rmst(df$Y[df$A==0], df$status[df$A==0],p[df$A==0],tau)
  delta <- mu0$mu - mu1$mu
  EW.mu <- c(mu0$mu, mu1$mu, delta)
  delta.sd <- sqrt(mu1$V + mu0$V)
  EW.sd <- c(sqrt(mu0$V), sqrt(mu1$V), delta.sd)
  return(list(mu = EW.mu, sd = EW.sd))
}
