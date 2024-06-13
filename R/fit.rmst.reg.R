# Fit the IPCW RMST regression with g(X) and its interaction with treatment
#@@df: data.frame
#@@varX: names of calibrated covariates
#@@gf: vector function of calibrated covariates (e.g, moment function)
#@@tau: truncation time t* for RMST
fit.rmst.reg <- function(df, varX, gf, tau){
  dfx <- t(apply(df[,varX], 1, gf))
  cov <- as.data.frame(cbind(df$A, dfx, df$A*dfx))
  rmst_fit <- my.rmst2reg(y = df$Y,
                          delta = df$status,
                          x = cov,
                          arm = df$A,
                          tau = tau)
  return(rmst_fit)
}
