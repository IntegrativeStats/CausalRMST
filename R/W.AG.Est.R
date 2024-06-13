# weighted Augmented RMST estimator
#@@df: data.frame
#@@p: weights (e.g., calibration weights, inverse propensity score weights)
#@@varX: names of calibrated covariates
#@@gf: vector function of calibrated covariates (e.g, moment function)
#@@tau: truncation time t* for RMST
W.AG.Est <- function(df, p, varX, gf, tau){
  require(geex)
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    mu1 <- data$mu1
    mu0 <- data$mu0
    function(theta){
      c(p*A*w*(Y-mu1-theta[1]),
        p*(1-A)*w*(Y-mu0-theta[2]),
        p*(mu1-theta[3]),
        p*(mu0-theta[4]))
    }
  }
  df$A <- as.numeric(df$A) - 1
  A <- df$A
  y <- pmin(df$Y, tau)
  d <- df$status
  d[y==tau]=1

  d1=d[A==1]; d0=d[A==0]
  y1=y[A==1]; y0=y[A==0]

  fit1=my.func_surv(y1, 1-d1)
  fit0=my.func_surv(y0, 1-d0)

  w1=d1/rep(pmax(fit1$surv,0.001), table(y1))
  w0=d0/rep(pmax(fit0$surv,0.001), table(y0))


  fit <- fit.rmst.reg(df, varX, gf, tau)
  gamma <- coef(fit)
  vgamma <- vcov(fit)
  dfx <- t(apply(df[,varX], 1, gf))
  mdf0 <- cbind(1,0,dfx, 0*dfx)
  mdf1 <- cbind(1,1,dfx, 1*dfx)
  #mdf0 <- cbind(1,0,df[,paste0("X",1:nX)], 0*df[,paste0("X",1:nX)])
  #mdf1 <- cbind(1,1,df[,paste0("X",1:nX)], 1*df[,paste0("X",1:nX)])
  m0 <- as.matrix(mdf0)%*%matrix(gamma,ncol = 1); m0 <- as.vector(m0)
  m1 <- as.matrix(mdf1)%*%matrix(gamma,ncol = 1); m1 <- as.vector(m1)


  dat_mod <- data.frame(p = c(p[A==1], p[A==0]),
                        w = c(w1, w0),
                        Y = c(y1, y0),
                        A = c(A[A==1], A[A==0]),
                        mu1 = c(m1[A==1], m1[A==0]),
                        mu0 = c(m0[A==1], m0[A==0]))


  res <- geex::m_estimate(
    estFUN = m_fun,
    data = dat_mod,
    root_control = setup_root_control(start = c(0.5,0.5,0.5,0.5))
  )

  coef <- coef(res)
  mu1 <- coef[1] + coef[3]
  mu0 <- coef[2] + coef[4]
  delta <- mu0 - mu1

  vcov <- vcov(res)
  c1 <- matrix(c(1,1), nrow = 1, ncol = 2)
  c2 <- matrix(c(-1,1,-1,1), nrow = 1, ncol = 4)
  sd1 <- as.numeric(sqrt(c1%*%vcov[c(1,3),c(1,3)]%*%t(c1)))
  sd0 <- as.numeric(sqrt(c1%*%vcov[c(2,4),c(2,4)]%*%t(c1)))
  delta.sd <- as.numeric(sqrt(c2%*%vcov%*%t(c2)))

  EW.mu <- c(mu0,mu1,delta)
  EW.sd <- c(sd0,sd1,delta.sd)

  return(list(mu = EW.mu, sd = EW.sd))
}

