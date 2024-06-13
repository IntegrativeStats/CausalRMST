# weighted Hajek RMST estimator
#@@df: data.frame
#@@p: weights (e.g., calibration weights, inverse propensity score weights)
#@@tau: truncation time t* for RMST
W.HJ.Est <- function(df, p, tau){
  require(geex)
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    function(theta){
      c(p*A*w*(Y-theta[1]),
        p*(1-A)*w*(Y-theta[2]))
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


  dat_mod <- data.frame(p = c(p[A==1], p[A==0]),
                        w = c(w1, w0),
                        Y = c(y1, y0),
                        A = c(A[A==1], A[A==0]))


  res <- geex::m_estimate(
    estFUN = m_fun,
    data = dat_mod,
    root_control = setup_root_control(start = c(0.5,0.5))
  )

  coef <- coef(res)
  mu1 <- coef[1]
  mu0 <- coef[2]
  delta <- mu0 - mu1

  vcov <- vcov(res)
  sd1 <- sqrt(vcov[1,1])
  sd0 <- sqrt(vcov[2,2])
  c1 = matrix(c(-1,1), nrow = 1, ncol = 2)
  delta.sd <- as.numeric(sqrt(c1%*%vcov%*%t(c1)))

  EW.mu <- c(mu0,mu1,delta)
  EW.sd <- c(sd0,sd1,delta.sd)

  return(list(mu = EW.mu, sd = EW.sd))
}
