# Calcualte the calibration weights
#@@df: data.frame
#@@varX: names of calibrated covariates
#@@gf: vector function of calibrated covariates (e.g, moment function)
#@@M: sample moments of g(X) in the target population
EW.p.M <- function(df,varX,gf,M){

  require(Rsolnp)
  na <- length(unique(df$A))
  fn <- function(x){
    sum(x*log(x))
  }
  constraints <- c(M,1)
  p <- NULL
  for (a in 1:na){
    dfa <- df[df$A == a,]
    dfa_x <- dfa[,varX]
    gx_a <- apply(dfa_x, 1, gf)
    eqn <- function(x){
      c <- NULL
      for (i in 1:nrow(gx_a)){
        c <- c(c, sum(x*gx_a[i,]))
      }
      c <- c(c,sum(x))
      return(c)
    }
    nr <- nrow(dfa)
    x0 <- rep(1/nr,nr)
    sol <- solnp(pars = x0, fun = fn, eqfun = eqn, eqB = constraints,
                 LB = rep(0,nr), UB = rep(1,nr), control = list(trace=0, tol=1e-6))
    p <- c(p, sol$pars)
  }

  return(p)
}
