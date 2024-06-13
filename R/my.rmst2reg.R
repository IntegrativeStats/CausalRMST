my.rmst2reg=function(y, delta, arm, x, tau, w=rep(1,length(y))){

  n=length(y)
  x=as.matrix(cbind(1, x))
  p=length(x[1,])

  y0=pmin(y, tau)
  d0=delta
  d0[y0==tau]=1

  d10=d0[arm==1]
  d00=d0[arm==0]
  y10=y0[arm==1]
  y00=y0[arm==0]
  x1=x[arm==1,]
  x0=x[arm==0,]
  n1=length(d10)
  n0=length(d00)

  id1=order(y10)
  y10=y10[id1]
  d10=d10[id1]
  x1=x1[id1,]

  id0=order(y00)
  y00=y00[id0]
  d00=d00[id0]
  x0=x0[id0,]

  fitc1=my.func_surv(y10, 1-d10)
  fitc0=my.func_surv(y00, 1-d00)

  weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
  weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))

  w1=w[arm==1]
  w0=w[arm==0]
  w1=w1[id1]
  w0=w0[id0]
  weights=c(weights1, weights0)*c(w1,w0)


  fitt=lm(c(y10,y00)~ rbind(x1, x0)-1, weights=weights)

  return(fitt)
}
