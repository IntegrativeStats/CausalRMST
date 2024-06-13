library(tidyverse)
library(nnet)
library(Rsolnp)
library(survRM2)
library(geex)
library(devtools)
#setwd("/Users/HuaKimi/Library/CloudStorage/Box-Box/ICSA2024 short course/Lab Session")
#source("source_CW_RMST.R")
install_github("https://github.com/IntegrativeStats/CausalRMST")


gen.dat <- function(n,beta,lambda,gamma){
  A <- rbinom(n, 1, 0.5)
  X1 <- X2 <- rep(NA, n)
  X1[A==0] <- runif(sum(A==0), 0, 0.9)
  X1[A==1] <- runif(sum(A==1), 0.1, 1)
  X2[A==0] <- rnorm(sum(A==0), 0.8, 1)
  X2[A==1] <- rnorm(sum(A==1), 1.2, 1)
  df <- data.frame(A = A, X1 = X1, X2 = X2)
  MM <- model.matrix(~ -1+A+X1+X2+A*X1+A*X2, df)
  # Exp
  u <- runif(n)
  time_e <- (-log(u)*exp(- MM %*% beta)/lambda)^(1/gamma)
  time_c <- rexp(sum(n))/0.1
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  df$A <- df$A + 1
  return(as.data.frame(df))
  #return(as.data.frame(df %>% group_by(A) %>% arrange(Y, .by_group = T)))
}





beta <- c(-0.5,-1,0.5,0.5,-0.5)
lambda <- 0.3
gamma <- 1
n <- 1000
set.seed(212)
df <- gen.dat(n, beta, lambda, gamma)
df <- df %>% group_by(A) %>% arrange(Y, .by_group = TRUE)
df <- as.data.frame(df)




## varX: selected covariates
nX <- 2
varX <- paste0("X",1:nX)
## gf: vector of moment functions for each calibrated covariates
## x & x^2 for continuous variable, and x for binary variable
gf <- function(x){
  res <- as.vector(c(x[1], x[1]^2, x[2], x[2]^2))
  return(as.numeric(res))
}

## M: vector of covariate moments in the target population (e.g., population average here)
M <- c(mean(df$X1), mean(df$X1^2),
       mean(df$X2), mean(df$X2^2))


# calibration weight
p_cw <- EW.p.M(df,varX,gf,M)





# CW-adjusted estimators of RMSTs
W.adj.RMST <- function(df_sample, p_cw, varX, gf, tau){
  adj_RMST <- adj_RMST_se <- adj_RMST_ci1 <- adj_RMST_ci2 <- matrix(NA, 5, 3)
  rownames(adj_RMST) <- rownames(adj_RMST_se) <-
    rownames(adj_RMST_ci1) <- rownames(adj_RMST_ci2) <-
    c("Naive","CW_KM","CW_GF","CW_HJ","CW_AG")
  colnames(adj_RMST) <- colnames(adj_RMST_se) <-
    colnames(adj_RMST_ci1) <- colnames(adj_RMST_ci2) <-
    c("A=1","A=2","Difference")

  ## Naive estimator
  p_naive <- rep(1,nrow(df_sample)) # weights all equal to 1 in the Naive estimator
  res_naive <- W.KM.Est(df_sample,p_naive,tau)
  adj_RMST[1,] <- res_naive$mu
  adj_RMST_se[1,] <- res_naive$sd
  adj_RMST_ci1[1,] <- res_naive$mu - qnorm(0.975)*res_naive$sd
  adj_RMST_ci2[1,] <- res_naive$mu + qnorm(0.975)*res_naive$sd

  ## CW-adjusted Kaplan-Meier (KM) estimator
  res_KM <- W.KM.Est(df_sample,p_cw,tau)
  adj_RMST[2,] <- res_KM$mu
  adj_RMST_se[2,] <- res_KM$sd
  adj_RMST_ci1[2,] <- res_KM$mu - qnorm(0.975)*res_KM$sd
  adj_RMST_ci2[2,] <- res_KM$mu + qnorm(0.975)*res_KM$sd

  ## CW-adjusted G-Formula (GF) estimator
  res_GF <- W.GF.Est(df_sample, p_cw, varX, gf, tau)
  #res_GF <- W.GF0.Est(df_sample, varX, gf, tau)
  adj_RMST[3,] <- res_GF$mu
  adj_RMST_se[3,] <- res_GF$sd
  adj_RMST_ci1[3,] <- res_GF$mu - qnorm(0.975)*res_GF$sd
  adj_RMST_ci2[3,] <- res_GF$mu + qnorm(0.975)*res_GF$sd

  ## CW-adjusted Hajek (HJ) estimator
  res_HJ <- W.HJ.Est(df_sample, p_cw, tau)
  adj_RMST[4,] <- res_HJ$mu
  adj_RMST_se[4,] <- res_HJ$sd
  adj_RMST_ci1[4,] <- res_HJ$mu - qnorm(0.975)*res_HJ$sd
  adj_RMST_ci2[4,] <- res_HJ$mu + qnorm(0.975)*res_HJ$sd

  ## CW-adjusted Augmented (AG) estimator
  res_AG <- W.AG.Est(df_sample, p_cw, varX, gf, tau)
  #res_AG <- W.AG0.Est(df_sample, p_cw, varX, gf, tau)
  adj_RMST[5,] <- res_AG$mu
  adj_RMST_se[5,] <- res_AG$sd
  adj_RMST_ci1[5,] <- res_AG$mu - qnorm(0.975)*res_AG$sd
  adj_RMST_ci2[5,] <- res_AG$mu + qnorm(0.975)*res_AG$sd

  return(list(RMST = as.data.frame(round(adj_RMST,3)),
              se = as.data.frame(round(adj_RMST_se,3)),
              ci1 = as.data.frame(round(adj_RMST_ci1,3)),
              ci2 = as.data.frame(round(adj_RMST_ci2,3))))

}

tau = 4
res <- W.adj.RMST(df, p_cw, varX, gf, tau)
res

