my_akm_rmst <- function(time, status, weight=NULL, tau=NULL){

  data <- data.frame(time, status, weight)

  #--- AKM ---
  # Based on 'adjusted.KM' function from {IPWsurvival} package
  # Author: F. Le Borgne and Y. Foucher
  tj <- c(0,sort(unique(data$time[data$status==1])))
  dj <- sapply(tj, function(x){sum(data$weight[data$time==x & data$status==1])})
  yj <- sapply(tj, function(x){sum(data$weight[data$time>=x])})
  st <- cumprod(1-(dj/yj))
  m <- sapply(tj, function(x){sum((data$weight[data$time>=x])^2)})
  mj <- ((yj^2)/m)
  #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
  ft <- data.frame(tj, yj, dj, st, mj)

  #--- RMST ---
  # Based on 'rmst1 function' from {survRM2} package
  # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
  rtime <- ft$tj<=tau
  tj_r <- sort(c(ft$tj[rtime],tau))
  st_r <- ft$st[rtime]
  yj_r <- ft$yj[rtime]
  dj_r <- ft$dj[rtime]
  time_diff <- diff(c(0, tj_r))
  areas <- time_diff * c(1, st_r)
  rmst <- sum(areas)

  #--- Variance ---
  mj_r <- ft$mj[rtime]
  var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
  #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
  var_r <- c(var_r,0)
  rmst_var <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])

  return(data.frame(mu = rmst, V = rmst_var))
}
