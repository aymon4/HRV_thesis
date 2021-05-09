robustvar = function(data, exogen = NULL, constant = TRUE, lags = lags, alpha = 0.01, ns = 500, nc = 10, delta = 0.01, cluster)
{
  T = dim(data)[1]
  nvar  = dim(data)[2]
  ydata = data[(lags+1):T,]
  if(constant){
    xdata = matrix(1, nrow = T - lags, ncol = 1)
    for(i in 1:lags){ xdata = cbind(xdata, data[((1+lags)-i):(T-i), ])}
  } else{
    xdata = NULL
    for(i in 1:lags){ xdata = cbind(xdata, data[((1+lags)-i):(T-i), ])}
  }
  if(!is.null(exogen)) xdata = cbind(xdata, tail(exogen, T - lags))
  if(!is.null(cluster)){
    datamlts = mlts.parallel(x = as.matrix(xdata), y = as.matrix(ydata), 
                             gamma = alpha, ns = ns, nc = nc, delta = delta, cluster = cluster)
  } else{
    datamlts = mlts(as.matrix(xdata), as.matrix(ydata), gamma = alpha, ns = ns, 
                    nc = nc, delta = delta)  
  }
  U = datamlts$sigmaR
  logdetSigma = log(det(U))
  cdelta = (1-delta)/pchisq(qchisq(1 - delta, nvar),nvar+2)
  ##### CHANGE IN THE CODE
  # from: m = sum(datamlts$nooutlier)
  # to: 
  m = sum(datamlts$nooutlier)
  phi = nvar*(lags*nvar+1)
  ll  = -(T-lags)/2 *(log(2*pi)*nvar+log(det(U)))- (m-nvar)*nvar/(2*cdelta)
  res  = datamlts
  res$AIC =-2/(T-lags)*(ll-phi)
  res$HQ =-2/(T-lags)*(ll-log(log(T-lags))*phi)
  res$SC =-2/(T-lags)*(ll-log(T-lags)*phi/2)
  return( res )
}

# VAR unconditional mean
.umeanvar = function(Bcoef, p){
  #Bcoef[,11]%*%t(solve(diag(10) - Bcoef[,-11]))
  m = dim(Bcoef)[1]
  Id = diag(m)
  C = Bcoef[, p * m]
  n = dim(Bcoef)[2]
  if(n > (m+1)){
    muEX = Bcoef[, ( p * m + 2 ):n, drop = FALSE] 
    muC  = Bcoef[, p * m + 1] + apply(muEX, 1, "sum")
  } else{
    muEX = NULL
    muC  = Bcoef[, p * m + 1]
  }
  idx = t(rugarch:::.embed(1:(m*p), k=m, by=m, ascending = TRUE))
  for(i in 1:p) Id = Id - C[,idx[,i]]
  muC %*% t(solve(Id - C))
}
