mlts = function(x, y, gamma, ns = 500, nc = 10, delta = 0.01)
{ 
  d = dim(x)
  n = d[1]
  p = d[2]
  q = ncol(y) 
  h = floor(n*(1-gamma))+1
  obj0 = 1e10 
  for (i in 1:ns)
  { 
    sorted = sort(runif(n), na.last = NA, index.return = TRUE)
    istart = sorted$ix[1:(p+q)]
    xstart = x[istart,]
    ystart = y[istart,]
    bstart = solve(t(xstart)%*%xstart,t(xstart)%*%ystart) 
    sigmastart = (t(ystart-xstart%*%bstart))%*%(ystart-xstart%*%bstart)/q
    for (j in 1:nc)
    {
      res  =  y - x %*% bstart
      tres = t(res)
      dist2 = colMeans(solve(sigmastart,tres)*tres)
      sdist2 = sort(dist2,na.last = NA,index.return = TRUE)
      idist2 = sdist2$ix[1:h]
      xstart = x[idist2,]
      ystart = y[idist2,]
      bstart = solve(t(xstart)%*%xstart,t(xstart)%*%ystart)
      sigmastart = (t(ystart-xstart%*%bstart))%*%(ystart-xstart%*%bstart)/(h-p)
    }
    obj = det(sigmastart)
    if(obj < obj0)
    { 
      result.beta = bstart
      result.sigma = sigmastart
      obj0 = obj
    }
  }
  cgamma = (1-gamma)/pchisq(qchisq(1-gamma,q),q+2)
  result.sigma = cgamma * result.sigma
  res = y - x %*% result.beta
  tres = t(res)
  result.dres = colSums(solve(result.sigma,tres)*tres)
  result.dres = sqrt(result.dres)
  
  qdelta = sqrt(qchisq(1-delta,q))
  good  = (result.dres <= qdelta)
  xgood = x[good,]
  ygood = y[good,]
  result.betaR = solve(t(xgood)%*%xgood,t(xgood)%*%ygood)
  result.sigmaR = (t(ygood-xgood%*%result.betaR)) %*% (ygood-xgood%*%result.betaR)/(sum(good)-p)
  cdelta = (1-delta)/pchisq(qdelta^2,q+2)
  result.sigmaR = cdelta*result.sigmaR
  resR = y - x%*%result.betaR
  tresR = t(resR)
  result.dresR = colSums(solve(result.sigmaR,tresR)*tresR)
  result.dresR = sqrt(result.dresR)
  return( list( beta = result.beta, sigma = result.sigma, dres = result.dres,
                betaR = result.betaR, sigmaR = result.sigmaR, dresR = result.dresR ) )
}