#### OUTLIERS  ####
library(rmgarch)
library(MTS)
install.packages("wesanderson")
library(wesanderson)
install.packages("mvtnorm")
library(mvtnorm)


# ADDITIVE OUTLIERS
 
set.seed(1995)
AR.sm <- list(order = c(1,0,0), ar = c(0.8)) 
inov <- rnorm(120,0,1)
yt1 <- arima.sim(n = 120, model = AR.sm, innov = inov)
yt2 <- yt1
yt2[80] <- yt2[80] + 10
ts.plot(as.ts(yt2[20:120]),as.ts(yt1[20:120]),col = c("red","black"))


# INNOVATION OUTLIERS

set.seed(1995)
AR.sm <- list(order = c(1,0,0), ar = c(0.8)) 
inov <- rnorm(120,0,1)
xt1 <- arima.sim(n = 120, model = AR.sm, innov = inov)

inov[80] <- 10
inov2 <- inov
xt2 <- arima.sim(n = 120, model = AR.sm,innov = inov2)

ts.plot(as.ts(xt2[20:120]),as.ts(xt1[20:120]), col = c("red","black"))



#####

install.packages("dae")
install.packages("mvtnorm")
library(mvtnorm)
p1=matrix(c(0.1,0.4,0.4,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)
nobs <- 600
arlags <- c(1)
malags <- NULL
cnst <- NULL
phi <- p1
sigma <- sig
theta <- NULL
skip <- 200


k = nrow(sigma)
nT = nobs + skip
at = rmvnorm(nT, rep(0, k), sigma)
at[700,1] <- at[700,1] + 10
at[700,2] <- at[700,2] + 10
nar = length(arlags)
p = 0
if (nar > 0) {
  arlags = sort(arlags)
  p = arlags[nar]
}
q = 0


ist = max(p, q) + 1
zt = matrix(0, nT, k)
cnst = rep(0, k)

for (it in ist:nT) {
  tmp = matrix(at[it, ], 1, k)
  if (nar > 0) {
    for (i in 1:nar) {
      idx = (i - 1) * k
      phj = phi[, (idx + 1):(idx + k)]
      ztm = matrix(zt[it - arlags[i], ], 1, k)
      tmp = tmp + ztm %*% t(phj)
    }
  }
  zt[it, ] = cnst + tmp
}
zt = zt[(1 + skip):nT, ]
at = at[(1 + skip):nT, ]
VARMAsimulation <- list(series = zt, noises = at)

yt <- VARMAsimulation$series
yt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))

ts.plot(yt, col = c("red","black"))

#####
#####
#####



set.seed(1995)
phi=matrix(c(0.1,0.4,0.8,0.1),2,2)
sigma=matrix(c(1,0.9,0.9,1),2,2)
sigma2= matrix(c(10,0.9,0.9,1),2,2)
nobs <- 600
arlags <- c(1)
malags <- NULL
cnst <- NULL
theta <- NULL
skip <- 200
epsilon <- 0.01
k = nrow(sigma)
nT = nobs + skip
mat5 <- matrix(nrow = 2, ncol = 1000)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)


for (j in 1:1000){
  at <- rmvnorm(nT, rep(0, k), sigma)
  ot <- rbinom(nT,1,epsilon)
  onorm <- runif(nT,5,15)
  onorm <- cbind(onorm,onorm)
  #onorm <- rmvnorm(nT, rep(0, k), sigma2)
  vo <- ot*onorm
  at = at + vo
  
  nar = length(arlags)
  p = 0
  if (nar > 0) {
    arlags = sort(arlags)
    p = arlags[nar]
  }
  q = 0
  
  
  
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  cnst = rep(0, k)
  
  for (it in ist:nT) {
    tmp = matrix(at[it, ], 1, k)
    if (nar > 0) {
      for (i in 1:nar) {
        idx = (i - 1) * k
        phj = phi[, (idx + 1):(idx + k)]
        ztm = matrix(zt[it - arlags[i], ], 1, k)
        tmp = tmp + ztm %*% t(phj)
      }
    }
    zt[it, ] = cnst + tmp
  }
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsimulation <- list(series = zt, noises = at)}


ts.plot(VARMAsimulation$series, col = c("black","red"))
