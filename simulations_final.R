#############################################
###                                       ###  
#  FINAL R SCRIPT FOR THE SIMULATION STUDY  #
###                                       ###  
#############################################


library(rmgarch)
library(MTS)


n <- 600
nfreq <- floor(n/2)
freq <- (0:nfreq)/(2*nfreq)


## FIRST PART: SIMULATION OF A VAR WITH HIGH CORRELATION ##

## ADDITIVE OUTLIERS ##

#MULTI#

# 1%

set.seed(1995)
p1=matrix(c(0.1,0.4,0.8,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)
mat <- matrix(nrow = 2,ncol = 1000)
epsilon <- 0.01
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

for (i in 1:1000){ 
   m1=VARMAsim(600,arlags=c(1),phi=p1,sigma=sig)
   zt <- m1$series
   zt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))
   v1 <- rbinom(600,1,epsilon)
   norm1 <- rnorm(600,0,10) # variance chosen to be 10
   vt1 <- v1*norm1
   v <- cbind(vt1,vt1)
   zt <- zt + v
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.01)
   mat[1,i] <- sum((t(r$coef)-p1)^2)
   mat[2,i] <- sum((rb$Bcoef-p1)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}

mat
boxplot(mat[1,],mat[2,],ylab = "Squared error", names = c("(1)","(2)"))


freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

# 0.5%

set.seed(1995)
p1=matrix(c(0.1,0.4,0.8,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)
mat2 <- matrix(nrow = 2,ncol = 1000)
epsilon2 <- 0.005
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

for (i in 1:1000){ 
   m1=VARMAsim(600,arlags=c(1),phi=p1,sigma=sig)
   zt <- m1$series
   zt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))
   v1 <- rbinom(600,1,epsilon2)
   norm1 <- rnorm(600,0,10) # variance chosen to be 10
   vt1 <- v1*norm1
   v <- cbind(vt1,vt1)
   zt <- zt + v
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.005)
   mat2[1,i] <- sum((t(r$coef)-p1)^2)
   mat2[2,i] <- sum((rb$Bcoef-p1)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}

mat2
boxplot(mat2[1,],mat2[2,],ylab = "Squared error", names = c("(1)","(2)"))

freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

# 5%

set.seed(1995)
p1=matrix(c(0.1,0.4,0.8,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)
mat3 <- matrix(nrow = 2,ncol = 1000)
epsilon3 <- 0.01
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

for (i in 1:1000){ 
   m1=VARMAsim(600,arlags=c(1),phi=p1,sigma=sig)
   zt <- m1$series
   zt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))
   v1 <- rbinom(600,1,epsilon3)
   norm1 <- rnorm(600,0,10) # variance chosen to be 10
   vt1 <- v1*norm1
   v <- cbind(vt1,vt1)
   zt <- zt + v
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.05)
   mat3[1,i] <- sum((t(r$coef)-p1)^2)
   mat3[2,i] <- sum((rb$Bcoef-p1)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}

mat3
boxplot(mat3[1,],mat3[2,],ylab = "Squared error", names = c("(1)","(2)"))

freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)




## UNIVARIATE 1% EACH ##

set.seed(1995)
p1=matrix(c(0.1,0.4,0.8,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)
mat4 <- matrix(nrow = 2,ncol = 1000)
epsilon <- 0.01
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

for (i in 1:1000){ 
   m1=VARMAsim(600,arlags=c(1),phi=p1,sigma=sig)
   zt <- m1$series
   zt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))
   v1 <- rbinom(600,1,epsilon)
   norm1 <- rnorm(600,0,10) # variance chosen to be 10
   vt1 <- v1*norm1
   v2 <- rbinom(600,1,epsilon)
   norm2 <- rnorm(600,0,10)
   vt2 <- v2*norm2
   v <- cbind(vt1,vt2)
   zt <- zt + v
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.05)
   mat4[1,i] <- sum((t(r$coef)-p1)^2)
   mat4[2,i] <- sum((rb$Bcoef-p1)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}

mat4
boxplot(mat4[1,],mat4[2,],ylab = "Squared error", names = c("(1)","(2)"))
ts.plot(zt)

freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)


### INNOVATION OUTLIERS ###

# MULTI

# 5%

install.packages("mvtnorm")
library(mvtnorm)

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
epsilon <- 0.05
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
   VARMAsimulation <- list(series = zt, noises = at)
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.05)
   mat5[1,j] <- sum((t(r$coef)-phi)^2)
   mat5[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat5
boxplot(mat5[1,],mat5[2,],ylab = "Squared error", names = c("(1)","(2)"))
freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

# 1%


 
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
epsilon2 <- 0.01
k = nrow(sigma)
nT = nobs + skip
mat6 <- matrix(nrow = 2, ncol = 1000)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)


for (j in 1:1000){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon2)
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
   VARMAsimulation <- list(series = zt, noises = at)
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.01)
   mat6[1,j] <- sum((t(r$coef)-phi)^2)
   mat6[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat6
boxplot(mat6[1,],mat6[2,],ylab = "Squared error", names = c("(1)","(2)"))
freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)


# 0.5%



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
epsilon3 <- 0.005
k = nrow(sigma)
nT = nobs + skip
mat7 <- matrix(nrow = 2, ncol = 1000)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)


for (j in 1:1000){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon3)
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
   VARMAsimulation <- list(series = zt, noises = at)
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.005)
   mat7[1,j] <- sum((t(r$coef)-phi)^2)
   mat7[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat7
boxplot(mat7[1,],mat7[2,],ylab = "Squared error", names = c("(1)","(2)"))
freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)



# UNIVARIATE

# 5%

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
epsilon <- 0.05
k = nrow(sigma)
nT = nobs + skip
mat8 <- matrix(nrow = 2, ncol = 100)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)


for (j in 1:1000){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon)
   onorm <- runif(nT,5,15)
   #onorm <- rmvnorm(nT, rep(0, k), sigma2)
   vo <- ot*onorm
   at[,1] = at[,1] + vo
   
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
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.05)
   mat8[1,j] <- sum((t(r$coef)-phi)^2)
   mat8[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat8
boxplot(mat8[1,],mat8[2,],ylab = "Squared error", names = c("(1)","(2)"))

freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/100,cova/100,freq)
fb <- ARspect(coefb/100,covb/100,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)



# 1%

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
mat9 <- matrix(nrow = 2, ncol = 100)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)


for (j in 1:100){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon)
   onorm <- runif(nT,5,15)
   #onorm <- rmvnorm(nT, rep(0, k), sigma2)
   vo <- ot*onorm
   at[,1] = at[,1] + vo
   
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
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.01)
   mat9[1,j] <- sum((t(r$coef)-phi)^2)
   mat9[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat9
boxplot(mat9[1,],mat9[2,],ylab = "Squared error", names = c("(1)","(2)"))
freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/100,cova/100,freq)
fb <- ARspect(coefb/100,covb/100,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)


# 0.5%

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
epsilon <- 0.005
k = nrow(sigma)
nT = nobs + skip
mat10 <- matrix(nrow = 2, ncol = 100)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

for (j in 1:100){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon)
   onorm <- runif(nT,5,15)
   #onorm <- rmvnorm(nT, rep(0, k), sigma2)
   vo <- ot*onorm
   at[,1] = at[,1] + vo
   
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
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.005)
   mat10[1,j] <- sum((t(r$coef)-phi)^2)
   mat10[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat10
boxplot(mat10[1,],mat10[2,],ylab = "Squared error", names = c("(1)","(2)"))

freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/100,cova/100,freq)
fb <- ARspect(coefb/100,covb/100,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)



### CORRELATION OUTLIERS ###

# 10%

set.seed(1995)
phi=matrix(c(0.1,0.4,0.8,0.1),2,2)
sigma=matrix(c(1,0.9,0.9,1),2,2)
sigma2= matrix(c(1,-0.7,-0.7,1),2,2)
nobs <- 600
arlags <- c(1)
malags <- NULL
cnst <- NULL
theta <- NULL
skip <- 200
epsilon <- 0.1
k = nrow(sigma)
nT = nobs + skip
mat11 <- matrix(nrow = 2, ncol = 1000)
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)


for (j in 1:1000){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon)
   ot <- cbind(ot,ot)
   onorm <- rmvnorm(nT, c(2,-2), sigma2)
   vo <- ot*onorm
   for (l in 1:nT){
      if (vo[l,] != c(0,0)){
         at[l,] <- vo[l,]
      }
   }
   
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
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.1)
   mat11[1,j] <- sum((t(r$coef)-phi)^2)
   mat11[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat11
boxplot(mat11[1,],mat11[2,],ylab = "Squared error", names = c("(1)","(2)"))
freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)




# 5%

set.seed(1995)
phi=matrix(c(0.1,0.4,0.8,0.1),2,2)
sigma=matrix(c(1,0.9,0.9,1),2,2)
sigma2= matrix(c(1,-0.7,-0.7,1),2,2)
nobs <- 600
arlags <- c(1)
malags <- NULL
cnst <- NULL
theta <- NULL
skip <- 200
epsilon <- 0.05
k = nrow(sigma)
nT = nobs + skip
mat12 <- matrix(nrow = 2, ncol = 1000)


for (j in 1:1000){
   at <- rmvnorm(nT, rep(0, k), sigma)
   ot <- rbinom(nT,1,epsilon)
   ot <- cbind(ot,ot)
   onorm <- rmvnorm(nT, c(2,-2), sigma2)
   vo <- ot*onorm
   for (l in 1:nT){
      if (vo[l,] != c(0,0)){
         at[l,] <- vo[l,]
      }
   }
   
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
   r <- VAR(zt, p = 1, include.mean = F)
   rb <- varxfit2(zt,p=1,constant = F,robust = T, gamma = 0.05)
   mat12[1,j] <- sum((t(r$coef)-phi)^2)
   mat12[2,j] <- sum((rb$Bcoef-phi)^2)
   coefa <- coefa + t(r$coef)
   coefb <- coefb + rb$Bcoef
   cova <- cova + r$Sigma
   covb <- covb + rb$SigmaR
}
mat12
boxplot(mat12[1,],mat12[2,],ylab = "Squared error", names = c("(1)","(2)"))

freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/1000,cova/1000,freq)
fb <- ARspect(coefb/1000,covb/1000,freq)

plot(freq,Re(freal[1,1,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[1,1,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[1,1,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)

plot(freq,Re(freal[2,2,]),type='l',lty=2,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(fa[2,2,]), type = 'l', lty = 1, col = "red")
lines(freq,Re(fb[2,2,]), type = 'l', lty = 1,col = "blue")
legend("topright", inset=.1,
       c("Real","RMLTS","ML"), fill=c("black","blue","red"), horiz=TRUE, cex=0.6)


