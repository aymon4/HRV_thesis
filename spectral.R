library(MTS)
library(astsa)
#=================
#VAR analysis
#=================
# Two necessary functions


# Obtain Spectrum from autoregressive model
ARspect <- function(phi,sigma,freq){
  dimen <- dim(phi)[2]
  len <- dim(phi)[1]
  bigphi <- array(0,dim=c(dimen,dimen,length(freq)))
  spect <- array(0,dim=c(dimen,dimen,length(freq)))
  for(k in 1:length(freq)){
    bigphi[,,k] <- diag(dimen)
    for(j in 1:(len/dimen)){
      if(j==1){
        bigmat <- phi[1:dimen,]*exp(-2*pi*(1i)*freq[k])
      }else{
        bigmat <- phi[((dimen*j-(dimen-1)):(dimen*j)),]*exp(-2*j*pi*(1i)*freq
                                                            [k])
      }
      bigphi[,,k] = bigphi[,,k] - bigmat
    }
    spect[,,k] = solve(bigphi[,,k])%*%sigma%*%solve(Conj(t(bigphi[,,k])))
  }
  return(spect/(2*pi))
}
#############



set.seed(1995)
p1=matrix(c(0.1,0.4,0.8,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)
mat <- matrix(nrow = 2,ncol = 10)
epsilon <- 0.05
coefa <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
coefb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
cova <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)
covb <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

for (i in 1:10){ 
  m1=VARMAsim(600,arlags=c(1),phi=p1,sigma=sig)
  zt <- m1$series
  zt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))
  v1 <- rbinom(600,1,epsilon)
  norm1 <- rnorm(600,0,10) # variance chosen to be 10
  vt1 <- v1*norm1
  v <- cbind(vt1,vt1)
  zt <- zt + v
  r <- VAR(zt, p = 1, include.mean = F)
  rb <- varxfit2(zt,p=1,constant = F,robust = T)
  mat[1,i] <- sum((t(r$coef)-p1)^2)
  mat[2,i] <- sum((rb$Bcoef-p1)^2)
  coefa <- coefa + t(r$coef)
  coefb <- coefb + rb$Bcoef
  cova <- cova + r$Sigma
  covb <- covb + rb$SigmaR
}




freal <- ARspect(p1,sig,freq)
fa <- ARspect(coefa/10,cova/10,freq)
fb <- ARspect(coefb/10,covb/10,freq)




plot(freq,Re(freal[1,1,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq,Re(f2[1,1,]), type = 'l', col = "red")
lines(freq,Re(f[1,1,]), type = 'l', col = "blue")




#############3
set.seed(1995)
p1=matrix(c(0.1,0.4,0.8,0.1),2,2)
sig=matrix(c(1,0.9,0.9,1),2,2)

epsilon <- 0.05


  m1=VARMAsim(600,arlags=c(1),phi=p1,sigma=sig)
  zt <- m1$series
  zt <- ts.union(as.ts(zt[,1]),as.ts(zt[,2]))



  model <- VAR(zt, p = 1, include.mean = F)
  model2 <- varxfit2(zt,p = 1,robust = T, constant = F)
n <- nrow(zt)
nfreq <- floor(n/2)
freq <- (0:nfreq)/(2*nfreq)
model$Sigma
model2$SigmaR
f <- ARspect(model$coef,model$Sigma,freq)
f2 <- ARspect(p1,sig,freq)
f

plot(freq,Re(f[1,1,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
plot(freq,Re(f2[2,2,]),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")

spec.pgram(zt,taper=0.1,pad = 0,fast=TRUE,demean=FALSE,detrend=TRUE)
k<- kernel("daniell", 30)

spec.pgram(zt,kernel=k,taper=0,pad=0,fast=TRUE,
           demean=FALSE,detrend=TRUE,plot=TRUE)

# Smoothed and tapered periodogram #
spec.pgram(zt,kernel=k,taper=0.1,pad=0,fast=TRUE,
           demean=FALSE,detrend=TRUE,plot=TRUE)
