#215
library(MTS)

setwd("~/Dropbox/UNIGE/MASTER_THESIS/215")

l1 <- read.table("215_0.txt")
l2 <- read.table("215_1.txt")

l1 <- as.ts(l1)
l2 <- as.ts(l2)


l1 <- as.ts(l1[201:1200,])
l2 <- as.ts(l2[201:1200,])
ts.plot(l1-mean(l1), l2-mean(l2), col = c("red","black"), ylab = "RR (in sec)", xlab = "Observation number")
legend("topleft", inset=.04, title="Lead",
       c("MLII","V1"), fill=c("red","black"), horiz=TRUE, cex=0.6)

mean(l1)*1000
mean(l2)*1000

final100 <- ts.union(l1,l2)
final10 <- cbind(l1-mean(l1),l2-mean(l2))




n <- length(final100)
sum(final100>0.9)
sum(final100<0.60)
200/1982
matt1 <- matrix(nrow = 3,ncol = 16)
for (i in 1:30){
  rmlts1 <- varxfit2(final10, p = i, constant = F, robust = TRUE, gamma = 0.1)
  matt1[1,i] <- rmlts1$AIC
  matt1[2,i] <- rmlts1$BIC
  matt1[3,i] <- rmlts1$SC
}
round(matt1[1:3,],2)
which.min(matt1[1,])
which.min(matt1[2,])
which.min(matt1[3,])
#14 for all 3 info criteria

rmlts_fin <- varxfit2(final100, p = 14, constant = F, robust = TRUE, gamma =0.1)




rrr <- robustvar(final100, constant = F, lags = 14, alpha = 0.1,cluster = NULL)

delta <- 0.01
plot(rrr$dres, type = 'l', xlab = "Observation number", ylab = "Robust residual distance")
abline(h = sqrt(qchisq(1-delta,2)), lty = 3)
sum(rrr$dres>sqrt(qchisq(1-0.01,2)))
# 120 observations are considered as outliers.

# Too slow to compute -> use the equivalent VAR function which does the LS estimation
matt2 <- matrix(nrow = 2,ncol = 12)
for (i in 1:12){
  mts <- VARMA(final100, p = i, q = 0, include.mean = F)
  matt2[1,i] <- mts$aic
  matt2[2,i] <- mts$bic
}

matt2
which.min(matt2[1,])
which.min(matt2[2,])

matt2 <- matrix(nrow = 3,ncol = 20)
for (i in 1:20){
  mts <- VAR(final100, p = i, include.mean = F)
  matt2[1,i] <- mts$aic
  matt2[2,i] <- mts$bic
  matt2[3,i] <- mts$hq
}

round(matt2,2)
which.min(matt2[1,])
which.min(matt2[2,])
which.min(matt2[3,])


finala <- VAR(final10, p = 10, include.mean = F)
finalb <- varxfit2(final10, p = 14, robust = TRUE, gamma = 0.25, constant = F)

n <- nrow(final10)
nfreq1 <- floor(n/2)
freq1 <- (0:nfreq1)/(2*nfreq1)

faa <- ARspect(finala$coef,finala$Sigma,freq1)
fbb <- ARspect(t(finalb$Bcoef),finalb$SigmaR,freq1)

faa


plot(freq1,Re(fbb[1,1,]*450),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq1,Re(faa[1,1,]*100), type = 'l', lty = 1, col = "red")
legend("topright", inset=.1,
       c("RMLTS","ML"), fill=c("black","red"), horiz=TRUE, cex=0.6)

plot(freq1,Re(fbb[2,2,]*450),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq1,Re(faa[2,2,]*100), type = 'l', lty = 1, col = "red")
legend("topright", inset=.1,
       c("RMLTS","ML"), fill=c("black","red"), horiz=TRUE, cex=0.6)



############################

#### RHRV Analysis ####

############################

setwd("~/Dropbox/UNIGE/MASTER_THESIS/215")

library(RHRV)
library(MTS)
rm("hrv.data.l1","hrv.data.l2")

hrv.data.l1 <- CreateHRVData()
hrv.data.l1 <- SetVerbose(hrv.data.l1, TRUE)

hrv.data.l2 <- CreateHRVData()
hrv.data.l2 <- SetVerbose(hrv.data.l2, TRUE)

hrv.data.l1 <- LoadBeatWFDB(hrv.data.l1,"215_0", annotator = "wqrs")
hrv.data.l2 <- LoadBeatWFDB(hrv.data.l2,"215_1", annotator = "wqrs")

hrv.data.l1 <- BuildNIHR(hrv.data.l1)
hrv.data.l2 <- BuildNIHR(hrv.data.l2)



# # Automatic filtering
# plot(hrv.data.l1$Beat[[3]][200:1200]-mean(hrv.data.l1$Beat[[3]][200:1200]), type ="l",
#      ylab = "RR in ms (demeaned)", xlab = "Observation number",ylim = c(-340,440) )
# hrv.data.l1 <- FilterNIHR(hrv.data.l1)
# plot(hrv.data.l1$Beat[[3]][200:1200]-mean(hrv.data.l1$Beat[[3]][200:1200]), type ="l")
# hrv.data.l1 <- FilterNIHR(hrv.data.l1)
# plot(hrv.data.l1$Beat[[3]][200:1200]-mean(hrv.data.l1$Beat[[3]][200:1200]), type ="l",
#      ylab = "RR in ms (demeaned)", xlab = "Observation number", ylim = c(-340,440))
# 
# # Manual filtering
# 
# hrv.data.l1 <- EditNIHR(hrv.data.l1)
# plot(hrv.data.l1$Beat[[3]][200:1200]-mean(hrv.data.l1$Beat[[3]][200:1200]), type ="l",
#      ylab = "RR in ms (demeaned)", xlab = "Observation number", ylim = c(-340,440))

hrv.data.l1 <- FilterNIHR(hrv.data.l1) #3234
hrv.data.l2 <- FilterNIHR(hrv.data.l2) #3256

hrv.data.l1 <- EditNIHR(hrv.data.l1) #3225 (14 to be removed)
hrv.data.l2 <- EditNIHR(hrv.data.l2) #3220 (36 to be removed)

final1 <- as.ts(hrv.data.l1$Beat[[3]][201:1200]-mean(hrv.data.l1$Beat[[3]][200:1200]))
final2 <- as.ts(hrv.data.l2$Beat[[3]][201:1200]-mean(hrv.data.l2$Beat[[3]][200:1200]))

ts.plot(final1,final2, col = c("red","black"), ylab = "RR (in ms)", xlab = "Observation number")
legend("bottomright", inset=.03, title="Lead",
      c("MLII","V1"), fill=c("red","black"), horiz=TRUE, cex=0.6)


final <- cbind(final1,final2)

matt3 <- matrix(ncol = 30,nrow = 3)

for (i in 1:30){
  mtss <- VAR(final,p = i, include.mean = F)
  matt3[1,i] <- mtss$aic
  matt3[2,i] <- mtss$bic
  matt3[3,i] <- mtss$hq
}

round(matt3,2)
which.min(matt3[1,])
which.min(matt3[2,])
which.min(matt3[3,])

mtssfinal <- VAR(final, p = 5, include.mean = F)
finalbbb <- varxfit2(final, p = 14, robust = TRUE, gamma = 0.1, constant = F)

n <- nrow(final)
nfreq2 <- floor(n/2)
freq2 <- (0:nfreq2)/(2*nfreq2)





faaa <- ARspect(mtssfinal$coef,mtssfinal$Sigma,freq2)
fbbb <- ARspect(t(finalbbb$Bcoef),finalbbb$SigmaR,freq2)

plot(freq2,Re(fbb[1,1,]*450),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq2,Re(faaa[1,1,]/1000), type = 'l', lty = 1, col = "red")
legend("topright", inset=.1,
       c("RMLTS","ML"), fill=c("black","red"), horiz=TRUE, cex=0.6)

plot(freq2,Re(fbb[2,2,]*450),type='l',lty=1,ylab="Density",xlab="Frequency (Hz)")
lines(freq2,Re(faaa[2,2,]/2000), type = 'l', lty = 1, col = "red")
legend("topright", inset=.1,
       c("RMLTS","ML"), fill=c("black","red"), horiz=TRUE, cex=0.6)

