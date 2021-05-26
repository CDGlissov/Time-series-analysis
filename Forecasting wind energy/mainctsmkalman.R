rm(list=ls())
library(msir)
library(mclust)
library(ctsmr)
library(verification)
library(stats)
library(rugarch)
library(xtable)
#### Load Data ####
setwd("C:/Users/Christian/Dropbox/7. Semester/Tidsrække Windpower/data")
X <- read.table("cex4WindDataInterpolated.csv", sep=",",
                header=TRUE, stringsAsFactors=FALSE)
X$t <- as.POSIXct(X$t, tz="UTC")
#### Tag kun nogle af dataen fordi jeg får errors ####
X = X[8500:17000,]
n<-nrow(X)
summary(X)
1/n*sum(abs(diff(X$p)))
sum(is.na(X))
#No nans
plot(complete.cases(X))


#### fit simpel ARMA(1,1) baseret på ACF og PACF ###
par(mfrow=c(2,2))
plot(X$t,X$p, type = 'l', main="Wind power over time", ylab="p (kW)", xlab="Time")
plot(X$t,X$Ws1, type = 'l', main="Wind speed over time", ylab="Ws (m/s)", xlab="Time")
plot(X$t,X$Wd1, type = 'l', main="Wind direction over time", ylab="Wd (deg)", xlab="Time")
plot(X$t,X$T1, type = 'l', main="Outdoor temperature over time", ylab="T (K)", xlab="Time")

par(mfrow=c(1,1))
plot(X$Ws1,X$p, type = 'p', main="Power curve for the wind speeds", ylab="p (kW)", xlab="Ws (m/s)")
points(X$Ws2,X$p,  col=2)
points(X$Ws3,X$p,  col=3)
legend("topright", c("1-Hour forecasted Ws", "2-Hour forecasted Ws", "3-Hour forecasted Ws"), col=c(1,2,3), pch=c(1,1,1))

Xtest<-X[7501:8501,]
X<-X[1:7500,]


#### Undersøg general powercurve for 1-hour windspeed prediction ####
# bar lav en simpel Loess uden optimal span og sådan #
plot(X$p ~ X$Ws1)
fit2 <- loess.sd(X$p ~ X$Ws1, span = bestspan)
points(X$Ws1, fit2$model$fitted, col = 2)
lines(fit2$x, fit2$upper, col = 2)
lines(fit2$x, fit2$lower, col = 2)
1/length(X$p)*sum(abs(fit2$model$residuals))

plot(X$p ~ X$t, type='l')
lines(X$t,fit2$model$fitted, col='red')

#### Gør det samme med 2-hour windspeed prediction ###
plot(X$p ~ X$Ws2)
fit3 <- loess.sd(X$p ~X$Ws2, span =0.25)
lines(fit3$x, fit3$y, col = 2)



#### undersøg om hvor vinden kommer fra har indflydelse.
#### Giv hver observation en direktion baseret på windderiktion for 1 time.
direction = rep(NA, length(X$p))
for (i in 1:length(X$p)) {
  if (X$Wd1[i] > 270) {
    direction[i] = 1
  }
  if (X$Wd1[i] <= 270 && X$Wd1[i]>180) {
    direction[i] = 2
  }
  if (X$Wd1[i] <= 180 && X$Wd1[i]>90) {
    direction[i] = 3
  }
  if (X$Wd1[i] <= 90) {
    direction[i] = 4
  }
}
X$direction = direction

# lav 4 subset for hver direction og fit loess
fst <- X[X$direction == 1,]
nd<- X[X$direction == 2,]
rd <- X[X$direction == 3,]
th <- X[X$direction == 4,]

fitst <- loess.sd(fst$p ~ fst$Ws1, span =0.3)
fitnd <- loess.sd(nd$p ~ nd$Ws1, span =0.3)
fitrd <- loess.sd(rd$p ~ rd$Ws1, span =0.3)
fitth <- loess.sd(th$p ~ th$Ws1, span =0.3)

## Nogle nice plots
par(mfrow=c(2,2))
plot(fst$p~fst$Ws1)
lines(fitst$x, fitst$y, col = 2)
plot(nd$p~nd$Ws1)
lines(fitnd$x, fitnd$y, col = 2)
plot(rd$p~rd$Ws1)
lines(fitrd$x, fitrd$y, col = 2)
plot(th$p~th$Ws1)
lines(fitth$x, fitth$y, col = 2)

#############################ARLS############################
###########################################################################
#CTSM 
options(digits=6)
X2 <- X[1:7500,]
ttest2 <- as.numeric(X2$t)/3600
X2$t <- ttest2
testme3<-windpow2(X2)

testme3$model$outputs

testme3


tmp<-predict(testme3)[[1]]
tmp <- tmp$output$pred$p

sum(X2$p-tmp)

plot(X2$t,X2$p/max(X2$p), type='l')
lines(X2$t,tmp, col='red')

##########################################KALMAN#######################
library(FKF)
## Fit Kalman filter ####
A <- matrix(c(1,0,0,1),nrow=2)
C <- matrix(c(1,0,0,1),nrow=2)
Sigma1 <- diag(0.01,nrow=2)
Sigma2 <- diag(0.01,nrow=2)

kftempt1 <- fkf(a0=c(0.1, 0.1), P0 = matrix(c(0.1,0,0,0.1),nrow=2), dt = matrix(c(0,0), nrow=2), Tt = A, ct=matrix(c(0,0), nrow=2), Zt = C, HHt = Sigma1, GGt = Sigma2,
                yt = rbind(X$p,X$Ws1))
par(mfrow=c(1,1))
plot(X$p, type = 'l', col = 1)
lines(kftempt1$att[1,], type = 'l', col = 2)

plot(X$Ws1, type = 'l', col = 1)
lines(kftempt1$att[2,], type = 'l', col = 2)

1/n*sum(abs(X$p-kftempt1$att[1,]))
plot(X$Ws1, X$p)
points(X$Ws1,kftempt1$att[1,], col='red')
plot(abs(X$p-kftempt1$att[1,]), type='l')

#######################################################################

