###GARCH

rm(list=ls())
library(msir)
library(mclust)
library(ctsmr)
library(verification)
library(stats)
library(rugarch)

results <- function(X,pred,pred_se){
  
  n <- length(pred)
  
  RMSE = sqrt(1/n * sum( (X$p-pred)^2 ) )
  MAE = 1/n * sum(abs( X$p-pred ))  
  
  crps <- c()
  for(j in 1:n){ #for each attribute
    mu <- pred[j]
    sigma <- pred_se[j]
    obs<-X$p[j]
    z <- (obs - mu)/sigma
    crps[j] <- sigma * (z * (2 * pnorm(z, 0, 1) - 1) +
                          2 * dnorm(z, 0, 1) - 1/sqrt(pi))
    
  }
  CRPS<-mean(crps)
  out = c(RMSE,MAE,CRPS)
  return(out)
}
#### Load Data ####
setwd("C:/Users/Christian/Dropbox/7. Semester/Tidsrække Windpower/data")
X <- read.table("cex4WindDataInterpolated.csv", sep=",",
                header=TRUE, stringsAsFactors=FALSE)
X$t <- as.POSIXct(X$t, tz="UTC")
#### Tag kun nogle af dataen fordi jeg får errors ####
X = X[8500:17000,]
n<-nrow(X)


acf(X$p)
pacf(X$p)

fit1 <- arima(X$p, order = c(2,0,0), xreg=X$Ws1)
fit1
par(mfrow=c(2,2))
acf(fit1$residuals, main="ACF of residuals")
pacf(fit1$residuals, main="PACF of residuals")
#still see lag in squared, need garch
acf(fit1$residuals^2, main="ACF of squared residuals")
pacf(fit1$residuals^2, main="PACF of squared residuals")


#GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1)),
                    mean.model     = list(armaOrder = c(2, 0),
                                          external.regressors = as.matrix(X$Ws1)),
                    distribution.model = "norm")

garch1 <- ugarchfit(spec = spec1, data = X$p, out.sample = 1000)

par(mfrow=c(2,2))
acf((residuals(garch1)/sigma(garch1))^2,main="ACF of squared standardized residuals")
pacf((residuals(garch1)/sigma(garch1))^2,main="ACF of squared standardized residuals")
cpgram((residuals(garch1)/sigma(garch1)), main="Cpgram of standardized residuals")
cpgram((residuals(garch1)/sigma(garch1))^2, main="Cpgram of squared standardized residuals")
par(mfrow=c(1,1))
plot(garch1@fit$residuals/garch1@fit$sigma, type='l', xlab="Time", ylab="Standardized residuals", main="ARX(2)-GARCH(1,1) model")

par(mfrow=c(1,1))
pred2<-garch1@fit$fitted.values
pred2[which(pred2<0)]=0
sum(abs(X$p-pred2))
plot(X$Ws1, X$p)
points(X$Ws1,pred2, col='red')

a=7200
b=7500
c=7500+501
plot(X$t[a:c], X$p[a:c], type='l', ylim=c(-10,40), lwd=2, xlab="Time", ylab='p (kW)', main="In and out of sample predictions for wind power using 1-hr forecasted wind speed")
lines(X$t[a:b],pred2[a:b], col='blue', lwd=2)
lines(X$t[a:b], pred2[a:b]+2*garch1@fit$sigma[a:b], col='blue', lty=2)
lines(X$t[a:b], pred2[a:b]-2*garch1@fit$sigma[a:b], col='blue', lty=2)

modelfor=ugarchforecast(garch1, data = NULL, n.ahead = 1, n.roll
                        = 1000, out.sample = 1000)
lines(X$t[b:c],as.numeric(fitted(modelfor))[1:502], col="red", lwd=2)
lines(X$t[b:c],(as.numeric(fitted(modelfor))+2*as.numeric(sigma(modelfor)))[1:502], col="red", lty=2)
lines(X$t[b:c],(as.numeric(fitted(modelfor))-2*as.numeric(sigma(modelfor)))[1:502], col="red", lty=2)
legend("topright", c("Fitted Values", "95% Confidence intervals", 
                     "Out-Of-Sample Predictions", "95% Prediction interval"), col=c("blue","blue","red","red"), lty=c(1,2,1,2))

npred <- length(fitted(modelfor))

1/npred*sum(abs(diff(X$p[7501:8501])))

#Evaluate
mae<-1/npred*sum(abs(X$p[7501:8501]-as.numeric(fitted(modelfor))))
rmse<-sqrt(1/npred*sum((X$p[7501:8501]-as.numeric(fitted(modelfor)))^2))

crps <- c()
for(j in 1:npred){ #for each attribute
  mu <- as.numeric(fitted(modelfor))[j]
  sigma <- as.numeric(sigma(modelfor))[j]
  obs<-X$p[7500+j]
  z <- (obs - mu)/sigma
  crps[j] <- sigma * (z * (2 * pnorm(z, 0, 1) - 1) +
                        2 * dnorm(z, 0, 1) - 1/sqrt(pi))
}
CRPS <- mean(crps)

mae
rmse
CRPS
plot(crps)

#MAKE FOR 3 step
modelfor=ugarchforecast(garch1, data = NULL, n.ahead = 3, n.roll
                        = 1000, out.sample = 1000)
pred2<-as.numeric(fitted(modelfor)[2,])
pred3<-as.numeric(fitted(modelfor)[3,])
pred2.sig<-as.numeric(sigma(modelfor)[2,])
pred3.sig<-as.numeric(sigma(modelfor)[3,])

mae<-1/npred*sum(abs(X$p[7501:8501]-pred2))
rmse<-sqrt(1/npred*sum((X$p[7501:8501]-pred2))^2)
mae
rmse

mae<-1/npred*sum(abs(X$p[7501:8501]-pred3))
rmse<-sqrt(1/npred*sum((X$p[7501:8501]-pred3))^2)
mae
rmse

crps <- c()
for(j in 1:npred){ #for each attribute
  mu <- pred2[j]
  sigma <- pred2.sig[j]
  obs<-X$p[7500+j]
  z <- (obs - mu)/sigma
  crps[j] <- sigma * (z * (2 * pnorm(z, 0, 1) - 1) +
                        2 * dnorm(z, 0, 1) - 1/sqrt(pi))
}
CRPS <- mean(crps)

crps <- c()
for(j in 1:npred){ #for each attribute
  mu <- pred3[j]
  sigma <- pred3.sig[j]
  obs<-X$p[7500+j]
  z <- (obs - mu)/sigma
  crps[j] <- sigma * (z * (2 * pnorm(z, 0, 1) - 1) +
                        2 * dnorm(z, 0, 1) - 1/sqrt(pi))
}
CRPS <- mean(crps)

par(mfrow=c(1,3))
plot(crps, type='l', xlab="Iteration", ylab="CRPS", main="Verification: CRPS")
plot(abs(X$p[7501:8501]-as.numeric(fitted(modelfor))), type='l', xlab="Iteration", ylab="Absolute Error", main="Verification: Absolute Error")
plot((X$p[7501:8501]-as.numeric(fitted(modelfor)))^2, type='l', xlab="Iteration", ylab="Squared Error", main="Verification: Squared Error")
      

spec1 <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1)),
                    mean.model     = list(armaOrder = c(2, 0),
                                          external.regressors = as.matrix(X$Ws2)),
                    distribution.model = "norm")

garch1 <- ugarchfit(spec = spec1, data = X$p, out.sample = 1000)
modelfor=ugarchforecast(garch1, data = NULL, n.ahead = 1, n.roll
                        = 1000, out.sample = 1000)
pred2<-as.numeric(fitted(modelfor)[1,])

pred2.sig<-as.numeric(sigma(modelfor)[1,])


results(X[7501:8501,],pred2, pred2.sig)


spec1 <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1)),
                    mean.model     = list(armaOrder = c(2, 0),
                                          external.regressors = as.matrix(X$Ws3)),
                    distribution.model = "norm")

garch1 <- ugarchfit(spec = spec1, data = X$p, out.sample = 1000)
modelfor=ugarchforecast(garch1, data = NULL, n.ahead = 1, n.roll
                        = 1000, out.sample = 1000)
pred3<-as.numeric(fitted(modelfor)[1,])
pred3.sig<-as.numeric(sigma(modelfor)[1,])
results(X[7501:8501,],pred3, pred3.sig)
