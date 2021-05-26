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
Xtest<-X[7501:8501,]
X<-X[1:8501,]

####ARLS


arls <- function(y, x, lambda){
  n<-length(y)
  np <- ncol(x)
  I <- diag(rep(1,np))
  theta1 <- matrix(0,ncol=1, nrow=np)
  store <- matrix(0, nrow=n, ncol=3+2*np)
  var1 <- 0
  colnames(store) <- c("e",paste("theta1",1:np,sep=""),paste("x",1:np,sep=""),"lambda", var1)
  store[1,] <- c(0,t(theta1), x[1,], lambda, var1)
  P<-I*0.1
  for(i in 2:n){
    xt<-t(x[i,,drop=FALSE])
    k <- P%*%xt / as.numeric((lambda + t(xt)%*%P%*%xt))
    e <- y[i] - t(xt)%*%theta1
    theta1 <- theta1 + k%*%e
    Ikx=I-k%*%t(xt)
    P=lambda^-1*Ikx%*%P;
    
    var1<-var(store[1:i,1])
    store[i,] <- c(e,t(theta1),t(xt),lambda, var1)
    
  }
  return(store)
}

obj.fun <- function(lambda, ynew, xnew){
  tmp <- arls(y=ynew, x=xnew, lambda=lambda)
  return(sum(tmp[-(1:100),1]^2))
}

opt<-optimise(obj.fun, interval = c(0.01,1), ynew=as.matrix(X$p), xnew=as.matrix(X$Ws1))
opt


test<-arls(as.matrix(X$p), as.matrix(X$Ws1), opt$minimum)
par(mfrow=c(1,1))
plot(X$Ws1,X$p, col='1', ylim=c(0,26), xlab="Wind Speed (m/s)", main="Power curve of predictions for the ARLS", ylab="p (kW)")
pred1 <- c()


for(i in 1:(length(test[,2])-1)){
  pred1[i] <-test[i,2]*X$Ws1[i+1]
}
pred1=c(0,pred1)
points(X$Ws1,pred1, col=2)
legend("topright", c("Observation", "Prediction"), col=c("1", "2"), pch=c(1,1))


par(mfrow=c(2,1))
plot(test[,1], type='l', main="Residuals of ARLS", ylab="Residual", xlab="Iteration")
plot(test[,2], type='l', main="Wind Speed parameter over time", ylab="Wind Speed Slope", xlab="Iteration")
par(mfrow=c(1,1))



plot(X$t[7501:8501],X$p[7501:8501], col='1', type='l', ylim=c(0,30), main="One step predictions of wind power", ylab="p (kW)", xlab="Time")
lines(X$t[7501:8501], pred1[7501:8501], col="red")

lower <-(pred1-sqrt(test[,5])*2)[7501:8501]
lower[which(lower<0)]=0

lines(X$t[7501:8501], (pred1+sqrt(test[,5])*2)[7501:8501], col="blue", lty=2)
lines(X$t[7501:8501], lower, col="blue", lty=2)
legend("topright", c("Fitted Values", "95% Confidence intervals", "Observation"), col=c("red","blue", "1"), lty=c(1,2,1,2))


#1-step:
results(Xtest,pred1[7501:8501], sqrt(test[7501:8501,5]))




opt<-optimise(obj.fun, interval = c(0.01,1), ynew=as.matrix(X$p), xnew=as.matrix(X$Ws2))
opt
test<-arls(as.matrix(X$p), as.matrix(X$Ws2), opt$minimum)
pred1 <- c()
for(i in 1:(length(test[,2])-1)){
  pred1[i] <-test[i,2]*X$Ws2[i+1]
}
pred1=c(0,pred1)
results(Xtest,pred1[7501:8501], sqrt(test[7501:8501,5]))


opt<-optimise(obj.fun, interval = c(0.01,1), ynew=as.matrix(X$p), xnew=as.matrix(X$Ws3))
opt
test<-arls(as.matrix(X$p), as.matrix(X$Ws3), opt$minimum)
pred1 <- c()
for(i in 1:(length(test[,2])-1)){
  pred1[i] <-test[i,2]*X$Ws3[i+1]
}
pred1=c(0,pred1)
results(Xtest,pred1[7501:8501], sqrt(test[7501:8501,5]))

