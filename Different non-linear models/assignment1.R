rm(list=ls())
par(mgp=c(2,0.8,0), mar=c(3,3,2,1),oma=c(0,0,0,0),las=0, pty="m", xpd=F)
par(mfrow=c(1,1))
setwd("C:/Users/Christian/Dropbox/7. Semester/Videregående Tidsrække/aflevering 1")

## Epanechnikov kernel for later use
kernelEp <- function(xall,x,h)
{
  ## Make the weights with an Epanechnikov kernel
  ## h has the same unit as x 
  #(i.e. it is on the same absolute scale, so if x is Watt, h is also given in Watt) 
  u <- abs(xall-x)
  u <- u / h
  w <- 3/(4) * (1 - u^2)
  ## Set values with |u|>1 to 0
  w[abs(u)>1] <- 0
  return(w)
}

####Q1.1#####

##Using seed(1) throughout the script
set.seed(1)
## Number of samplepoints
n <- 300
## Errors
r <- rnorm(n)
## Make a time series y with a regime model
y <- rep(NA,n)
y[1] <- r[1]
#######################SETAR(2,1,1)################
for(t in 2:n){
  if(y[t-1] < 0)
  {
    y[t] <- 4+0.5 * y[t-1] + r[t]
  }
  else
  {
    y[t] <- -4-0.5 * y[t-1] + r[t]
  }
}

## Plot it
plot(y[-1], type='l', xlab="t", ylab="x[t]", 
     main="Time series of SETAR(2,1,1)")
plot(y=y[-1], x=y[-n], ylab = "X[t]", xlab="X[t-1]",
     main="Simulation using 300 sample points of SETAR(2,1,1)")
lines(x=rep(0,151), y=seq(-50, 100), lty=2)
legend("topright", "Regime Shift", lty=2)

###################IGAR(2,1)#######################
set.seed(1)
n <- 300
r <- rnorm(n)
y <- rep(NA,n)
y[1] <- r[1]
x <- runif(n,0,1)

for(t in 2:n){
  if(x[t-1] < 0.5)
  {
    y[t] <- 4 + 0.5 * y[t-1] + r[t]
  }
  else
  {
    y[t] <- -4 - 0.5 * y[t-1] + r[t]
  }
}

## Plot it

plot(y[-n], y[-1],ylab = "X[t]", xlab="X[t-1]", 
     main="Simulation using 300 sample points of IGAR(2,1)")
plot(y[-1], type='l', xlab="t", ylab="x[t]", 
     main="Time series of IGAR(2,1)")

########################MMAR(2,1)######################
set.seed(1)
#set current state
state = 1
states = c()
n <- 400
r <- rnorm(n)
y <- rep(NA,n)
y[1] <- r[1]
P<-matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2)


for(t in 2:n){
  if(state == -1){
    state == -1
    if(runif(1,0,1) > P[1,1]){
      state = 1
    }
  } else if(state == 1){
    state=1
    if(runif(1,0,1) > P[2,2]){
      state = -1
    }
  }
  states <- c(states,state)
  if(state == 1)
  {
    y[t] <- 4 + 0.5 * y[t-1] + r[t]
  }
  else
  {
    y[t] <- -4 - 0.5 * y[t-1] + r[t]
  }
}


## Plot simulations and states
plot(y[-n], y[-1],,ylab = "X[t]", xlab="X[t-1]", 
     main="Simulation using 400 sample points of MMAR(2,1)")
par(mfrow=c(2,1))
plot(y=states, x=1:399, type='l', xlab="t", 
     ylab="Current State", main="States Over Time", ylim=c(-2,2))
plot(y[-1], type='l', xlab="t", ylab="x[t]", 
     main="Time series of MMAR(2,1)")
par(mfrow=c(1,1))


######1.2########
set.seed(1)
n <- 1000
r <- rnorm(n)
y <- rep(NA,n)
y[1] <- r[1]

##Using setar(2,1,1) model
for(t in 2:n)
{
  if(y[t-1] < 0)
  {
    y[t] <- 4 + 0.5 * y[t-1] + r[t]
  }
  else
  {
    y[t] <- -4 - 0.5 * y[t-1] + r[t]
  }
}

##Using kernelep and estimating the conditional mean (2.34)
m1<-function(x, h){
  k<-kernelEp(D$y, x, h)
  
  tsum = 0
  N<-length(D$y)
  for( i in 2:N){
    tsum = tsum + D$y[i] * k[i-1]
  }
  
  bsum = 0
  for (i in 2:(N+1)){
    bsum = bsum + k[i-1]
  }
  
  return((1/(N-1) * tsum) / (1/N * bsum))
}

#plot data
D <- data.frame(y=y[-1], y1=y[-n])
plot(D$y1, D$y, xlab="X[t-1]", ylab="X[t]", 
     main="Different bandwidth (h) using Epanechnikov Kernel")

#Estimate mean
mean1<-c()
lin<-seq(min(D$y),max(D$y), length.out=n)
spans <- c(0.1,0.5,1,1.5,2,3)
for(k in 1:length(spans)){
  for(i in 2:n){
    mean1[i-1]<-m1(lin[i-1], spans[k])
  }
  lines(lin[2:1000],mean1, col=k+1, lwd=1.5)
}
lines(seq(min(D$y),0, length.out = 50),4 + 0.5 * 
        seq(min(D$y),0, length.out = 50), lwd=2, col=8)
lines(seq(0,max(D$y), length.out = 50),-4 - 0.5 * 
        seq(0,max(D$y), length.out = 50), lwd=2, col=8)
legend("topright",lty=1, c(paste0(rep("h=",5),spans), 
                           "Theoretical Cond. Mean"),col=c(2,3,4,5,6,7,8))

##Finding optimal bandwidth
leaveOneOut <- function(D, span)
{
  ## Find the bandwidth giving the best balance between bias and variance
  ## of the model by leave one out.
  ## Matrix for keeping the residuals
  R <- matrix(nrow=nrow(D), ncol=length(span))
  ## Do it for each bandwidth
  for(ii in 1:length(span))
  {
    print(paste("  Fitting for bandwidth",ii,"of",length(span)))
    ## Do the local 2-order polynomial regression one time for each point
    ## leaving the point out while fitting, and then predicting the point.
    for(i in 1:nrow(D))
    {
      wx <- kernelEp(D$y, D$y1[i], h=span[ii])
      w <- wx
      ok <- w>0
      #R[i,ii] <- D[i,"y"] - predict(loess(y ~ y1, dat=D[-i,], span=span[ii]), D[i,])
      R[i,ii] <- D[i,"y"] - predict(lm(y ~ y1, 
                                       weights=w[ok], data=D[-i,][ok,]), D[i,])
      
      }
  }
  ## Find the best bandwidth
  RSSkAll <- apply(R, 2, function(x){sum(x^2,na.rm=TRUE)})
  return(RSSkAll)
}
spans <- c(seq(0.4,1.5, length.out=100))
rss.test<-leaveOneOut(D,spans)
spanBest <- spans[which.min(rss.test)]

# Plot CV
plot(y=rss.test,x=spans, xlab="Bandwidth, h", ylab="RSS", 
     main="Cross-validation of the bandwidth, h" )
points(y=min(rss.test),x=spanBest, col="red", lwd=2, pch=16)

#Find mean with optimal bandwidth
mean1<-c()
lin<-seq(min(D$y),max(D$y), length.out=n)
for(i in 2:n){
  mean1[i-1]<-m1(lin[i-1], spanBest)
}

# find conditional variance (2.36)
v1<-function(x, h){
  k<-kernelEp(D$y, x, h)
  
  tsum = 0
  N<-length(D$y)
  for( i in 2:N){
    tsum = tsum + (D$y[i]^2) * k[i-1]
  }
  
  bsum = 0
  for (i in 2:(N+1)){
    bsum = bsum + k[i-1]
  }
  
  return((1/(N-1) * tsum) / (1/N * bsum)-m1(x,h)^2)
}

var1<-c()
lin<-seq(min(D$y),max(D$y), length.out=n)
for(i in 2:n){
  var1[i-1]<-v1(lin[i-1], spanBest)
}

#Plot
plot(D$y1, D$y, ylim=c(-10,10), xlab="X[t-1]", ylab="X[t]", 
     main=paste0("Epanechnikov smoothing with optimal bandwidth, h="
                 ,round(spanBest,digits=3)))
lines(seq(min(D$y),0, length.out = 50),4 + 0.5 * 
        seq(min(D$y),0, length.out = 50), lwd=2, col=8)
lines(seq(0,max(D$y), length.out = 50),-4 - 0.5 * 
        seq(0,max(D$y), length.out = 50), lwd=2, col=8)
lines(lin[2:1000],mean1, col=2, lwd=2)
lines(lin[2:1000],mean1+2*sqrt(var1), col=2, lty=2,lwd=2)
lines(lin[2:1000],mean1-2*sqrt(var1), col=2, lty=2,lwd=2)
legend("topright", col=c(2,2,8), lty=c(1,2,1), c("Estimation",
                                     "95% confidence interval", "Theoretical Cond. Mean"), lwd=2 )

#######1.3########
#calculate cumulative mean theoretical and estimate, do it for 4 different number of sample points
par(mfrow=c(2,2))
for(k in 1:4){
set.seed(1)
choose.n = c(500,1000,3000,12000)
n <- choose.n[k]
r <- rnorm(n)
y <- rep(NA,n)
y[1] <- r[1]

for(t in 2:n)
{
  if(y[t-1] < 0)
  {
    y[t] <- 4 + 0.5 * y[t-1] + r[t]
  }
  else
  {
    y[t] <- -4 - 0.5 * y[t-1] + r[t]
  }
}

D <- data.frame(y=y[-1], y1=y[-n])


## Script for calculation of a confidence bands using the cumulative means technique.
x<-D$y
## Parameters for the histogram regression
## Number of intervals 
n.bin <- 10
## The breaks between the intervals 
breaks <- seq(-5,5,len=n.bin+1)
## Initialize
h <- diff(breaks)[1]
lambda <- gamma <- f.hat <- h.hat <- numeric(n.bin)
##----------------------------------------------------------------


##----------------------------------------------------------------
## Cut into intervals conditioned on x_{t-1}
L <- split(x[-1], cut(x[-length(x)],breaks))
## Check if there are at least 5 points in each interval
if(!all(sapply(L,length)>=5)){ 
  print('Stopped: There are less than 5 points in one of the intervals'); break;}
## Calc the hist regressogram, i.e. for each interval
for(i in 1:n.bin){
  x.bin <- L[[i]]
  lambda[i] <- mean(x.bin)
  f.hat[i] <- (n.bin*h)^(-1) * length(x.bin)
  gamma[i] <- sum((x.bin - lambda[i])^2) / length(x.bin)
}
## Make confidence bands for the cumulated function. Def. (3.10).
## 95% confidence band, c is found in table 3.1
c.alpha <- 1.273
##
Lambda <- cumsum(lambda*h)
for(i in 1:n.bin)
{
  h.hat[i] <- gamma[i]/f.hat[i];
}
H.hat <- cumsum(h.hat*h);
##
H.hat.b <- H.hat[n.bin];
Lambda.lower <- Lambda - c.alpha * n.bin^(-0.5) * H.hat.b^(0.5) * (1 + H.hat/H.hat.b);
Lambda.upper <- Lambda + c.alpha * n.bin^(-0.5) * H.hat.b^(0.5) * (1 + H.hat/H.hat.b);
#############################################################
integral1<-function(a,b){
  return(b^2*0.25+4*b-(a^2*0.25+4*a))
}

integral2<-function(a,b){
  return(-b^2*0.25-4*b+a^2*0.25+4*a)
}


plot(y=Lambda, x=breaks[-1],lty=2, ylim=c(-15,20), xlim=c(-4.2,5), 
     main="", xlab="Interval (t)", ylab="Estimated Cumulative Mean (Lambda)")
if(k==1){
  title("Cumulative Conditional Mean",outer=T, line=-1.5)
}
lines(x=breaks[-1], Lambda.lower, lty=2)
lines(x=breaks[-1], Lambda.upper, lty=2)
tx<-seq(min(breaks[-1]),0, 0.01)
tx1<-seq(0,max(breaks[-1]), 0.01)
lines(x=tx, y=integral1(-5,tx), col=2)
lines(x=tx1, y=integral2(0,tx1)+integral1(-5,0), col=2)
legend("bottomleft", col=c(1,2,1,1), lty=c(NA,1,2,NA),pch=c(1,NA,NA,NA), 
       c("Estimated", "Theoretical","95% Confidence Bands",
         paste0("n=", choose.n[k])))
}
par(mfrow=c(1,1))


#####Q1.4####
library(msir)

dat<-read.csv("DataPart4.CSV", header=T, skip=0, sep=",", 
              colClasses = c(rep("numeric",4)))

#phi is heat load
#Ua is heat loss coef
#Ti is internal temp
#Te is external temp 
#W is the wind speed

dat$Ua<-dat$Ph/(dat$Ti-dat$Te)
plot(dat$W,dat$Ua, ylab="Heat Loss, Ua", xlab="Wind Speed, 
     W", main="Plot of the Heat loss vs Wind speed")


weights<-abs((dat$Ti-dat$Te))/max(dat$Ti-dat$Te)

leaveOneOut2 <- function(D,span)
{
  ## Find the bandwidth giving the best balance between bias and variance
  ## of the model by leave one out.
 
  ## Matrix for keeping the residuals
  R <- matrix(nrow=nrow(D), ncol=length(span))
  ## Do it for each bandwidth
  for(ii in 1:length(span))
  {
    print(paste("  Fitting for bandwidth",ii,"of",length(span)))
    ## Do the local 2-order polynomial regression one time for each point
    ## leaving the point out while fitting, and then predicting the point.
    for(i in 1:nrow(D))
    {
      R[i,ii] <-as.numeric(D[i,"x"]-predict(loess(x ~ xk, 
                                                  dat=D[-i,], 
                                                  span=span[ii], 
                                                  weights = weights[-i]), 
                                            D[i,]))
    }
  }
  ## Find the best bandwidth
  RSSkAll <- apply(R, 2, function(x){sum(x^2,na.rm=TRUE)})
  
  return(RSSkAll)
}

D<-dat
D$x<-D$Ua
D$xk<- D$W

spans2 <- c(seq(0.1,2, length.out=100))
rss.test2<-leaveOneOut2(D, spans2)
spanBest2 <- spans2[which.min(rss.test2)]

#Plot CV
plot(y=rss.test2,x=spans2, xlab="Bandwidth, h", 
     ylab="RSS", main="Cross-validation of the bandwidth, h" )
points(y=min(rss.test2),x=spanBest2, col="red", 
       lwd=2, pch=16)


#Prepare data for fit
U<-dat$Ua[order(dat$W)]
weights2 <- weights[order(dat$W)]
W<-sort(dat$W)

#Make LOESS fit
fit1<-loess(U~W, span=spanBest2, weights=weights2)
pred1<-predict(fit1, se=T)
fit1.sd<-loess.sd(U~W, nsigma=1.96,span=spanBest2, 
                  weights=weights2)

#Plot it
plot(dat$W,dat$Ua, type='p',col=1, ylab="Heat Loss Coefficient (Ua)", 
     xlab="Wind speed (W)"
     , main="LOWESS of Ua(W)",ylim=c(140,240))
lines(fit1.sd$x,fit1.sd$y,col=2, lwd=2)
lines(fit1.sd$x, fit1.sd$upper, col=2, lty=3,lwd=2)
lines(fit1.sd$x, fit1.sd$lower, col=2, lty=3,lwd=2)
lines(W, pred1$fit +qt(0.975,pred1$df)*pred1$se.fit, col=2, lty=2,lwd=2)
lines(W, pred1$fit -qt(0.975,pred1$df)*pred1$se.fit, col=2, lty=2,lwd=2)
legend("topleft", c("Fitted Values","95% Prediction Interval", "95% Confidence Interval", 
                    "Observations"), 
       lty=c(1,3,2,NA), pch=c(NA,NA,NA,1), col=c(2,2,2,1), lwd=c(2,2,2,1))


######Q1.5#######
dat<-read.csv("DataPart5.CSV", header=T, skip=0, sep=",", 
              colClasses = c(rep("numeric",1)))
plot(dat$x)

par(mfrow=c(1,2))
acf(dat$x, main="", lag.max=80)
title("ACF")
pacf(dat$x, main="", lag.max=80)
title("PACF")
par(mfrow=c(1,1))

model1<-arima(dat$x, order=c(0,0,2), fixed=c(0,NA,NA))
par(mfrow=c(1,2))
acf(model1$residuals, main="", lag.max=80)
title("ACF of ARIMA(0,0,2)")
pacf(model1$residuals, main="", lag.max=80)
title("PACF of ARIMA(0,0,2)")
par(mfrow=c(1,1))
#1 borderline significant acf, however seems fine!

plot(model1$residuals[c(-799,-800)],model1$residuals[c(-1, -2)], 
     xlab="e[t-2]", ylab="e[t]", main="Residuals vs residuals plot")

#Plot LDF
lagsldf<-ldf(model1$residuals, lags=seq(1,30))


series = data.frame(x=dat$x[c(-1,-2)], xk=dat$x[c(-799,-800)])

#Make Regressions
plot(series$x, x=series$xk, main="Linear regression fitted to each regime.", 
     xlab="X[t-2]", ylab="X[t]")
region1<-series[series$xk< (-2),]
model1<-lm(region1$x~region1$xk)
arima(region1$x, order=c(0,0,2), fixed=c(0,NA,NA))

region2<-series[(series$xk > (-2) & series$xk < (2)),]
model2<-lm(region2$x~region2$xk)
region3<-series[(series$xk < (5) & series$xk > (2)),]
model3<-lm(region3$x~region3$xk)
region4<-series[(series$xk > (5) & series$xk <= (8)),]
model4<-lm(region4$x~region4$xk)

lines(x=seq((-6),(-2)), y=2.281-0.68*seq((-6),(-2)), col=2,lwd=2)
lines(x=seq((-2),(2)), y=0.9009+2.7342*seq((-2),(2)), col=2, lwd=2)
lines(x=seq((2),(5)), y=-1.002+0.804*seq((2),(5)), col=2, lwd=2)
lines(x=seq((5),(8)), y=-5.996+1.02*seq((5),(8)), col=2, lwd=2)
legend("topright", c("Fitted Values", "Observation"), pch=c(NA,1), lty=c(1,NA), col=c(2,1))

#Fit LOESS

par(mfrow=c(1,1))
test<-series$x[order(series$xk)]
test2<-sort(series$xk)
fit = loess(test ~ test2, span = 0.11)
pred2<-predict(fit,newdata=data.frame(x=seq(1,10,length.out=798)), se=T)
fit = loess.sd(test ~ test2, span = 0.11, nsigma=1.96)

#Plot
plot(y=series$x, x=series$xk, ylim=c(-6,9), ylab="X[t]", xlab="X[t-2]",
     main = "LOESS fitted to observations, h=0.1")
lines(test2,pred2$fit, col=2, lwd=2)

lines(test2,pred2$fit+qt(0.975, pred2$df)*pred2$se, col=2,lty=2, lwd=2)
lines(test2,pred2$fit-qt(0.975, pred2$df)*pred2$se, col=2,lty=2, lwd=2)
lines(test2,fit$upper, col=2, lty=3, lwd=2)
lines(test2,fit$lower, col=2, lty=3, lwd=2)
legend("topright", c("Fitted Values","95% Prediction Interval", "95% Confidence Interval", 
                    "Observations"), 
       lty=c(1,3,2,NA), pch=c(NA,NA,NA,1), col=c(2,2,2,1), lwd=c(2,2,2,1))

#Look at residuals of LOESS
fit2<-loess(series$x ~ series$xk, span = 0.11)
par(mfrow=c(1,2))
acf(fit2$residuals, main="", lag.max=80)
title("ACF of LOESS")
pacf(fit2$residuals, main="", lag.max=80)
title("PACF of LOESS")
par(mfrow=c(1,1))
ldf(fit$residuals, lags=seq(1,30))

