setwd("C:/Users/Jonas/SkyDrive/Dokumenter/DTU/02427 Advanceret Tidsrække/Excercise 4")
dat <- read.table("data/cex4WindDataInterpolated.csv", sep=",",
                  header=TRUE, stringsAsFactors=FALSE)
dat$t <- as.POSIXct(dat$t, tz="UTC")

dat = dat[8500:17000,]

X_test = dat[7501:8501,]
X = dat[1:7500,]

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
  CRPS<-mean(crps,na.rm = 1)
  out = c(RMSE,MAE,CRPS)
  return(out)
}



##################################EPANECHNIKOV KERNEL#############

kernelEp <- function(xall,x,h)
{
  ## Make the weights with an Epanechnikov kernel
  ## h has the same unit as x (i.e. it is on the same absolute scale, so if x is Watt, h is also given in Watt)
  u <- abs(xall-x)
  u <- u / h
  w <- 3/4 * (1 - u^2)
  ## Set values with |u|>1 to 0
  w[abs(u)>1] <- 0
  return(w)
}

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
      wx <- kernelEp(D$Ws1, D$Ws1[i], h=span[ii])
      w <- wx
      ok <- w>0
      R[i,ii] <- D[i,"p"] - predict(loess(p ~ Ws1, dat=D[-i,], span=span[ii]), D[i,])
      #R[i,ii] <- D[i,"p"] - predict(lm(p ~ Ws1, 
      #                                 weights=w[ok], data=D[-i,][ok,]), D[i,])
      
    }
  }
  ## Find the best bandwidth
  RSSkAll <- apply(R, 2, function(x){sum(x^2,na.rm=TRUE)})
  return(RSSkAll)
}
spans <- c(seq(0.1,1.5, length.out=50))
rss.test<-leaveOneOut(X[600:1800,3:4],spans)
bestspan<-spans[which.min(rss.test)]
#bestspan<-0.3






epi_kern_pred = rep(NA,length(X_test$p))
epi_kern_se = rep(NA,length(X_test$p))
pred_upper = rep(NA,length(X_test$p))
pred_lower = rep(NA,length(X_test$p))
for (i  in 1:length(X_test$p)) {
  wy <- kernelEp(X$Ws3, X_test$Ws3[i], h=bestspan)
  w <- wy
  ## Do it only with positive weights
  ok <- w>0
  ## Note that this is local first order polynomial regression, but can easily be made 2'nd
  fit <- lm(p ~ Ws3, weights=w[ok], data=X[ok,])
  epi_pred = predict(fit, data.frame(Ws3=X_test$Ws3[i]),se.fit = TRUE,interval = "prediction")
  epi_kern_pred[i] = epi_pred$fit[1]
  epi_kern_se[i] = epi_pred$se.fit
  pred_upper[i] = epi_pred$fit[3]
  pred_lower[i] = epi_pred$fit[2]
}




epi_kern_fit = rep(NA,length(X$p))
epi_kern_se_fit = rep(NA,length(X$p))
conf_upper = rep(NA,length(X$p))
conf_lower = rep(NA,length(X$p))
for (i  in 1:length(X$p)) {
  wy <- kernelEp(X$Ws1, X$Ws1[i], h=bestspan)
  w <- wy
  ## Do it only with positive weights
  ok <- w>0
  ## Note that this is local first order polynomial regression, but can easily be made 2'nd
  fit <- lm(p ~ Ws1, weights=w[ok], data=X[ok,])
  epi_pred = predict(fit, data.frame(Ws1=X$Ws1[i]),se.fit = TRUE,interval = "confidence")
  epi_kern_fit[i] = epi_pred$fit[1]
  epi_kern_se_fit[i] = epi_pred$se.fit
  conf_upper[i] = epi_pred$fit[3]
  conf_lower[i] = epi_pred$fit[2]
}


plot(X_test$p~X_test$Ws1,xlab = "Windspeed 1 Hour prediction", ylab = "Power Production")
points(epi_kern_pred~X_test$Ws1, col = 2)
legend("bottomright",legend = c("Test data","Predicted"),col=c("black","red"),pch=c(1,1))

plot(X$p~X$Ws1,xlab = "Windspeed 1 Hour prediction", ylab = "Power Production")
points(epi_kern_fit~X$Ws1, col = 2)
legend("bottomright",legend = c("Model data","fitted"),col=c("black","red"),pch=c(1,1))


results(X_test,epi_kern_pred,epi_kern_se)

plot(spans,rss.test,xlab = "Span",ylab = "RSS")

par(mfrow=c(1,1))
res = epi_kern_pred-X_test$p
plot(X_test$t, res,xlab="Time",ylab="Residuals",type="l")




plot(dat$t[7400:8501],dat$p[7400:8501],type="l",xlab="Time",ylab="Power",ylim=c(0,30))

lines(X_test$t,epi_kern_pred,col=2)
lines(X_test$t,pred_lower,col=2,lty=3)
lines(X_test$t,pred_upper,col=2,lty=3)

lines(X$t[7200:7500],epi_kern_fit[7200:7500],col=4)
lines(X$t[7200:7500],conf_lower[7200:7500],col=4,lty=3)
lines(X$t[7200:7500],conf_upper[7200:7500],col=4,lty=3)

legend("topleft",legend = c("Observations","Fit","Confidence Interval","Prediction","Prediction Interval"),col = c(1,4,4,2,2),lty=c(1,1,3,1,3),cex = 0.8)
