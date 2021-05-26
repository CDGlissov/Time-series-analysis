##----------------------------------------------------------------
## EKF algorithm for use in Part 4 of computer exercise 2 in
## Advanced Time Series Analysis
##----------------------------------------------------------------
par(mgp=c(2,0.8,0), mar=c(3,3,2,1),oma=c(0,0,0,0),las=0, pty="m", xpd=F)
##----------------------------------------------------------------
## Do the simulation here and keep the values in y
##----------------------------------------------------------------
zn<-20
n=3000
aZ<-vector(mode="list", zn)
for(i in 1:zn){
  time<-seq(1:n)
  a=0.4
  ev <- rnorm(n,sd=1)
  ee <- rnorm(n,sd=1)
  xt<- rep(NA,n)
  yt <-rep(NA,n)
  xt[1]<-ev[1]
  yt[1] <- ee[1]
  for(t in 2:n){
    xt[t] <- a*xt[t-1]+ev[t-1]
    yt[t] <- xt[t]+ee[t]
  }
  
  
  y=yt
  ##----------------------------------------------------------------
  ## Estimation with the EKF
  ##----------------------------------------------------------------
  ## aInit : The starting guess of the AR coefficient estimate
  aInit <- 0.5
  ## aVarInit : The initial variance for estimation of the AR coefficient
  aVarInit <- 1
  ## sigma.v : Standard deviation of the system noise of x in the filter
  sigma.v <- sqrt(10)
  
  ## Initialize
  ## Init the state vector estimate
  zt <- c(0,aInit)
  ## Init the variance matrices
  Rv <- matrix(c(sigma.v^2,0,0,0), ncol=2)
  ## sigma.e : Standard deviation of the measurement noise in the filter
  Re <- sqrt(1) 
  
  ## Init the P matrix, that is the estimate of the state variance
  Pt <- matrix(c(Re,0,0,aVarInit), nrow=2, ncol=2)
  ## The state is [X a] so the differentiated observation function is
  Ht <- t(c(1,0))
  ## Init a vector for keeping the parameter a variance estimates
  aVar <- rep(NA,length(y))
  ## and keeping the states
  Z <- matrix(NA, nrow=length(y), ncol=2)
  Z[,1] <- zt
  
  ## The Kalman filtering
  for(t in 1:(length(y)-1))
  {
    ## Derivatives (Jacobians)
    Ft <- matrix(c(zt[2],0,zt[1],1), ncol=2)  # F_t-1
    # Ht does not change 
    
    ## Prediction step
    zt = c(zt[2]*zt[1],zt[2]) #z_t|t-1 f(z_t-1|t-1)
    Pt = Ft %*% Pt %*% t(Ft) + Rv #P_t|t-1
    
    ## Update step
    res = y[t] - zt[1] # the residual at time t
    St =  Ht %*% Pt %*% t(Ht) + Re # innovation covariance
    Kt = Pt %*% t(Ht) %*% St^-1 # Kalman gain
    zt = zt + Kt * res # z_t|t
    Pt = (diag(2) - Kt%*%Ht)%*%Pt #P_t|t
    
    ## Keep the state estimate
    Z[t+1,] <- zt
    ## Keep the P[2,2], which is the variance of the estimate of a
    aVar[t+1] <- Pt[2,2]
    
  }
  aZ[[i]] <- Z[,2]
}
test<-matrix(unlist(aZ), ncol = 20, byrow=F)


plot(test[,1],type='l', ylim=c(min(test, na.rm=T),max(test, na.rm=T)), 
     main=paste(20, "simulations using EKF:", paste0("a[init]=", aInit, ","),paste0("Var[a[init]]=", aVarInit,","),
                paste0("sigma_v^2=", sigma.v^2)), xlab="Observations",
     ylab="Estimated parameter, a")
for(i in 2:zn){
  lines(test[,i], col=1)
}
abline(h=a, lwd=2, col=2)
#abline(h=mean(rowMeans(simplify2array(aZ)), na.rm=T), lwd=2, col=4)
lines(1:3000, rowMeans(simplify2array(aZ)), col=3, lwd=2)
legend("topright", c("True parameter", "Mean over time of estimated parameter"), text.col=c(2,3))
