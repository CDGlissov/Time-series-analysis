rm(list=ls())
set.seed(200)
## Use rgl
library(rgl)
par(mfrow=c(1,1))
par(mgp=c(2,0.8,0), mar=c(3,3,2,1),oma=c(0,0,0,0),las=0, pty="m", xpd=F)
##------------------------------------------------


##------------------------------------------------
## Make some data
## Number of samplepoints

setar<-function(n){
x <- seq(-n/2+1,n/2)
## Errors
r <- rnorm(n)
## Make a time series y with a regime model
y <- rep(NA,n)
y[1] <- r[1]
for(t in 2:n)
{
  if(y[t-1] < 0)
  {
    y[t] <-4 + 0.5 * y[t-1] + r[t]
  }
  else
  {
    y[t] <- -4 - 0.5 * y[t-1] + r[t]
  }
}
## Put it into a data.frame, and make x1 and y1 which are lagged one step
return(data.frame(y=y[-1], x1=x[-n], y1=y[-n]))
}

loss <- function(theta,D1, n){
  r1 <- theta[1]
  r2 <- theta[2]
  r3 <- theta[3]
  r4 <- theta[4]
  RSSE <- 0
  for(i in 1:(n-1))
    if(D1$y1[i] < 0)
    {
      RSSE <- RSSE + (D1$y[i]-(r3+r1*D1$y1[i]))^2
    }
  else
  {
    RSSE <- RSSE + (D1$y[i]-(r4+r2*D1$y1[i]))^2
  }
  
  return (RSSE)
}
D<-setar(1000)
opt<-optim(par=c(1,1,1,1), loss, D1=D, n=1000)
opt

opt2 <- matrix(0, 9, 4)
k=1
for(i in c(10,20, 40,80, 160, 320, 640, 1280, 1280*2)){
  D<-setar(i)
  opt2[k,]<-optim(par=c(1,1,1,1), loss, D1=D, n=i)$par
  k=k+1
}



matplot(abs(sweep(opt2, 2, c(0.5,-0.5,4,-4))),type='l', lwd=2, xaxt='n', ylab="Parameter Deviation", xlab="Number of Sample Points", main="Convergence of parameter estimates")
axis(1, 1:9,c(10,20, 40,80, 160, 320, 640, 1280, 1280*2) )
legend("topright", c(expression(a[0]^{(1)}), expression(a[0]^{(2)}),expression(a[1]^{(1)}),expression(a[1]^{(2)})), text.col=c(3,4,1,2))



plot(D$y1, D$y,ylab = "X[t]", xlab="X[t-1]",
     main="Simulation using 1000 sample points of SETAR(2,1,1)")
lines(x=seq(min(D$y),0, 0.01), y=4+0.5*seq(min(D$y),(0), 0.01), col=2,lwd=2)
lines(x=seq(0, max(D$y), 0.01), y=-4-0.5*seq(0, max(D$y), 0.01), col=2,lwd=2)
lines(x=seq(min(D$y),0, 0.01), y=3.9604616+0.5140395*seq(min(D$y),(0), 0.01), col=3,lwd=2)
lines(x=seq(0, max(D$y), 0.01), y=-3.8778064-0.5177555*seq(0, max(D$y), 0.01), col=3,lwd=2)
legend("topright", c("Theoretical Mean", "Estimated Mean"), text.col=c(2,3))


#####Opgave 2######
n<-3000
set.seed(200)
D<-setar(n)
loss2 <- function(p1,p2,N1,N2,D){
  RSSE <- 0
  dat <- D$y[N1:N2]
  dat1 <- D$y1[N1:N2]
  for(i in 1:(length(dat)-1))
    if(dat1[i] < 0)
    {
      RSSE <- RSSE + (dat[i]-(4+p1*dat1[i]))^2
    }
  else
  {
    RSSE <- RSSE + (dat[i]-(-4+p2*dat1[i]))^2
  }
  
  return (RSSE)
}

n1 = 200
p1 <- seq(-1, 1,length.out=n1)
p2<-seq(-1, 1, length.out=n1)

l2<-matrix(0,n1,n1)

subset1<-c(1, 1,1,1001,1001,3000,300,30,1300,1030)
dim(subset1)=c(5,2)

for(j in 1:5){
  N1<-subset1[j,1]
  N2<-subset1[j,2]
  for(i in 1:n1){
    for(k in 1:n1){
      l2[i,k]<-loss2(p1[i],p2[k], N1, N2,D)
    }
  }
  m1<-which(l2==min(l2), arr.ind=T)
  par(mgp=c(2,0.8,0), mar=c(3,4,2,1),oma=c(0,0,0,0),las=0, pty="m", xpd=F)
  filled.contour(p1,p2,l2, 
                 color.palette=colorRampPalette(c("#00007F","blue","#007FFF",
                                                  "cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
                 ,plot.axes = { axis(1); axis(2); points(0.5, -0.5, col=2, lwd=1); 
                   points(p1[m1[1]], p2[m1[2]], col=3, lwd=1);})
  title(xlab=expression(a[1]^{(1)}), ylab=expression(a[1]^{(2)}), line=2)
  title(paste("Contour curves of loss function, interval from",N1, "to",N2))
  legend(x=0.12, y=1, c("True Parameter", "Estimated Parameter"), pch=1,col=c(2,3))
}


################## OPGAVE 3
n <- 500
mu<-0.8
ry <- rnorm(n,sd=4)
rphi <- rnorm(n, sd=0.1)
phi <- rep(NA,n)
yt <- rep(NA,n)
theta<-0.4
yt[1] <- ry[1]
phi[1] <- rphi[1]
delta <- mu*(1-theta)
for(t in 2:n){
  phi[t]=theta*phi[t-1] +delta+rphi[t]
  delta=delta
  yt[t]=phi[t]*yt[t-1]+ ry[t]
}
plot(phi, type='l')
plot(yt, type='l')


#AR(2)-ARMA(2)
n <- 500
mu<-0.5
reps <- rnorm(n,sd=0.4)
rpsi <- rnorm(n, sd=0.1)
phi <- rep(NA,n)
yt <- rep(NA,n)
phi1<- 0.85
phi2<- -0.6
delta <- mu*(1-phi1 - phi2)

yt[1] <- reps[1]
yt[2] <- reps[2]
phi[1] <- rpsi[1]
phi[2]<-rpsi[2]
theta <- 0.4
for(t in 3:n){
  phi[t]=phi1*phi[t-1] + phi2*phi[t-2]+delta+rpsi[t]
  phi[t-1] = phi[t-1] + theta*rpsi[t]
  delta=delta
  yt[t]=phi[t]*yt[t-1]+phi[t-1]*yt[t-2] + reps[t]
}
plot(phi, type='l', xlab = "Time", ylab=expression(Phi[t]), main="Simulation of Phi[t]")
plot(yt, type='l', xlab = "Time", ylab=expression(y[t]), main="Simulation of y[t]")

par(mfrow=c(2,1))
acf(yt, main="")
title("ACF of y[t]")
pacf(yt, main="")
title("PACF of y[t]")
acf(phi, main="")
title("ACF of Phi[t]")
pacf(phi, main="")
title("PACF of Phi[t]")

##### OPGAVE 4

n=20
a=0.4
ev <- rnorm(n,sd=1)
ee <- rnorm(n,sd=1)
xt<- rep(NA,n)
yt <-rep(NA,n)
xt[1]<-ev[1]
yt[1] <- ee[1]
for(t in 2:n){
  xt[t] <- a*xt[t-1]+ev[t-1]
  yt[t] <- xt[t-1]+ee[t-1]
}

par(mfrow=c(2,1))
plot(xt, type='l')
plot(yt, type='l')
par(mfrow=c(1,1))



