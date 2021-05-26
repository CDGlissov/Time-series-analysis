rm(list=ls())
set.seed(200)
## Use rgl
library(ctsmr)
library(rgl)
library(plot3D)
par(mfrow=c(1,1))
par(mgp=c(2,0.8,0), mar=c(3,3,2,1),oma=c(0,0,0,0),las=0, pty="m", xpd=F)
##------------------------------------------------

Delta = 2^(-9)
theta1 = 0.7
theta2 = 0.8
theta3 = 3
theta4 = -0.34
n=2500


sim1<-function(sigma){
  DW <- rep(0,n+1)
  for (i in 2:length(DW)) {
    DW[i] <-  DW[ i - 1] + rnorm(1, 0, Delta)
  }
  y1 <- rep(NA,n)
  y2 <- rep(NA,n)
  y1[1] = -1.9
  y2[1] = 1.2
  for (t in seq(2, n)){
    y1[t] = y1[t-1]+theta3*(y1[t-1]+y2[t-1]- (1/3)*y1[t-1]^3+theta4)*Delta+sigma*DW[t]
    y2[t] = y2[t-1]- (1/theta3)*(y1[t-1]+theta2*y2[t-1]-theta1)*Delta
  }
  return(data.frame(y1=y1, y2=y2))
}

par(mfrow=c(1,3))
D<-sim1(0)
plot(1:n,D$y1, type='l', xlab="Time", ylab="Y_k^1", main=paste0("Realization of Y_k^1, sigma=", 0))
plot(1:n,D$y2, type='l', xlab="Time", ylab="Y_k^2", main=paste0("Realization of Y_k^2, sigma=", 0))
plot(D$y1,D$y2, xlab="Y_k^1", ylab="Y_k^2", main=paste0("Phaseplot of (Y_k^1, Y_k^2), sigma=", 0))

par(mfrow=c(2,3))
D<-list()
sig1<-c(0.1,0.2,0.3,0.4)
for(i in 1:4){
D[[i]]<-sim1(sig1[i])
plot(1:n,D[[i]]$y1, type='l', xlab="Time", ylab="Y_k^1", main=paste0("Realization of Y_k^1, sigma=", sig1[i]))
plot(1:n,D[[i]]$y2, type='l', xlab="Time", ylab="Y_k^2", main=paste0("Realization of Y_k^2, sigma=", sig1[i]))
plot(D[[i]]$y1,D[[i]]$y2, xlab="Y_k^1", ylab="Y_k^2", main=paste0("Phaseplot of (Y_k^1, Y_k^2), sigma=", sig1[i]))
}
par(mfrow=c(1,1))





#1bpar(mfrow=c(2,3))
sig2<-c(0.1,0.2,0.3,0.4)
for(i in 1:4){
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

library(MASS)
h1 <- hist(D[[i]]$y1, breaks=100, plot=F)
h2 <- hist(D[[i]]$y2, breaks=100, plot=F)
top <- max(h1$counts, h2$counts)
k <- kde2d(D[[i]]$y1, D[[i]]$y2, n=100)

# margins
oldpar <- par()
par(mar=c(3,3,1,1))
layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
image(k, col=r) #plot the image
lines(D[[i]]$y1, D[[i]]$y2, lwd=2)
par(mar=c(0,2,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col='red')
par(mar=c(2,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col='red', horiz=T)
title(paste0("3D Realization of (Y_k^1, Y_k^2), sigma=", sig2[i]), outer=T, line=-1)
}