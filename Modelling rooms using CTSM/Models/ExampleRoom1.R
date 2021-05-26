rm(list=ls())
library(ctsmr)
options(scipen=999)
library(splines)
setwd("C:/Users/Christian/Dropbox/7. Semester/Videregående Tidsrække/Aflevering 3")
source("sdeTiTm.R")
files <- dir("functions", full.names=TRUE)
for(i in 1:length(files)) source(files[i])
load("Exercise3.RData")
par(mgp=c(2,0.8,0), mar=c(3,3,2,1),oma=c(0,0,0,0),las=0, pty="m", xpd=F)
reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                       cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                       col.axis = "black", col.lab = "black", col.main = "black",
                       col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                       fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                       ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                       lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                       lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                         4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                        1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                         0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                               0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                       tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                           1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                       yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                "yaxt", "ylbias"))
  par(op)
}
library(extRemes)

str(AllDat)
#Summary statistics#####################################################
summary(AllDat)

#Room temp
matplot(AllDat[,3:6], type='l', xaxt='n', lty=c(1,1,1,1), ylab="Temperature, C", xlab="Date",
        main="Temperature in each of the 4 rooms")
axis(1, round(seq(100,2910, length.out=3)), AllDat$date[round(seq(100,2910, length.out=3))])

#Ambient Temp
plot(AllDat$date,AllDat$Ta, type='l', ylab="Temperature, C", xlab="Date", main="Ambient Temperature")

#GHSR
plot(AllDat$date,AllDat$Gv, type='l', ylab="GHSR, W/m^2", xlab="Date", main="Global Horizontal Solar Radiation")

# Heating power north
plot(AllDat$date,AllDat$Ph1, type='l', ylab="Heating power, W", xlab="Date", main="Heating power in northern circuit")

#Heating power south
plot(AllDat$date,AllDat$Ph2, type='l', ylab="Heating power, W", xlab="Date", main="Heating power in southern circuit")
lines(AllDat$date,AllDat$Ph1, col='red')

#More heating in southern part
plot(AllDat$date,AllDat$Ph1-AllDat$Ph2, ylab="Heating power, W", xlab="Date", main="Difference of heating power in northern vs southern circuit")

####################################################################


#fit1 <- sdeTiTm(AllDat,AllDat$yTi1,AllDat$Ph1)

#summary(fit1,extended=TRUE)

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

#Pred <- predict(fit1)
#plot(Pred[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour)
# What is going on 10 AM?
# People are going to meeting up.
# Try to fit a varying effective window area

# Looking at the radiation is varying
plot(AllDat$Gv ~ Hour, xlab="Hours", ylab="Gv, W/m^2", main="Global Horizontal Solar Radiation")

# It is impossible to fit a window area for the hours without any sun, 
# so we limit the window area estimation to the hours with sun.
idx <- (Hour>8 & Hour < 23) 
bs = bs(Hour[idx],df=7,intercept=TRUE, degree=3) 
#df = 6 extra spline

# What does the splines look like?
plot(bs[14:27,1],type='l')
lines(bs[ 14:27,2])
lines(bs[ 14:27,3])
lines(bs[ 14:27,4])
lines(bs[ 14:27,5])
lines(bs[ 14:27,6])
lines(bs[ 14:27,7])
#bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- bs6 <- numeric(dim(AllDat)[1])
bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- bs6 <- bs7  <- numeric(dim(AllDat)[1])

bs1[idx] = bs[,1]
bs2[idx] = bs[,2]
bs3[idx] = bs[,3]
bs4[idx] = bs[,4]
bs5[idx] = bs[,5]
bs6[idx] = bs[,6]
bs7[idx] = bs[,7]

AllDat$bs1 = bs1
AllDat$bs2 = bs2
AllDat$bs3 = bs3
AllDat$bs4 = bs4
AllDat$bs5 = bs5
AllDat$bs6 = bs6
AllDat$bs7 = bs7

###########################SUBSET OF DATA
SubDat <- AllDat

#########################################
###SIMPLEST MODEL:
source("sdeTitm.R")
fit1 <- sdeTiTm(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit1,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
par(.pardefault)
summary(fit1)

###############
### SPLINE MODEL:
source("sdeTitm2.R")
fit2 <- sdeTiTm2(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit2,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
par(.pardefault)
summary(fit2, extended=T)
fit2$loglik


###############################ITERATION 1

#Add heater
source("sdeTitm4.R")
fit4 <- sdeTiTm4(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit4,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
par(.pardefault)
summary(fit4, extended=T)
fit4$loglik
summary(fit4)

#add envelope
source("sdeTitm5.R")
fit5 <- sdeTiTm5(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit5,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
reset_par()
summary(fit5)
fit5$loglik
summary(fit5, extended=T)






source("sdeTitm6.R")
fit6<- sdeTiTm6(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit6,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
par(.pardefault)
summary(fit6)
fit6$loglik
summary(fit6, extended=T)


source("sdeTitm8.R")
fit8 <- sdeTiTm5(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit8,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
summary(fit8)
fit8$loglik
summary(fit8, extended=T)

lr.test(127,108, df=2)

source("sdeTitm9.R")
fit9 <- sdeTiTm9(SubDat,SubDat$yTi1,SubDat$Ph1)
analyzeFit(fit9,tPer=c(min(SubDat$date), max(SubDat$date)), newdev=F)
fit9$loglik
summary(fit9, extended=T)



tmp <- predict(fit2)[[1]]
SubDat$residuals <- SubDat$yTi1 - tmp$output$pred$yTi
SubDat$yTiHat <- tmp$output$pred$yTi

acf(SubDat$residuals, lag=100)
cpgram(SubDat$residuals)

plot(bs[14:27,1]*fit2$xm[3]+bs[14:27,2]*fit2$xm[4]+bs[14:27,3]*fit2$xm[5]+bs[14:27,4]*fit2$xm[6]+bs[14:27,5]*fit2$xm[7],type='l')


source("sdeTiTm4Model.R")
SubDat<-AllDat[1:200,]
fit10 <- sdeTiTm4room(SubDat)
reset_par()
tmp <- predict(fit10)[[1]]
par(mfrow=c(2,2))
for(i in 1:4){
SubDat$residuals <- SubDat[i+2,] - tmp$output$pred[i]
acf(SubDat$residuals, lag=24*8)
}
par(mfrow=c(2,2))
for(i in 1:4){
  SubDat$residuals <- SubDat[i+2,] - tmp$output$pred[i]
  spec.pgram(SubDat$residuals, main="Raw periodogram")
}
par(mfrow=c(2,2))
for(i in 1:4){
  SubDat$residuals <- SubDat[i+2,] - tmp$output$pred[i]
  cpgram(SubDat$residuals, main="Cumulated periodogram")
}


source("sdeTiTm4Model.R")
SubDat<-AllDat[1:1000,]
fit11 <- sdeTiTm4room(SubDat)
reset_par()
tmp <- predict(fit11)[[1]]
par(mfrow=c(2,2))
for(i in 1:4){
  SubDat$residuals <- SubDat[,i+2] - tmp$output$pred[i]
  acf(SubDat$residuals, lag=24*8)
}
par(mfrow=c(2,2))
for(i in 1:4){
  SubDat$residuals <- SubDat[,i+2] - tmp$output$pred[i]
  spec.pgram(SubDat$residuals, main="Raw periodogram")
}
par(mfrow=c(2,2))
for(i in 1:4){
  SubDat$residuals <- SubDat[,i+2] - tmp$output$pred[i]
  cpgram(SubDat$residuals, main="Cumulated periodogram")
}

fit11$loglik
