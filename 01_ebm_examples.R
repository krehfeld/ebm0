# EBM_examples

library(deSolve)
library(zoo)
library(multitaper)
library(PaleoSpec)
library(nest)
source("ebm_functions.R")
source("var_from_spec.R")


## solve with step function flux
equilyears<-100
tmax<-500
steplength<-50
Flux<-cbind(c(0,seq(steplength,tmax,by=steplength)),c(0,rep(c(10,0),tmax/steplength/2)))
#Flux<-matrix(ncol=2,byrow=TRUE,data=c(1,Q_0,10,1.5*Q_0,20,Q_0,30,1.5*Q_0,40,Q_0,50,1.5*Q_0,60,Q_0))
plot(Flux)
E<-getebm.0d(Flux,seq(-equilyears,tmax),equilyears = equilyears)
plot(E$tsout)
par(new=TRUE)
plot(Flux,col="blue")


Fyears<-c(seq(-equilyears,2000,by=1))
Fval<-rep(0,length(Fyears))

Fval[(equilyears+1):length(Fyears)]<-rnorm(length(Fyears)-equilyears,mean = 0,sd = 5) #dQ
Flux<-cbind(Fyears,Fval)




plot(Flux)
E<-getebm.0d(Flux,seq(-equilyears,1000))
plot(E$tsout)
acf(E$tsout)


## Matrix with forcing

FMat<-replicate(20,{Fval<-rnorm(length(Fyears),mean = 0,sd = 5)})
FMat[1:equilyears,]<-0

TMat<-apply(FMat,2,function(Fval,inttime){getebm.0d(Flux=cbind(inttime,Fval),equilyears=100,inttime)$tsout},inttime=Fyears)

matplot(TMat,type="l")

TL<-apply(TMat,2,function(x){sp<-spec.pgram(ts(data=x,frequency=1),kernel="daniell",c(10,20),plot=FALSE,detrend = TRUE);sp$dof<-rep(sp$df,length(sp$spec));return(sp)})



## bivariate example
#SP.sm<-LogSmooth(SP)
tsc.in<-c(10,50)
vtsc<-mean(sapply(TL,function(x,tsc.in){var_from_spec(x,tsc.in)$var.tsc},tsc.in=tsc.in))
mean(sapply(TL,function(x,tsc.in){var_from_spec(x,tsc.in)$var.tsc},tsc.in=c(10,50))+sapply(TL,function(x,tsc.in){var_from_spec(x,tsc.in)$var.tsc},tsc.in=c(50,1000))+sapply(TL,function(x,tsc.in){var_from_spec(x,tsc.in)$var.tsc},tsc.in=c(2,10)))
SP.m<-MeanSpectrum(specList = TL)
LPlot(AddConfInterval(SP.m$spec),col = "blue",conf=TRUE,axes=FALSE,xlim=1/c(200,2))
abline(h=vtsc,col="green")

abline(h=mean(SP.m$spec$spec),col="red")
axis(1,at=axTicks(1),labels = 1/axTicks(1))
axis(2)
box()
abline(v=1/tsc.in)


spec.win<-SP.m$spec
plot(1/spec.win$freq,spec.win$spec,log="xy",xlim=c(1000,1),type="l")

tsc.in<-c(10,50)
abline(v=tsc.in)
fcut=1/tsc.in; #omega.upper omega.lower

ind<-c(which.min(abs(fcut[1]-spec.win$freq)):which.min(abs(fcut[2]-spec.win$freq)))
points(1/spec.win$freq[ind],spec.win$spec[ind])

var.tsc<-mean(spec.win$spec[ind])*(max(spec.win$f[ind])-min(spec.win$f[ind]))*2
abline(h=mean(spec.win$spec[ind]))
# read forcing
# get spectra for high/low/original volcanic forcing
# compare mean and variance