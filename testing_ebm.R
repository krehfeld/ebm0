## 0DEBM - implementation and testing

library(deSolve)
library(zoo)



ebm0<-function(t,y,Q,C=210241850,alpha=0.3,eps=0.6){
sigma<-5.67*10^(-8) 
  list(1/C*((1-alpha)*Q-sigma*eps*y^4))
}

alpha<-0.3

# effective heat capacity of a 1 square meter area of the earth
# C = f*rho*c_water*h
# f fraction of Earth covered by water (.7)
# rho density of sea water 1025kg/m^3
# c specific heat of water (4186 J/kgK)
# h depth of the mixed layer (70m)
f<-0.7
rho<-1025
c<-4186
h<-100
C<-f*rho*c*h
# to change the temperature of 1m^2 of the Earth by 1K will take on average 2.08*10^8 J/(m^2K)

S<-1365 #solar constant in W/m^2
Q_0<-1/4*S


eps<-0.6 # emissivity

# ebm0(0,15,C,Q,alpha,eps)

x_initial<-285
yout<-ode(y=x_initial,func=ebm0,times=seq(0,60*60*24*365*30,by=60*60),parms=Q_0)
plot(yout[,1],yout[,2]-273.15)


times<-seq(0,60*60*24*365*30,by=60*60*24)
Q<-rnorm(length(times),mean=Q_0)

yout<-ode(y=x_initial,func=ebm0,times=times,parms=Q_0)
plot(yout[,1],yout[,2]-273.15)


Flux<-matrix(ncol=2,byrow=TRUE,data=c(1,Q_0,24*60*60*365*10,1.5*Q_0,24*60*60*365*20,Q_0))

Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "constant",rule=2)
Frc2<-approxfun(x=Flux[,1],y=Flux[,2],method = "linear",rule=2)

parms<-c(C=C,alpha=0.3,eps=0.62,Q_0=1365/4)
#parms<-c(210241850,0.3,0.6)
ebm1<-function(t,y,parms){
  C=parms[1]
  alpha=parms[2]
  eps=parms[3]
  Q<-Frc(t)+Q_0
  sigma<-5.67*10^(-8) 
  list(1/C*((1-alpha)*Q-sigma*eps*y^4))
}

yout<-ode(y=255,func=ebm1,times=times,parms=parms)
plot(yout[,1]/(60*60*24*365),yout[,2]-273.15)
grid()

## Utility function changes

# in: time (in years), forcing (in W/m^2), params: C, alpha, eps, adjust 
# inside: adapt time,get equilibrium solution,initialize with equilibrium solution, calculate effective dQ, call ebm; adjust returned time axis
# return: time/temperature

equilyears<-500

Fyears<-c(seq(-equilyears,2000,by=1))
Fval<-rep(0,length(Fyears))
Fval[(equilyears+1):length(Fyears)]<-rnorm(length(Fyears)-equilyears,mean = 0,sd = 5) #dQ
tunit<-"years"
if (tunit=="years") Ftime<-Fyears*60*60*24*365
Frc2<-approxfun(x=Ftime,y=Fval,method = "linear",rule=2)

ebm2<-function(t,y,parms){
  C=parms[1]
  alpha=parms[2]
  eps=parms[3]
  Q<-Frc2(t)+Q_0
  sigma<-5.67*10^(-8) 
  list(1/C*((1-alpha)*Q-sigma*eps*y^4))
}
yout<-ode(y=288,func=ebm2,times=Ftime,parms=parms)
plot(yout[,1]/(60*60*24*365),yout[,2]-273.15)

yout.cut<-yout[-c(1:equilyears),]
tsout<-zoo(yout.cut[,2],order.by=yout.cut[,1]/(60*60*24*365))
plot(tsout)
acf(tsout)

forcing.ts<-sapply(index(tsout),Frc2)


#tsout<-yout[-c(1:500),]
plot(sapply(Ftime,Frc2),yout[,2])



