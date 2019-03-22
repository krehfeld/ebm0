## EBM functions
## 0DEBM - implementation and testing


ebm0<-function(t,y,Q,C=210241850,alpha=0.3,eps=0.6){
  sigma<-5.67*10^(-8) 
  list(1/C*((1-alpha)*Q-sigma*eps*y^4))
}
# 
# alpha<-0.3
# 
# 
# f<-0.7
# h<-100
# C<-f*rho*c*h
# to change the temperature of 1m^2 of the Earth by 1K will take on average 2.08*10^8 J/(m^2K)
# S<-1365 #solar constant in W/m^2

#Flux<-cbind((seq(10,100,by=10)),(rep(c(Q_0,1.5*Q_0),5)))
#Flux<-matrix(ncol=2,byrow=TRUE,data=c(1,Q_0,10,1.5*Q_0,20,Q_0,30,1.5*Q_0,40,Q_0,50,1.5*Q_0,60,Q_0))
#Flux[,1]<-Flux[,1]*60*60*24*365
heatcap<-function(h=70,f=0.7,rho=1025,c=4186){
  # effective heat capacity of a 1 square meter area of the earth
  # C = f*rho*c_water*h
  # f fraction of Earth covered by water (.7)
  # rho density of sea water 1025kg/m^3
  # c specific heat of water (4186 J/kgK)
  # h depth of the mixed layer (70m)
  C<-f*rho*c*h
  print(C/(60*60*24*365))
  return(C)
}

heatcap.inv<-function(tau,f=0.7,rho=1025,c=4186){
  # effective heat capacity of a 1 square meter area of the earth
  # C = f*rho*c_water*h
  # f fraction of Earth covered by water (.7)
  # rho density of sea water 1025kg/m^3
  # c specific heat of water (4186 J/kgK)
  # h depth of the mixed layer (70m)
  h<-(tau*60*60*24*365)/(f*rho*c)
  C=f*rho*c*h
  return(list(C=C,h=h))
}


ebm.0d<-function(time.out,Flux,t.initial=288,interp="constant",plot=FALSE,parms=list(f=0.7,h=70,S=1365,eps=0.6,alpha=0.3)){
#  attach(parms)
# effective heat capacity of a 1 square meter area of the earth
# C = f*rho*c_water*h
# f fraction of Earth covered by water (.7)
# rho density of sea water 1025kg/m^3
# c specific heat of water (4186 J/kgK)
# h depth of the mixed layer (70m)
Q_0<-1/4*parms$S
rho<-1025
c<-4186
C<-parms$f*rho*c*parms$h

if (interp=="constant") Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "constant",rule=2)
if (interp=="linear") Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "linear",rule=2)

#params<-c(C=C,alpha=parms$alpha,eps=parms$eps,Q_0=Q_0)
ebm1<-function(t,y,parms){
  Q<-Frc(t)+Q_0
  sigma<-5.67*10^(-8) 
  list(1/C*((1-parms$alpha)*Q-sigma*parms$eps*y^4))
}

yout<-ode(y=t.initial,func=ebm1,times=time.out,parms=parms)
if (plot) {plot(yout[,1],yout[,2])
  grid()
}
forcing.ts<-sapply(time.out,Frc)
rout<-list(yout=yout,forcing=forcing.ts)
return(rout)
}


ebm.0d.icealb<-function(time.out,Flux,t.initial=288,interp="constant",plot=FALSE,parms=list(f=0.7,h=70,S=1365,eps=0.6,alphamax=0.3,kappa=1)){
  #  attach(parms)
  # effective heat capacity of a 1 square meter area of the earth
  # C = f*rho*c_water*h
  # f fraction of Earth covered by water (.7)
  # rho density of sea water 1025kg/m^3
  # c specific heat of water (4186 J/kgK)
  # h depth of the mixed layer (70m)
  Q_0<-1/4*parms$S
  rho<-1025
  c<-4186
  C<-parms$f*rho*c*parms$h
  
  if (interp=="constant") Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "constant",rule=2)
  if (interp=="linear") Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "linear",rule=2)
  
  #params<-c(C=C,alpha=parms$alpha,eps=parms$eps,Q_0=Q_0)
  ebm2<-function(t,y,parms){
    Q<-Frc(t)+Q_0
    Tc=273
    c2=0.7
    c1=0.15
#    alphat<-parms$alphamax+b(1-tanh(parms$kappa(y-272)))/2
       alphat<-c1+c2*(1-tanh(parms$kappa(y-Tc)))/2
    sigma<-5.67*10^(-8)
    list(1/C*((1-alphat)*Q-sigma*parms$eps*y^4))
  }
  
  yout<-ode(y=t.initial,func=ebm2,times=time.out,parms=parms)
  if (plot) {plot(yout[,1],yout[,2])
    grid()
  }
  forcing.ts<-sapply(time.out,Frc)
  rout<-list(yout=yout,forcing=forcing.ts)
  return(rout)
}
# yout<-ebm.0d(seq(0,60*60*24*365*100,by=60*60*24),Flux=cbind(0,10),plot=TRUE)
## Utility function changes

# in: time (in years), forcing (in W/m^2), params: C, alpha, eps, adjust 
# inside: adapt time,get equilibrium solution,initialize with equilibrium solution, calculate effective dQ, call ebm; adjust returned time axis
# return: time/temperature



# Fyears<-c(seq(-equilyears,2000,by=1))
# Fval<-rep(0,length(Fyears))


# random forcing
#Fval[(equilyears+1):length(Fyears)]<-rnorm(length(Fyears)-equilyears,mean = 0,sd = 5) #dQ
# Fval[(equilyears+1):length(Fyears)]<-rnorm(length(Fyears)-equilyears,mean = 0,sd = 5) #dQ
# Flux<-cbind(Fyears,Fval)
# step forcing
#Flux<-matrix(ncol=2,byrow=TRUE,data=c(1,Q_0,30,1.5*Q_0,80,Q_0,100,1.5*Q_0))

getebm.0d<-function(Flux,inttime,tunit="years",equilyears=500,plot=FALSE,interp="constant",parms=list(f=0.7,h=70,S=1365,eps=0.6,alpha=0.3)){
if (tunit=="years") {
  Flux[,1]<-Flux[,1]*60*60*24*365
  inttime<-inttime*60*60*24*365} else {tunit<-"seconds"}
#Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "constant",rule=2)
#inttime<-Flux[,1]
rout<-ebm.0d(inttime,Flux=Flux,t.initial=288,parms=parms,interp=interp,plot=FALSE)
yout<-rout$yout
yout.cut<-yout[-c(1:equilyears),]
forcing.ts<-rout$forcing[-c(1:equilyears)]
tsout<-zoo(yout.cut[,2],order.by=yout.cut[,1]/(60*60*24*365))

if (plot) plot(tsout)

rlist<-list(tsout=tsout,forcing.ts=forcing.ts,misc=list(tunit=tunit,yout=yout,equilyears=equilyears,parms=parms))
return(rlist)
}
#E<-getebm.0d(seq(0,1000),Flux = cbind(seq(0,1000),rnorm(1001)))
#plot(E$tsout)



#' drawfrom
#'
#' @param SP Segmented Forcing data (list of data frames with at least one column "Volcanic")
#' @param yearsout number of years of forcing to be generated (integer)
#' @param plot if true, plot (logical)
#' @param repzero number of times zero eruption segments are added (integer)
#' @param byname draw b
#'
#' @return Concatenated forcing table
#' @export
#'
#' @examples
drawfrom<-function(SP,yearsout,plot=FALSE,repzero=0,byname="dR_Volc",fixvarout=FALSE,targetvarout=1){
  # repzero<-number of times SP is added as a zero/novolcanic segment list, i.e. 1/repzero is the fraction of no volcanic eruption segments
  Ftab<-do.call(rbind,SP)
  qvar<-0  
  
    
  SP.zero<-lapply(SP,function(x,val){x[[byname]]<-rep(val,length(x[,1]));return(x)},val=max(Ftab[[byname]]))
  SPmix<-c(SP,rep(SP.zero,repzero))
  SPmixtab<-do.call(rbind,SPmix)
  
  
  # change to SPmix below
  nseg<-length(SPmix)

  # get variance in original
  addnoise<-function(seg,qvar){
    z<-seg[[byname]];
    znew<-z+rnorm(length(seg[[byname]]),sd=sqrt(qvar));
    # print(var(z))
    # print(var(znew))
    seg[[byname]]<-znew;
    return(seg)
  }
  
  if (fixvarout){
    var.in<-var(SPmixtab[[byname]])
    qvar<-targetvarout-var.in
    # print(qvar)
    if (qvar<0) warning("target variance smaller than original variance! Check.")
    # add white noise variance 
    SPmix<-lapply(SPmix,addnoise,qvar)
  }
    # print(var(do.call(rbind,SPmix)[[byname]]))
    # print(var(do.call(rbind,SP)[[byname]]))
  
  
  nyears<-length(do.call(rbind,SPmix)[[byname]])
  #2*yearsout/nyears
  # find out whether the ratio of yearsout/length training set is appropriate
  #indices
  ind<-round(runif(nseg*(2*yearsout/nyears),min=0,max=nseg))
  #reassemble
  permuted<-(do.call(rbind,SPmix[ind]))[1:yearsout,]
  permuted$YearOld<-permuted$Year
  permuted$Year<-seq(1,yearsout)
  # now cut to length
  # length(permuted$Volcanic)
  
  if (plot) plot(permuted[[byname]],type="l")
  return(permuted)
}


plotxy<-function(tsx1,tsy2,col1="indianred",col2="darkblue",timunit="year CE",ax2unit=expression(bold(paste("GMT "," ", group("[",K,"]")))), ax4unit=expression(bold(paste(bold(Forcing)," ", group("[",W*m^{-2},"]")))),ylimy1=range(tsx1),ylimy2=range(tsy2),yax2=FALSE){
  plot(tsx1-median(tsx1),col=col1,xlim=trange,axes=FALSE,xlab="",ylab="",ylim=ylimy1)
  axis(1,at=axTicks(1),col="black",labels = FALSE)
  mtext(side=1,at=axTicks(1),col="black",line=0.5,text = axTicks(1))
  mtext(side=1,line=1.5,timunit,font=2,col="black")
  axis(2,col=col1,col.ticks = col1,at=axTicks(2),labels=FALSE)
  mtext(side=2,at=axTicks(2),cex=1,las=1,line=0.6,text=round(axTicks(2),digits = 1),col=col1)
  # mtext(side=2,"GMT response",font=2,col=col1,line=2)
  mtext(side=2,col=col1,ax2unit,font=2,line=2,las=0)
  #lines(ts2,col=col2)

    par(new=TRUE)

  if (ylimy2==range(tsy2)){
  addfac<-diff(range(tsy2))
  ylimy<-c(min(tsy2)-1.1*addfac,max(tsy2)) # to set at the bottom, for a negative scale
  # response
  } else {ylimy=ylimy2}
  plot(tsy2,col=col2,axes=FALSE,xlim=trange,ylim=rev(ylimy),xlab="",ylab="")
  if (!yax2){
  yax<-seq(round(min(tsy2)/10)*10,round(max(tsy2)/10)*10,length.out = 5)
  } else {yax=yax2}
  axis(4,at=yax,col.ticks = col2,col.axis=col2,col=col2,labels=FALSE)
  mtext(side = 4,col=col2,at=yax,text = yax,cex = 1,las=1,line=0.8)
  mtext(side=4,col=col2,ax4unit,font=2,line=2.5,at=mean(range(tsy2)),las=0)
}



plotxy_inset<-function(tsx1,tsy2,col1="indianred",col2="darkblue",timunit="year CE",ax2unit=expression(bold(paste("GMT response"," ", group("[",K,"]")))), ax4unit=expression(bold(paste(bold(Forcing)," ", group("[",W*m^{-2},"]")))),ylimy1=range(tsx1),ylimy2=range(tsy2),yax2=FALSE){
  plot(tsx1-median(tsx1),col=col1,xlim=trange,axes=FALSE,xlab="",ylab="",ylim=ylimy1)
  axis(1,at=axTicks(1),col="black",labels = FALSE)
  mtext(side=1,at=axTicks(1),col="black",line=0.5,text = axTicks(1))
  mtext(side=1,line=1.5,timunit,font=1,col="black")
  axis(2,col=col1,col.ticks = col1,at=axTicks(2),labels=FALSE)
  mtext(side=2,at=axTicks(2),las=1,line=0.6,text=round(axTicks(2),digits = 1),col=col1)
  # mtext(side=2,"GMT response",font=2,col=col1,line=2)
  mtext(side=2,col=col1,ax2unit,font=2,line=2,las=0)
  #lines(ts2,col=col2)
  
  par(new=TRUE)
  
  if (ylimy2==range(tsy2)){
    addfac<-diff(range(tsy2))
    ylimy<-c(min(tsy2)-1.1*addfac,max(tsy2)) # to set at the bottom, for a negative scale
    # response
  } else {ylimy=ylimy2}
  plot(tsy2,col=col2,axes=FALSE,xlim=trange,ylim=rev(ylimy),xlab="",ylab="")
  if (!yax2){
    yax<-seq(round(min(tsy2)/10)*10,round(max(tsy2)/10)*10,length.out = 5)
  } else {yax=yax2}
  axis(4,at=yax,col.ticks = col2,col.axis=col2,col=col2,labels=FALSE)
  mtext(side = 4,col=col2,at=yax,text = yax,las=1,line=0.8)
  mtext(side=4,col=col2,ax4unit,font=2,line=2.5,at=mean(range(tsy2)),las=0)
}



#' Get the 0d Energy Balance including a temperature-dependent albedo
#'
#' @param Flux 2D vector with time and external radiative flux in W/m^2
#' @param inttime integration time
#' @param tunit time unit (default "years")
#' @param equilyears number of equilibration years that will be removed from inttime for the initial model equilibration
#' @param plot logical TRUE/FALSE
#' @param interp interpolation scheme for the radiative flux ("constant", "linear")
#' @param parms additional parameters for the EBM
#'
#' @return list 
#' @export
#'
#' @examples
#' inttime<-seq(0,100)
#' Flux<-rep(0,length(inttime))
#' getebm.0d.tempdep(Flux,inttime,plot=TRUE)
getebm.0d.tempdep<-function(Flux,inttime,tunit="years",equilyears=0,plot=FALSE,interp="constant",parms=list(f=0.7,h=70,S=1365,eps=0.6,alphamax=0.3,kappa=1,Tinit=288)){
  if (tunit=="years") {
    Flux[,1]<-Flux[,1]*60*60*24*365
    inttime<-inttime*60*60*24*365} else {tunit<-"seconds"}
  #Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "constant",rule=2)
  #inttime<-Flux[,1]
  rout<-ebm.0d.icealb(inttime,Flux=Flux,parms=parms,interp=interp,plot=FALSE)
  yout<-rout$yout
  yout.cut<-yout[-c(1:equilyears),]
  forcing.ts<-rout$forcing[-c(1:equilyears)]
  tsout<-zoo(yout.cut[,2],order.by=yout.cut[,1]/(60*60*24*365))
  c1<-0.15 #min albedo
  c2<-0.7 # max albedo = c1+c2
  albedo=zoo(c1+c2*(1-tanh(parms$kappa*(coredata(tsout)-273)))/2,order.by=index(tsout))
  #GHG
  T06<-1.9*10^{-15} # [K^{-6}]sellers 1969, ghil 1976, T_0^{-6]}
  m<-0.4
  gofT<-1-m*tanh((coredata(tsout)^6*T06))
  sigma<-5.67*10^(-8)
  
  Rout<-zoo(sigma*gofT*coredata(tsout)^4,order.by=index(tsout))
  Rin<-((1-coredata(albedo))*(parms$S/4)+coredata(forcing.ts))
  
  if (plot) {
    par(mfcol=c(3,1),mar=c(4,4,0,1))
    plot(tsout,axes=FALSE,xlab="",ylab="")
    axis(2);mtext(side=2,line=2,"GMST [K]")
    abline(h=273,col="indianred")
    plot(albedo,axes=FALSE,xlab="",ylab="")
    abline(h=c(c1,c1+c2),lty=2,col="indianred")
    axis(2);mtext(side=2,line=2,"albedo [0-1]")
    plot(Rout,xlab="",ylab="")
    mtext(side=2,line=2,"R_out [W/m^2]")
  }
  
  rlist<-list(tsout=tsout,forcing.ts=forcing.ts,misc=list(tunit=tunit,yout=yout,equilyears=equilyears,parms=parms),albedo=albedo,Rout=Rout,Rin=Rin)
  return(rlist)
}


#' Energy balance model with temperature-dependent albedo (internal function)
#'
#' @param time.out vector of integration time
#' @param Flux vector with radiative forcing 
#' @param interp interpolation type for Flux ("constant","linear")
#' @param plot logical
#' @param parms additional parameters to pass to EBM
#'
#' @return list
#' @export
#'
#' @examples
ebm.0d.icealb<-function(time.out,Flux,interp="constant",plot=FALSE,parms=list(f=0.7,h=70,S=1365,eps=0.6,alphamax=0.3,kappa=1,Tinit=288)){
  #  attach(parms)
  # effective heat capacity of a 1 square meter area of the earth
  # C = f*rho*c_water*h
  # f fraction of Earth covered by water (.7)
  # rho density of sea water 1025kg/m^3
  # c specific heat of water (4186 J/kgK)
  # h depth of the mixed layer (70m)
  Q_0<-1/4*parms$S
  rho<-1025
  c<-4186
  C<-parms$f*rho*c*parms$h
  
  if (interp=="constant") Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "constant",rule=2)
  if (interp=="linear") Frc<-approxfun(x=Flux[,1],y=Flux[,2],method = "linear",rule=2)
  
  #params<-c(C=C,alpha=parms$alpha,eps=parms$eps,Q_0=Q_0)
  ebm2<-function(t,y,parms){
    Q<-Frc(t)+Q_0
    Tc=273
    c2=0.7
    c1=0.15
    T06<-1.9*10^{-15} # [K^{-6}]sellers 1969, ghil 1976, T_0^{-6]}
    m<-0.4
    gofT<-1-m*tanh((y^6*T06))
    #    alphat<-parms$alphamax+b(1-tanh(parms$kappa(y-272)))/2
    alphat<-c1+c2*(1-tanh(parms$kappa*(y-Tc)))/2
    sigma<-5.67*10^(-8)
    #    eb<-1/C*((1-alphat)*Q)-sigma*parms$eps*y^4 # use gofT instead of parms$eps
    #    list(eb)
    list(1/C*((1-alphat)*Q-sigma*gofT*y^4))
  }
  #print(parms$Tinit)
  yout<-ode(y=parms$Tinit,func=ebm2,times=time.out,parms=parms)
  head(yout)
  if (plot) {plot(yout[,1],yout[,2])
    grid()
  }
  
  
  forcing.ts<-sapply(time.out,Frc)
  rout<-list(yout=yout,forcing=forcing.ts)
  return(rout)
}

