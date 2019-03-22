var_from_spec<-function(spec.win,tsc.in,pval=0.1){


### Definitions: 
results<-list(std=NA,var.tsc=NA,var.tot<-NA,dof=NA,var.ci=list(up=NA,lo=NA))
std<-var<-dof<-var.tsc<-var.tot<-ts.used<-c(NA);

### Checks:
# print(any(is.na(timser)))

# Checking the order of the supplied timescales
if (tsc.in[2]<tsc.in[1]) tsc.in<-tsc.in[2:1]
fcut=1/tsc.in; #omega.upper omega.lower





#spec.win<-spec.pgram(tmp,pad=0,taper=FALSE,plot=FALSE,detrend=detrend,fast=FALSE)
#contains p [power] f [frequency] df [degrees of freedom] and 

# find the indices of the spectral estimates in the tsc.in-window 
ind<-c(which.min(abs(fcut[1]-spec.win$f)):which.min(abs(fcut[2]-spec.win$f)))

if (length(ind)< 3) warning("Less than three spectral estimates in the frequency window")

# take the mean over these spectral estimates and multiply by 2* bandwidth to get a correct variance
var.tsc<-mean(spec.win$spec[ind])*(max(spec.win$f[ind])-min(spec.win$f[ind]))*2
var.tot<-mean(spec.win$spec)

# get the degrees of freedom of this variance estimate
dof<-length(ind)*spec.win$df

std<-sqrt(var.tsc)




if (pval>0.1) warning("p-value larger than .1 - is this on purpose or misspecified significance level?")
var.ci<-var_ci(var.tsc,dof,pval=pval)

results<-list(var.tsc=var.tsc,var.tot=var.tot,dof=dof,var.ci=var.ci,spec.win=spec.win)
return(results)
}



