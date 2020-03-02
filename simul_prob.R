# params = (lambda=1,q=.2,h=.1,x0=1,disp=1)
library(readr)
# download the data


#data_ts<-read.table("ncov_ts_2802.csv",sep=",",header=T)
#data_rec<-read.table("ncov_rec_2802.csv",sep=",",header=T)
#data_death<-read.table("ncov_death_2802.csv",sep=",",header=T)

cat("Download data\n")

file_conf = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv"
rawdata_conf<-read_csv(file_conf)

file_rec = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv"
rawdata_rec<-read_csv(file_rec)

file_deaths = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv"
rawdata_deaths<-read_csv(file_deaths)

country=rawdata_conf[,"Country/Region"]

rawdata_conf=as.matrix(rawdata_conf[,5:ncol(rawdata_conf)])
rawdata_rec=as.matrix(rawdata_rec[,5:ncol(rawdata_rec)])
rawdata_deaths=as.matrix(rawdata_deaths[,5:ncol(rawdata_deaths)])

data_conf<-t(apply(rawdata_conf,1,function(x) diff(as.numeric(x))))
data_rec<-t(apply(rawdata_rec,1,function(x) diff(as.numeric(x))))
data_conf[data_conf<0]=0
data_rec[data_rec<0]=0

cov.sim <- function(params,tf,dt=.1){
	names(params)<-c("lambda","q","h","x0")

	t=seq(0,tf,dt)
	
	x=rep(0,tf/dt)
	y=rep(0,tf/dt)
	
	Rx=rep(0,tf)
	Ry=rep(0,tf)
	
	x[1]=params['x0']
	y[1]=0

	Rx[1]=params['x0']
	Ry[1]=0

	counter=0
	for(i in 2:length(t)){
		x[i] = x[i-1]+dt*(params['lambda']*(x[i-1])-params['q']*y[i-1]*x[i-1])
		y[i] = y[i-1]+dt*(params['h']*x[i-1])
		if(t[i]>counter) {
			counter=counter+1
			Rx[counter] = x[i-1]
			Ry[counter] = y[i-1]			
		}
	}



	return(list(Rx,Ry))
}

cov.LogLikelihood <- function(X,params,A,B,tf){
	res=cov.sim(params,tf)
	LogLik = sum(dpois(X,lambda=res[[1]],log=T)) + sum(dgamma(params,A,B,log=T))
	return(LogLik)
}

cov.GlobalLogLikelihood <- function(lam,paramlist,Al,Bl,A,B,tf,print=F){
	LogLik=dgamma(lam,Al,Bl,log=T)
	for(i in 1:nrow(data)){
		p=c(lam,paramlist[i,])
		res=cov.sim(p,tf)
		LogLik = LogLik+sum(dpois(data[i,],lambda=res[[1]],log=T)) + sum(dgamma(paramlist[i,],A,B,log=T))
		if(print){
			cat(i,' ',LogLik,'\n')
		}
	}

	return(LogLik)
}

cov.MCMC <- function(time_series,NITER){
	A=c(1,1,1,1,1)
	B=c(10,10,10,1,1)
	nvar=5
	tf=length(time_series)

	prop_size=100

	params=rgamma(nvar,A,B)
	params_samples=matrix(NA,NITER,nvar)
	params_samples[1,]=params
	
	LogLik = cov.LogLikelihood(time_series,params,A,B,tf) 

	for(i in 1:NITER){
		
		var=sample(1:nvar,1)
		fac=rgamma(1,prop_size,prop_size)
		params_new=params
		params_new[var]=fac*params_new[var]
		LogLik_new = cov.LogLikelihood(time_series,params_new,A,B,tf) 
		logalpha = LogLik_new- LogLik + dgamma(1/fac,prop_size,prop_size,log=T) - dgamma(fac,prop_size,prop_size,log=T) 
		u=runif(1)
		if(log(u)<logalpha){
			params=params_new
			LogLik=LogLik_new
			params_samples[i,]=params_new
		} else {
			params_samples[i,]=params
		}
	}

	return(params_samples)
}

cov.GlobalMCMC <- function(NITER){
	Al=1
	Bl=10
	A=c(1,1,1)
	B=c(10,10,1)
	nvar=4
	nc=nrow(data)
	tf=ncol(data)

	prop_size_lam=1000
	prop_size=100

	params=matrix(NA,nc,nvar-1)
	lambda=rgamma(1,Al,Bl)

	for(i in 1:nc){ 
		params[i,]=rgamma(nvar-1,A,B)
	}

	params_samples=vector("list",NITER)
	params_samples[[1]]=list(param=params,lambda=lambda)
	
	LogLik = cov.GlobalLogLikelihood(lambda,params,Al,Bl,A,B,tf) 

	for(i in 1:NITER){
		
		var=sample(1:nvar,1)
		if(var==1){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			params_new=params
			lambda_new=fac*lambda
			LogLik_new = cov.GlobalLogLikelihood(lambda_new,params_new,Al,Bl,A,B,tf) 
			logalpha = LogLik_new- LogLik + dgamma(1/fac,prop_size_lam,prop_size_lam,log=T) - dgamma(fac,prop_size_lam,prop_size_lam,log=T) 
		} else {
			fac=rgamma(nc,prop_size,prop_size)
			params_new=params
			params_new[,var-1]=fac*params[,var-1]
			lambda_new=lambda
			LogLik_new = cov.GlobalLogLikelihood(lambda_new,params_new,Al,Bl,A,B,tf) 
			logalpha = LogLik_new- LogLik + sum(dgamma(1/fac,prop_size,prop_size,log=T) - dgamma(fac,prop_size,prop_size,log=T))
		}

		u=runif(1)
		if(log(u)<logalpha){
			params=params_new
			lambda=lambda_new
			LogLik=LogLik_new
			params_samples[[i]]=list(param=params_new,lambda=lambda_new)
		} else {
			params_samples[[i]]=list(param=params,lambda=lambda)
		}
	}

	return(params_samples)
}

cov.plot_param <- function(X,param,tf){
	plot(X)

	res=cov.sim(param,tf)
	lines(cumsum(res[[1]]),col="red")
}

cov.plot <- function(X,samples,tf){
	plot(1:tf,X)

	par_mean=colMeans(samples)
	res=cov.sim(par_mean,tf)
	lines(cumsum(res[[1]]),col="red")
}
		
cov.plot_last <- function(X,samples,tf,n=10,ymax=max(X)){
	plot(X,xlim=c(1,tf),ylim=c(0,ymax))
	Nsam=nrow(samples)
	for(i in 1:n){
		res=cov.sim(samples[Nsam-i+1,],tf)
	       	lines(res[[1]],col=rgb(1,0,0,.5),lty=2)
	}
}

cov.plotFromGlobal <- function(ind,samples,tf,n=10,ymax=max(data[ind,])){
	plot(data[ind,],xlim=c(1,tf),main=country[ind,],ylim=c(0,ymax))
	Nsam=length(samples)
	for(i in 1:n){
		ranind=sample((Nsam/2):Nsam,1)
		res=cov.sim(c(samples[[ranind]]$lambda,samples[[ranind]]$param[ind,]),tf)
	       	lines(res[[1]],col=rgb(1,0,0,.5),lty=2)
	}
}


cov.GlobalPlot<- function(samples,n=10,tf=ncol(data)){
	layout(matrix(1:120,12,10))
	par(mar=c(0,0,0,0))
	nc=nrow(data)
	Nsam=length(samples)
	for(i in 1:nc){
		plot(data[i,],xlim=c(0,tf))
		for(k in 1:n){
			p=c(samples[[Nsam-k+1]]$lambda,
			    samples[[Nsam-k+1]]$param[i,])
			res=cov.sim(p,tf)
	       		lines(res[[1]],col=rgb(1,0,0,.5),lty=2)
		}
	}
}

