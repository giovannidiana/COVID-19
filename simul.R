# params = (lambda=1,q=.2,h=.1,x0=1)

data<-read.table("ncov.csv",sep=",",header=T)

TF=ncol(data)-4

cov.sim_vec <- function(params,tf){
	t=seq(0,tf,1)
	dt=1
	
	x=matrix(0,nrow(params),tf)
	y=matrix(0,nrow(params),tf)
	
	x[,1]=params[4]
	y[,1]=0

	for(i in 2:length(t)){
		x[,i] = x[,i-1]+dt*(params[,1]*x[,i-1]-params[,2]*y[,i-1]*x[,i-1])
		y[,i] = y[,i-1]+dt*(params[,3]*x[,i-1])
	}

	return(list(x,y))
}

cov.sim <- function(params,tf){
	t=seq(0,tf,1)
	dt=1
	
	x=rep(0,tf)
	y=rep(0,tf)
	
	x[1]=params[4]
	y[1]=0

	for(i in 2:length(t)){
		x[i] = x[i-1]+dt*(params[1]*x[i-1]-params[2]*y[i-1]*x[i-1])
		y[i] = y[i-1]+dt*(params[3]*x[i-1])
	}

	return(list(x,y))
}

cov.LogLikelihood <- function(X,params,A,B,tf){
	res=cov.sim(params,tf)
	LogLik = -tf/2*log(2*pi*params[5]^2)-0.5*sum((X-cumsum(res[[1]]))^2)/params[5]^2 + sum(dgamma(params,A,B,log=T))
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
	       	lines(cumsum(res[[1]]),col=rgb(1,0,0,.5),lty=2)
	}
}

		
