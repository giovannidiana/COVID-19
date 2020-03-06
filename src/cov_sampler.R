# params
#    lambda=1,
#    h=.1,
#    t0=1,
#    x0=1
#    k=1


library(readr)
# download the data


#data_ts<-read.table("ncov_ts_2802.csv",sep=",",header=T)
#data_rec<-read.table("ncov_rec_2802.csv",sep=",",header=T)
#data_death<-read.table("ncov_death_2802.csv",sep=",",header=T)

cat("Download data\n")

file_conf = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv"
rawdata_conf<-read_csv(file_conf)

file_rec = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv"
rawdata_rec<-read_csv(file_rec)

file_deaths = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv"
rawdata_deaths<-read_csv(file_deaths)

country=rawdata_conf[,"Country/Region"]
province=rawdata_conf[,"Province/State"]

rawdata_conf=as.matrix(rawdata_conf[,5:ncol(rawdata_conf)])
rawdata_rec=as.matrix(rawdata_rec[,5:ncol(rawdata_rec)])
rawdata_deaths=as.matrix(rawdata_deaths[,5:ncol(rawdata_deaths)])

data_conf<-t(apply(rawdata_conf,1,function(x) diff(as.numeric(x))))
data_rec<-t(apply(rawdata_rec,1,function(x) diff(as.numeric(x))))
data_conf[data_conf<0]=0
data_rec[data_rec<0]=0
data=data_conf[rowSums(data_conf)>0,]
country=country[rowSums(data_conf)>0,]
province=province[rowSums(data_conf)>0,]

cov.MixedGammaVec <- function(x1,a1,b1,x2,a2,b2){
	
	l1=dgamma(x1,a1,b1,log=T)
	l2=dgamma(x2,a2,b2,log=T)
	sel=which(l2>l1)
	nsel=which(l1>=l2)
	r=rep(NA,length(x1))
	r[sel] = -log(2)+l2[sel]+log(1+dgamma(x1[sel],a1[sel],b1[sel])/dgamma(x2[sel],a2[sel],b2[sel]))
	r[nsel] = -log(2)+l1[nsel]+log(1+dgamma(x2[nsel],a2[nsel],b2[nsel])/dgamma(x1[nsel],a1[nsel],b1[nsel]))

	if(any(is.na(r))) cat("there is a problem\n")
	return(r)
}

cov.MixedGamma <- function(x1,a1,b1,x2,a2,b2){
	
	l1=dgamma(x1,a1,b1,log=T)
	l2=dgamma(x2,a2,b2,log=T)
	
	if(l2>l1){
		r = -log(2)+l2+log(1+dgamma(x1,a1,b1)/dgamma(x2,a2,b2))
	} else {
		r = -log(2)+l1+log(1+dgamma(x2,a2,b2)/dgamma(x1,a1,b1))
	}

	return(r)
}

cov.sim <- function(params,tf,dt=.01){
	names(params)<-c("lambda","h","x0","t0","k")

	t=seq(0,tf,dt)
	
	x=rep(0,tf/dt)
	x[1]=params['x0']
	y=rep(0,tf/dt)
	
	Rx=rep(0,tf)
	Ry=rep(0,tf)
	
	ind0=1;while(t[ind0]<params['t0'] && ind0<length(t)) ind0=ind0+1 
	
	counter=1
	for(i in 2:length(t)){
		x[i] = x[i-1]+dt*(params['lambda']*(x[i-1])-y[i-1]*x[i-1])
		y[i] = (i>ind0)*(y[i-1]+dt*(params['h']*x[i-1]))
		if(t[i]+dt>counter) {
			Rx[counter] = x[i]
			Ry[counter] = y[i]			
			counter=counter+1
		}
	}

	return(list(Rx,Ry))
}

cov.SIM <- function(params,tf,dt=0.1){
	nr=nrow(params)
	t=seq(0,tf,l=ceiling(tf/dt))
	
	x=matrix(0,nr,ceiling(tf/dt))
	y=matrix(0,nr,ceiling(tf/dt))
	
	Rx=matrix(0,nr,tf)
	Ry=matrix(0,nr,tf)
	
	ind0=apply(params,1,
		   function(x){
			   ind=which(t>=x['t0'])
			   if(length(ind)>0){
				   return(ind[1])
			   } else {
				   return(length(t))
			   }
		   })

	x[,1]=params[,'x0']

	counter=1
	for(i in 2:length(t)){
		x[,i] = (x[,i-1]+dt*(params[,'lambda']*(x[,i-1])-y[,i-1]*x[,i-1]))
		y[,i] = (i>ind0)*(y[,i-1]+dt*(params[,'h']*x[,i-1]))
		if(t[i]+dt>counter) {
			Rx[,counter] = x[,i]
			Ry[,counter] = y[,i]
			counter=counter+1
						
		}
	}

        
	if(any(is.na(c(Rx,Ry)))) {
		return(cov.SIM(params,tf,dt/10))
	} else if(any(c(Rx,Ry)<0)){ 
		#warning('Negative values')
		#cat('Negative values\n')
		#cat("recursive calc",'\n')
		return(cov.SIM(params,tf,dt/10))
	} else {
		#cat("done\n")
		return(list(Rx,Ry))
	}
}




cov.LogLikelihood <- function(X,params,A,B,tf){
	names(params)<-c("lambda","h","x0","t0","k")
	res=cov.sim(params,tf)
	pp=res[[1]]/(params['k']+res[[1]])
	LogLik = sum(dpois(X,lambda=res[[1]]*pp,log=T)) + sum(dgamma(params,A,B,log=T))
	return(LogLik)
}

cov.GlobalLogLikelihood <- function(lam,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(paramlist)
	p=data.frame(lambda=lam,
		     h=paramlist[,1],
		     x0=paramlist[,2],
		     t0=paramlist[,3],
		     k=paramlist[,4]
		     )
	res=cov.SIM(p,tf)
	pp=res[[1]]/(p[,'k']+res[[1]])
	LogLik = dgamma(lam,Al,Bl,log=T)+sum(dpois(data,lambda=res[[1]]*pp,log=T)) + sum(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	return(LogLik)
}

cov.GlobalLogLikelihood_vec <- function(lam,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(data)
	p=data.frame(lambda=lam,
		     h=paramlist[,1],
		     x0=paramlist[,2],
		     t0=paramlist[,3],
		     k=paramlist[,4]
		     )
	res=cov.SIM(p,tf)
	pp=res[[1]]/(p[,'k']+res[[1]])
	if(ncol(data)!=ncol(res[[1]])) warning("data and solutions differ in dimensions\n")
	LogLikVec = rowSums(dpois(data,lambda=res[[1]]*pp,log=T)) + rowSums(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	LogLik=sum(LogLikVec)+dgamma(lam,Al,Bl,log=T)
	return(list(LogLik=temperature*LogLik,LogLikVec=temperature*LogLikVec))
}

cov.MCMC <- function(time_series,NITER){
	A=c(1,1,1,1)
	B=c(1000,1,1,1)
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
	        if(!is.na(logalpha)){	
		if(log(u)<logalpha){
			params=params_new
			LogLik=LogLik_new
			params_samples[i,]=params_new
		} else {
			params_samples[i,]=params
		}}

	}

	return(params_samples)
}
cov.PriorSample <- function(NITER){
	Al=1
	Bl=2
	A=c(1,1,1)
	B=c(10,10,.2)
	nvar=4
	nc=nrow(data)
	tf=ncol(data)

	params_samples=vector("list",NITER)
	
	for(i in 1:NITER){
		params=matrix(NA,nc,nvar-1)
		lambda=rgamma(1,Al,Bl)
	
		for(k in 1:nc){ 
			params[k,]=rgamma(nvar-1,A,B)
		}

		params_samples[[i]]=list(param=params,lambda=lambda)
	}

	return(params_samples)
}

cov.GlobalMCMC <- function(NITER,id_show=1,adaptive=FALSE,init=NA){
	Al=1
	Bl=2
	A=c(1,1,1,1)
	B=c(10,.1,10,10)
	nvar=5
	nc=nrow(data)
	tf=ncol(data)

	temp=1
	prop_size_lam=50
	prop_size=matrix(40,nc,nvar-1) #* rowSums(data)

	params=matrix(NA,nc,nvar-1)

	## init using single MCMC

	lambda=rgamma(1,Al,Bl)
	for(i in 1:nc){ 
		params[i,]=rgamma(nvar-1,A,B)
		#params[i,]=cov.MCMC(data[i,],100)[100,2:5]
	}

	params_samples=vector("list",NITER)
	prop_choice=rep("black",NITER)
	params_samples[[1]]=list(param=params,lambda=lambda)

	if(!is.na(init)){
		params_samples[[1]]=init
		params=params_samples[[1]]$param
		lambda=params_samples[[1]]$lambda
	}
	
	GLL=cov.GlobalLogLikelihood_vec(lambda,params,Al,Bl,A,B,tf,temperature=temp)

	for(i in 1:NITER){
		## first update lambda
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambda_new=fac*lambda
			prop_ind=1
		} else {
			lambda_new=rgamma(1,Al,Bl)
			prop_ind=2
		}

		cat(i,"     \r")	
		GLL_new=cov.GlobalLogLikelihood_vec(lambda_new,params,Al,Bl,A,B,tf,temperature=temp)
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			   cov.MixedGamma(lambda,Al,Bl,lambda/lambda_new,prop_size_lam,prop_size_lam)-
			   cov.MixedGamma(lambda_new,Al,Bl,lambda_new/lambda,prop_size_lam,prop_size_lam)


		## and update lambda
		u=runif(1)
	        if(!is.na(logalpha)){	
		if(log(u)<logalpha){
			lambda=lambda_new
			GLL=GLL_new
			prop_choice[i]=c("black","red")[prop_ind]
		}
		}

		## then update the other parameters independently

		for(var in 1:ncol(params)){
			which_from_prior=runif(nc)<0.5
			params_new=params


			fac=rgamma(nc,prop_size[,var],prop_size[,var])
			rprior=rgamma(nc,A[var],B[var])			

			params_new[!which_from_prior,var]=fac[!which_from_prior]*params[!which_from_prior,var]
			params_new[which_from_prior,var]=rprior[which_from_prior]

			GLL_new = cov.GlobalLogLikelihood_vec(lambda,params_new,Al,Bl,A,B,tf,temperature=temp) 
			logalpha = GLL_new$LogLikVec - GLL$LogLikVec + 
				   cov.MixedGammaVec(params[,var],rep(A[var],nc),rep(B[var],nc),params[,var]/params_new[,var],prop_size[,var],prop_size[,var]) - 
				   cov.MixedGammaVec(params_new[,var],rep(A[var],nc),rep(B[var],nc),params_new[,var]/params[,var],prop_size[,var],prop_size[,var])

			u=runif(nc)
	        	if(all(!is.na(logalpha))){	
			selec=log(u)<logalpha
			params[selec,var]=params_new[selec,var]
			GLL$LogLikVec[selec]=GLL_new$LogLikVec[selec]
			}
			if(adaptive){
			###################################################
			## Adaptive proposal
			ar.par[,var,]=cbind(ar.par[,var,-1],selec)
			ar.par.means=apply(ar.par,c(1,2),mean)
			nonzero=ar.timer.par[,var]>0
			ar.timer.par[nonzero,var]=ar.timer.par[nonzero,var]-1

			zeros_and_low = ar.par.means[,var]<0.25 & ar.timer.par[,var]==0
			prop_size[zeros_and_low,var] = prop_size[zeros_and_low,var]*2

			zeros_and_high = ar.par.means[,var]>0.8 & ar.timer.par[,var]==0
			prop_size[zeros_and_high,var] = prop_size[zeros_and_high,var]*2

			ar.timer.par[zeros_and_low | zeros_and_high,var] = 100
			###################################################
			}
		}

		params_samples[[i]]=list(param=params,lambda=lambda)
		if(i>20 && i%%10==0){
			par(mar=c(0,3,0,0))
			layout(matrix(1:5,5))
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$lambda)),col=prop_choice[max(1,i-200):i])
		       	for(l in 1:4) 
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$param[43,l])),col=prop_choice[max(1,i-200):i])
		}
	}

	return(params_samples)
}

cov.plot_param <- function(X,param,tf){
	plot(X)

	res=cov.sim(param,tf)
	lines(cumsum(res[[1]]),col="red")
}

cov.plot_last <- function(X,samples,tf,n=10,ymax=max(X)){
	plot(X,xlim=c(1,tf),ylim=c(0,ymax))
	Nsam=nrow(samples)
	randsam=sample((Nsam/2):Nsam,n)
	for(i in 1:n){
		params=samples[randsam[i],]
		names(params)<-c("lambda","h","x0","t0","k")
		res=cov.sim(params,tf)
	        pp=res[[1]]/(params['k']+res[[1]])
	       	lines(res[[1]]*pp,col=rgb(1,0,0,.5),lty=2)
	}
}

cov.plotFromGlobal <- function(ind,samples,tf,n=10,ymax=max(data[ind,])){
	plot(data[ind,],xlim=c(1,tf),main=country[ind,],ylim=c(0,ymax))
	Nsam=length(samples)
	randsam=sample((Nsam/2):Nsam,n)
	for(i in 1:n){
		ranind=sample((Nsam/2):Nsam,1)
		params=c(samples[[randsam[i]]]$lambda,samples[[randsam[i]]]$param[ind,])
		names(params)<-c("lambda","h","x0","t0","k")
		res=cov.sim(params,tf)
		pp=res[[1]]/(params['k']+res[[1]])
	       	lines(res[[1]]*pp,col=rgb(1,0,0,.5),lty=2)
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
			names(params)<-c("lambda","h","x0","t0","k")
			res=cov.sim(p,tf)
	                pp=res[[1]]/(p['k']+res[[1]])
	       		lines(res[[1]]*pp,col=rgb(1,0,0,.5),lty=2)
		}
	}
}


