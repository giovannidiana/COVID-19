# params
#    lambda=1,
#    h=.1,
#    x0=1
#    k=1
#    g


library(readr)
# download the data

cat("Download data\n")

file_conf = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
rawdata_conf<-read_csv(file_conf)

file_rec = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
rawdata_rec<-read_csv(file_rec)

file_deaths = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
rawdata_deaths<-read_csv(file_deaths)

country.conf=rawdata_conf[,"Country/Region"]
province.conf=rawdata_conf[,"Province/State"]
country.rec=rawdata_rec[,"Country/Region"]
province.rec=rawdata_rec[,"Province/State"]
country.deaths=rawdata_deaths[,"Country/Region"]
province.deaths=rawdata_deaths[,"Province/State"]

nc=ncol(rawdata_conf)
rawdata_conf=as.matrix(rawdata_conf[,5:nc])
rawdata_rec=as.matrix(rawdata_rec[,5:nc])
rawdata_deaths=as.matrix(rawdata_deaths[,5:nc])

# Clean confirmed
data.conf=rawdata_conf
data.rec=rawdata_rec
data.deaths=rawdata_deaths

#for(i in 1:nrow(rawdata_conf)){
#	for(k in ncol(rawdata_conf):2){
#		if(data.C[i,k]<data.C[i,k-1]) data.C[i,k-1]=data.C[i,k]
#		if(data.R[i,k]<data.R[i,k-1]) data.R[i,k-1]=data.R[i,k]
#		if(data.D[i,k]<data.D[i,k-1]) data.D[i,k-1]=data.D[i,k]
#	}
#}


## check that deaths = rec are always less than confirmed after correction
#if(any(data.C-data.R-data.D<0)) stop("confirmed less than deaths plus recovered")
#
# if everything went fine then we can calculate the number of infected individuals as

## set labels 

labels.conf=rep("",nrow(data.conf))
labels.rec=rep("",nrow(data.rec))
labels.deaths=rep("",nrow(data.deaths))

with_province=!is.na(province.conf[,1])
for(i in 1:nrow(data.conf)){
    if(with_province[i]){
        labels.conf[i]=paste(country.conf[i,1],",",province.conf[i,1],sep="")
    } else {        
        labels.conf[i]=paste(country.conf[i,1],sep="")
    } 
}

with_province=!is.na(province.rec[,1])
for(i in 1:nrow(data.rec)){
    if(with_province[i]){
        labels.rec[i]=paste(country.rec[i,1],",",province.rec[i,1],sep="")
    } else {        
        labels.rec[i]=paste(country.rec[i,1],sep="")
    } 
}

with_province=!is.na(province.deaths[,1])
for(i in 1:nrow(data.deaths)){
    if(with_province[i]){
        labels.deaths[i]=paste(country.deaths[i,1],",",province.deaths[i,1],sep="")
    } else {        
        labels.deaths[i]=paste(country.deaths[i,1],sep="")
    } 
}

## read labels for which the population size is available
pop<-read.table("../data/countries_population2.csv",header=T,sep=';')
labels=pop[ pop[,1] %in% labels.conf &
            pop[,1] %in% labels.rec &
            pop[,1] %in% labels.deaths,1]
data.conf=data.conf[match(labels,labels.conf),]
data.rec=data.rec[match(labels,labels.rec),]
data.deaths=data.deaths[match(labels,labels.deaths),]

datar=data.rec+data.deaths
data=data.conf-datar
data[data<0]=0
data.I=data
data.C=data.conf
data.R=data.rec
data.D=data.deaths

country=country.conf[match(labels,labels.conf),]
province=province.conf[match(labels,labels.conf),]
popsize=pop[ pop[,1] %in% labels.conf &
             pop[,1] %in% labels.rec &
             pop[,1] %in% labels.deaths,2]


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
	names(params)<-c("lambda","h","x0","k","g")

	t=seq(0,tf,dt)
	
	x=rep(0,tf/dt)
	x[1]=params['x0']
	y=rep(0,tf/dt)
	
	Rx=rep(0,tf)
	Ry=rep(0,tf)
	
	counter=1
	for(i in 2:length(t)){
		x[i] = x[i-1]+dt*(params['lambda']-y[i-1])*x[i-1]
		pp=1/(1+(params['k']/x[i-1])^params['g'])
		y[i] = y[i-1]+dt*(params['h']*x[i-1]*pp)
		if(t[i]+dt>counter) {
			Rx[counter] = x[i]
			Ry[counter] = y[i]			
			counter=counter+1
		}
	}

	return(list(Rx=Rx,Ry=Ry))
}

cov.SIM <- function(params,tf,dt=0.1){
	nr=nrow(params)
	t=seq(0,tf,l=ceiling(tf/dt))
	
	x=matrix(0,nr,ceiling(tf/dt))
	y=matrix(0,nr,ceiling(tf/dt))
	
	Rx=matrix(0,nr,tf)
	Ry=matrix(0,nr,tf)

	x[,1]=params[,'x0']

	counter=1
	for(i in 2:length(t)){
		x[,i] = (x[,i-1]+dt*(params[,'lambda']-y[,i-1])*x[,i-1])
		pp=1/(1+(params[,'k']/x[,i-1])^params[,'g'])
		y[,i] = y[,i-1]+dt*(params[,'h']*x[,i-1]*pp)
		if(t[i]+dt>counter) {
			Rx[,counter] = x[,i]
			Ry[,counter] = y[,i]
			counter=counter+1
						
		}
	}

	return(list(Rx=Rx,Ry=Ry))
	
	#if(any(is.na(c(Rx,Ry)))) {
	#	return(cov.SIM(params,tf,dt/10))
	#} else if(any(c(Rx,Ry)<0)){ 
	#	return(cov.SIM(params,tf,dt/10))
	#} else {
	#	return(list(Rx=Rx,Ry=Ry))
	#}
}

cov.LogLikelihood <- function(ind,params,A,B,tf){
	names(params)<-c("lambda","h","x0","k","g")
	res=cov.sim(params,tf)
	pp=1/(1+(params['k']/res$Rx)^params['g'])
	
	LogLik = sum(dpois(data[ind,],lambda=res$Rx*pp,log=T)) + 
		 sum(dgamma(params,A,B,log=T))
	return(LogLik)
}

cov.GlobalLogLikelihood <- function(gparlist,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(paramlist)
	p=data.frame(lambda=gparlist,
		     h=paramlist[,1],
		     x0=paramlist[,2],
		     k=paramlist[,3],
		     g=paramlist[,4]
		     )
	res=cov.SIM(p,tf)
	pp=1/(1+(p[,'k']/res$Rx)^p[,'g'])
	LogLik = dgamma(gparlist,Al,Bl,log=T)+sum(dpois(data,lambda=res$Rx*pp,log=T)) + sum(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	return(LogLik)
}

cov.GlobalLogLikelihood_vec <- function(gparlist,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(data)
	p=data.frame(lambda=gparlist,
		     h=paramlist[,1],
		     x0=paramlist[,2],
		     k=paramlist[,3],
		     g=paramlist[,4]
		     )
	res=cov.SIM(p,tf)
	pp=1/(1+(p[,'k']/res$Rx)^p[,'g'])
	if(ncol(data)!=ncol(res$Rx)) warning("data and solutions differ in dimensions\n")
	LogLikVec = rowSums(dpois(data,lambda=res$Rx*pp,log=T)) + rowSums(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	LogLik=sum(LogLikVec)+dgamma(gparlist,Al,Bl,log=T)
	return(list(LogLik=temperature*LogLik,LogLikVec=temperature*LogLikVec))
}

cov.MCMC <- function(ind,NITER){
	A=c(1,1,3,1,1)
	B=c(10,10,50,0.1,0.10)
	nvar=5
	tf=ncol(data)

	prop_size=100

	params=rgamma(nvar,A,B)
	params_samples=matrix(NA,NITER,nvar)
	LL_samples=rep(NA,NITER)
	params_samples[1,]=params
	
	LogLik = cov.LogLikelihood(ind,params,A,B,tf) 
	LL_samples[1]=LogLik

	for(i in 1:NITER){
		
		var=sample(1:nvar,1)
		fac=rgamma(1,prop_size,prop_size)
		params_new=params
		params_new[var]=fac*params_new[var]
		LogLik_new = cov.LogLikelihood(ind,params_new,A,B,tf) 
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
		
		LL_samples[i]=LogLik
		
		if(i>20 && i%%10==0){
			par(mar=c(0,3,0,0))
			layout(matrix(1:6,6))
			plot(LL_samples[max(1,i-200):i],type='l')
		       	for(l in 1:5) 
				plot(params_samples[max(1,i-200):i,l],type='l')
		}


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

cov.GlobalMCMC <- function(NITER,show=1,adaptive=FALSE,init=NA,init_labels=NA,verbose=F){
	Al=1
	Bl=10
	A=c(3,1,1,2)
	B=c(300,0.01,0.001,2)
	nvar=5
	nc=nrow(data)
	tf=ncol(data)

	temp=1
	prop_size_lam=100
	prop_size=matrix(10,nc,nvar-1) #* rowSums(data)


	params=matrix(NA,nc,nvar-1)

	## init using single MCMC

	lambda=rgamma(1,Al,Bl)
	for(i in 1:nc){ 
		params[i,]=rgamma(nvar-1,A,B)
	}

	params_samples=vector("list",NITER)
	LL_samples=rep(NA,NITER)
	prop_choice=rep("black",NITER)
	params_samples[[1]]=list(param=params,lambda=lambda)

	if(!is.na(init)){		
		existing_labels=init_labels[init_labels %in% labels]
		params[match(existing_labels,labels),]=init$param[which(init_labels %in% labels),]
		lambda=init$lambda
		params_samples[[1]]=list(param=params,lambda=lambda)
	}
	
	GLL=cov.GlobalLogLikelihood_vec(lambda,params,Al,Bl,A,B,tf,temperature=temp)

	while(any(is.infinite(GLL$LogLikVec))){
		print(labels[which(is.infinite(GLL$LogLikVec))])
		lambda=rgamma(1,Al,Bl)
		for(i in 1:nc){ 
			params[i,]=rgamma(nvar-1,A,B)
		}
	        	
		GLL=cov.GlobalLogLikelihood_vec(lambda,params,Al,Bl,A,B,tf,temperature=temp)
	}

	params_samples[[1]]=list(param=params,lambda=lambda)
	LL_samples[1]=GLL$LogLik

	for(i in 1:NITER){
		
		cat(i,"     \r")	
		
		## first update lambda
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambda_new=fac*lambda
			prop_ind=1
		} else {
			lambda_new=rgamma(1,Al,Bl)
			prop_ind=2
		}

		## recalculate the LogLikelihood
		GLL$LogLik=dgamma(lambda,Al,Bl,log=T)+sum(GLL$LogLikVec)
		LL_samples[i]=GLL$LogLik


		
		if(is.infinite(GLL$LogLik)) stop("log likelihood is infinite")
		GLL_new=cov.GlobalLogLikelihood_vec(lambda_new,params,Al,Bl,A,B,tf,temperature=temp)
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			   cov.MixedGamma(lambda,Al,Bl,lambda/lambda_new,prop_size_lam,prop_size_lam)-
			   cov.MixedGamma(lambda_new,Al,Bl,lambda_new/lambda,prop_size_lam,prop_size_lam)


		## and update lambda
		u=runif(1)
		if(verbose){
			cat("proposal 1 = ", cov.MixedGamma(lambda,Al,Bl,lambda/lambda_new,prop_size_lam,prop_size_lam),'\n',
		    "proposal 2 = ", cov.MixedGamma(lambda_new,Al,Bl,lambda_new/lambda,prop_size_lam,prop_size_lam),'\n',
		    "newlik = ", GLL_new$LogLik,'\n',
		    "oldlik = ", GLL$LogLik, '\n'
		    )
		}

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
			
		}

		params_samples[[i]]=list(param=params,lambda=lambda)
		if(i>20 && i%%20==0 && show){
			par(mar=c(0,3,0,0))
			layout(matrix(1:6,6))
				plot(LL_samples[max(1,i-500):i],type='l')
				plot(unlist(lapply(params_samples[max(1,i-500):i],function(x) x$lambda)),col=prop_choice[max(1,i-500):i],type='l')
		       	for(l in 1:4) 
				plot(unlist(lapply(params_samples[max(1,i-500):i],function(x) x$param[show,l])),col=prop_choice[max(1,i-500):i],type='l')
		}
	}

	return(params_samples)
}

cov.plot_param <- function(ind,param,tf,ymax){
	
	res=cov.sim(param,tf)
	pp=1/(1+(param['k']/res$Rx)^param['g'])
	plot(data[ind,],xlim=c(1,tf),ylim=c(0,ymax))
	lines(res$Rx*pp)
}

cov.plot_last <- function(X,samples,tf,n=10,ymax=max(X)){
	plot(X,xlim=c(1,tf),ylim=c(1e-3,ymax))
	Nsam=nrow(samples)
	randsam=sample((Nsam/2):Nsam,n)
	for(i in 1:n){
		params=samples[randsam[i],]
		names(params)<-c("lambda","h","x0","k","g")
		res=cov.sim(params,tf)
		pp=1/(1+(params['k']/res$Rx)^params['g'])
	       	lines(res$Rx*pp,col="grey",lty=2)
	}
}

cov.plotFromGlobal <- function(ind,samples,tf,n=10,ymax=max(data[ind,])){
	plot(data[ind,],xlim=c(1,tf),main=country[ind,],ylim=c(1e-2,ymax))
	Nsam=length(samples)
	randsam=sample(1:Nsam,n)
	for(i in 1:n){
		params=c(samples[[randsam[i]]]$lambda,samples[[randsam[i]]]$param[ind,])
		names(params)<-c("lambda","h","x0","k","g")
		res=cov.sim(params,tf)
		pp=1/(1+(params['k']/res$Rx)^params['g'])
	       	lines(res$Rx*pp,col=rgb(1,0,0,.5),lty=2)
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
	                pp=1/(1+(p['k']/res$Rx)^p['g'])
	       		lines(res$Rx*pp,col=rgb(1,0,0,.5),lty=2)
		}
	}
}

cov.plotData <- function(country_index){
	barplot(rbind((data.I[country_index,]),
	      (data.R[country_index,]), 
	      (data.D[country_index,])),
	col=c("red","green","purple" ),las=2,names.arg=names(rawdata_conf))
legend("topleft",legend=c("infected","recovered","deaths"),fill=c("red","green","purple" ))
}


