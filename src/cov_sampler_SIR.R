# params
#    lambda=1,
#    h=.1,
#    x0=1
#    k=1
#    g


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

nc=ncol(rawdata_conf)-1
rawdata_conf=as.matrix(rawdata_conf[,5:nc])
rawdata_rec=as.matrix(rawdata_rec[,5:nc])
rawdata_deaths=as.matrix(rawdata_deaths[,5:nc])

# Clean confirmed
data.C=rawdata_conf
data.R=rawdata_rec
data.D=rawdata_deaths

for(i in 1:nrow(rawdata_conf)){
	for(k in ncol(rawdata_conf):2){
		if(data.C[i,k]<data.C[i,k-1]) data.C[i,k-1]=data.C[i,k]
		if(data.R[i,k]<data.R[i,k-1]) data.R[i,k-1]=data.R[i,k]
		if(data.D[i,k]<data.D[i,k-1]) data.D[i,k-1]=data.D[i,k]
	}
}


# check that deaths = rec are always less than confirmed after correction
if(any(data.C-data.R-data.D<0)) stop("confirmed less than deaths plus recovered")

# if everything went fine then we can calculate the number of infected individuals as
data.I=data.C-data.R-data.D
data.I.diff=t(apply(data.I,1,diff))
data.I.diff[data.I.diff<0]=0
data.R.diff=t(apply(data.R,1,diff))
data.R.diff[data.R.diff<0]=0
data.D.diff=t(apply(data.D,1,diff))
data.D.diff[data.D.diff<0]=0

data=data.I
datar=data.R+data.D


## set labels 

labels=rep("",nrow(data))

with_province=!is.na(province[,1])
for(i in 1:nrow(data)){
    if(with_province[i]){
        labels[i]=paste(country[i,1],",",province[i,1],sep="")
    } else {        
        labels[i]=paste(country[i,1],sep="")
    } 
}


## read labels for which the population size is available
pop<-read.table("../data/countries_population2.csv",header=T,sep=';')
selected_labels=match(pop[,1],labels)
labels=labels[selected_labels]
data=data[selected_labels,]
datar=datar[selected_labels,]
country=country[selected_labels,]
province=province[selected_labels,]
popsize=pop[,2]

cov.MixedGammaVec <- function(x1,a1,b1,x2,a2,b2){
	
	l1=dgamma(x1,a1,b1,log=T)
	l2=dgamma(x2,a2,b2,log=T)
	sel=which(l2>l1)
	nsel=which(l1>=l2)
	r=rep(NA,length(x1))
	r[sel] = -log(2)+l2[sel]+log(1+dgamma(x1[sel],a1[sel],b1[sel])/dgamma(x2[sel],a2[sel],b2[sel]))
	r[nsel] = -log(2)+l1[nsel]+log(1+dgamma(x2[nsel],a2[nsel],b2[nsel])/dgamma(x1[nsel],a1[nsel],b1[nsel]))

	if(any(is.na(r))){
		w=which(is.na(r))
		print(l1[w],l2[w])
	}
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

cov.sim <- function(params,tf,dt=.01,N=10000){
	names(params)<-c("lambda","lambdaR",
			 "hR","kR","gR",
			 "hA","kA","gA",
			 "hT","kT","gT",
			 "x0")

	t=seq(0,tf,dt)
	
	x0=rep(0,tf/dt)
	r0=rep(0,tf/dt)
	x=rep(0,tf/dt)
	r=rep(0,tf/dt)
	s=rep(0,tf/dt)
	fA=rep(0,tf/dt)
	fR=rep(0,tf/dt)
	fT=rep(0,tf/dt)
	
	x0[1]=params['x0']
	fA[1]=1
	fR[1]=1
	fT[1]=0
	s[1]=N-x0[1]

	Rx=rep(0,tf)
	Rx0=rep(0,tf)
	Rr=rep(0,tf)
	Rr0=rep(0,tf)
	
	counter=1
	
	for(i in 2:length(t)){
		fA[i] = 1+params['hA']/(1+(params['kA']/t[i])^params['gA'])
		fR[i] = 1+params['hR']/(1+(params['kR']/t[i])^params['gR'])
		fT[i] =   params['hT']/(1+(params['kT']/t[i])^params['gT'])
		x0[i] = x0[i-1]+dt*((params['lambda']/fA[i-1])*s[i-1]/N-params['lambdaR']*fR[i-1])*x0[i-1]
		r0[i] = r0[i-1]+dt*(                                    params['lambdaR']*fR[i-1])*x0[i-1]
		s[i]  = s[i-1] -dt*((params['lambda']/fA[i-1])*s[i-1]/N                          )*x0[i-1]
		x[i]  = x0[i]*fT[i]
		r[i]  = r0[i]*fT[i]
		if(t[i]+dt>counter) {
			Rx[counter] = x[i]
			Rx0[counter] = x0[i]
			Rr[counter] = r[i]
			Rr0[counter] = r0[i]
			counter=counter+1
		}
	}

	return(list(Rx=Rx,Rr=Rr,Rx0=Rx0,Rr0=Rr0))
}

cov.SIM <- function(p,tf,dt=0.1){
	nr=nrow(p)
	t=seq(0,tf,l=ceiling(tf/dt))
	
	x=matrix(0,nr,ceiling(tf/dt))
	x0=matrix(0,nr,ceiling(tf/dt))
	s=matrix(0,nr,ceiling(tf/dt))
	r=matrix(0,nr,ceiling(tf/dt))
	r0=matrix(0,nr,ceiling(tf/dt))
	fA=matrix(0,nr,ceiling(tf/dt))
	fR=matrix(0,nr,ceiling(tf/dt))
	fT=matrix(0,nr,ceiling(tf/dt))
	
	Rx=matrix(0,nr,tf)
	Rr=matrix(0,nr,tf)
	
	x0[,1]=p[,'x0']
	fR[,1]=1
	fA[,1]=1
	fT[,1]=0
	s[,1]=popsize-p[,'x0']

	counter=1
	for(i in 2:length(t)){
		fA[,i] = 1+p[,'hA']/(1+(p[,'kA']/t[i])^p[,'gA'])
		fR[,i] = 1+p[,'hR']/(1+(p[,'kR']/t[i])^p[,'gR'])
		fT[,i] =   p[,'hT']/(1+(p[,'kT']/t[i])^p[,'gT'])
		
		x0[,i] = x0[,i-1]+dt*(p[,'lambda']/fA[,i-1]*s[,i-1]/popsize-p[,'lambdaR']*fR[,i-1])*x0[,i-1]
		s[,i]  = s[,i-1] -dt*(p[,'lambda']/fA[,i-1]*s[,i-1]/popsize                       )*x0[,i-1]
		r0[,i] = r0[,i-1]+dt*(                                      p[,'lambdaR']*fR[,i-1])*x0[,i-1]
		x[,i] = x0[,i]*fT[,i] 
		r[,i] = r0[,i]*fT[,i] 
		
		if(t[i]+dt>counter) {
			Rx[,counter] = x[,i]
			Rr[,counter] = r[,i]
			counter=counter+1
		}
	}

	return(list(Rx=Rx,Rr=Rr))
	
	#if(any(is.na(c(Rx,Ry)))) {
	#	return(cov.SIM(params,tf,dt/10))
	#} else if(any(c(Rx,Ry)<0)){ 
	#	return(cov.SIM(params,tf,dt/10))
	#} else {
	#	return(list(Rx=Rx,Ry=Ry))
	#}
}

cov.LogLikelihood <- function(ind,params,A,B,tf){
	names(params)<-c("lambda","lambdaR",
			 "hR","kR","gR",
			 "hA","kA","gA",
			 "hT","kT","gT",
			 "x0")
	res=cov.sim(params,tf,N=popsize[ind])
	
	LogLik = sum(dpois(data[ind,],lambda=res$Rx,log=T)) + 
		 sum(dpois(datar[ind,],lambda=res$Rr,log=T)) +
		 sum(dgamma(params,A,B,log=T))
	return(LogLik)
}

cov.GlobalLogLikelihood <- function(gparlist,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(paramlist)
	p=cbind(lambda=gparlist[1],
		lambdaR=gparlist[2],
		paramlist)
	
	res=cov.SIM(p,tf)
	LogLik = sum(dgamma(gparlist,Al,Bl,log=T))+
		 sum(dpois(data,lambda=res$Rx,log=T)) + 
		 sum(dpois(datar,lambda=res$Rr,log=T)) + 
		 sum(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	return(LogLik)
}

cov.GlobalLogLikelihood_vec <- function(gparlist,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(data)
	p=cbind(lambda=gparlist[1],
		lambdaR=gparlist[2],
		paramlist)
	res=cov.SIM(p,tf)
	
	LogLikVec = rowSums(dpois(data,lambda=res$Rx,log=T)) + 
		    rowSums(dpois(datar,lambda=res$Rr,log=T)) +
		    rowSums(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	LogLik=sum(LogLikVec)+sum(dgamma(gparlist,Al,Bl,log=T))
	return(list(LogLik=temperature*LogLik,LogLikVec=temperature*LogLikVec))
}

cov.MCMC <- function(ind,NITER){
	A=c(1,1,
	    1,1,2,
	    1,1,2,
	    1,1,2,
	    1)
	B=c(10,20,
	    1,.01,2,
	    1,.01,2,
	    1,.1,2,
	    1)
	nvar=12
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
			#par(mar=c(0,3,0,0))
			#layout(matrix(1:6,6))
			#plot(LL_samples[max(1,i-200):i],type='l')
		       	#for(l in 1:5) 
			#	plot(params_samples[max(1,i-200):i,l],type='l')
			cov.plot_last(ind,params_samples[max(1,i-200):i,],100,2)
		}


	}

	return(params_samples)
}

cov.GlobalMCMC <- function(NITER,show=1,adaptive=FALSE,init=NA,init_labels=NA,verbose=F){
	Al=c(1,1)
	Bl=c(1,2)

	A=c(1,1,2, # R 
	    1,1,2, # A
	    1,1,2, # T
	    1)     # x0
	B=c(.01,.01,2,
	    .01,.01,2,
	    100,.1,2,
	    1)


	nvar=12
	nc=nrow(data)
	tf=ncol(data)

	temp=1
	prop_size_lam=100
	prop_size=matrix(10,nc,nvar-2) #* rowSums(data)

	params=matrix(NA,nc,nvar-2)
	colnames(params)=c('hR','kR','gR',
                           'hA','kA','gA',
                           'hT','kT','gT',
                           'x0')

	## init using single MCMC

	lambda=rgamma(1,Al[1],Bl[1])
	lambdaR=rgamma(1,Al[2],Bl[2])
	for(i in 1:nc){ 
		params[i,]=rgamma(nvar-2,A,B)
	}

	params_samples=vector("list",NITER)
	LL_samples=rep(NA,NITER)
	prop_choice=rep("black",NITER)
	params_samples[[1]]=list(param=params,lambda=lambda,lambdaR=lambdaR)

	if(!is.na(init)){		
		existing_labels=init_labels[init_labels %in% labels]
		params[match(existing_labels,labels),]=init$param[which(init_labels %in% labels),]
		lambda=init$lambda
		lambdaR=init$lambdaR
		params_samples[[1]]=list(param=params,lambda=lambda,lambdaR=lambdaR)
	}
	
	GLL=cov.GlobalLogLikelihood_vec(c(lambda,lambdaR),params,Al,Bl,A,B,tf,temperature=temp)

	LL_samples[1]=GLL$LogLik

	while(any(is.na(GLL$LogLikVec))){
		print(labels[which(is.na(GLL$LogLikVec))])
		lambda=rgamma(1,Al[1],Bl[1])
		lambdaR=rgamma(1,Al[2],Bl[2])
		for(i in 1:nc){ 
			params[i,]=rgamma(nvar-2,A,B)
		}
	        	
		GLL=cov.GlobalLogLikelihood_vec(c(lambda,lambdaR),params,Al,Bl,A,B,tf,temperature=temp)
	}

	for(i in 1:NITER){
		
		cat(i,' ',GLL$LogLik,"     \r")	
		
		## first update lambda
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambda_new=fac*lambda
			prop_ind=1
		} else {
			lambda_new=rgamma(1,Al,Bl)
			prop_ind=2
		}

		## recalculate the LogLikelihood from previous iteration
		GLL$LogLik=sum(dgamma(c(lambda,lambdaR),Al,Bl,log=T))+sum(GLL$LogLikVec)
		LL_samples[i]=GLL$LogLik
		
		GLL_new=cov.GlobalLogLikelihood_vec(c(lambda_new,lambdaR),params,Al,Bl,A,B,tf,temperature=temp)
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			   cov.MixedGamma(lambda,Al[1],Bl[1],lambda/lambda_new,prop_size_lam,prop_size_lam)-
			   cov.MixedGamma(lambda_new,Al[1],Bl[1],lambda_new/lambda,prop_size_lam,prop_size_lam)


		## and update lambda
		u=runif(1)

	        if(!is.na(logalpha)){	
		if(log(u)<logalpha){
			lambda=lambda_new
			GLL=GLL_new
			prop_choice[i]=c("black","red")[prop_ind]
		}
		}
	
		##  update lambdaR
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambdaR_new=fac*lambdaR
		} else {
			lambdaR_new=rgamma(1,Al[2],Bl[2])
		}

		GLL_new=cov.GlobalLogLikelihood_vec(c(lambda,lambdaR_new),params,Al,Bl,A,B,tf,temperature=temp)
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			   cov.MixedGamma(lambdaR,Al[2],Bl[2],lambdaR/lambdaR_new,prop_size_lam,prop_size_lam)-
			   cov.MixedGamma(lambdaR_new,Al[2],Bl[2],lambdaR_new/lambdaR,prop_size_lam,prop_size_lam)


		## and update lambdaR
		u=runif(1)

	        if(!is.na(logalpha)){	
		if(log(u)<logalpha){
			lambdaR=lambdaR_new
			GLL=GLL_new
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

			GLL_new = cov.GlobalLogLikelihood_vec(c(lambda,lambdaR),params_new,Al,Bl,A,B,tf,temperature=temp) 
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

		params_samples[[i]]=list(param=params,lambda=lambda,lambdaR=lambdaR)
		if(i>20 && i%%20==0 && show){
			par(mar=c(0,3,0,0))
			layout(matrix(1:7,7))
				plot(LL_samples[max(1,i-200):i],type='l',xaxt='n',xlab="")
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$lambda)),type='l',xaxt='n',xlab="")
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$lambdaR)),type='l',xaxt='n',xlab="")
		       	for(l in 1:4) 
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$param[show,l])),type='l',xaxt='n',xlab="")
		}
	}

	return(params_samples)
}

cov.plot_param <- function(ind,param,tf,ymax){

        nr=nrow(data)	
	res=cov.sim(param,tf,N=popsize[ind])
	plot(data[ind,],xlim=c(1,tf),ylim=c(0,ymax))
	lines(res$Rx,col="red")
	lines(res$Rr,col="green")
}

cov.plot_last <- function(ind,samples,tf,n=10,ymax=max(data[ind,])){
	plot(data[ind,],xlim=c(1,tf),ylim=c(1e-3,ymax))
	points(datar[ind,],col="green")
	Nsam=nrow(samples)
	ransam=sample((Nsam/2):Nsam,n)
        nr=nrow(data)	
	for(i in 1:n){
		params=samples[ransam[i],]
		names(params)<-c("lambda","lambdaR",
				 "hR","kR","gR",
				 "hA","kA","gA",
			 	 "hT","kT","gT",
			 	 "x0")
		res=cov.sim(params,tf,N=popsize[ind])
	       	lines(res$Rx,col="grey",lty=2)
	       	lines(res$Rr,col="green",lty=2)
	}
}

cov.plotFromGlobal <- function(ind,samples,tf,n=10,ymax=NA,log=""){
	Nsam=length(samples)
	randsam=sample((Nsam-1000+1):Nsam,n)
	params=c(samples[[Nsam]]$lambda,samples[[Nsam]]$lambdaR,samples[[Nsam]]$param[ind,])
	names(params)<-c("lambda","lambdaR",
			 "hR","kR","gR",
			 "hA","kA","gA",
		 	 "hT","kT","gT",
		 	 "x0")
	res=cov.sim(params,tf,N=popsize[ind])
	if(is.na(ymax)) ymax=max(c(res$Rx,res$Rr))

	plot(data[ind,],xlim=c(1,tf),main=country[ind,],ylim=c(1,ymax),pch=19,log=log)
	points(datar[ind,],xlim=c(1,tf),col="green",pch=19)
	for(i in 1:n){
		params=c(samples[[randsam[i]]]$lambda,samples[[randsam[i]]]$lambdaR,samples[[randsam[i]]]$param[ind,])
		names(params)<-c("lambda","lambdaR",
				 "hR","kR","gR",
				 "hA","kA","gA",
			 	 "hT","kT","gT",
			 	 "x0")
		res=cov.sim(params,tf,N=popsize[ind])
	       	lines(res$Rx,col=rgb(1,0,0,.5),lty=2)
	       	lines(res$Rr,col=rgb(0,1,0,.5),lty=2)
	       	
		#lines(res$Rx0,col=rgb(1,0,0,.5),lty=3)
	       	#lines(res$Rr0,col=rgb(0,1,0,.5),lty=3)
	}
}

cov.plotData <- function(country_index){
	barplot(rbind((data.I[country_index,]),
	      (data.R[country_index,]), 
	      (data.D[country_index,])),
	col=c("red","green","purple" ),las=2,names.arg=names(rawdata_conf))
legend("topleft",legend=c("infected","recovered","deaths"),fill=c("red","green","purple" ))
}

cov.plotHill <- function(h,k,g,tf,one=T){
	t=1:tf
	plot(1:tf,one+h/(1+(k/1:tf)^g))
}



