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

data=data.conf
datar=data.rec+data.deaths
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

	if(any(is.na(r))){
		w=which(is.na(r))
		cat(l1[w],' ',l2[w])
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

cov.evol <- function(params,t,Y,N){
	fun=rep(NA,3)
	fA = 1+params['hA']/(1+(params['kA']/t)^params['gA'])
	fR = 1+params['hR']/(1+(params['kR']/t)^params['gR'])
	fun[1] = ((params['lambda']/fA)*Y[3]/N-params['lambdaR']*fR)*Y[1]
	fun[2] = (                             params['lambdaR']*fR)*Y[1]
	fun[3] =-((params['lambda']/fA)*Y[3]/N                     )*Y[1]
	
	return(fun)
}

cov.EVOL <- function(params,t,Y,N){
	fun=matrix(NA,nrow(params),3)
	fA = 1+params[,'hA']/(1+(params[,'kA']/t)^params[,'gA'])
	fR = 1+params[,'hR']/(1+(params[,'kR']/t)^params[,'gR'])
	fun[,1] = ((params[,'lambda']/fA)*Y[,3]/N-params[,'lambdaR']*fR)*Y[,1]
	fun[,2] = (                               params[,'lambdaR']*fR)*Y[,1]
	fun[,3] =-((params[,'lambda']/fA)*Y[,3]/N                      )*Y[,1]
	
	return(fun)
}

cov.sim <- function(params,tf,dt=.01,N=10000){
	names(params)<-c("lambda","lambdaR",
			 "hR","kR","gR",
			 "hA","kA","gA",
			 "hT","kT","gT",
			 "x0")

	t=seq(dt,tf,dt)
	
	x0=rep(0,tf/dt)
	r0=rep(0,tf/dt)
	x=rep(0,tf/dt)
	r=rep(0,tf/dt)
	s=rep(0,tf/dt)
	fA = 1+params['hA']/(1+(params['kA']/t)^params['gA'])
	#fR = 1+params['hR']/(1+(params['kR']/t)^params['gR'])
	fR = 1+params['hR']/(1+(params['kA']/t)^params['gR'])
	#fR = rep(1,tf/dt)
	fT =   params['hT']/(1+(params['kT']/t)^params['gT'])
	
	x0[1]=params['x0']
	s[1]=N-x0[1]

	Rx=rep(0,tf)
	Rr=rep(0,tf)
	Rs=rep(0,tf)
	
	counter=1
	
	for(i in 2:length(t)){
		x0[i] = x0[i-1]+dt*((params['lambda']/fA[i-1])*s[i-1]/N-params['lambdaR']*fR[i-1])*x0[i-1]
		r0[i] = r0[i-1]+dt*(                                    params['lambdaR']*fR[i-1])*x0[i-1]
		s[i]  = s[i-1] -dt*((params['lambda']/fA[i-1])*s[i-1]/N                          )*x0[i-1]
		x[i]  = x0[i]*fT[i]
		r[i]  = r0[i]*fT[i]
		if(t[i]+dt>counter) {
			Rx[counter] = x[i]
			Rr[counter] = r[i]
			Rs[counter] = s[i]
			counter=counter+1
		}
	}

	return(list(Rx=Rx,Rr=Rr,Rs=Rs))
}

cov.sim_rk <- function(params,tf,dt=.01,N=10000){
	names(params)<-c("lambda","lambdaR",
			 "hR","kR","gR",
			 "hA","kA","gA",
			 "hT","kT","gT",
			 "x0")

	t=seq(dt,tf,dt)

        Y=matrix(0,3,tf/dt)	
	x=rep(0,tf/dt)
	r=rep(0,tf/dt)
	fT=params['hT']/(1+(params['kT']/t)^params['gT'])
	
	Y[1,1]=params['x0']
	Y[3,1]=N-Y[1,1]

	Rx=rep(0,tf)
	Rr=rep(0,tf)
	
	counter=1
	
	for(i in 2:length(t)){
		k1=cov.evol(params,t[i]-dt,Y[,i-1],N)
		k2=cov.evol(params,t[i]-dt/2,Y[,i-1]+dt/2*k1,N)
		k3=cov.evol(params,t[i]-dt/2,Y[,i-1]+dt/2*k2,N)
		k4=cov.evol(params,t[i],Y[,i-1]+dt*k3,N)
		Y[,i] = Y[,i-1]+dt/6*(k1+2*k2+2*k3+k4)
		if(t[i]+dt>counter) {
			Rx[counter] = Y[1,i]*fT[i]
			Rr[counter] = Y[2,i]*fT[i]
			counter=counter+1
		}
	}

	return(list(Rx=Rx,Rr=Rr))
}

cov.SIM <- function(p,tf,dt=0.1){
	nr=nrow(p)
	t=seq(dt,tf,l=ceiling(tf/dt))
	
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
		#fR[,i] = 1+p[,'hR']/(1+(p[,'kR']/t[i])^p[,'gR'])
		fR[,i] = 1+p[,'hR']/(1+(p[,'kA']/t[i])^p[,'gR'])
		fT[,i] =   p[,'hT']/(1+(p[,'kT']/t[i])^p[,'gT'])
		
		x0[,i] = x0[,i-1]+dt*(p[,'lambda']/fA[,i-1]*s[,i-1]/popsize-p[,'lambdaR']*fR[,i-1])*x0[,i-1]
		s[,i]  = s[,i-1] -dt*(p[,'lambda']/fA[,i-1]*s[,i-1]/popsize                       )*x0[,i-1]
		r0[,i] = r0[,i-1]+dt*(                                      p[,'lambdaR']*fR[,i-1])*x0[,i-1]
		x[,i] = x0[,i]*fT[,i] 
		r[,i] = r0[,i]*fT[,i] 
		
		if(t[i]+dt>counter) {
			Rx[,counter] = pmax(1e-300,x[,i])
			Rr[,counter] = pmax(1e-300,r[,i])
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

cov.SIMP <- function(p,tf,dt=0.1){
	ncolPart=ncol(p)
	nrowData=nrow(data)
	npart=ncolPart/nrowData

	t=seq(dt,tf,l=ceiling(tf/dt))
	
	popsizeRep=rep(popsize,npart)

    x<-matrix(0,ncolPart,tf/dt)
	x0<-matrix(0,ncolPart,tf/dt)
	s<-matrix(0,ncolPart,tf/dt)
	r<-matrix(0,ncolPart,tf/dt)
	r0<-matrix(0,ncolPart,tf/dt)
	fA<-matrix(0,ncolPart,tf/dt)
	fR<-matrix(0,ncolPart,tf/dt)
	fT<-matrix(0,ncolPart,tf/dt)
	
	Rx<-matrix(0,ncolPart,tf)
	Rr<-matrix(0,ncolPart,tf)

	x0[,1]<-p['x0',]
	s[,1]<-popsizeRep-p['x0',]

	fA <- 1+p['hA',]/(1+outer(p['kA',],t,"/")^p['gA',])
	fR <- 1+p['hR',]/(1+outer(p['kA',],t,"/")^p['gR',])
	fT <-   p['hT',]/(1+outer(p['kT',],t,"/")^p['gT',])

	counter=1
	for(i in 2:length(t)){
	    U1=p['lambda',]/fA[,i-1]*s[,i-1]/popsizeRep*x0[,i-1]	
		U2=p['lambdaR',]*fR[,i-1]*x0[,i-1]

		x0[,i] = pmax(0,x0[,i-1]+dt*(U1-U2))
		s[,i]  = s[,i-1] -dt*U1
		r0[,i] = r0[,i-1]+dt*U2
		
		if(t[i]+dt>counter) {
			Rx[,counter] <- pmax(0,x0[,i]*fT[,i])
			Rr[,counter] <- pmax(0,r0[,i]*fT[,i])
			counter=counter+1
		}

		if(any(is.na(x0[,i]))) {
				ind=which(is.na(x0[,i]))[1]
				cat("na detected at time ",i,"\n")
						print(x0[ind,])
						print(p[,ind])
						exit()
		}
	}

	return(list(Rx=Rx,Rr=Rr))

}

cov.SIM_rk <- function(p,tf,dt=0.1){
	nr=nrow(p)
	t=seq(dt,tf,l=tf/dt)
	
	Y=array(0,dim=c(nr,3,tf/dt))
	fT=matrix(0,nr,tf/dt)
	
	Rx=matrix(0,nr,tf)
	Rr=matrix(0,nr,tf)
	
	Y[,1,1]=p[,'x0']
	fT = p[,'hT']/(1+(p[,'kT']/matrix(t,nr,tf/dt,byrow=1))^p[,'gT'])
	Y[,3,1]=popsize-p[,'x0']

	counter=1
	for(i in 2:length(t)){		
		
		k1=cov.EVOL(p,t[i]-dt,Y[,,i-1],popsize)
		k2=cov.EVOL(p,t[i]-dt/2,Y[,,i-1]+dt/2*k1,popsize)
		k3=cov.EVOL(p,t[i]-dt/2,Y[,,i-1]+dt/2*k2,popsize)
		k4=cov.EVOL(p,t[i],Y[,,i-1]+dt*k3,popsize)
		Y[,,i] = Y[,,i-1]+dt/6*(k1+2*k2+2*k3+k4)
		
		if(t[i]+dt>counter) {
			Rx[,counter] = Y[,1,i]*fT[,i]
			Rr[,counter] = Y[,2,i]*fT[,i]
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

cov.GlobalLogLikelihood_vec <- function(gparlist,paramlist,Al,Bl,A,B,tf,temperature=1,print=F){
	
	nr=nrow(data)
	p=cbind(lambda=gparlist[1],
		lambdaR=gparlist[2],
		paramlist)
	res=cov.SIM(p,tf,dt=.1)
	
	LogLikVec = temperature*rowSums(dpois(data[,1:tf],lambda=res$Rx,log=T)) + 
		    temperature*rowSums(dpois(datar[,1:tf],lambda=res$Rr,log=T)) +
		    rowSums(dgamma(paramlist,matrix(A,nr,length(A),byrow=1),matrix(B,nr,length(B),byrow=1),log=T))

	LogLik=sum(LogLikVec)+sum(dgamma(gparlist,Al,Bl,log=T))
	return(list(LogLik=LogLik,LogLikVec=LogLikVec))
}

cov.GlobalLogLikelihood_par <- function(particles,Al,Bl,A,B,tf,temperature=1){
	
	ncolPart=ncol(particles)
	nrowData=nrow(data)
	npart=ncolPart/nrowData

	res=cov.SIMP(particles,tf,dt=.1)

	if(any(is.na(res$Rx))){
		exit()
	}
	
	data_rep = data[rep(1:nrow(data),npart),]
	datar_rep = datar[rep(1:nrow(datar),npart),]

	LogLikVecMap = temperature*rowSums(dpois(matrix(data_rep[,1:tf],ncolPart,tf),lambda=matrix(res$Rx,ncolPart,tf),log=T)) + 
				   temperature*rowSums(dpois(matrix(datar_rep[,1:tf],ncolPart,tf),lambda=matrix(res$Rr,ncolPart,tf),log=T)) +
		           colSums(dgamma(particles[3:12,],matrix(A,length(A),ncolPart),matrix(B,length(B),ncolPart),log=T))

	if(any(is.na(datar_rep))){
			print(dim(datar_rep[,1:tf]))
			print(dim(res$Rr))
			cat("datar is NA\n")
	 		exit()
	}
	LogLikVec = matrix(LogLikVecMap,nrowData,npart)
	LogLik=colSums(LogLikVec)+
	       colSums(dgamma(particles[1:2,seq(1,ncolPart,nrowData)],
					  matrix(Al,2,npart),
					  matrix(Bl,2,npart),
					  log=T))
	
		   return(list(LogLik=LogLik,LogLikVec=LogLikVec))
}

cov.GlobalLogPrior_par <- function(particles,Al,Bl,A,B,temperature=1,print=F){
	
	ncolPart=ncol(particles)
	nrowData=nrow(data)
	npart=ncolPart/nrowData

	LogLikVecMap = colSums(dgamma(particles[3:12,],
								  matrix(A,length(A),ncolPart),
								  matrix(B,length(B),ncolPart),log=T))

	LogLikVec = matrix(LogLikVecMap,nrowData,npart)
	LogLik=colSums(LogLikVec)+
	       colSums(dgamma(particles[1:2,seq(1,ncolPart,nrowData)],
					  matrix(Al,2,npart),
					  matrix(Bl,2,npart),
					  log=T))
	
		   return(list(LogLik=LogLik,LogLikVec=LogLikVec))
}

cov.MCMC <- function(ind,NITER){
	A=c(1,1,
	    1,2,2, # R 
	    1,2,2, # A
	    1,1,2, # T
	    1) 
	B=c(1,10,
	    .01,.1,.2,
	    .01,.1,.2,
	    10,.1,.2,
	    .1)

	nvar=12
	tf=ncol(data)

	prop_size=10

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

cov.GlobalMCMC <- function(NITER,show=1,init=NA,init_labels=NA,init_temp=0,verbose=F){
	Al=c(1,1)
	Bl=c(1,10)

	A=c(1,2,2, # R 
	    1,2,2, # A
	    1,1,2, # T
	    1)     # x0
	B=c(1,.1,.2,
	    0.01,.1,.2,
	    10,  .1,.2,
	    .1)


	nvar=12
	nc=nrow(data)
	tf=ncol(data)

	temp=init_temp
	prop_size_lam=100
	prop_size=matrix(50,nc,nvar-2) #* rowSums(data)

	params=matrix(NA,nc,nvar-2)
	colnames(params)=c('hR','kR','gR',
                           'hA','kA','gA',
                           'hT','kT','gT',
                           'x0')

	## init using single MCMC

	lambda=0.1
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

	while(any(is.na(GLL$LogLikVec) | is.infinite(GLL$LogLikVec))){
		lambda=0.1
		lambdaR=rgamma(1,Al[2],Bl[2])
		for(i in 1:nc){ 
			params[i,]=rgamma(nvar-2,A,B)
		}
	        	
		GLL=cov.GlobalLogLikelihood_vec(c(lambda,lambdaR),params,Al,Bl,A,B,tf,temperature=temp)
	}

	for(i in 1:NITER){
		
		
		cat(i,' ',GLL$LogLik,' ',temp,"     \r")	
		
		## first update lambda
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambda_new=fac*lambda
			prop_ind=1
		} else {
			lambda_new=rgamma(1,Al,Bl)
			prop_ind=2
		}

		temp=init_temp+(1-init_temp)/(1+((NITER/2)/i)^7)
		## recalculate the LogLikelihood from previous iteration
		#GLL$LogLik=sum(dgamma(c(lambda,lambdaR),Al,Bl,log=T))+sum(GLL$LogLikVec)
		GLL=cov.GlobalLogLikelihood_vec(c(lambda,lambdaR),params,Al,Bl,A,B,tf,temperature=temp)
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
		if(i>0 && i%%10==0 && show){
			if(1){
				cov.plotFromGlobal(show,params_samples[max(1,i-200):i],ncol(data)+100,10,ymax=max(data[show,]+50))
			} else {


			par(mar=c(0,3,0,0))
			layout(matrix(1:7,7))
				plot(LL_samples[max(1,i-200):i],type='l',xaxt='n',xlab="")
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$lambda)),type='l',xaxt='n',xlab="")
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$lambdaR)),type='l',xaxt='n',xlab="")
		       	for(l in 1:4) 
				plot(unlist(lapply(params_samples[max(1,i-200):i],function(x) x$param[show,l])),type='l',xaxt='n',xlab="")
			}
		}
	}

	return(params_samples)
}

cov.GlobalSMC.moveParticle <- function(p,GlobalLogLikelihood,Al,Bl,A,B,prop_size_lam,prop_size){

	    nc=nrow(data)
		tf=ncol(data)
		particle=p
		GLL=GlobalLogLikelihood

		## first update lambda
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambda_new=fac*particle$lambda
		} else {
			lambda_new=rgamma(1,Al,Bl)
		}

		GLL_new=cov.GlobalLogLikelihood_vec(c(lambda_new,particle$lambdaR),particle$params,Al,Bl,A,B,tf,temperature=1)
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			   cov.MixedGamma(particle$lambda,Al[1],Bl[1],particle$lambda/lambda_new,prop_size_lam,prop_size_lam)-
			   cov.MixedGamma(lambda_new,Al[1],Bl[1],lambda_new/particle$lambda,prop_size_lam,prop_size_lam)


		## and update lambda
		u=runif(1)

	    if(!is.na(logalpha)){	
		if(log(u)<logalpha){
			particle$lambda=lambda_new
			GLL=GLL_new
		}
		}
	
		##  update lambdaR
		if(runif(1)<0.5){
			fac=rgamma(1,prop_size_lam,prop_size_lam)
			lambdaR_new=fac*particle$lambdaR
		} else {
			lambdaR_new=rgamma(1,Al[2],Bl[2])
		}

		GLL_new=cov.GlobalLogLikelihood_vec(c(particle$lambda,lambdaR_new),particle$params,Al,Bl,A,B,tf,temperature=1)
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			   cov.MixedGamma(particle$lambdaR,Al[2],Bl[2],particle$lambdaR/lambdaR_new,prop_size_lam,prop_size_lam)-
			   cov.MixedGamma(lambdaR_new,Al[2],Bl[2],lambdaR_new/particle$lambdaR,prop_size_lam,prop_size_lam)


		## and update lambdaR
		u=runif(1)

	    if(!is.na(logalpha)){	
		if(log(u)<logalpha){
			particle$lambdaR=lambdaR_new
			GLL=GLL_new
		}
		}

		## then update the other parameters independently

		for(var in 1:ncol(particle$params)){
			which_from_prior=runif(nc)<0.5
			params_new=particle$params


			fac=rgamma(nc,prop_size[,var],prop_size[,var])
			rprior=rgamma(nc,A[var],B[var])			

			params_new[!which_from_prior,var]=fac[!which_from_prior]*particle$params[!which_from_prior,var]
			params_new[which_from_prior,var]=rprior[which_from_prior]

			GLL_new = cov.GlobalLogLikelihood_vec(c(particle$lambda,particle$lambdaR),params_new,Al,Bl,A,B,tf,temperature=1) 
			logalpha = GLL_new$LogLikVec - GLL$LogLikVec + 
				   cov.MixedGammaVec(particle$params[,var],rep(A[var],nc),rep(B[var],nc),particle$params[,var]/params_new[,var],prop_size[,var],prop_size[,var]) - 
				   cov.MixedGammaVec(params_new[,var],rep(A[var],nc),rep(B[var],nc),params_new[,var]/particle$params[,var],prop_size[,var],prop_size[,var])

			u=runif(nc)
	        if(all(!is.na(logalpha))){	
				selec=log(u)<logalpha
				particle$params[selec,var]=params_new[selec,var]
				GLL$LogLikVec[selec]=GLL_new$LogLikVec[selec]
			}
			
		}

		GLL=cov.GlobalLogLikelihood_vec(c(particle$lambda,particle$lambdaR),particle$params,Al,Bl,A,B,tf,temperature=1)
		
		return(list(p=particle,GLL=GLL))
}

cov.GlobalSMC.moveAllParticles <- function(p,GlobalLogLikelihood,Al,Bl,A,B,prop_size_lam,prop_size,tf,temp=1){

	    nrowData=nrow(data)
		ncolPart=ncol(p)
		npart=ncolPart/nrowData
		particles=p
		GLL=GlobalLogLikelihood

		## first update lambda
		fac=rgamma(npart,prop_size_lam,prop_size_lam)
		lambda_new=rep(fac,each=nrowData)*particles['lambda',]
		pnew=particles
		pnew['lambda',]=lambda_new

		GLL_new=cov.GlobalLogLikelihood_par(pnew,Al,Bl,A,B,tf,temperature=temp)
		if(any(is.na(GLL_new$LogLik))) {
				print(GLL_new$LogLik)
				cat("problem while sampling lambda\n")
		}
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			       dgamma(1/fac,prop_size_lam,prop_size_lam,log=T)-
			       dgamma(fac,prop_size_lam,prop_size_lam,log=T)


		## and update lambda
		u=runif(npart)

		selection=(log(u)<logalpha & !is.na(logalpha))
		particles['lambda',selection]=pnew['lambda',selection]
		GLL$LogLik[selection]=GLL_new$LogLik[selection]
		GLL$LogLikVec[,selection]=GLL_new$LogLikVec[,selection]

        ## second update lambdaR
		fac=rgamma(npart,prop_size_lam,prop_size_lam)
		lambdaR_new=rep(fac,each=nrowData)*particles['lambdaR',]
		pnew=particles
		pnew['lambdaR',]=lambdaR_new

		GLL_new=cov.GlobalLogLikelihood_par(pnew,Al,Bl,A,B,tf,temperature=temp)
		if(any(is.na(GLL_new$LogLik))) {
				print(GLL_new$LogLik)
				cat("problem while sampling lambdaR\n")
		}
		
		logalpha = GLL_new$LogLik- GLL$LogLik + 
			       dgamma(1/fac,prop_size_lam,prop_size_lam,log=T)-
			       dgamma(fac,prop_size_lam,prop_size_lam,log=T)


		## and update lambda
		u=runif(npart)

		selection=(log(u)<logalpha & !is.na(logalpha))
		particles['lambdaR',selection]=pnew['lambdaR',selection]
		GLL$LogLik[selection]=GLL_new$LogLik[selection]
		GLL$LogLikVec[,selection]=GLL_new$LogLikVec[,selection]

		## then update the other parameters independently

		for(var in 1:10){
			pnew=particles

			fac=rgamma(ncolPart,prop_size,prop_size)
			fac.asMat=matrix(fac,nrowData,npart)

			pnew[2+var,]=fac*particles[2+var,]

			GLL_new = cov.GlobalLogLikelihood_par(pnew,Al,Bl,A,B,tf,temperature=temp) 
		if(any(is.na(GLL_new$LogLik))) {
				print(GLL_new$LogLik)
				cat("problem while sampling var ",var,"\n")
		}
			logalpha = GLL_new$LogLikVec - GLL$LogLikVec + 
				       dgamma(1/fac.asMat,prop_size,prop_size,log=T) - 
				       dgamma(fac.asMat,prop_size,prop_size,log=T)  

			u=matrix(runif(ncolPart),nrowData,npart)
			selection=(log(u)<logalpha & !is.na(logalpha))
			selection.flat=as.vector(log(u)<logalpha & !is.na(logalpha))
			particles[2+var,selection.flat]=pnew[2+var,selection.flat]
			GLL$LogLikVec[selection]=GLL_new$LogLikVec[selection]
		}
			
	
	GLL=cov.GlobalLogLikelihood_par(particles,Al,Bl,A,B,tf,temperature=temp)
		
	return(list(p=particles,GLL=GLL))
}

cov.GlobalSMC <- function(NPART,STEPS){
	Al=c(1,1)
	Bl=c(1,10)

	A=c(1,2,2, # R 
	    1,2,2, # A
	    1,1,2, # T
	    1)     # x0
	B=c(0.1,.1,.2,
	    0.01,.1,.2,
	    10,  .1,.2,
	    .1)

	nvar=12
	nrowData=nrow(data)
	tf=ncol(data)
	dt=0.1
	ncolPart=nrowData*NPART

	temp=1
	prop_size_lam=100
	prop_size=10 #* rowSums(data)

	## init using single MCMC

	allParticles=matrix(NA,nvar,ncolPart)
	
## =======================================
## Initialize memory list
## =======================
#	memList$x<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$x0<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$s<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$r<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$r0<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$fA<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$fR<<-matrix(0,ncolPart,ceiling(tf/dt))
#	memList$fT<<-matrix(0,ncolPart,ceiling(tf/dt))
#	
#	memList$Rx<<-matrix(0,ncolPart,tf)
#	memList$Rr<<-matrix(0,ncolPart,tf)
## =======================================

	cat("Initializing particles...\n")
	for(k in 1:NPART){

		allParticles[1:2,1:nrowData+(k-1)*nrowData]=matrix(rgamma(2,Al,Bl),nrowData,2,byrow=1)

		for(i in 1:nrowData){
				allParticles[3:12,i+(k-1)*nrowData]=rgamma(10,A,B)
		}
	}
	cat("done\n")

    LogML=0
	protocol=seq(0,1,l=STEPS)
	
	W = rep(1.0/NPART,NPART)
	ESS = rep(NA,STEPS)
    logW=-log(NPART)

    log_w_inc = rep(NA,NPART)

	for(i in 1:STEPS){
			cat(i,"     \r");
			ESS[i]=1.0/sum(W^2)
	
            #if ESS<npart/2.0:
			if(1){
				allParticles.sample=sample(1:NPART,NPART,replace=T,prob=W)
                dim(allParticles)=c(nvar,nrowData,NPART)
				allParticles=allParticles[,,allParticles.sample]
				dim(allParticles)=c(nvar,ncolPart)

				rownames(allParticles)=c('lambda','lambdaR',
					   'hR','kR','gR',
                       'hA','kA','gA',
                       'hT','kT','gT',
                       'x0')

            	W=rep(1/NPART,NPART)
            	logW=-log(NPART)
			}

        	delta=protocol[i+1]-protocol[i]
            
			GLL=cov.GlobalLogLikelihood_par(
						        allParticles,
								Al,Bl,A,B,tf,temperature=temp)

#=============================================================
# Recalculate weights
#====================
			log_w_inc = delta*GLL$LogLik
			
        	log_w_un = log_w_inc + logW
			
			lwun_max=max(log_w_un)
			
			W = exp(log_w_un-lwun_max)
			W = W/sum(W)

            logML_inc = lwun_max+log(sum(exp(log_w_un-lwun_max)))
			logW=log_w_un-logML_inc
			
			LogML=LogML+logML_inc
#=============================================================


            SMCmove=cov.GlobalSMC.moveAllParticles(allParticles,GLL,Al,Bl,A,B,prop_size_lam,prop_size)
			
			allParticles=SMCmove$p

	}

    return(list(ESS,allParticles,W,LogML))

}

cov.Global.IBIS <- function(NPART,tf){
	Al=c(1,1)
	Bl=c(10,10)

	A=c(1,2,2, # R 
	    1,2,2, # A
	    1,1,2, # T
	    1)     # x0
	B=c(0.1,.1,.2,
	    0.01,.1,.2,
	    10,  .1,.2,
	    .1)

	nvar=12
	nrowData=nrow(data)
	dt=0.1
	ncolPart=nrowData*NPART

	temp=1
	prop_size_lam=100
	prop_size=100 #* rowSums(data)

	## init using single MCMC

	allParticles=matrix(NA,nvar,ncolPart)

	cat("Initializing particles...\n")
	for(k in 1:NPART){

		allParticles[1:2,1:nrowData+(k-1)*nrowData]=matrix(rgamma(2,Al,Bl),2,nrowData)

		for(i in 1:nrowData){
				allParticles[3:12,i+(k-1)*nrowData]=rgamma(10,A,B)
		}
	}
	cat("done\n")

	GLL=cov.GlobalLogPrior_par(
				        allParticles,
						Al,Bl,A,B,temperature=temp)

    LogML=0
	
	W = rep(1.0/NPART,NPART)
	ESS = rep(NA,tf)
	ESS[1]=1.0/sum(W^2)
    logW=-log(NPART)

    log_w_inc = rep(NA,NPART)

	for(time in 1:tf){
			cat(time,"     \r");
			ESS[time]=1.0/sum(W^2)
	
            #if ESS<npart/2.0:
			if(1){
				allParticles.sample=sample(1:NPART,NPART,replace=T,prob=W)
                dim(allParticles)=c(nvar,nrowData,NPART)
				allParticles=allParticles[,,allParticles.sample]
				dim(allParticles)=c(nvar,ncolPart)
				GLL$LogLikVec=GLL$LogLikVec[,allParticles.sample]
				GLL$LogLik=GLL$LogLik[allParticles.sample]

				rownames(allParticles)=c('lambda','lambdaR',
					   'hR','kR','gR',
                       'hA','kA','gA',
                       'hT','kT','gT',
                       'x0')

            	W=rep(1/NPART,NPART)
            	logW=-log(NPART)
			}

            SMCmove=cov.GlobalSMC.moveAllParticles(allParticles,GLL,Al,Bl,A,B,prop_size_lam,prop_size,time,temp)
			
			allParticles=SMCmove$p
			GLL_new=SMCmove$GLL

#=============================================================
# Recalculate weights
#====================
			log_w_inc = GLL_new$LogLik-GLL$LogLik
			if(any(is.na(log_w_inc))) {
					print(cbind(GLL_new$LogLik,GLL$LogLik))
							stop()
			}
			
        	log_w_un = log_w_inc + logW
			
			lwun_max=max(log_w_un)
			
			W = exp(log_w_un-lwun_max)
			W = W/sum(W)

            logML_inc = lwun_max+log(sum(exp(log_w_un-lwun_max)))
			logW=log_w_un-logML_inc
			
			LogML=LogML+logML_inc
#=============================================================

			GLL=GLL_new

	}

    return(list(ESS,allParticles,W,LogML))

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
	       	lines(res$Rx,col="red",lty=2)
	       	lines(res$Rr,col="green",lty=2)
	}
}

cov.plotFromGlobal <- function(ind,samples,tf,n=10,ymax=NA,log=""){
	Nsam=length(samples)
	randsam=sample((Nsam-10+1):Nsam,n)
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
		abline(v=params[c('kA','kT')],col=c("red","grey"),lty=2)
	}
}

cov.plotFromSMC <- function(ind,samples,tf,n=10,ymax=NA,log=""){
	nrowData=nrow(data)
	npart=ncol(samples)/nrowData
	randsam=1:n

	params=samples[,ind]
	
	res=cov.sim(params,tf,N=popsize[ind])
	if(is.na(ymax)) ymax=max(c(res$Rx,res$Rr))
     
	plot(data[ind,],xlim=c(1,tf),main=country[ind,],ylim=c(0,ymax),pch=19)
	points(datar[ind,],xlim=c(1,tf),col="green",pch=19)
	for(i in 1:n){
		params=samples[,ind+(n-1)*nrowData]
		res=cov.sim(params,tf,N=popsize[ind])
	       	lines(res$Rx,col=rgb(1,0,0,.5),lty=2)
	       	lines(res$Rr,col=rgb(0,1,0,.5),lty=2)
	       	
		abline(v=params[c('kA','kT')],col=c("red","grey"),lty=2)
	}
}

cov.plotHill <- function(h,k,g,tf,one=T){
	t=1:tf
	plot(1:tf,one+h/(1+(k/1:tf)^g))
}



