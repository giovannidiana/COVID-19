#load("samples_sir_17mar.RData")

png("../Figures/figpost_sir_henan.png",600,600)

ind=grep("Henan",labels)[1]
p=c(lambda=sam_glob_sir[[10000]]$lambda,
    lambdaR=sam_glob_sir[[10000]]$lambdaR,
    sam_glob_sir[[10000]]$param[ind,])

tf=100
plot(data[ind,],xlim=c(0,tf),ylim=c(0,ymax),pch=19,xlab="days",log="")

res=cov.sim(p,tf,0.1,N=popsize[ind])
lines(res$Rx,col='orange')

counter=0
NL=6
line_name=rep("",NL)
for(k in seq(0,2,l=NL)){
	counter=counter+1
	line_name[counter]=paste("factor =",format(k,digits=3))
	p_no_int=p
	p_no_int['hA']=k*p['hA']
	p_no_int['hR']=k*p['hR']
	res=cov.sim(p_no_int,tf,0.1,N=popsize[ind])
	lines(res$Rx,col='red',lty=counter)
}

legend("topright",legend=c(line_name,"estimated","data"),
       col=c(rep("red",NL),"orange","black"),lty=c(1:(NL),1,NA),
       pch=c(rep(NA,NL+1),19))

dev.off()
