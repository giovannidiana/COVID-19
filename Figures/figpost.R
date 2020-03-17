load("samples_13mar.RData")

p=rep(NA,5)
ind=grep("Henan",labels)[1]
p[1]=mean(unlist(lapply(sam_globX[c(1000:2000)],"[[",'lambda')))
for(i in 1:4) p[i+1]=mean(unlist(lapply(sam_globX[c(1000:2000)],function(x) x$param[ind,i])))
names(p)=c('lambda','h','x0','k','g')

tf=60
#ymax=2000
ymax=2000
plot(data[ind,],xlim=c(0,tf),ylim=c(0,ymax),pch=19,xlab="days")

res=cov.sim(p,tf,0.01)
pp=1/(1+(p['k']/res$Rx)^p['g'])
lines(res$Rx*pp,col='orange')

counter=0
NL=6
line_name=rep("",NL)
for(k in seq(0,2*p['h'],l=NL)){
	counter=counter+1
	line_name[counter]=paste("h =",format(k,digits=4))
	p_no_int=p
	p_no_int['h']=k
	res=cov.sim(p_no_int,tf,0.01)
	pp=1/(1+(p['k']/res$Rx)^p['g'])
	lines(res$Rx*pp,col='red',lty=counter)
}

legend("topleft",legend=c(line_name,"estimated","data"),
       col=c(rep("red",NL),"orange","black"),lty=c(1:(NL),1,NA),
       pch=c(rep(NA,NL+1),19))



