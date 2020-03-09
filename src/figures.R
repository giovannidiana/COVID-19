load(file="samples_9mar.RData")
samples=sam_globX


# Figure 1

if(FALSE){
png("../Figures/Figure_1.png",width=400,height=400)
plot(cov.sim(c(.1,.1,.01,1,.01),100,.1)[[1]],type='l',lwd=2,xlab="time (days)", ylab="Infected cases")
grid()
dev.off()

# Figure stat

png("../Figures/Figure_stat_1.png",width=600,height=400)
country_index=which(province[,1]=="Henan")
barplot(rbind((data.I[country_index,]),
	      (data.R[country_index,]), 
	      (data.D[country_index,])),
	col=c("red","green","purple" ),las=2,names.arg=names(rawdata_conf))
legend("topright",legend=c("infected","recovered","deaths"),fill=c("red","green","purple" ))

dev.off()
}
# Figure param lambda

par_comb=vector("list",5)
labels=rep("",nrow(data))

with_province=!is.na(province[,1])
for(i in 1:nrow(data)){
    if(with_province[i]){
        labels[i]=paste(province[i,],",",country[i,],sep="")
    } else {        
        labels[i]=paste(country[i,],sep="")
    } 
}


par_comb[[1]]=par_comb[[1]]=unlist(lapply(samples,"[[",'lambda'))
for(i in 2:5){
       	par_comb[[i]]=t(do.call(cbind,lapply(samples,function(x) x$param[,i-1])))
	colnames(par_comb[[i]])=labels
}

png("../Figures/Figure_stat_lambda.png",width=600,height=400)
hist(unlist(lapply(samples,function(x) x$lambda)),
     main="",
     xlab="worldwide daily infection rate",
     ylab="frequency",freq=0)
dev.off()


for(i in 2:5){
	png(paste("../Figures/Figure_stat_",i,".png",sep=""),width=700,height=3000)
	means=colMeans(par_comb[[i]])
	ord=order(means,decreasing=T)
	par(mar=c(3,15,3,3))
	boxplot(par_comb[[i]][,ord], horizontal=T,las=2,outline=F)
	dev.off()
}
