# Figure 1

png("../Figures/Figure_1.png",width=400,height=400)
plot(cov.sim(c(.1,.1,.01,1,.01),100,.1)[[1]],type='l',lwd=2,xlab="time (days)", ylab="Infected cases")
grid()
dev.off()

# Figure stat

png("../Figures/Figure_stat_1.png",width=600,height=400)
country_index=12
barplot(rbind(diff(rawdata_conf[country_index,]),
	      diff(rawdata_rec[country_index,]), 
	      diff(rawdata_deaths[country_index,])),
	col=c("red","green","purple" ),las=2,names.arg=names(rawdata_conf)[-1])
legend("topright",legend=c("confirmed","recovered","deaths"),fill=c("red","green","purple" ))

dev.off()

# Figure param lambda

par_comb=vector("list",5)
labels=rep("",nrow(data))
for(j in 1:nrow(data)) labels[j]<-paste(as.character(country[j,1]),as.character(province[j,1]))

par_comb[[1]]=par_comb[[1]]=unlist(lapply(sam_glob3,"[[",'lambda'))
for(i in 2:5){
       	par_comb[[i]]=t(do.call(cbind,lapply(sam_glob3,function(x) x$param[,i-1])))
	colnames(par_comb[[i]])=labels
}

png("../Figures/Figure_stat_lambda.png",width=600,height=400)
hist(unlist(lapply(sam_glob3,function(x) x$lambda)),
     main="",
     xlab="worldwide infection rate per day",
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
