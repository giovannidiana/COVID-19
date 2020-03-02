# Figure 1

png("Figures/Figure_1.png",width=400,height=400)
plot(cov.sim(c(.1,.1,.01,1),100,.1)[[1]],type='l',lwd=2,xlab="time (days)", ylab="Infected cases")
grid()
dev.off()

# Figure stat

png("Figures/Figure_stat_1.png",width=600,height=400)
country_index=12
barplot(rbind(diff(data_conf[country_index,]),
	      diff(data_rec[country_index,]), 
	      diff(data_death[country_index,])),
	col=c("red","green","purple" ),las=2,names.arg=substring(names(data_ts)[5:ncol(data_ts)],2)[-1])
legend("topright",legend=c("confirmed","recovered","deaths"),fill=c("red","green","purple" ))

dev.off()
