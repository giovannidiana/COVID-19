library(plotly)

cov.getArea <- function(ind,samples,tf,n=10){
	Nsam=length(samples)
	traces=matrix(NA,n,tf)
	counter=1
	for(i in (Nsam-n+1):Nsam){
		params=c(samples[[i]]$lambda,samples[[i]]$param[ind,])
                names(params)<-c("lambda","q","h","x0","k")
                res=cov.sim(params,tf)
                pp=res[[1]]/(res[[1]]+params['k'])
                traces[counter,]=res[[1]]*pp
		counter=counter+1
        }
	return(list(data=data[ind,],
	    	    mins=apply(traces,2,min),
	    	    maxes=apply(traces,2,max)))
}

cov.getAllAreas <- function(samples,tf,n){
	AreaList<-vector("list",nrow(data))
	for(i in 1:nrow(data)) AreaList[[i]] = cov.getArea(i,samples,tf,n)
	return(AreaList)
}

AreaList=cov.getAllAreas(sam_glob,ncol(data)+50,10000)
fig<-plot_ly()

# function to set the visible
setvisible <- function(i,n){
	bv=rep(FALSE,n)
	bv[i]=TRUE

	return(rep(bv,each=3))
}

for(i in 1:nrow(data)){
	fig <- fig %>% add_trace(x=1:nrow(data),y=data[i,],type="scatter", mode="markers",color="black") %>%
		add_lines(x=1:nrow(data),y=AreaList[[i]]$mins,fill='tonexty',fillcolor="green") %>%
		add_lines(x=1:nrow(data),y=AreaList[[i]]$maxes)
}

buttons = vector("list",nc)
for(i in 1:nrow(data)){
	buttons[[i]]=list(method = "restyle",
      			  args = list("visible",setvisible(i)),
			  label = i)
}

fig <- fig %>% layout(
    title = "Drop down menus - Plot type",
    xaxis = list(domain = c(1, 10)),
    updatemenus = list(
      list(
        y = 0.8,
	buttons=buttons       
	)
    ))

fig

#htmlwidgets::saveWidget(fig, "/home/giovanni/plotly_chart.html")

