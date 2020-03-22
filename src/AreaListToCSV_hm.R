load("AreaList_HM-2203.RData")

tf=ncol(AreaList$data)+300
time_axis=format(as.Date(seq(1,tf,1),origin="21 Jan 2020",format="%d %b %Y"),"%d %b %Y")

rownames(AreaList$q1)=paste(labels,"_qi1",sep="")
rownames(AreaList$q2)=paste(labels,"_qi2",sep="")
rownames(AreaList$q3)=paste(labels,"_qi3",sep="")

dataPad=matrix(NA,length(labels),ncol(AreaList$q1))
dataPad[1:nrow(AreaList$data),1:ncol(AreaList$data)]=AreaList$data

AreaList2=AreaList
AreaList2$data=dataPad

rownames(AreaList2$data)=paste(labels,"_INF",sep="")

AreaListCsv = cbind(t=time_axis,t(do.call(rbind,AreaList2[c('data','q1','q2','q3')])))

write.table(AreaListCsv,file="AreaListCsv_hm_2203.csv",row.names=F,sep=",")

