load("AreaList_sir-2103.RData")

AreaList$datar=datar

tf=ncol(AreaList$data)+300
time_axis=format(as.Date(seq(1,tf,1),origin="21 Jan 2020",format="%d %b %Y"),"%d %b %Y")

rownames(AreaList$q1)=paste(labels,"_qi1",sep="")
rownames(AreaList$q2)=paste(labels,"_qi2",sep="")
rownames(AreaList$q3)=paste(labels,"_qi3",sep="")
rownames(AreaList$qr1)=paste(labels,"_qr1",sep="")
rownames(AreaList$qr2)=paste(labels,"_qr2",sep="")
rownames(AreaList$qr3)=paste(labels,"_qr3",sep="")

dataPad=matrix(NA,length(labels),ncol(AreaList$q1))
dataPadR=matrix(NA,length(labels),ncol(AreaList$q1))
dataPad[1:nrow(AreaList$data),1:ncol(AreaList$data)]=AreaList$data
dataPadR[1:nrow(AreaList$datar),1:ncol(AreaList$datar)]=AreaList$datar


AreaList2=AreaList
AreaList2$data=dataPad
AreaList2$datar=dataPadR

rownames(AreaList2$data)=paste(labels,"_INF",sep="")
rownames(AreaList2$datar)=paste(labels,"_REC",sep="")

AreaListCsv = cbind(t=time_axis,t(do.call(rbind,AreaList2[c('data','datar','q1','q2','q3','qr1','qr2','qr3')])))


