## This script should be ran right after generating the samples so that by loading the data through cov_sampler_SIR.R we have the same labels,data and datar that have been used to generate the posterior samples.

#source("../src/cov_sampler_SIR.R")
#load("../src/samples_sir_23mar.RData")
samples=sam_glob_sir[9000:10000]

cov.getAllAreas <- function(samples,tf,n){
        Nsam=length(samples)
        traces=array(NA,dim=c(nrow(data),tf,n))
        tracesR=array(NA,dim=c(nrow(data),tf,n))
        counter=1
        for(i in (Nsam-n+1):Nsam){
            s=samples[[i]]
            p=data.frame(lambda=s$lambda,
                         lambdaR=s$lambdaR,
                         s$param
                     )
              
                res=cov.SIM(p,tf)
                traces[,,counter]=res$Rx
                tracesR[,,counter]=res$Rr
                counter=counter+1
        }
        return(list(data=data,
                    q1=apply(traces,c(1,2),quantile,prob=0.05),
                    q2=apply(traces,c(1,2),quantile,prob=0.5),
                    q3=apply(traces,c(1,2),quantile,prob=0.95),
                    
                    qr1=apply(tracesR,c(1,2),quantile,prob=0.05),
                    qr2=apply(tracesR,c(1,2),quantile,prob=0.5),
                    qr3=apply(tracesR,c(1,2),quantile,prob=0.95),
                    
                    peaks=apply(
                              apply(traces,c(1,3),max),1,quantile,probs=c(0.05,.5,0.95) 
                           ),
                    
                    peaksR=apply(
                              apply(tracesR,c(1,3),max),1,quantile,probs=c(0.05,.5,0.95) 
                           )
                    
                    
                   ))
                    
}


## Collecting a posterior samples. This can take a few minutes
AreaList=cov.getAllAreas(samples,ncol(data)+300,1000)

#load(file="AreaList_sir-2303.RData")

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

write.table(AreaListCsv,file="AreaListCsv_sir.csv",row.names=F,sep=",")
