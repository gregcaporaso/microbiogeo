
rm(list = ls())
setwd("/Users/antoniog/Desktop/Biogeo/FullTreeUniFracVariograms")
require(geoR)
require(ncf)
distmat=read.table('distmat.txt',header=F,sep="\t")
ufdistmat=read.table('unifracdistmat.txt',header=F,sep="\t")
cols=colnames(datamod[1:ncol(datamod)])

distmatv=distmat[1,1]
for(i in 2:ncol(distmat)){
	distmatv=c(distmatv,distmat[1:i,i])
}

ufdistmatv=0 #ufdistmat[1,1]
for(j in 2:ncol(ufdistmat)){
	ufdistmatv=c(ufdistmatv,ufdistmat[1:j,j])
}
j=length(ufdistmatv)
ufdistmatv=ufdistmatv[2:j-1]
j=length(ufdistmatv)
k=1
ufdistm=ufdistmatv[1]
for(k in 2:j)
	if(ufdistmatv[k]!=0)
		ufdistm=c(ufdistm,ufdistmatv[k])

lags=c(seq(0,10,2.5),seq(20,400,10),seq(450,1000,50))
D=distmatv

gamma=list(0)
for (c in 2:length(lags))
	gamma=c(gamma,list(0))


weights=rep(0,length(lags));
weights[1]=ncol(distmat)

for (k in 2:length(lags)){
    minD = lags[k-1];
    maxD = lags[k];
    for(a in 1:(length(D))){
	 if(D[a]>minD && D[a]<=maxD)
		if(gamma[[k]][1]>0){
			gamma[[k]] =c(gamma[[k]],ufdistmatv[a])
			weights[k]=weights[k]+1
		}else {gamma[[k]]=ufdistmatv[a]
			weights[k]=weights[k]+1}
	}
}
v=mean(gamma[[1]])
std=0
for(b in 2:length(lags)){
	v=c(v,mean(gamma[[b]]))
	std=c(std,sd(gamma[[b]]))
}

xlags=c(0,seq(1.25,9.75, 2.5),seq(15,400,10),seq(425,975,50))
rnum=nrow(data)

hvarmodel=variofit(hvar,cov.model="spherical",ini.cov.pars=c(0.03,240),nugget=0.2101)

plot(xlags,v,ylab="Community Dissimilarity",xlab="Distance(m)",xlim=c(0,400),ylim=c(0.2,0.26),lwd=2,font.axis=2,font.lab=2)
#plot(xlags,v,ylab="Community Dissimilarity",xlab="Distance(m)",lwd=2,font.axis=2,font.lab=2)

#lines.variomodel.variofit(hvarmodel,lwd=2)
#range=array(as.numeric(hvarmodel$cov.pars),dim=c(1,2))
#print(abline(v=range[1,2],lty=2,lwd=2))
#print(abline(h=hvarmodel$nugget,lty=2,lwd=2))
title(main="Weighted Unifrac Distance Semivariogram")

##
## weights
## [1] 85  0  8  9 11 35  8 24 52 49 44 39 58 71 53 64 42 46 53 43 36 45 46 39 45 43 42 26 29 26 25 25 17 30 22 22 30 30 18 22 32 31 21 25 81 58 39 21 32 31 39 34 44 52 48 26
##
##