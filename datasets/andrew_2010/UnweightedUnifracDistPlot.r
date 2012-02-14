require(geoR)
distmat=read.table('distmat.txt',header=F,sep="\t")
ufdistmat=read.table('Uunifracdistmat.txt',header=F,sep="\t")
distmatv=distmat[1,1]
for(i in 2:ncol(distmat)){
	distmatv=c(distmatv,distmat[1:i,i])
}
ufdistmatv=ufdistmat[1,1]
for(j in 2:ncol(ufdistmat)){
	ufdistmatv=c(ufdistmatv,ufdistmat[1:j,j])
}

#win.graph()
#plot(distmatv,ufdistmatv,xlim=c(0,500))

#lags=c(seq(0,50,5),seq(60,300,10))   
#lags=c(seq(0,50,5),seq(60,300,10),seq(350,1000,50))
lags=c(seq(0,1000,5))
D=distmatv

#gamma = rep(0,length(lags));
gamma=list(0)
for (c in 2:length(lags))
	gamma=c(gamma,list(0))


weights=rep(0,length(lags));
weights[1]=85

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

#xlags=c(0,seq(2.5,47.5, 5),seq(55,295,10))
#xlags=c(0,seq(2.5,47.5, 5),seq(55,295,10),seq(325,975,50))
xlags=c(0,seq(2.5,997.5, 5))
data=read.table('Extract_FullTal1.txt',header=T,sep="\t")

rnum=nrow(data)

#H2O

rm(hmat)
hmat=c(data[1,2:3],data[1,15])

for(f in 2:rnum){
	if(!is.na(data[f,15]))
		hmat=rbind(hmat,c(data[f,2:3],data[f,15]))

	}
r=nrow(hmat)
hmat=array(as.numeric(hmat),dim=c(r,3))
hcoord=array(hmat[1:r,1:2])

hmat=as.geodata(hmat,coords.col=1:2,data.col=3)
bins=lags

hvar=variog(hmat,breaks=bins,option="bin")

hvar$v=v[2:length(v)]
hvar$n=weights[2:length(weights)]
hvar$sd=std[2:length(std)]
hvar$var.mark=var(ufdistmatv)
hvar$beta.ols=10
hvar$n.data=85


hvarmodel=variofit(hvar,cov.model="spherical",ini.cov.pars=c(0.10,300),fix.nugget=TRUE,nugget=0.74)
plot(xlags,v,ylab="Unifrac Distance",xlab="Distance(m)",ylim=c(0.73,0.83))
lines.variomodel.variofit(hvarmodel)
range=array(as.numeric(hvarmodel$cov.pars),dim=c(1,2))
print(abline(v=range[1,2],lty=2))
print(abline(h=0.74,lty=2))
title(main="Unweighted Unifrac Distance Semivariogram")




