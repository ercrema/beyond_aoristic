library(baorista)
library(nimbleCarbon)
library(emdbook)
library(here)
library(dplyr)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))


# Fixed Exponential Experiment (Experiment 1) ----
load(here('new_sim','expFixed.RData'))


# setEPS()
# postscript(here('new_sim','figure1.eps'),width=8,height=10)

tiff(file=here('new_sim','figure1.tiff'),width=6,height=7,units='in',res=800)

layout(cbind(1:4,c(5,7,9,11),c(6,8,10,12)),width=c(0.1,1,1),height=c(0.1,1,1,1))
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F)
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
text('n=100',x=0.5,y=0.5,las=2,adj=0.5,srt=90,cex=1.5)
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
text('n=250',x=0.5,y=0.5,las=2,adj=0.5,srt=90,cex=1.5)
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
text('n=500',x=0.5,y=0.5,las=2,adj=0.5,srt=90,cex=1.5)
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
text('Regression on Aoristic Sum',x=0.5,y=0.5,las=2,adj=0.5,cex=1.5)
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
text('Bayesian Exponential Fit',x=0.5,y=0.5,las=2,adj=0.5,cex=1.5)



par(mar=c(1,4,0,0))
for (k in 1:length(sample.sizes))
{
	plot(NULL,pch=20,ylim=c(exp.fixed[[k]]$r[1]-0.005,exp.fixed[[k]]$r[1]+0.005),main='',xlim=c(1,20),ylab='',xlab='',axes=F)
	arrows(x0=1:20,x1=1:20,y0=exp.fixed[[k]]$c95lo,y1=exp.fixed[[k]]$c95hi,length=0.04,code=3,angle=90,col=ifelse(exp.fixed[[k]]$c95acc,'darkblue','darkred'))
	abline(h=exp.fixed[[k]]$r[1],lty=2)
	axis(2)
	mtext(side=2,'r',las=2,line=3)


	plot(NULL,pch=20,ylim=c(exp.fixed[[k]]$r[1]-0.005,exp.fixed[[k]]$r[1]+0.005),main='',xlim=c(1,20),ylab='',xlab='',axes=F)
	arrows(x0=1:20,x1=1:20,y0=exp.fixed[[k]]$hpdilo,y1=exp.fixed[[k]]$hpdihi,length=0.04,code=3,angle=90,col=ifelse(exp.fixed[[k]]$hpdiacc,'darkblue','darkred'))
	abline(h=exp.fixed[[k]]$r[1],lty=2)
	axis(2)
	mtext(side=2,'r',las=2,line=3)
}

dev.off()

### Randomised Model (Experiment 2) ----
load(here('new_sim','res_coarse.RData'))
load(here('new_sim','res_fine.RData'))

# setEPS()
# postscript(here('new_sim','figure2.eps'),width=9,height=8)
tiff(file=here('new_sim','figure2.tiff'),width=8,height=7,res=800,units='in')

par(mfrow=c(2,2),mar=c(5,4,3,1))
plot(res.coarse$n,res.coarse$r,pch=20,col=ifelse(res.coarse$acc.lm,'darkgrey','indianred'),xlab='n',ylab='r',main=paste('Aoristic Sum (100yrs res); ',round(sum(res.coarse$acc.lm)/nrow(res.coarse)*100,1),'% acccuracy'))
plot(res.coarse$n,res.coarse$r,pch=20,col=ifelse(res.coarse$acc.bayes,'darkgrey','indianred'),xlab='n',ylab='r',main=paste('Bayesian (100yrs res); ',round(sum(res.coarse$acc.bayes)/nrow(res.coarse)*100,1),'% acccuracy'))
plot(res.fine$n,res.fine$r,pch=20,col=ifelse(res.fine$acc.lm,'darkgrey','indianred'),xlab='n',ylab='r',main=paste('Aoristic Sum (1yr res); ',round(sum(res.fine$acc.lm)/nrow(res.fine)*100,1),'% acccuracy'),cex=1)
plot(res.fine$n,res.fine$r,pch=20,col=ifelse(res.fine$acc.bayes,'darkgrey','indianred'),xlab='n',ylab='r',main=paste('Bayesian (1yr res); ',round(sum(res.fine$acc.bayes)/nrow(res.fine)*100,1),'% acccuracy'),cex=1)
dev.off()

### Shape Aoristic SUM ----
load(here('new_sim','res_icar_shape.RData'))

phaserect <- function(x,w=0.15)
{
	ybottom=quantile(par('usr')[3:4],1-w)
	ytop=par('usr')[4]
	for (i in 1:nrow(x))
	{
		rect(ybottom=ybottom,ytop=ytop,xleft=x[i,2],xright=x[i,3],border=ifelse(i%%2,'lightgrey','darkgrey'),col=ifelse(i%%2,'lightgrey','darkgrey'))
	}
}


# setEPS()
# postscript(here('new_sim','figure3.eps'),width=5,height=9)
tiff(file=here('new_sim','figure3.tiff'),width=5,height=7,res=800,units='in')
par(mfrow=c(3,2),mar=c(4,4,2,1))
for (i in 1:nrow(params))
{
	ymax=max(apply(probMats[[i]]$pmat,2,sum))*1.25
	plot(probMats[[i]],ylim=c(0,ymax),main=letters[i])
	phaserect(get(params$phases[i]))
	lines(5000:3001,P2*probMats[[i]]$resolution*probMats[[i]]$n,col='darkred',lty=2,lwd=2)
	box()
}
dev.off()


### Shape ICAR ----
# setEPS()
# postscript(here('new_sim','figure4.eps'),width=5,height=9)
tiff(file=here('new_sim','figure4.tiff'),width=5,height=7,res=800,units='in')
par(mfrow=c(3,2),mar=c(4,4,2,1))
for (i in 1:nrow(params))
{
	plot(icarfitted[[i]],plot.legend=ifelse(i==6,TRUE,FALSE),legend.arg=list(x='topright',bty='n'),hpd=c(0.5,0.95))
	title(letters[i])
	lines(5000:3001,P2*probMats[[i]]$resolution,col='darkred',lty=2,lwd=2)
	box()
}
dev.off()


### Equifinality plot ----
# setEPS()
# postscript(here('new_sim','figure5.eps'),width=8,height=8)
tiff(file=here('new_sim','figure5.tiff'),width=8,height=8,res=800,units='in')
load(here('new_sim','equifinality.RData'))
par(mfrow=c(2,2),mar=c(4,3,2,1))
# panel a
c1  <- dLogisticGrowth2(800:301,r=rr[1],a=800,b=301,m=mm[1],log=F)
c2  <- dLogisticGrowth2(800:301,r=rr[2],a=800,b=301,m=mm[2],log=F)
c3  <- dLogisticGrowth2(800:301,r=rr[3],a=800,b=301,m=mm[3],log=F)
plot(800:301,c1,type='l',xlim=c(800,301),ylim=c(0,0.009),lwd=2,axes=F,xlab='BP',ylab='Probability',main='a')
phaserect(x=data.frame(phasenames=letters[1:2],starts=c(800,500),ends=c(501,301)),w=0.05)
lines(800:301,c2,lty=2,lwd=2)
lines(800:301,c3,lty=3,lwd=2)
axis(1);axis(2)
legend(x=800,y=0.008,legend=c(paste0('r=',round(rr[1],5),' & m=',mm[1]),paste0('r=',round(rr[2],5),' & m=',mm[2]),paste0('r=',round(rr[3],5),' & m=',mm[3])),lty=1:3,lwd=2,bty='n',cex=0.9)

# panel b
plot(fit,hpd=c(0.95,0.95),col1='lightblue',col2='lightblue',plot.legend=FALSE,pch=NA,main='b',ylim=c(0,0.009),xlab='BP')

# panel c
plot(params$m,params$r,pch=20,col='lightgrey',xlab='m',ylab='r',main='c')
points(sparam$m,sparam$r,pch=20,col='indianred',xlab='m',ylab='r')

# panel d
plot(fit$posterior.m,fit$posterior.r,pch=1,col=rgb(0.67,0.84,0.9,0.1),xlim=c(300,800),ylim=c(0,0.03),main='d',xlab='m',ylab='r')
dev.off()

