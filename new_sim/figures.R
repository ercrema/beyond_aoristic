library(baorista)
library(nimbleCarbon)
library(here)
library(dplyr)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))
load(here('new_sim','expFixed.RData'))


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
