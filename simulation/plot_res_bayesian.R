# Load Library and Data ----
library(here)
library(nimbleCarbon)
source(here('src','aoristicsum.R'))
source(here('src','unifdisc.R'))
load(here('simulation','sim_analyses_1.RData'))
load(here('simulation','sim_analyses_2.RData'))
load(here('simulation','sim_analyses_3.RData'))
load(here('simulation','sim_analyses_4.RData'))
load(here('simulation','simulated_data.RData'))

post.sample.1  <- do.call(rbind.data.frame,post.sample.1)
post.sample.2  <- do.call(rbind.data.frame,post.sample.2)
post.sample.3  <- do.call(rbind.data.frame,post.sample.3)
post.sample.4  <- do.call(rbind.data.frame,post.sample.4)


pdf(width=7,height=7,file=here('figures','figure_sim_bayes.pdf')) 
expected.dens  <- dDoubleExponentialGrowth(x=timeRange[1]:timeRange[2],a=timeRange[1],b=timeRange[2],r1=r1,r2=r2,mu=changepoint,log=FALSE)
par(mfrow=c(2,2),lend=2)

modelPlot(model=dDoubleExponentialGrowth,a=timeRange[1],b=timeRange[2],params=list(r1=post.sample.1$r1,r2=post.sample.1$r2,mu=post.sample.1$mu),nsample = 1000,type='envelope',main='a',interval = 0.95)
lines(timeRange[1]:timeRange[2],expected.dens,lty=2,lwd=2,col='steelblue')
asum1 = aoristicsum(x=x1,timeRange=timeRange,blockSize=50)
lines(apply(asum1[,-3],1,median),asum1[,3]/nrow(x1)/50,pch=20,col=1,type='b')

modelPlot(model=dDoubleExponentialGrowth,a=timeRange[1],b=timeRange[2],params=list(r1=post.sample.2$r1,r2=post.sample.2$r2,mu=post.sample.2$mu),nsample = 1000,type='envelope',main='b',interval = 0.95)
lines(timeRange[1]:timeRange[2],expected.dens,lty=2,lwd=2,col='steelblue')
asum2 = aoristicsum(x=x2,timeRange=timeRange,blockSize=50)
lines(apply(asum2[,-3],1,median),asum2[,3]/nrow(x1)/50,pch=20,col=1,type='b')

modelPlot(model=dDoubleExponentialGrowth,a=timeRange[1],b=timeRange[2],params=list(r1=post.sample.3$r1,r2=post.sample.3$r2,mu=post.sample.3$mu),nsample = 1000,type='envelope',main='c',interval = 0.95)
lines(timeRange[1]:timeRange[2],expected.dens,lty=2,lwd=2,col='steelblue')
asum3 = aoristicsum(x=x3,timeRange=timeRange,blockSize=50)
lines(apply(asum3[,-3],1,median),asum3[,3]/nrow(x1)/50,pch=20,col=1,type='b')

modelPlot(model=dDoubleExponentialGrowth,a=timeRange[1],b=timeRange[2],params=list(r1=post.sample.4$r1,r2=post.sample.4$r2,mu=post.sample.4$mu),nsample = 1000,type='envelope',main='d',interval=0.95)
lines(timeRange[1]:timeRange[2],expected.dens,lty=2,lwd=2,col='steelblue')
asum4 = aoristicsum(x=x4,timeRange=timeRange,blockSize=50)
lines(apply(asum4[,-3],1,median),asum4[,3]/nrow(x1)/50,pch=20,col=1,type='b')

legend('topright',legend=c('Aoristic Sum','True Density','Posterior Mean','Posterior 95% HPDI'),lty=c(1,1,2,2),col=c(1,'lightblue',1,'lightgrey'),pch=c(20,NA,NA,NA),lwd=c(1,1,1,7),bty='n',cex=0.7)

dev.off()

