library(nimbleCarbon)
library(baorista)
library(here)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))

timeRange  <- c(5000,3001)
BP  <- c(5000,4750,4500,4250,4000,3750,3500,3350,3300)
logP  <- c(5,2,1,0,1,0.8,0.6,0.2,1.5) 
fit  <- lm(logP~poly(BP,3))
Pseq <- predict(fit,data.frame(BP=timeRange[1]:timeRange[2]))
P <- (Pseq-min(Pseq)) / sum(Pseq-min(Pseq))
P2 <- P
plot(timeRange[1]:timeRange[2],P,type='l',xlim=timeRange)
set.seed(123245323)
true.dates.50  <- sample(timeRange[1]:timeRange[2],size=50,replace=T,prob=P)
true.dates.500  <- sample(timeRange[1]:timeRange[2],size=500,replace=T,prob=P)
phases.1 <- diristick(alpha=rep(2,3),timeRange=c(5000,3001))
phases.2 <- diristick(alpha=rep(2,5),timeRange=c(5000,3001))
phases.3 <- diristick(alpha=rep(2,8),timeRange=c(5000,3001))
params <- expand.grid(samples=c('true.dates.50','true.dates.500'),phases=c('phases.1','phases.2','phases.3'),stringsAsFactors = F)
icarfitted <- probMats <- vector('list',length=nrow(params))

for (i in 1:nrow(params))
{
	x <- get(params$samples[i])
	y <- get(params$phases[i])
	probMats[[i]] <- time2phase(x,y)[,2:3] |> createProbMat(timeRange=c(5000,3001),resolution=50)
}

par(mfrow=c(3,2))
for (i in 1:length(probMats))
{
	plot(probMats[[i]])
	lines(timeRange[1]:timeRange[2],P*probMats[[i]]$resolution*probMats[[i]]$n,col='red')
}

# hist(true.dates.50,xlim=c(5000,3001),breaks=30)
# hist(true.dates.500,xlim=c(5000,3001),breaks=30)

for (i in 1:length(probMats))
{

	icarfitted[[i]] <- icarfit(probMats[[i]],niter=4000000,nburnin=3000000,thin=100,parallel=T)
}


save.image(file='eaa_shape.RData')
load('eaa_shape.RData')
par(mfrow=c(3,1),mar=c(4,4,2,1))
# plot(scenario.2.50.fit,plot.legend=F,hpd=c(0.5,0.95))
# lines(5000:3001,P2*50,col=2,lty=2)
# title('n=50')
plot(scenario.2.100.pm,type='dens')
plot(scenario.2.100.fit,plot.legend=F,hpd=c(0.5,0.95))
lines(5000:3001,P2*50,col=2,lty=2,lwd=2)
title('n=100')
plot(scenario.2.500.fit,plot.legend=F,hpd=c(0.5,0.95))
lines(5000:3001,P2*50,col=2,lty=2,lwd=2)
title('n=500')
plot(scneario.2.2000.fit,plot.legend=T,legend.arg=list(x='topright',bty='n'),hpd=c(0.5,0.95))
lines(5000:3001,P2*50,col=2,lty=2,lwd=2)
title('n=2000')

# plot(scneario.1.200.fit,plot.legend=F,ylim=c(0,0.2))
# lines(5000:3001,P1*100,col=2)

plot(scneario.1.5000.fit,plot.legend=F,ylim=c(0,0.2))
lines(5000:3001,P1*100,col=2)

# plot(scneario.2.500.fit,plot.legend=F,ylim=c(0,0.2))
# lines(5000:3001,P2*100,col=2)

plot(scneario.2.5000.fit,plot.legend=F,ylim=c(0,0.2))
lines(5000:3001,P2*100,col=2)
