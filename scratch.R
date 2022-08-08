library(nimbleCarbon)
library(here)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','aoristicsum.R'))
source(here('src','randtimespan.R'))
source(here('src','mcsim.R'))
source(here('src','ribbon.R'))
source(here('src','unif2.R'))

set.seed(123)
nsample = 1000
x = replicate(n=nsample, rDoubleExponentialGrowth(a=5000,b=3001,r1=0.002,r2=-0.002,mu=4300))

# abutting, perfect info
phases1 <- data.frame(phasename=letters[1:7])
phases1$starts = c(5000,4700,4500,3800,3600,3500,3400)
phases1$ends   = c(4701,4501,3801,3601,3501,3401,3001)
x1  <- time2phase(x,phases1)

# abutting, regular chronology of 200 years
phases2 <- data.frame(phasename=letters[1:10])
phases2$starts = seq(5000,3200,-200)
phases2$ends   = seq(4801,3001,-200)
x2  <- time2phase(x,phases2)

# phase uncertainty
phases3 <- data.frame(phasename=letters[1:14])
phases3$starts = c(5000,5000,4700,4700,4700,4700,4700,4500,4500,4500,4000,4000,3800,3500)
phases3$ends   = c(4001,4701,4501,4001,3801,3501,3001,4001,3801,3001,3801,3001,3001,3001)
x3  <- time2phase(x,phases3)

# true time-span
x4  <- randtimespan(x,seed=123,limMin=10,limMax=300,edge=FALSE)



# Aoristic Sum ----

pdf(width=7,height=7,file=here('figures','figure1.pdf')) 
par(mfrow=c(2,2))
resolution = 50
tt  <- 5000:3001
dd  <- dDoubleExponentialGrowth(x=tt,a=5000,b=3001,r1=0.002,r2=-0.002,mu=4300,log=FALSE)*resolution*nsample
timeRange = c(5000,3001)
res1 = aoristicsum(x=x1,timeRange=timeRange,blockSize=resolution)
plot(apply(res1[,-3],1,median),res1[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(dd)),main='a')
# abline(v=4300,lty=2,col=2)
lines(tt,dd,col='steelblue',lty=2,lwd=2)

res2 = aoristicsum(x=x2,timeRange=timeRange,blockSize=resolution)
plot(apply(res2[,-3],1,median),res2[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(dd)),main='b')
# abline(v=4300,lty=2,col=2)
lines(tt,dd,col='steelblue',lty=2,lwd=2)

res3 = aoristicsum(x=x3,timeRange=timeRange,blockSize=resolution)
plot(apply(res3[,-3],1,median),res3[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(dd)),main='c')
# abline(v=4300,lty=2,col=2)
lines(tt,dd,col='steelblue',lty=2,lwd=2)

res4 = aoristicsum(x=x4,timeRange=timeRange,blockSize=resolution)
plot(apply(res4[,-3],1,median),res4[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(dd)),main='d')
# abline(v=4300,lty=2,col=2)
lines(tt,dd,col='steelblue',lty=2,lwd=2)
legend('topright',legend=c('Aoristic Sum','True Density'),lty=c(1,2),col=c(1,'steelblue'),pch=c(20,NA),lwd=c(1,2))
dev.off()



# Monte Carlo Simulation ----

res1 = mcsim(x=x1,nsim=1000,timeRange=timeRange,blockSize=resolution)
res2 = mcsim(x=x2,nsim=1000,timeRange=timeRange,blockSize=resolution)
res3 = mcsim(x=x3,nsim=1000,timeRange=timeRange,blockSize=resolution)
res4 = mcsim(x=x4,nsim=1000,timeRange=timeRange,blockSize=resolution)


pdf(width=7,height=7,file=here('figures','figure2.pdf')) 
mdp = apply(res1[[1]],1,median)
par(mfrow=c(2,2),lend=2)
plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(c(res1[[2]],dd))),xlab='cal BP',main='a')
ribbon(x=mdp,y=res1[[2]])
lines(mdp,apply(res1[[2]],1,median),type='b',pch=20)
lines(tt,dd,col='steelblue',lty=2,lwd=2)


plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(c(res1[[2]],dd))),xlab='cal BP',main='b')
ribbon(x=mdp,y=res2[[2]])
lines(mdp,apply(res2[[2]],1,median),type='b',pch=20)
lines(tt,dd,col='steelblue',lty=2,lwd=2)


plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(c(res1[[2]],dd))),xlab='cal BP',main='c')
ribbon(x=mdp,y=res3[[2]])
lines(mdp,apply(res3[[2]],1,median),type='b',pch=20)
lines(tt,dd,col='steelblue',lty=2,lwd=2)


plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(c(res1[[2]],dd))),xlab='cal BP',main='d')
ribbon(x=mdp,y=res4[[2]])
lines(mdp,apply(res4[[2]],1,median),type='b',pch=20)
lines(tt,dd,col='steelblue',lty=2,lwd=2)
legend('topright',legend=c('Aoristic Sum','MC Simulation','True Density'),lty=c(1,1,2),col=c(1,'lightgrey','steelblue'),pch=c(20,NA,NA),lwd=c(1,7,2),bty='n')
dev.off()


### Bayesian Analyses ###
odel <- nimbleCode({
      for (i in 1:N){
        # Growth Model Likelihood
        theta[i] ~ dDoubleExponentialGrowth(a=start,b=end,r1=r1,r2=r2,mu=mu);
        X[i] ~ dunif2(m=theta[i],s=s[i]);
      }
      # Prior
      r1 ~ dnorm(0,0.1); # Prior
      r2 ~ dnorm(0,0.1); # Prior
      mu ~ dunif(end,start); # Prior
    })  

d = list(X=round(apply(x1[,-1],1,median)),s=apply(x1[,-1],1,function(x){abs(diff(x))/2}))
inits = list(r1=0.01,r2=0.01,mu=4500,theta=d$X)
constants = list(N=nrow(x1),start=5000,end=3001)
modelCode  <- nimbleModel(code=model,constants=constants,inits=inits)
mcmc.samples<- nimbleMCMC(code = model,constants = constants,data = d,niter = 10000, nchains = 2, thin=1, nburnin = 2000, monitors=c('r1','r2','mu'), inits=inits, samplesAsCodaMCMC=TRUE,setSeed=c(123,456))
mcmc.samples <- do.call(rbind,mcmc.samples)
save(mcmc.samples,file='temp.RData')







