library(nimbleCarbon)
library(here)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','aoristicsum.R'))
source(here('src','randtimespan.R'))
source(here('src','mcsim.R'))
source(here('src','ribbon.R'))
source(here('src','unif2.R'))


x = replicate(n=1000, rDoubleExponentialGrowth(a=5000,b=3001,r1=0.002,r2=-0.002,mu=4300))

# abutting, perfect info
phases1 <- data.frame(phasename=letters[1:7])
phases1$starts = c(5000,4700,4500,4000,3600,3500,3400)
phases1$ends   = c(4701,4501,4001,3601,3501,3401,3001)
x1  <- time2phase(x,phases1)
# phase uncertainty
phases2 <- data.frame(phasename=letters[1:10])
phases2$starts = c(5000,4700,4700,4700,4500,4000,4000,3500,3300,3500)
phases2$ends   = c(4201,4501,4101,4001,4001,3501,3001,3300,3001,3001)
x2  <- time2phase(x,phases2)
# true time-span
x3  <- randtimespan(x,seed=123,limMin=10,limMax=200,edge=FALSE)


par(mfrow=c(1,3))
timeRange = c(5000,3001)
resolution = 50
res1 = aoristicsum(x=x1,timeRange=timeRange,blockSize=resolution)
plot(apply(res1[,-3],1,median),res1[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum')
abline(v=4300,lty=2,col=2)

res2 = aoristicsum(x=x2,timeRange=timeRange,blockSize=resolution)
plot(apply(res2[,-3],1,median),res2[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum')
abline(v=4300,lty=2,col=2)


res3 = aoristicsum(x=x3,timeRange=timeRange,blockSize=resolution)
plot(apply(res3[,-3],1,median),res3[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum')
abline(v=4300,lty=2,col=2)




res1 = mcsim(x=x1,nsim=1000,timeRange=timeRange,blockSize=resolution)
res2 = mcsim(x=x2,nsim=1000,timeRange=timeRange,blockSize=resolution)
res3 = mcsim(x=x3,nsim=1000,timeRange=timeRange,blockSize=resolution)


mdp = apply(res1[[1]],1,median)
par(mfrow=c(1,3))
plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(res1[[2]])))
ribbon(x=mdp,y=res1[[2]])
lines(mdp,apply(res1[[2]],1,median),type='b',pch=20)
abline(v=4300,lty=2,col=2)

plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(res2[[2]])))
ribbon(x=mdp,y=res2[[2]])
lines(mdp,apply(res2[[2]],1,median),type='b',pch=20)
abline(v=4300,lty=2,col=2)

plot(NULL,xlim=rev(range(mdp)),ylim=c(0,max(res3[[2]])))
ribbon(x=mdp,y=res3[[2]])
lines(mdp,apply(res3[[2]],1,median),type='b',pch=20)
abline(v=4300,lty=2,col=2)


### Bayesian Analyses ###
model <- nimbleCode({
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







