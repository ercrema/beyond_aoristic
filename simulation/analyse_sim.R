library(nimbleCarbon)
library(here)
source(here('src','unif2.R'))
load(here('simulation','simulated_data.RData'))


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







