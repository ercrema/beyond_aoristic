# Load library and custom R functions ----
library(here)
library(nimbleCarbon)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','aoristicsum.R'))
source(here('src','randtimespan.R'))
source(here('src','mcsim.R'))
source(here('src','ribbon.R'))
# source(here('src','unif2.R'))


# Simulate observed data in calendar time ----
set.seed(123)
nsample = 1000
timeRange = c(5000,3001)
changepoint = 4300
r1 = 0.002
r2 = -0.002
x = replicate(n=nsample, rDoubleExponentialGrowth(a=timeRange[1],b=timeRange[2],r1=r1,r2=r2,mu=changepoint))

# Convert observed data in to ordinal variable (i.e. periodisation) ----

## Scenario 1 - abutting phases with different lengths
phases1 <- data.frame(phasename=letters[1:7])
phases1$starts = c(5000,4700,4500,3800,3600,3500,3400)
phases1$ends   = c(4701,4501,3801,3601,3501,3401,3001)
x1  <- time2phase(x,phases1)

## Scenario 2 - abutting phases with same temporal duration
phases2 <- data.frame(phasename=letters[1:10])
phases2$starts = seq(5000,3200,-200)
phases2$ends   = seq(4801,3001,-200)
x2  <- time2phase(x,phases2)

## Scenario 3 - overlapping phses (emulating phase assignment uncertainty) 
phases3 <- data.frame(phasename=letters[1:14])
phases3$starts = c(5000,5000,4700,4700,4700,4700,4700,4500,4500,4500,4000,4000,3800,3500)
phases3$ends   = c(4001,4701,4501,4001,3801,3501,3001,4001,3801,3001,3801,3001,3001,3001)
x3  <- time2phase(x,phases3)

## Scenario 4 - random time-span error 
x4  <- randtimespan(x,seed=123,limMin=10,limMax=300,edge=FALSE)


# Aoristic Sum and MC simulation ----

# Set resolution of time-blocks:
resolution = 50 

# Aoristic Sum
asum1 = aoristicsum(x=x1,timeRange=timeRange,blockSize=resolution)
asum2 = aoristicsum(x=x2,timeRange=timeRange,blockSize=resolution)
asum3 = aoristicsum(x=x3,timeRange=timeRange,blockSize=resolution)
asum4 = aoristicsum(x=x4,timeRange=timeRange,blockSize=resolution)

# MC simulation
mcs1 = mcsim(x=x1,nsim=1000,timeRange=timeRange,blockSize=resolution)
mcs2 = mcsim(x=x2,nsim=1000,timeRange=timeRange,blockSize=resolution)
mcs3 = mcsim(x=x3,nsim=1000,timeRange=timeRange,blockSize=resolution)
mcs4 = mcsim(x=x4,nsim=1000,timeRange=timeRange,blockSize=resolution)

# Plot Aoristic Sum / MC simulation ----
# Generate expected density 
expected.dens  <- dDoubleExponentialGrowth(x=timeRange[1]:timeRange[2],a=timeRange[1],b=timeRange[2],r1=r1,r2=r2,mu=changepoint,log=FALSE)*resolution*nsample

# Plot Data
pdf(width=7,height=7,file=here('figures','figure_sim.pdf')) 
par(mfrow=c(2,2),lend=2)
plot(apply(asum1[,-3],1,median),asum1[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(c(expected.dens,asum1[,3],mcs1[[2]]))),main='a')
mdp = apply(mcs1[[1]],1,median)
ribbon(x=mdp,y=mcs1[[2]])
lines(mdp,apply(mcs1[[2]],1,median),lty=2,col='indianred',lwd=1.5)
lines(timeRange[1]:timeRange[2],expected.dens,col='steelblue',lty=2,lwd=2)


plot(apply(asum2[,-3],1,median),asum2[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(c(expected.dens,asum2[,3],mcs2[[2]]))),main='b')
mdp = apply(mcs2[[1]],1,median)
ribbon(x=mdp,y=mcs2[[2]])
lines(mdp,apply(mcs2[[2]],1,median),lty=2,col='indianred',lwd=1.5)
lines(timeRange[1]:timeRange[2],expected.dens,col='steelblue',lty=2,lwd=2)


plot(apply(asum3[,-3],1,median),asum3[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(c(expected.dens,asum3[,3],mcs3[[2]]))),main='a')
mdp = apply(mcs3[[1]],1,median)
ribbon(x=mdp,y=mcs3[[2]])
lines(mdp,apply(mcs3[[2]],1,median),lty=2,col='indianred',lwd=1.5)
lines(timeRange[1]:timeRange[2],expected.dens,col='steelblue',lty=2,lwd=2)


plot(apply(asum4[,-3],1,median),asum4[,3],xlim=timeRange,type='b',pch=20,xlab='calBP',ylab='Aoristic Sum',ylim=c(0,max(c(expected.dens,asum4[,3],mcs4[[2]]))),main='a')
mdp = apply(mcs4[[1]],1,median)
ribbon(x=mdp,y=mcs4[[2]])
lines(mdp,apply(mcs4[[2]],1,median),lty=2,col='indianred',lwd=1.5)
lines(timeRange[1]:timeRange[2],expected.dens,col='steelblue',lty=2,lwd=2)

legend('topright',legend=c('Aoristic Sum','MC Simulation','MC Simulation Median','True Density'),lty=c(1,1,2,2),col=c(1,'lightgrey','indianred','steelblue'),pch=c(20,NA,NA,NA),lwd=c(1,7,1,2),bty='n')
dev.off()

# Store simulated data ----
save(x1,x2,x3,x4,r1,r2,changepoint,timeRange,file=here('simulation','simulated_data.RData'))
