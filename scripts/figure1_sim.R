library(nimbleCarbon)
library(here)
source(here('src','time2phase.R'))
source(here('src','randtimespan.R'))
source(here('src','aoristicsum.R'))
source(here('src','unifdisc.R'))
source(here('src','mcsim.R'))

set.seed(123)
nsample = 1000
timeRange = c(5000,3001)
changepoint = 4300
r1 = 0.002
r2 = -0.002
x = replicate(n=nsample, rDoubleExponentialGrowth(a=timeRange[1],b=timeRange[2],r1=r1,r2=r2,mu=changepoint))

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
## Scenario 1 - abutting phases with different lengths
phases1 <- data.frame(phasename=letters[1:7])
phases1$starts = c(5000,4700,4500,3800,3600,3500,3400)
phases1$ends   = c(4701,4501,3801,3601,3501,3401,3001)
x1  <- time2phase(x,phases1)

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


expected.dens  <- dDoubleExponentialGrowth(x=timeRange[1]:timeRange[2],a=timeRange[1],b=timeRange[2],r1=r1,r2=r2,mu=changepoint,log=FALSE)*resolution*nsample

save(phases1,phases2,phases3,mcs1,mcs2,mcs3,mcs4,asum1,asum2,asum3,asum4,expected.dens,timeRange,file=here('results','figure1_res.RData'))

