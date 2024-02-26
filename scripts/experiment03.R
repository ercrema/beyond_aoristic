library(here)
library(nimbleCarbon)
library(baorista)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))


# Define parameter ranges
params = expand.grid(r=seq(0.001,0.03,length.out=100),m=seq(302,799,by=4))
params$s <- NA

# Compute cumulative probability between 800 and 501 for each parameter combination
for (i in 1:nrow(params))
{
	params$s[i] <- round(sum(dLogisticGrowth2(800:501,r=params$r[i],a=800,b=301,m=params$m[i],log=F)),2)
}

# Subset to all runs where the time-interval of the first phase (800~501) has a cumulative probability mass of 0.13:
sparam <- subset(params,s==0.13)

# Consider three example cases
rr <- sparam$r[c(1,70,138)]
mm <- sparam$m[c(1,70,138)]
c1  <- dLogisticGrowth2(800:301,r=rr[1],a=800,b=301,m=mm[1],log=F)
c2  <- dLogisticGrowth2(800:301,r=rr[2],a=800,b=301,m=mm[2],log=F)
c3  <- dLogisticGrowth2(800:301,r=rr[3],a=800,b=301,m=mm[3],log=F)


# Fit Logistic model via baorista
caldates  <- replicate(500,rLogisticGrowth2(1,r=rr[1],a=800,b=301,m=mm[1]))
probMatEqui <- time2phase(caldates,data.frame(names=letters[1:2],starts=c(800,500),ends=c(501,301)))[,2:3] |> createProbMat(timeRange=c(800,301),resolution=1)
fit  <- logisticfit(probMatEqui,rPrior='dunif(0.001,0.03)',mPrior='dunif(1,z)')

#Store results
save(rr,mm,params,sparam,fit,file=here('results','exp03_res.RData'))


