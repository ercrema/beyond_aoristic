library(here)
library(nimbleCarbon)
library(baorista)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))


# Parameter Space Plot ----
params = expand.grid(r=seq(0.001,0.03,length.out=1000),m=799:302)
params$s <- NA
for (i in 1:nrow(params))
{
	params$s[i] <- round(sum(dLogisticGrowth2(800:501,r=params$r[i],a=800,b=301,m=params$m[i],log=F)),2)
}

# Subset to all runs where the time-interval of the first phase (800~501) has a cumulative probability mass of 0.13:

sparam <- subset(params,s==0.13)
# plot(params$r,params$m,col='lightgrey',pch=20)
# points(sparam$r,sparam$m,col='indianred',pch=20)
rr <- sparam$r[c(1,5483,3000)]
mm <- sparam$m[c(1,5483,3000)]
# c1  <- dLogisticGrowth2(800:301,r=rr[1],a=800,b=301,m=mm[1],log=F)
# c2  <- dLogisticGrowth2(800:301,r=rr[5483],a=800,b=301,m=mm[5483],log=F)
# c3  <- dLogisticGrowth2(800:301,r=rr[3000],a=800,b=301,m=mm[3000],log=F)
# plot(800:301,c1,type='l',xlim=c(800,301),ylim=c(0,0.01))
# lines(800:301,c2,col=2)
# lines(800:301,c3,col=3)


# Fit Logistic model
# rr[1];mm[1] #0.02988388 and 509
caldates  <- replicate(500,rLogisticGrowth2(1,r=rr[1],a=800,b=301,m=mm[1]))
probMatEqui <- time2phase(caldates,data.frame(names=letters[1:2],starts=c(800,500),ends=c(501,301)))[,2:3] |> createProbMat(timeRange=c(800,301),resolution=1)
fit  <- logisticfit(probMatEqui,rPrior='dunif(0.001,0.03)',mPrior='dunif(1,z)')
save(rr,mm,sparam,fit,here('results','exp03_res.RData'))


