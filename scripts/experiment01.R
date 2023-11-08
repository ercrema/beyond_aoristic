library(baorista)
library(nimbleCarbon)
library(here)
library(dplyr)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))


nsim=20 #number of replicates
sample.sizes  <- c(100,250,500) #sample size
r <- 0.002 #'true' growth rate
resolution <- 10 #timeblock resolution
exp.fixed  <- vector('list',length=length(sample.sizes))

for (k in 1:length(sample.sizes))
{
	exp.fixed[[k]] <- data.frame(id=1:nsim,r=NA,c95lo=NA,c95hi=NA,c95prec=NA,c95acc=NA,hpdilo=NA,hpdihi=NA,hpdiprec=NA,hpdiacc=NA)
}

for (i in 1:nsim)
{
	for (k in 1:length(sample.sizes))
	{
		set.seed(i)
		
		#generate samples from exponential growth
		cal.dates  <-  replicate(n=sample.sizes[k], rExponentialGrowth(a=4999,b=3002,r=r))
		
		# Generate random phases and convert dates to timespans
		nphases <- sample(3:10,replace=T,size=1)
		phases <- diristick(alpha=rep(0.5,nphases),timeRange=c(5000,3001)) 
		adata <- time2phase(cal.dates,phases)
		
		# Data setup
		x <- createProbMat(adata[,2:3],timeRange=c(5000,3001),resolution=resolution)
		yy <- apply(x$pmat,2,sum)
		dd  <- data.frame(xx=-apply(x$tblocks,1,median),yy=yy)

		# Estimate growth rate via regression
		fit  <- lm(log(yy+0.001)~xx,data=dd)
		c95 <- confint(fit)[2,]

		# Estimate growth rate via baorista
		fitB <- expfit(x,rPrior='dunif(-1,1)',parallel=TRUE)
		hpdi <- HPDinterval(mcmc(fitB$posterior.r)) |> as.numeric()

		# Store output
		exp.fixed[[k]]$r[i] <- r
		exp.fixed[[k]]$c95lo[i] <- c95[1]
		exp.fixed[[k]]$c95hi[i] <- c95[2]
		exp.fixed[[k]]$c95prec[i] <- abs(diff(c95))
		exp.fixed[[k]]$c95acc[i] <- (c95[1]<=r & c95[2]>=r)
		exp.fixed[[k]]$hpdilo[i] <- hpdi[1]
		exp.fixed[[k]]$hpdihi[i] <- hpdi[2]
		exp.fixed[[k]]$hpdiprec[i] <- abs(diff(hpdi))
		exp.fixed[[k]]$hpdiacc[i] <- (hpdi[1]<=r & hpdi[2]>=r)
	}
}

save(exp.fixed,sample.sizes,r,file=here('results','exp01_res.RData'))
