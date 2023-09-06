library(baorista)
library(nimbleCarbon)
library(dplyr)
library(nimble)
library(coda)
library(here)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))


dAoristicExponentialGrowth_vector2=nimbleFunction(
  run = function(x = double(2),z=integer(0),r=double(0), log = integer(0)) {
    returnType(double(0))
    t = 1:z
    n = numeric(z)
    for (i in 1:z)
    {
      n[i] = (1+r)^t[i]
    }
    p = n/sum(n)
    pg = x %*% p
    k = which(pg>0)
    pg2 = pg[k,]
    logProb = sum(log(pg2))
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

suppressMessages(registerDistributions(list(
  dAoristicExponentialGrowth_vector2 = list(
    BUGSdist = "dAoristicExponentialGrowth_vector2(z,r)",
    pqAvail = FALSE,
    types = c('value = double(2)', 'z = integer(0)','r=double(0)')
  ))))


nsim  <- 1000
max.sample  <- 500
resolution <- 10
res <- data.frame(nsim.id=1:nsim,r=NA,n=NA,ci.lo=NA,ci.hi=NA,hpdi.lo=NA,hpdi.hi=NA,prec.lm=NA,prec.bayes=NA,acc.lm=NA,acc.bayes=NA)

for (i in 1:nsim)
{
	print(i)
	set.seed(i)
# 	res$r[i]  <- rnorm(1,mean=0,sd=0.002)
	res$r[i]  <- runif(1,-0.002,0.002)
	res$n[i] <- round(runif(1,min=100,max=max.sample))
	cal.dates  <-  replicate(n=res$n[i], rExponentialGrowth(a=4999,b=3002,r=res$r[i]))
	nphases <- sample(3:10,replace=T,size=1)
	phases <- diristick(alpha=rep(5,nphases),timeRange=c(5000,3001))
	adata <- time2phase(cal.dates,phases)
	x <- createProbMat(adata[,2:3],timeRange=c(5000,3001),resolution=resolution)
	yy <- apply(x$pmat,2,sum)
	dd  <- data.frame(xx=-apply(x$tblocks,1,median),yy=yy)
	fit  <- lm(log(yy+0.001)~xx,data=dd)
	ci <- confint(fit)[2,]
	if (i == 1)
	{

		d  <- list(theta=rbind(x$pmat,matrix(0,nrow=max.sample-res$n[i],ncol=ncol(x$pmat))))
		constants  <- list()
		constants$n.tblocks  <- ncol(d$theta) 
		expmodel  <- nimbleCode({
			theta[,] ~ dAoristicExponentialGrowth_vector2(r=r,z=n.tblocks)
			r ~ dunif(-1,1)
		})

		assign('dAoristicExponentialGrowth_vector2',dAoristicExponentialGrowth_vector2,envir=.GlobalEnv)
		inits  <- vector('list',length=4)
		for (k in 1:4)
		{
			inits[[k]]  <- list(r=rnorm(1,0,0.05))
		}
		model  <- nimbleModel(expmodel,constants=constants,data=d,inits=inits[[1]])
		assign('rAoristicExponentialGrowth_vector2',rAoristicExponentialGrowth_vector2,envir=.GlobalEnv)
		cModel <- compileNimble(model)
		conf  <- configureMCMC(model) 
		MCMC <- buildMCMC(conf)
		cMCMC <- compileNimble(MCMC,project=model)
		out <- runMCMC(cMCMC, niter = 100000, thin=10,nburnin =50000 ,inits=inits,samplesAsCodaMCMC = T,nchains=4,progressBar=F)
		post <- do.call(rbind.data.frame,out)/x$resolution
	}

	if (i>1)
	{

		cModel$resetData()
		cModel$setData(theta=rbind(x$pmat,matrix(0,nrow=max.sample-res$n[i],ncol=ncol(x$pmat))))
# 		MCMC <- buildMCMC(conf)
# 		cMCMC <- compileNimble(MCMC)
		out <- runMCMC(cMCMC, niter = 100000, thin=10,nburnin =50000 ,inits=inits,samplesAsCodaMCMC = T,nchains=4,progressBar=F)
		post <- do.call(rbind.data.frame,out)/x$resolution
	}

# 	fitB <- expfit(x,rPrior='dunif(-1,1)')
# 	hpdi <- HPDinterval(mcmc(fitB$posterior.r)) |> as.numeric()
	hpdi <- HPDinterval(mcmc(post)) |> as.numeric()
	res$ci.lo[i] <- ci[1]
	res$ci.hi[i] <- ci[2]
	res$prec.lm[i] <- abs(diff(ci))
	res$hpdi.lo[i] <- hpdi[1]
	res$hpdi.hi[i] <- hpdi[2]
	res$prec.bayes[i] <- abs(diff(hpdi))
	res$acc.lm[i] <- between(res$r[i],ci[1],ci[2])
	res$acc.bayes[i] <- between(res$r[i],hpdi[2],hpdi[2])
}
res.fine <- res
save(res.fine,file='res_fine.RData')

