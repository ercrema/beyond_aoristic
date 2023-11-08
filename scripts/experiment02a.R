library(baorista)
library(nimbleCarbon)
library(dplyr)
library(nimble)
library(coda)
library(here)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))


# Register nimble function with for inference
# This allows the use of the same compiled model saving processing time 
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


nsim  <- 1000 #Number of replicates
max.sample  <- 500 #Maximum sample size
resolution <- 100 #Time-block resolution
res <- data.frame(nsim.id=1:nsim,r=NA,n=NA,ci.lo=NA,ci.hi=NA,hpdi.lo=NA,hpdi.hi=NA,prec.lm=NA,prec.bayes=NA,acc.lm=NA,acc.bayes=NA)

for (i in 1:nsim)
{
	print(i)
	set.seed(i)
	#Randomly sample growth rate and sample size
	res$r[i]  <- runif(1,-0.002,0.002)
	res$n[i] <- round(runif(1,min=100,max=max.sample))
	#Randomly sample dates via nimbleCarbon
	cal.dates  <-  replicate(n=res$n[i], rExponentialGrowth(a=4999,b=3002,r=res$r[i]))
	#Dates to time-span conversion (and random generation of archaeological periodisations)
	nphases <- sample(3:10,replace=T,size=1)
	phases <- diristick(alpha=rep(5,nphases),timeRange=c(5000,3001))
	adata <- time2phase(cal.dates,phases)
	#Data setup
	x <- createProbMat(adata[,2:3],timeRange=c(5000,3001),resolution=resolution)
	yy <- apply(x$pmat,2,sum)
	dd  <- data.frame(xx=-apply(x$tblocks,1,median),yy=yy)
	#Regression based estimate
	fit  <- lm(log(yy+0.001)~xx,data=dd)
	ci <- confint(fit)[2,]

	#The lines below replaces baorista::expfit() using the same compile C object
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
		out <- runMCMC(cMCMC, niter = 100000, thin=10,nburnin =50000 ,inits=inits,samplesAsCodaMCMC = T,nchains=4,progressBar=F)
		post <- do.call(rbind.data.frame,out)/x$resolution
	}

	hpdi <- HPDinterval(mcmc(post)) |> as.numeric()

	# Store output
	res$ci.lo[i] <- ci[1]
	res$ci.hi[i] <- ci[2]
	res$prec.lm[i] <- abs(diff(ci))
	res$hpdi.lo[i] <- hpdi[1]
	res$hpdi.hi[i] <- hpdi[2]
	res$prec.bayes[i] <- abs(diff(hpdi))
	res$acc.lm[i] <- between(res$r[i],ci[1],ci[2])
	res$acc.bayes[i] <- between(res$r[i],hpdi[1],hpdi[2])
}

res.coarse <- res
save(res.coarse,file=here('results','exp02a_res.RData'))

