#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load library and Data ----
library(nimbleCarbon)
library(parallel)
library(here)
load(here('simulation','simulated_data.RData'))

# Setup Data ----
d  <- list()
d$X  <- x4$m
d$s  <- x4$s
constants = list(N=nrow(x4),a=5000,b=3001)

# Setup MCMC ----
ncores  <- as.numeric(args[[1]])
cl  <- makeCluster(ncores)
seeds  <- c(12,34,56,78)[1:ncores]
niter  <- as.numeric(args[[2]])
nburnin  <- round(ceiling(niter/2))
thin  <- ceiling((niter-nburnin)/10000)

# Create Run Script ----
### Bayesian Analyses ###
run  <- function(seed,d,constants,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','unif2.R'))
	growthmodel <- nimbleCode({
		for (i in 1:N){
			# Growth Model Likelihood
			theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=r1,r2=r2,mu=mu);
			X[i] ~ dunif2(m=theta[i],s=s[i]);
		}
		# Prior
		r1 ~ dnorm(mean=0,sd=0.1); # Prior
		r2 ~ dnorm(mean=0,sd=0.1); # Prior
		mu ~ dunif(b,a); # Prior
	})  
	set.seed(seed)
	inits  <- list()
	inits$theta  <- d$X
	inits$r1  <- rnorm(1,0,0.1)
	inits$r2  <- rnorm(1,0,0.1)
	inits$mu  <- round(runif(1,constants$b,constants$a))


	model <- nimbleModel(growthmodel,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('r1')
	conf$addMonitors('r2')
	conf$addMonitors('mu')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(results)
}

# Run MCMC ----
print(paste0('Running with nchains=', ncores, '; niter=', niter,'; nburnin=', nburnin, ',; thin=',thin))
out.4  <- parLapply(cl=cl,X=seeds,fun=run,d=d,constants=constants,niter=niter,nburnin=nburnin,thin=thin)
stopCluster(cl)

# Diagnostic and Posterior Processing
post.sample.4  <- coda::mcmc.list(out.4)
diag.4  <- coda::gelman.diag(post.sample.4)
print('Diagnostics:\n')
print(diag.4$psrf)
save(post.sample.4,diag.4,file=here('simulation','sim_analyses_4.RData'))





