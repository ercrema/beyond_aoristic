# Load Library, functions, and R images
library(here)
library(nimbleCarbon)
source(here('src','unif2.R'))
load(here('data','sitedata.RData'))

# Define total number of sites
N <- nrow(sitedb)

# Create Region ID and add to c14db
sitedb$RegionID = match(sitedb$Prefecture,win$name_vi)#check whether is actually win$Prefecture
Nregions  <- nrow(win)

# Define Data
d <- list(X=round(sitedb$m),s=round(sitedb$s))

# Define Constants
constants  <-  list(Nregions=Nregions,N=N,adj=nbInfo$adj,weights=nbInfo$weights,num=nbInfo$num,L=length(nbInfo$adj))
constants$a <- 2900
constants$b <- 1700
constants$id.region <-sitedb$RegionID

# Define theta init
theta.init <- d$X


# Core Analysis ----
runFun <- function(seed, d, constants, theta.init, nburnin, niter, thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','unif2.R'))
	## ICAR  Model
	icarmodel <- nimbleCode({
		for (i in 1:N)
		{
			theta[i] ~ dExponentialGrowth(a=a,b=b,r=r[id.region[i]])
			X[i] ~ dunif2(m=theta[i],s=s[i])
		}

		# ICAR Model Prior
		r[1:Nregions] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Nregions], tau, zero_mean =0)
		tau <- 1/sigma^2
		sigma ~ dunif(0,100)
	})

	## Setup Init
	set.seed(seed)
	inits <- list(sigma = runif(1,0,100), r = rnorm(constants$Nregions,sd=0.001),theta=theta.init)

	#MCMC
	model <- nimbleModel(icarmodel, constants = constants, data = d, inits = inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model, monitors = c('r','sigma'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC, project = cModel)
	samples <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(samples)
}

# Run MCMC ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds  <-  c(12,34,56,78)
niter  <- 250000
nburnin  <- 125000
thin  <- 5

chain_output  <- parLapply(cl = cl, X = seeds, fun = runFun, d = d,constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
icar.samples <- coda::mcmc.list(chain_output)
rhats.yayoi  <- coda::gelman.diag(icar.samples)
icar.yayoi  <- do.call(rbind.data.frame,icar.samples)
save(rhats.yayoi,icar.yayoi,file=here('results','icar_yayoi.RData'))
