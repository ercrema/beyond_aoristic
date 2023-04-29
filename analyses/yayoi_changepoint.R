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
constants  <-  list()
constants$Nregions  <- Nregions
constants$N  <- N
constants$adj <- nbInfo$adj
constants$weights  <- nbInfo$weights
constants$num  <- nbInfo$num
constants$L  <- length(nbInfo$adj)
constants$a <- 2900
constants$b <- 1700
constants$id.region <-sitedb$RegionID
contants$const_sequence  <- rep(1,constants$Nregions)

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
			theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=r1[id.region[i]],r2=r2[id.region[i]],mu=chp[id.region[i]])
			X[i] ~ dunif2(m=theta[i],s=s[i])
		}

		# ICAR Model Prior
		r1[1:Nregions] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Nregions], tau1, zero_mean =0)
		r2[1:Nregions] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Nregions], tau2, zero_mean =0)
		w[1:Nregions] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Nregions], tau3, zero_mean =0)

		tau1 <- 1/sigma1^2
		tau2 <- 1/sigma2^2
		tau3 <- 1/sigma3^2
		sigma1 ~ dunif(0,100)
		sigma2 ~ dunif(0,100)
		sigma3 ~ dunif(0,100)
		mu ~ dunif(b+1,a-1)
		for (k in 1:Nregions)
		{
			chp[k] <- round(w[k] + mu)
			const_sequence[k] ~ dconstraint(chp[k]>a & chp[k]<b)
		}

	})

	## Setup Init
	set.seed(seed)
	inits <- list()
	inits$sigma1  <- runif(1,0,100)
	inits$sigma2  <- runif(1,0,100)
	inits$sigma3  <- runif(1,0,100)
	inits$r1  <- rnorm(constants$Nregions,mean=0,sd=0.001)
	inits$r2  <- rnorm(constants$Nregions,mean=0,sd=0.001)
	inits$theta  <- theta.init
	cond  <- TRUE
	while(cond)
	{
		inits$w  <- rnorm(constants$Nregions,mean=0,sd=10)
		inits$mu  <- runif(1,constants$b+1,constants$a-1)
		if (all((inits$w+inits$mu)>constants$b & (inits$w+inits$mu)<constants$a)){cond <- FALSE}
	}

	#MCMC
	model <- nimbleModel(icarmodel, constants = constants, data = d, inits = inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model, monitors = c('r1','r2','chp'))
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
rhats.yayoi2  <- coda::gelman.diag(icar.samples)
icar.yayoi2  <- do.call(rbind.data.frame,icar.samples)
save(rhats.yayoi2,icar.yayoi2,file=here('results','icar_yayoi2.RData'))
