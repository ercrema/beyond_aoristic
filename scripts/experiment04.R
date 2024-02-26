library(nimbleCarbon)
library(baorista)
library(here)
source(here('src','time2phase.R'))
source(here('src','unifdisc.R'))
source(here('src','diristick.R'))


# Generate some arbitrary underlying population
timeRange  <- c(5000,3001)
BP  <- c(5000,4750,4500,4250,4000,3750,3500,3350,3300)
logP  <- c(5,2,1,0,1,0.8,0.6,0.2,1.5) 
fit  <- lm(logP~poly(BP,3))
Pseq <- predict(fit,data.frame(BP=timeRange[1]:timeRange[2]))
P <- (Pseq-min(Pseq)) / sum(Pseq-min(Pseq))
P2 <- P
# plot(timeRange[1]:timeRange[2],P,type='l',xlim=timeRange)


# Sample Dates from the distribution
set.seed(123245323)
true.dates.50  <- sample(timeRange[1]:timeRange[2],size=50,replace=T,prob=P)
true.dates.500  <- sample(timeRange[1]:timeRange[2],size=500,replace=T,prob=P)
# Assign to random phases
phases.1 <- diristick(alpha=rep(2,3),timeRange=c(5000,3001))
phases.2 <- diristick(alpha=rep(2,5),timeRange=c(5000,3001))
phases.3 <- diristick(alpha=rep(2,8),timeRange=c(5000,3001))
params <- expand.grid(samples=c('true.dates.50','true.dates.500'),phases=c('phases.1','phases.2','phases.3'),stringsAsFactors = F)
icarfitted <- probMats <- vector('list',length=nrow(params))


# Setup data
for (i in 1:nrow(params))
{
	x <- get(params$samples[i])
	y <- get(params$phases[i])
	probMats[[i]] <- time2phase(x,y)[,2:3] |> createProbMat(timeRange=c(5000,3001),resolution=50)
}

# ICAR fit
for (i in 1:length(probMats))
{

	icarfitted[[i]] <- icarfit(probMats[[i]],niter=4000000,nburnin=3000000,thin=100,parallel=T)
}

# Store results
save.image(file=here('results','exp04_res.RData'))
