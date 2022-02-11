# x = round(runif(100,min=4001,max=5000))
# library(here)
# source(here('src','unifdisc.R'))
# timeRange=c(5000,4001)
# limMin=10
# limMax=200

randtimespan  <- function(x,seed,limMin=10,limMax,timeRange,edge=TRUE)
{
	true.lim.max = rep(limMax,length(x))
	true.lim.min = rep(limMin,length(x))

	if (edge)
	{
	distToEdge = cbind(timeRange[1] - x, x - timeRange[2])
	true.lim.min = apply(cbind(distToEdge,true.lim.min),1,min)
	true.lim.max = apply(cbind(distToEdge,true.lim.max),1,min)
	} 

	limMatrix = cbind(true.lim.min,true.lim.max)
# 	set.seed(seed)
	actualLim = apply(limMatrix,1,function(x){runifdisc(1,x[1],x[2])})
	starts = x+actualLim
	ends = x-actualLim
	return(data.frame(phases=1:length(starts),starts=starts,ends=ends)) 
}
