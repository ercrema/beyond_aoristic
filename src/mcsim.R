mcsim  <- function(x,nsim,timeRange,blockSize)
{
	tblocks = data.frame(starts =  seq(from=timeRange[1],to=timeRange[2],-blockSize), ends = seq(from=timeRange[1]-blockSize + 1, to=timeRange[2],-blockSize))
	tblock.mat = matrix(NA,nrow=nrow(tblocks),ncol=nsim)

	for (i in 1:nsim)
	{
		tmp = apply(x[,-1],1,function(x){runifdisc(1,x[2],x[1])})
		tmp2 = unlist(sapply(tmp,function(x,tblocks){which(x<=tblocks$starts & x>=tblocks$ends)},tblocks=tblocks))
		tblock.mat[,i] = as.numeric(table(c(1:nrow(tblocks),tmp2))-1)
	}
	return(list(tblocks,tblock.mat))
}


