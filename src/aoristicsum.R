aoristicsum  <- function(x,timeRange,blockSize)
{

	aoristic = data.frame(starts =  seq(from=timeRange[1],to=timeRange[2],-blockSize), ends = seq(from=timeRange[1]-blockSize + 1, to=timeRange[2],-blockSize))
	aoristic$sum = NA

	for (i in 1:nrow(aoristic))
	{
		aoristic$sum[i] = sum(apply(x[,-1],1,function(x,a,b){sum(dunifdisc(a:b,max=x[1],min=x[2]))},a=aoristic[i,2],b=aoristic[i,1]))
	}

	return(aoristic)
}
