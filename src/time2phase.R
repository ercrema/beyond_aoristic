# x ... calendar time of each event (in BP)
# phases ... data.frame with three columns, phasename, start, and end date (in BP)
# examples:
# x  <- runif(50,min=4001,max=5000)
# phases1  <- data.frame(phasename=letters[1:10])
# phases1$starts = c(5000,4900,4800,4700,4600,4500,4400,4300,4200,4100)
# phases1$ends   = c(4901,4801,4701,4601,4501,4401,4301,4201,4101,4001)
# 
# phases2  <- data.frame(phasename=letters[1:10])
# phases2$starts = c(5000,4800,4600,4400,4200,4900,4500,4900,4200,4150)
# phases2$ends   = c(4801,4601,4401,4201,4001,4501,4301,4201,4101,4001)
# time2phase(x,phases2)

time2phase  <- function(x,phases)
{
	res.index = numeric(length(x))
	for (i in 1:length(x))
	{
		ii = which(x[i] <= phases$starts & x[i] >= phases$ends)
		if (length(ii)==1){res.index[i]=ii} else {res.index[i]=sample(ii,size=1)}
	}
	res = phases[res.index,]
	row.names(res)=NULL
	return(res)
}


