phaserect <- function(x,w=0.15)
{
	ybottom=quantile(par('usr')[3:4],1-w)
	ytop=par('usr')[4]
	for (i in 1:nrow(x))
	{
		rect(ybottom=ybottom,ytop=ytop,xleft=x[i,2],xright=x[i,3],border=ifelse(i%%2,'lightgrey','darkgrey'),col=ifelse(i%%2,'lightgrey','darkgrey'))
	}
}
