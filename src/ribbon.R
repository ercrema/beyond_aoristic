ribbon  <- function(x,y,index=1,interval=0.95,col=rgb(0,0,0,0.2))
{
	lo = apply(y,index,quantile,prob=(1-interval)/2)
	hi = apply(y,index,quantile,prob=interval+(1-interval)/2)
	polygon(c(x,rev(x)),c(lo,rev(hi)),border=NA,col=col)
}
