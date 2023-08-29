diristick<-function(alpha,timeRange=c(5000,3001))
{
  n <- length(alpha)
  x <- matrix(rgamma(n, alpha), ncol = n , byrow = TRUE)
  res <- x %*% rep(1, n)
  res <-  x/as.vector(res)
  res <- cumsum(res)
  dur <- timeRange[1]-timeRange[2]
  st  <- c(timeRange[1],ceiling(timeRange[1]-res[1:(length(res)-1)]*dur))
  en  <- c(st[-1]+1,timeRange[2])
  return(data.frame(phasename=letters[1:length(alpha)],starts=st,ends=en))
}
