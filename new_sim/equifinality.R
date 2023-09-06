library(nimbleCarbon)


r1  <-  0.01
r2  <-  0.027
a  <- dExponentialGrowth(800:301,r=r1,a=800,b=301,log=F)
b  <- dLogisticGrowth2(800:301,r=r2,a=800,b=301,m=500,log=F)
plot(800:301,a,type='l',xlim=c(800,301))
lines(800:301,b,col=2)
round(sum(a[1:300]),2)
round(sum(b[1:300]),2)


# Parameter Space Plot
r = seq(0.001,0.03,length.out=1000)
m = 801:302
params = expand.grid(r=r,m=m)
params$s <- NA
for (i in 1:nrow(params))
{
print(i)
params$s[i] <- round(sum(dLogisticGrowth2(800:500,r=params$r[i],a=800,b=301,m=params$m[i],log=F)),2)
}

sparam <- subset(params,s==0.13)
plot(params$r,params$m,col='lightgrey',pch=20)
points(sparam$r,sparam$m,col='indianred',pch=20)


c1  <- dLogisticGrowth2(800:301,r=sparam$r[1],a=800,b=301,m=sparam$m[1],log=F)
c2  <- dLogisticGrowth2(800:301,r=sparam$r[5483],a=800,b=301,m=sparam$m[5483],log=F)
c3  <- dLogisticGrowth2(800:301,r=sparam$r[3000],a=800,b=301,m=sparam$m[3000],log=F)
plot(800:301,c1,type='l',xlim=c(800,301),ylim=c(0,0.01))
lines(800:301,c2,col=2)
lines(800:301,c3,col=3)

