# Load Libraries, Functions, and Data from Crema and Kobayashi 2020 ----
library(here)
library(trapezoid)
# scripts
source("https://raw.githubusercontent.com/ercrema/jomonPhasesAndPopulation/master/R/mcsim.R")
source("https://raw.githubusercontent.com/ercrema/jomonPhasesAndPopulation/master/R/utilities.R")
# R images
load(url("https://github.com/ercrema/jomonPhasesAndPopulation/raw/master/R_images/pithouseData.RData"))
load(url("https://github.com/ercrema/jomonPhasesAndPopulation/raw/master/R_images/posteriorSamples.RData"))


# Sample dates via MCMC ----
nsim  <- 10000
simTrapezoid  <- mcsim(pithouseData[,-c(3,4)],nsim=nsim,posterior=postTrapezoid,weights="variance") 

# Generate probability matrix ----
resolution  <- 10
minMax  <- c(ceiling(max(simTrapezoid)/100)*100,floor(min(simTrapezoid)/100)*100) 
midPoints  <- seq(minMax[2]+resolution/2,minMax[1]-resolution/2,resolution) 
mat  <- t(apply(simTrapezoid,1,function(x,minMax,resolution,nsim){(hist(x,breaks=seq(minMax[1],minMax[2],-resolution),plot=F)$counts)/nsim},minMax=minMax,resolution=resolution,nsim=nsim))
colnames(mat) <- midPoints
# barplot(rev(apply(mat,2,sum)))

# Subsample to 7000 to 3000
ii <- which(midPoints>=7000)[1]:which(midPoints<=3000)[sum(midPoints<=3000)]
cumultaive.prob  <- apply(mat[,ii],1,sum)
mat  <- mat[which(cumultaive.prob>0.8),ii]
midPoints  <- as.numeric(colnames(mat))
housedb  <- prop.table(mat,margin = 1)
# barplot(apply(mat,2,sum))

save(housedb,midPoints,file=here('data','housedata.RData'))
