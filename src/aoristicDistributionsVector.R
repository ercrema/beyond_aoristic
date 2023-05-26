library(nimble)

dAoristicExponentialGrowth_vector=nimbleFunction(
  run = function(x = double(2),z=integer(0),r=double(0), log = integer(0)) {
    returnType(double(0))
    t = 1:z
    n = numeric(z)
    for (i in 1:z)
    {
      n[i] = (1+r)^t[i]
    }
    p = n/sum(n)
    pg = x %*% p
    logProb = sum(log(pg))
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

registerDistributions(list(
  dAoristicExponentialGrowth_vector = list(
    BUGSdist = "dAoristicExponentialGrowth_vector(z,r)",
    pqAvail = FALSE,
    types = c('value = double(2)', 'z = integer(0)','r=double(0)')
  )))


dAoristicDoubleExponentialGrowth_vector=nimbleFunction(
  run = function(x = double(2),z=integer(0),r1=double(0),r2=double(0),mu=double(0), log = integer(0)) {
    returnType(double(0))
    mu = round(mu)
    t1 = 1:(mu-1)
    t2 = 1:(z-mu+1)
    t1final = mu-1
    t2final = z-mu+1
    n1 = numeric(mu-1)
    n2 = numeric(z-mu+1)
    for (i in 1:t1final)
    {
      n1[i] = (1+r1)^t1[i]
    }
    for (i in 1:t2final)
    {
      n2[i] = ((1+r1)^(mu-1))  * (1+r2)^t2[i]
    }
    n = c(n1,n2)
    p = n/sum(n)
    pg = x %*% p
    logProb = sum(log(pg))
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   

registerDistributions(list(
  dAoristicDoubleExponentialGrowth_vector = list(
    BUGSdist = "dAoristicDoubleExponentialGrowth_vector(z,r1,r2,mu)",
    pqAvail = FALSE,
    types = c('value = double(2)', 'z = integer(0)','r1=double(0)','r2=double(0)','mu=double(0)')
  )))

dAoristicLogisticGrowth_vector=nimbleFunction(
  run = function(x = double(2),z=integer(0),k=double(0),r=double(0), log = integer(0)) {
    returnType(double(0))
    t = 1:z
    n = 1/(1+((1-k)/k)*exp(-r*t))
    p = n/sum(n)
    pg = x %*% p
    logProb = sum(log(pg))
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
  })   


registerDistributions(list(
  dAoristicLogisticGrowth_vector = list(
    BUGSdist = "dAoristicLogisticGrowth_vector(z,k,r)",
    pqAvail = FALSE,
    types = c('value = double(2)', 'z = integer(0)','k=double(0)','r=double(0)')
  )))



