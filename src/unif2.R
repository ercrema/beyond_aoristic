#' @title Uniform Distribution
#' @aliases dunif2 runif2
#' @name  dunif2
#' @description  Density and random generation of a Uniform distribution with altenrative parametrisation.
#' @param x value to be computed.
#' @param m location parameter.
#' @param s scale parameter.
#' @param log TRUE or 1 to return log probability. FALSE or 0 to return probability.
#' @param n number of random draws. Currently only n = 1 is supported, but the argument exists for standardization of "r" functions.
#' @details The function provides an alternative parametrisation for the uniform distribution, so that the upper and lower limits of the distribution is equal to \code{m} - \code{s} and \code{m} + \code{s}.
#' @author Enrico Crema

NULL

#' @rdname dunif2
#' @import nimble
#' @export

dunif2=nimbleFunction(
  run = function(x = integer(0),m=double(0),s=double(0), log = integer(0)) {
    returnType(double(0))
    a = m - s
    b = m + s
    dens = dunif(x,a,b,log=log)
    return(dens)
  })   

#' @rdname dTrapezoidal
#' @export
#' 
runif2=nimbleFunction(
  run = function(n = integer(0),m=double(0),s=double(0)) {
    returnType(double(0))
    a = m - s
    b = m + s
    res = runif(1,a,b)
    return(res)
  })

registerDistributions(list(
  dunif2 = list(
    BUGSdist = "dunif2(m,s)",
    Rdist = "dunif2(m,s)",
    types = c('m=double(0)','s=double(0)'),
    pqAvail = FALSE
  )))
