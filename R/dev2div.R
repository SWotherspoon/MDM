##' Converts deviances to diversities
##'
##' Takes either (1) a mdm or (2) a scalar, vector or matrix of deviances, and
##' extracts or converts then to diverities. The relationship between deviance
##' (dev) and diversity (d) is given by div = exp(dev/2/n) where n is the number
##' of units (typically rows of a matrix) over which deviance is being averaged.
##'
##'
##' @param x a mdm or a scaler, vector or matrix of deviances.
##' @param n if x is not a mdm, then the divisor in the conversion as defined as
##' above.
##' @return The diversity of x.
##' @seealso \code{\link{ed}}, \code{\link{ed1}}, \code{\link{eds}},
##' \code{\link{eds1}}
##' @examples
##'
##' x <- c(5,10,15)
##' dev2div(x,n=10)
##'
##' @export
dev2div <- function(x, n)  {
  if (class(x)[[1]]=="mdm") exp(x$deviance/2/nrow(x$fitted.values))
  else exp(x/2/n)
}

