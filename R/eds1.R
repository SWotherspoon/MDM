##' Calculate alpha, beta and gamma parametric entropies and diversities
##'
##' Calculates parametric diversities of individual sites of a data set
##'
##'
##' @param x the input matrix or data frame.
##' @param q the order of diversity; typically 0, 1 or 2.
##' @param retq if TRUE then parametric diversities are returned; if FALSE the
##' entropies for alpha and gamma are returned.
##' @return a vector of entropies or diversities
##' @seealso \code{\link{dev2div}}, \code{\link{ed1}}, \code{\link{eds}},
##' \code{\link{eds1}}
##' @examples
##'
##' data(spider6)
##' eds1(spider6[,1:6])
##' eds1(spider6[,1:6],q=0)
##' eds1(spider6[,1:6],q=2)
##' eds1(spider6[,1:6],retq=FALSE)
##' data(spider6)
##' eds1(spider6[,1:6])
##'
##' @export
eds1 <- function(x, q=1, retq = TRUE) {
  bsums <- function(x) {
    if (all(x==0)) return(0)
    sum(-x*log(x),na.rm=TRUE)
  }
  x[x!=0] <- x[x!=0]^q
  if (!is.null(dim(x))) a <- apply(x,1,eds1,retq=retq)
  else {
    x <- matrix(x,nrow=1)
    rs <- sum(x)
    if (rs==0) stop("Zero data")
    x <- x/rs
    a <- bsums(x)
    if (retq) {
      a <- exp(a)
    }
  }
  a
}

