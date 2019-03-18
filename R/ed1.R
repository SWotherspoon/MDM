##' Calculate alpha, beta and gamma true entropies and diversities
##'
##' Calculates true diversities of individual sites of a data set
##'
##'
##' @param x the input matrix or data frame.
##' @param q the order of diversity; typically 0, 1 or 2.
##' @param retq if TRUE then diversities are returned; if FALSE the entropies
##' for alpha and gamma are returned.
##' @return a vector of entropies or diversities
##' @seealso \code{\link{dev2div}}, \code{\link{ed1}}, \code{\link{eds}},
##' \code{\link{eds1}}
##' @examples
##'
##' data(spider6)
##' ed1(spider6[,1:6])
##' ed1(spider6[,1:6],q=0)
##' ed1(spider6[,1:6],q=2)
##' ed1(spider6[,1:6],retq=FALSE)
##' data(spider6)
##' ed1(spider6[,1:6])
##'
##' @export
ed1 <- function(x, q=1, retq = TRUE) {
  bsums <- function(x, q=1) {
    if (all(x==0)) return(0)
    if (q==0) sum(x>0)
    else if (q==1) sum(-x*log(x),na.rm=TRUE)
    else sum(x^q)
  }
  if (!is.null(dim(x))) a <- apply(x,1,ed1,q=q,retq=retq)
  else {
    x <- matrix(x,nrow=1)
    rs <- sum(x)
    if (rs==0) stop("Zero data")
    x <- x/rs
    a <- bsums(x,q=q)
    if (retq) {
      if (q==1) {
        a <- exp(a)
      }
      else if (q!=1) {
        a <- a^(1/(1-q))
      }
    }
  }
  a
}

