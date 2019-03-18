##' Calculate alpha, beta and gamma parametric entropies and diversities
##'
##' Calculates alpha, beta and gamma parametric diversities of a data set. The
##' matrix or data frame is automatically scaled to row sums of one. Unlike true
##' diversities where weighting is done through the function ed, weighting for
##' parameterized diversity is done within the MDM, or more genrelaly by using
##' case weights.
##'
##'
##' @param x the input matrix or data frame.
##' @param w weights if required.
##' @param q the order of diversity; typically 0, 1 or 2.
##' @param retq if TRUE then parametric diversities are returned; if FALSE the
##' entropies for alpha and gamma are returned.
##' @return a vector of entropies or diversities
##' @seealso \code{\link{dev2div}}, \code{\link{ed1}}, \code{\link{ed}},
##' \code{\link{eds1}}
##' @examples
##'
##' data(spider6)
##' eds(spider6[,1:6])
##' eds(spider6[,1:6],q=0)
##' eds(spider6[,1:6],q=2)
##' eds(spider6[,1:6],retq=FALSE)
##' data(spider6)
##' eds(spider6[,1:6])
##'
##' @export
eds <- function(x, q = 1, w = 1, retq = TRUE) {
  bsums <- function(x, w) {
    if (all(x == 0)) return(0)
    sum(-w*x*log(x), na.rm=TRUE)
  }
  if (is.null(dim(x))) x <- matrix(x,nrow=1)
  x[x!=0] <- x[x!=0]^q
  rs <- rowSums(x)
  if (length(w) != 1 & length(w) != nrow(x))
    cat("Warning: n weights NE n rows of x !!")
  if (any(rs==0)) {
    drops <- (1:nrow(x))[rs==0]
    rs <- rs[rs!=0]
    w <- w[rs != 0]
    x <- x[rs!=0,]
    cat("Dropping zero sum rows: ",drops,"\n")
  }
  x <- x/rs
  wa <- w^q/mean(w^q)
  a <- bsums(x, w = wa)/nrow(x)
  wg <- w/mean(w)
  g <- bsums(colMeans(wg * x, na.rm = TRUE), w = 1)
  if (retq) {
    a <- exp(a)
    g <- exp(g)
    c(alpha=a, beta=g/a, gamma=g)
  }
  else c(absums=a, gbsums=g)
}

