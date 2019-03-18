##' Calculate alpha, beta and gamma true entropies and diversities
##'
##' Calculates alpha, beta and gamma true diversities of a data set. The matrix
##' or data frame is automatically scaled to row sums of one.
##'
##'
##' @param x the input matrix or data frame.
##' @param q the order of diversity; typically 0, 1 or 2.
##' @param w weights if required.
##' @param retq if TRUE then diversities are returned; if FALSE the entropies
##' for alpha and gamma are returned.
##' @return a vector of entropies or diversities
##' @seealso \code{\link{dev2div}}, \code{\link{ed1}}, \code{\link{eds}},
##' \code{\link{eds1}}
##' @examples
##'
##' data(spider6)
##' ed(spider6[,1:6])
##' ed(spider6[,1:6],q=0)
##' ed(spider6[,1:6],q=2)
##' ed(spider6[,1:6],retq=FALSE)
##' data(spider6)
##' ed(spider6[,1:6])
##'
##' @export
ed <- function (x, q = 1, w = 1, retq = TRUE) {
  bsums <- function(x, q = 1, w = 1) {
    if (all(x == 0))
      return(0)
    if (q == 0)
      sum(x > 0)
    else if (q == 1)
      sum(-w * x * log(x), na.rm = TRUE)
    else sum(w * x^q)
  }
  rs <- rowSums(x)
  if (length(w) != 1 & length(w) != nrow(x))
    cat("Warning: n weights NE n rows of x !!")
  if (any(rs == 0)) {
    drops <- (1:nrow(x))[rs == 0]
    rs <- rs[rs != 0]
    x <- x[rs != 0,]
    w <- w[rs != 0]
    cat("Dropping zero sum rows: ", drops, "\n")
  }
  x <- x/rs
  wa <- w^q/mean(w^q)
  a <- bsums(x, q = q, w = wa)/nrow(x)
  wg <- w/mean(w)
  g <- bsums(colMeans(wg * x, na.rm = TRUE), q = q, w = 1)
  if (retq) {
    if (q == 1) {
      a <- exp(a)
      g <- exp(g)
    }
    else if (q!= 1) {
      a <- a^(1/(1 - q))
      g <- g^(1/(1 - q))
    }
    c(alpha = a, beta = g/a, gamma = g)
  }
  else c(absums=a, gbsums=g)
}

