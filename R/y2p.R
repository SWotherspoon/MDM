##' Scale rows of a species-sites data set
##'
##' The response matrix for mdm diversity analyses requires that species
##' (columns of the data matrix) are scaled to proportions i.e. that they sum to
##' one for each site (row of the data matrix). Two matrices are particularly
##' useful; these correspond to alpha and gamma diversities. The former is
##' particularly useful for \code{\link{mdm}} anovas that include the alpha
##' model which is computationally expensive using the function
##' \code{\link{mdm}}.
##'
##'
##' @param y matrix or data frame of numeric values to be transformed
##' @param mean if mean = TRUE then each row is replaced by the species means
##' scaled to row sums of one
##' @return A matrix of the same dimensions as the input matrix. Each row of the
##' matrix will sum to one.
##' @examples
##'
##' mydata <- matrix(0:8,nrow=3,ncol=3)
##' mydata
##' y2p(mydata)
##'
##' @export
y2p <- function(y, mean = FALSE) {
  dn <- dimnames(y)
  y  <- as.matrix(y)/apply(y,1,sum)
  if (mean) y <- matrix(rep(apply(y,2,mean),each=nrow(y)),nrow=nrow(y),ncol=ncol(y),dimnames=dimnames(y))
  dimnames(y) <- dn
  y
}
