##' Entropy plot
##'
##' Plots components of entropy by species or sites.
##'
##'
##' @param lst a list of \code{\link{mdm}} fitted models
##' @param y optional -- locations of the components on the vertical axis
##' @param x optional -- used to provide the range of the horizontal axis
##' @param ord if TRUE orders the components by total entropy of each fitted
##' model.
##' @param type should the plot break down entropy by species (columns =
##' default) or sites (rows).
##' @param labs an (optional) vector of labels for species (or sites if type ==
##' "sites").
##' @param segs use segments to plot end of bars as vertical lines
##' @param wid.seg widths of bar ends
##' @param pchs some nice plotting characters
##' @param \dots other arguements passed to the plot function (e.g. cex).
##' @return a list of plotted values of the length of lst, each component of
##' which is of length equal to the number of species (sites).
##' @seealso \code{\link{mdm}}
##' @examples
##'
##' data(spider6)
##' fit0 <- mdm(y2p(spider6[,1:6])~1,data=spider6)
##' fit1 <- mdm(y2p(spider6[,1:6])~Water+Herbs,data=spider6)
##' fit2 <- mdm(y2p(spider6[,1:6])~Site,data=spider6,alpha=TRUE)
##' anova(fit0,fit1,fit2)
##' entropy.plot(list(fit0,fit2,fit1))
##'
##' @importFrom graphics axis box plot points segments
##' @export
entropy.plot <- function(lst, y, x, ord=TRUE, type=c("species","sites")[1], labs,
                         segs = TRUE, wid.seg=0.25, pchs = c(15,19,0,2,6,2,15,17,21,17), ...) {
  plnp <- function(x) {
    x <- -x*log(x)
    x[is.na(x)] <- 0
    x
  }

  if (class(lst)!="list") stop("'lst' must be a list of fitted values from mdms !!")
  n <- length(lst)
  if (n<2) stop("needs at least two matrices")
  types <- sapply(lst,function(x) class(x)[1])
  for (i in 1:length(types)) {
    if (types[i]=="mdm") lst[[i]] <- lst[[i]]$fitted.values
  }
  rc <- (type=="species")+1
  nplt <- dim(lst[[1]])[rc]
  if (missing(labs)) labs <- dimnames(lst[[1]])[[rc]]

  for (i in 1:n) lst[[i]] <- plnp(lst[[i]])

  if(ord) {
    ords <- order(apply(lst[[1]],rc,mean))
    labs <- labs[ords]
  }
  else ords <- 1:nplt

  for (i in 1:n) lst[[i]] <- apply(lst[[i]],rc,mean)[ords]
  if (missing(y)) y <- 1:length(lst[[1]])
  if (missing(x)) xlim <- range(unlist(lst))
  plot(lst[[1]],y,type="n",xlab="Entropy",ylab="",axes=FALSE, xlim=xlim, ...)
  axis(2,at=y,labels=labs, ...)
  axis(1)
  box()
  # could use arrows below but I want to line thickness of bars #
  if (segs) {
    segments(lst[[1]],y+wid.seg,lst[[1]],y-wid.seg,lty=1,lwd=1)
    segments(lst[[2]],y+wid.seg,lst[[2]],y-wid.seg,lty=1,lwd=1)
    segments(lst[[1]],y,lst[[2]],y,lty=1)
    if (n>2) for (i in 3:n) points(lst[[i]],y,pch=pchs[i-2],...)
  }
  else {
    segments(lst[[1]],y,lst[[2]],y,lty=1, col="grey75")
    for (i in 1:n) points(lst[[i]],y,pch=pchs[i],...)
  }
  invisible(lst)
}

