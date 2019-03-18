##' Species abundance data simulator
##'
##' Simulates species abundance data along a one-dimensional gradient
##'
##'
##' @param d the (optional) locations of the species along the 1-D gradient. If
##' d is given then it will define both the number of species and also the
##' locations on the gradient e.g. d = rep(1:10,each=3) will generate species at
##' locations 1,1,1,2,2,2,...,10,10,10. If d is not specified then d.rand = TRUE
##' will randomly allocate the species modes along a gradient on [0, 1], but if
##' d.rand = FALSE will uniformly distribute the species modes along a gradient.
##' @param p number of species.
##' @param n number of sites.
##' @param strip0 if TRUE the sites with zero total abundance are omitted.
##' @param extreme number typically in the range -1 to +1 with larger numbers
##' reducing the range of species.
##' @param ret if TRUE the generated data are returned
##' @param k.rand should the be random (TRUE) or fixed
##' @param d.rand should the be random (TRUE) or fixed
##' @param mu.rand should the be random (TRUE) or fixed
##' @param s the spans of the species response curves; s is the standard
##' deviation of the spread
##' @param amp the amplitudes of the species response curves
##' @param skew the skewness of the distribution; range (>0 to 5), 1 =
##' symmetric.
##' @param ampfun any function to modify the amplitude
##' @param lst if lst == TRUE then both the systematic and random values are
##' returned
##' @param err if err == 0 then the values are systematic with no random
##' variation
##' @param err.type type of error; p = poison, g = gaussian
##' @param as.df if return returns a data frame, otherwise a matrix
##' @param plotit if TRUE then the data are plotted
##' @param ptype species plot types e.g. "l" gives lines
##' @param plty species plot line types
##' @param pcols species plot colours
##' @param add.rug should a rug be added?
##' @param ... other arguments passed to plot.
##' @return if lst == FALSE then a data frame with variables "Locations",
##' "Taxa.1" -- "Taxa.N" where N is number of species. if lst == TRUE then two
##' data frames "x" and "xs" with variables "Locations", "Taxa.1" -- "Taxa.N"
##' and additionally, components "sigma", "amp" and "mu" that represent the
##' spans, amplitudes and locations of the N species along the 1-D gradient.
##' @examples
##'
##' mydata <- simdata()
##' summary(mydata)
##' mydata <- simdata(p=5, n=50, amp=1, err=0, d.rand=FALSE,
##' mu.rand=FALSE, plotit = TRUE)
##' summary(mydata)
##'
##' @importFrom grDevices rainbow
##' @importFrom graphics rug matplot
##' @importFrom stats runif rpois rnorm
##' @export
simdata <- function (d, p = 10, n = 100, strip0 = TRUE, extreme = 0, ret = TRUE, k.rand = FALSE,
                     d.rand = TRUE, mu.rand = TRUE, s = rep((4:8)/10, length = p), amp = c(sample(1:5),rep=TRUE,length=p),
                     skew = 1, ampfun, lst = FALSE, err = 1, err.type = c("p","n")[1], as.df = TRUE,
                     plotit = TRUE, ptype = "l", plty = 1, pcols = rainbow(p),  add.rug = FALSE, ...) {
  if (!missing(d)) {
    n <- length(d)
    d <- (d-min(d))/(max(d)-min(d))
  }
  else if (d.rand)
    d <- runif(n, 0, 1)
  else d <- (0:(n - 1))/(n - 1)
  amp <- rep(amp, length = p)
  skew <- rep(skew, length = p)
  s <- rep(s, length = p)
  plty <- rep(plty, length = p)
  pcols <- rep(pcols, length = p)
  x <- matrix(NA, nrow = n, ncol = p + 1)
  if (err) xs <- matrix(NA, nrow = n, ncol = p + 1)
  else xs <- NULL
  rowins <- rep(TRUE, p + 1)
  colins <- rep(TRUE, p + 1)
  if (mu.rand) {
    if (!extreme)
      mu <- runif(p, 0, 1)
    else mu <- runif(p, -0.1, 1.1)
  }
  else if (!extreme)
    mu <- (0:(p - 1))/(p - 1)
  else mu <- (0:(p - 1))/(p - 1) * extreme - extreme/2 + 0.5
  if (k.rand) {
    k1 <- sample(1:p, p, replace = TRUE)
    k2 <- sample(1:p, p, replace = TRUE)
    k3 <- sample(1:p, p, replace = TRUE)
  }
  else k1 <- k2 <- k3 <- rep(1,p)
  if (!missing(ampfun)) amp <- amp*ampfun(d)
  for (i in 1:p) {
    if (!err) {
      x[, i + 1] <- amp[k2[i]] * exp(-0.5 * ((d^(skew[k1[i]]) -
                                                mu[i])/s[k3[i]])^2)
    }
    else {
      xs[, i + 1] <- amp[k2[i]] * exp(-0.5 * ((d^(skew[k1[i]]) -
                                                 mu[i])/s[k3[i]])^2)
      if (err.type == "p")
        x[, i + 1] <- rpois(n, xs[, i + 1])
      else if (err.type == "n")
        x[, i + 1] <- xs[, i + 1] + err * rnorm(n, 0,
                                                1)
    }
  }
  x[, 1] <- d
  if (err) xs[, 1] <- d
  if (strip0) {
    rowins <- apply(x[, -1], 1, sum) != 0
    x <- x[rowins, ]
    if (err) xs <- xs[rowins, ]
    colins <- c(TRUE, (apply(x[, -1], 2, sum) != 0))
    x <- x[, colins]
    if (err) xs <- xs[, colins]
  }
  ord <- order(x[, 1])
  x <- x[ord, ]
  if (err) xs <- xs[ord, ]
  dimnames(x) <- list(as.character(1:length(x[, 1])), c("Locations",
                                                        paste("Taxa", 1:(dim(x)[2] - 1), sep = ".")))
  if (err) dimnames(xs) <- dimnames(x)
  if (plotit) {
    z <- rev(order(x[, 1]))
    matplot(x[z, 1], as.matrix(x[z, -1]), type = ptype, lty = plty,
            col = pcols, xlab = "Locations", ylab = "Response",
            main = "Simulated Data")
    if (add.rug)
      rug(x[z, 1], side = 3)
  }
  if (ret) {
    if (lst)
      list(x = x, xs = xs, sigma = s[k1[1:p]][colins[-1]],
           amp = amp[k2[1:p]][colins[-1]], mu = mu[k3[1:p]][colins[-1]])
    else if (as.df)
      data.frame(x)
    else x
  }
}

