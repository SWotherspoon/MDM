##' Fits the parameteric diversity model
##'
##' The parametric diversity model (mdm) is a new method for directly relating
##' diversity to environmental predictors. It is based on three components: (1)
##' parametric diversity, a new parametric model of diversity that can represent
##' any configuration of species proportional abundances, (2) the multinomial
##' logit model (MLM) that can relate species proportional abundances to complex
##' predictors and (3) the link between parametric diversity and the likelihood
##' function of the MLM.
##'
##' The \code{mdm} is fitted using the \code{multinom} package from \code{nnet}.
##' Parametric diversities and true diversities can also be calculated for data
##' sets using the functions \code{\link{eds}}, \code{\link{eds1}},
##' \code{\link{ed}}, \code{\link{ed1}}.
##'
##' \code{\link{mdm}} calls \code{\link[nnet]{multinom}} which calls
##' \code{\link{nnet}}. The variables on the rhs of the formula should be
##' roughly scaled to [0,1] or the fit will be slow or may not converge at all.
##'
##' @param formula a formula expression as for regression models, of the form
##' response ~ predictors. The response should be a matrix with K columns
##' comprising proportions for each of K classes. A log-linear model is fitted,
##' with coefficients zero for the first class. An offset can be included: it
##' should be a numeric matrix with K columns. See the documentation of
##' formula() for other details.
##' @param data an optional data frame, list or environment (or object coercible
##' by as.data.frame to a data frame) containing the variables in the model. If
##' not found in data, the variables are taken from environment(formula),
##' typically the environment from which \code{\link{mdm}} is called.
##' @param weights optional case weights in fitting.
##' @param subset expression saying which subset of the rows of the data should
##' be used in the fit. All observations are included by default.
##' @param na.action a function to filter missing data.
##' @param MaxNWts The maximum allowable number of weights. There is no limit in
##' the code, but MaxNWts is set to the exact number of required as specified in
##' the formula. Thus it should not need to be changed when fitting
##' \code{\link{mdm}}.
##' @param maxit maximum number of iterations. Default 1000.
##' @param contrasts a list of contrasts to be used for some or all of the
##' factors appearing as variables in the model formula.
##' @param Hess logical for whether the Hessian (the O/E information matrix)
##' should be returned.
##' @param censored If Y is a matrix with K > 2 columns, interpret the entries
##' as one for possible classes, zero for impossible classes, rather than as
##' counts.
##' @param model logical. If true, the model frame is saved as component model
##' of the returned object.
##' @param use.shortcut logical. If true, and the model is ~1 (a constant) or
##' ~sites (a factor with one level for each site) then the model is not fitted
##' since the fitted values are known in each case. The first (alpha) model has
##' fitted values equal to the input data and the second (gamma) model fits the
##' row means. Fitting the alpha-model using nnet can be prohibitively expensive
##' in computational time and is unneccessary. The returned models when
##' use.shortcut == TRUE has the same components as when use.shortcut == FALSE,
##' and hence can be used in anova tables and plotting. Using use.shortcut ==
##' TRUE can result in saving > 99\% of computational time for a collection of
##' models.
##' @param \dots additional arguments for nnet.
##' @return A nnet object with additional components:
##'
##' \item{deviance}{the residual deviance, compared to the full saturated model
##' (that explains individual observations exactly). Also, minus twice
##' log-likelihood} \item{edf}{the (effective) number of degrees of freedom used
##' by the model} \item{AIC}{the AIC for this fit} \item{Hessian}{if Hess is
##' true} \item{model}{if model is true} \item{entropy}{the entropy of the
##' fitted values} \item{diversity}{the diversity of the fitted values}
##' @note \code{\link{mdm}} is a modifed version of \code{\link[nnet]{multinom}}
##' in the \code{\link{nnet}} package.
##' @seealso \code{\link[nnet]{multinom}}, \code{\link{nnet}}
##' @references De'ath, G. (2011) \emph{The Multinomial Diversity Model: Linking
##' Shannon Diversity To Multiple Predictors}.
##'
##' Ripley, B. D. (1996) \emph{Pattern Recognition and Neural Networks}.
##' Cambridge.
##'
##' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
##' with S.} Fourth edition. Springer.
##' @keywords mdm
##' @examples
##'
##' library(MDM)
##' data(spider6)
##' fit0 <- mdm(y2p(spider6[,1:6])~1,data=spider6)
##' fit1 <- mdm(y2p(spider6[,1:6])~Water,data=spider6)
##' fit2 <- mdm(y2p(spider6[,1:6])~Water+Herbs,data=spider6)
##' fit3 <- mdm(y2p(spider6[,1:6])~Site,data=spider6)
##' anova(fit0,fit1,fit2,fit3)
##'
##' @importFrom stats model.matrix model.response model.weights model.offset
##' @importFrom stats .getXlevels weighted.mean
##' @importFrom nnet nnet.default
##' @export
mdm <- function (formula, data, weights, subset, na.action,
                 MaxNWts, maxit = 1000, contrasts = NULL, Hess = FALSE,
                 censored = FALSE, model = TRUE, use.shortcut = TRUE, ...) {

  # file nnet/multinom.R
  # copyright (C) 1994-2006 W. N. Venables and B. D. Ripley
  #
  #  This program is free software; you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation; either version 2 or 3 of the License
  #  (at your option).
  #
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.
  #
  #  A copy of the GNU General Public License is available at
  #  http://www.r-project.org/Licenses/
  #

  # file MDM/mdm.R
  # mdm is modified version of multinom from nnet package by B.R.Ripley
  # and Bill Venables

  class.ind <- function(cl) {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1L:n) + n * (as.vector(unclass(cl)) - 1L)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$MaxNWts <- m$maxit <- m$summ <- m$Hess <- m$contrasts <- m$censored <- m$q <- m$model <- m$use.shortcut <- m$... <- NULL
  m[[1L]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  X <- model.matrix(Terms, m, contrasts)
  cons <- attr(X, "contrasts")
  Xr <- qr(X)$rank
  Y <- model.response(m)
  w <- model.weights(m)
  if (length(w) == 0L)
    if (is.matrix(Y))
      w <- rep(1, dim(Y)[1L])
  else w <- rep(1, length(Y))
  offset <- model.offset(m)
  r <- ncol(X)
  if (is.matrix(Y)) {
    use.shortcut <-  use.shortcut & ((Xr==1)|(Xr==nrow(X)))
    p <- ncol(Y)
    sY <- Y %*% rep(1, p)
    if (any(sY == 0))
      stop("some case has no observations")
    Y <- Y/matrix(sY, nrow(Y), p)
    w <- w * sY
    if (use.shortcut) maxit <- 0
    if (missing(MaxNWts)) MaxNWts <- (r + 1) * p
    if (length(offset) > 1L) {
      if (ncol(offset) != p)
        stop("ncol(offset) is wrong")
      mask <- c(rep(FALSE, r + 1L + p), rep(c(FALSE, rep(TRUE,
                                                         r), rep(FALSE, p)), p - 1L))
      X <- cbind(X, offset)
      Wts <- as.vector(rbind(matrix(0, r + 1L, p), diag(p)))
      fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask,
                          size = 0, skip = TRUE, softmax = TRUE, censored = FALSE,
                          rang = 0, maxit = maxit, MaxNWts = MaxNWts, ..., PACKAGE="nnet")
    }
    else {
      mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE,
                                                     r)), p - 1L))
      fit <- nnet.default(X, Y, w, mask = mask, size = 0,
                          skip = TRUE, softmax = TRUE, censored = FALSE,
                          rang = 0, maxit = maxit, MaxNWts = MaxNWts, ..., PACKAGE="nnet")
    }
    if (use.shortcut) {
      if (Xr==1) {
        if (missing(weights)) fit$fitted.values  <- matrix(apply(Y, 2, mean), nrow=nrow(Y),
                                                           ncol=ncol(Y), byrow=TRUE, dimnames=dimnames(Y))
        else fit$fitted.values  <- matrix(apply(Y, 2, weighted.mean, w = c(w)), nrow=nrow(Y),
                                          ncol=ncol(Y), byrow=TRUE, dimnames=dimnames(Y))
      }
      else if (Xr==nrow(Y)) fit$fitted.values  <-  Y
      if (missing(weights)) fit$value <- sum(-fit$fitted.values * log(fit$fitted.values), na.rm = TRUE)
      else fit$value <- sum(-fit$fitted.values * log(fit$fitted.values) * c(w), na.rm = TRUE)
    }
  }
  else {
    if (length(offset) <= 1L) {
      mask <- c(FALSE, rep(TRUE, r))
      fit <- nnet.default(X, Y, w, mask = mask, size = 0,
                          skip = TRUE, entropy = TRUE, rang = 0, maxit = maxit,
                          MaxNWts = MaxNWts, ..., PACKAGE="nnet")
    }
    else {
      mask <- c(FALSE, rep(TRUE, r), FALSE)
      Wts <- c(rep(0, r + 1L), 1)
      X <- cbind(X, offset)
      fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask,
                          size = 0, skip = TRUE, entropy = TRUE, rang = 0,
                          maxit = maxit, MaxNWts = MaxNWts, ..., PACKAGE="nnet")
    }
  }
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$terms <- Terms
  fit$call <- call
  fit$weights <- w
  fit$deviance <- 2 * fit$value
  fit$entropy <-  fit$deviance/2/sum(c(w))
  fit$diversity <- exp(fit$entropy)
  fit$rank <- Xr
  edf <- (ncol(Y) - 1) * Xr
  if (length(dn <- colnames(Y)) > 0)
    fit$lab <- dn
  else fit$lab <- 1L:ncol(Y)
  fit$coefnames <- colnames(X)
  fit$vcoefnames <- fit$coefnames[1L:r]
  fit$na.action <- attr(m, "na.action")
  fit$contrasts <- cons
  fit$xlevels <- .getXlevels(Terms, m)
  fit$edf <- edf
  fit$AIC <- fit$deviance + 2 * edf
  if (model) fit$model <- m
  class(fit) <- c("mdm", "multinom", "nnet")
  fit
}

