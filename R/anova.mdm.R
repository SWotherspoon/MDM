##' Analysis of Deviance, Entropy and Diversity Tables
##'
##' Provides an analysis of deviance, entropy and diversity for a collection of
##' diversity models (outputs from \code{mdm}).
##'
##' Specifying a single object gives a sequential analysis of deviance table for
##' that fit. That is, the reductions in the residual sum of deviances, plus the
##' residual deviance The deviances are converted into entropies and
##' diversities, and differences in deviances are converted into differnces in
##' entropies and diversities.
##'
##' If more than one object is specified, the table has a row for the residual
##' degrees of freedom and sum of deviances for each model.  For all but the
##' first model, the change in degrees of freedom and deviances. This only makes
##' statistical sense if the models are nested.  It is conventional to list the
##' models from smallest to largest, but this is up to the user.
##'
##' @param object, objects of class \code{mdm}, usually, a result of a call to
##' \code{\link{mdm}}.
##' @param ... objects of class \code{mdm}, usually, a result of a call to
##' \code{\link{mdm}}.
##' @param topnote If TRUE then model descriptions appear above the anova table,
##' if FALSE they appear as the first column of the table
##' @param cols The list of colums to print out. Defaults to all columns.
##' @return The analysis of deviance, entropy and diversity for a collection of
##' diversity models.
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
##' @importFrom stats deviance formula pchisq
##' @export
anova.mdm <- function(object, ..., topnote = TRUE, cols = c("df","dev","pr","ent","div")) {
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
             rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named))
    warning("the following arguments to 'anova.mdm' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.mdm <- unlist(lapply(dotargs, function(x) inherits(x,
                                                        "mdm")))
  dotargs <- dotargs[is.mdm]
  object <- c(list(object), dotargs)
  nt <- length(object)
  dflis <- sapply(object, "[[",'edf')
  s <- order(dflis)
  dflis <- nrow(object[[1]]$residuals) * (ncol(object[[1]]$residuals) - 1) - dflis
  object <- object[s]
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  rsp <- unique(sapply(object, function(x) paste(formula(x)[2L])))
  mds <- sapply(object, function(x) paste(formula(x)[3L]))
  nr <- nrow(object[[1]]$residuals)
  dfs <- dflis[s]
  lls <- sapply(object, function(x) deviance(x))
  tss <- c("", paste(1L:(nt - 1), 2L:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1L], df[-1L]))
  ent <- sapply(object, "[[",'entropy')
  dent <- c(NA, -diff(ent))
  div <- exp(ent)
  ddiv <- exp(dent)
  ins <- !is.na(match(rep(c("df","dev","pr","ent","div"),times=c(2,2,1,2,2)),cols))
  variables <- lapply(object, function(x) paste(deparse(formula(x)),
                                                collapse = "\n"))
  top <- paste("Model ", format(1L:nt), ": ", variables,
               sep = "", collapse = "\n")
  out <- data.frame(Resid.df = dfs, df = df, Deviance = lls, ddev = x2, pr = pr,
                    ent = ent, dent = dent, div = div, ddiv = ddiv)[,ins]
  names(out) <- c("DF", "DF-Diff", "Dev", "Dev-Diff", "Pr",
                  "Ent", "Ent-Diff","Div", "Div-Ratio")[ins]
  rownames(out) <- as.character(1:nt)
  if (!topnote) {
    rownames(out) <- paste("Model ", format(1L:nt), ": ", variables,
                           sep = "")
    attr(out, "heading") <- c("Deviances, Entropies and Diversities of Parametric Diversity Models\n")
  }
  else
    attr(out, "heading") <- c("Deviances, Entropies and Diversities of Parametric Diversity Models\n",
                              paste("Response:", rsp,"\n"),paste(top,"\n"))
  class(out) <- c("anova", "data.frame")
  out
}

