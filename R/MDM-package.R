

##' Multinomial Diversity Model
##'
##' The multinomial diversity model is a toolbox for relating diversity to
##' complex predictors. It is based on (1) Shannon diversity; (2) the
##' multinomial logit model, and (3) the link between Shannon diversity and the
##' log-likelihood of the MLM.
##'
##' \tabular{ll}{ Package: \tab MDM\cr Type: \tab Package\cr Version: \tab
##' 1.0\cr Date: \tab 2011-09-08\cr License: \tab GPL (version 2 or newer)\cr
##' LazyLoad: \tab yes\cr }
##'
##' @name MDM-package
##' @aliases MDM-package MDM
##' @docType package
##' @author Glenn De'ath: \email{g.death@@aims.gov.au}
##' @references De'ath, G. (2011) \emph{The Multinomial Diversity Model: Linking
##' Shannon Diversity To Multiple Predictors}
##' @keywords diversity
##' @examples
##'
##' library(MDM)
##' data(spider6)
##' fit0 <- mdm(y2p(spider6[,1:6])~1,data=spider6)
##' fit1 <- mdm(y2p(spider6[,1:6])~Water,data=spider6)
##' fit2 <- mdm(y2p(spider6[,1:6])~Water+Herbs,data=spider6)
##' fit3 <- mdm(y2p(spider6[,1:6])~Site,data=spider6,alpha=TRUE)
##' anova(fit0,fit1,fit2,fit3)
##'
NULL





##' The spider data set
##'
##' Data set on abundances of spiders and environmental predictors. This is a
##' subset of a larger data comprising 12 species and 6 environmental
##' predictors. All variables are rated on a 0-9 scale.
##'
##'
##' @name spider6
##' @docType data
##' @format A data frame with 28 observations on the following 9 variables.
##' \describe{
##'   \item{Pard.lugu}{a numeric vector}
##'   \item{Pard.pull}{a numeric vector}
##'   \item{Troc.terr}{a numeric vector}
##'   \item{Pard.mont}{a numeric vector}
##'   \item{Alop.acce}{a numeric vector}
##'   \item{Alop.fabr}{a numeric vector}
##'   \item{Water}{a numeric vector}
##'   \item{Herbs}{a numeric vector}
##'   \item{Site}{a factor with 28 levels}
##' }
##' @references De'ath G. (2002) \emph{Multivariate regression trees: A new
##' technique for modeling species-environment relationships}. Ecology, 2002,
##' 83:1105--1117.
##' @source package mvpart
##' @keywords datasets
##' @examples
##'
##' data(spider6)
##' summary(spider6)
##' fit0 <- mdm(y2p(spider6[,1:6])~1,data=spider6)
##' fit1 <- mdm(y2p(spider6[,1:6])~Water+Herbs,data=spider6)
##' fit2 <- mdm(y2p(spider6[,1:6])~Site,data=spider6,alpha=TRUE)
##' anova(fit0,fit1,fit2)
##'
"spider6"



