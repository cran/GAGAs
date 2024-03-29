#' compute C index for a Cox model
#'
#' Computes Harrel's C index for predictions from a \code{"cox"} object.
#'
#' Computes the concordance index, taking into account censoring.
#' This file fully references the Cindex.R file in glmnet package.

#' @param pred Predictions from a \code{"cox"} object
#' @param y a survival response object - a matrix with two columns "time" and
#' "status"; see documentation for "glmnet" or see documentation for "GAGA"
#' @param weights optional observation weights
#' @return Harrel's C index
#'
#' @author Trevor Hastie  <hastie@@stanford.edu>
#' @seealso \code{cv.glmnet}
#' @references Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996)
#' \emph{Tutorial in biostatistics: multivariable prognostic models: issues in
#' developing models, evaluating assumptions and adequacy, and measuring and
#' reducing error}, Statistics in Medicine, 15, pages 361--387.
#' @examples
#'
#' set.seed(10101)
#' N = 1000
#' p = 30
#' nzc = p/3
#' x = matrix(rnorm(N * p), N, p)
#' beta = rnorm(nzc)
#' fx = x[, seq(nzc)] %*% beta/3
#' hx = exp(fx)
#' ty = rexp(N, hx)
#' tcens = rbinom(n = N, prob = 0.3, size = 1)  # censoring indicator
#' y = cbind(time = ty, status = 1 - tcens)  # y=Surv(ty,1-tcens) with library(survival)
#' fit = GAGAs(x, y, family = "cox")
#' pred = predict(fit, newx = x)
#' cat("\n Cindex:", cal.cindex(pred, y))
#'
#'
#' @export cal.cindex
cal.cindex=function(pred,y,weights=rep(1,nrow(y))){
  ###  This function links to the concordance function in the survival package
  if(!survival::is.Surv(y))y=survival::Surv(y[,"time"],y[,"status"])
  f=-pred
  if(missing(weights))
    survival::concordance(y~f)$concordance
  else
    survival::concordance(y~f,weights=weights)$concordance
}
