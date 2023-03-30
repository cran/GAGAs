#' Print a summary of GAGA object
#'
#' @param object Fitted "GAGA" object.
#' @param ... some other params
#'
#' @exportS3Method summary GAGA
#' @export summary.GAGA
#'
#' @examples
#' set.seed(2022)
#' p_size = 30
#' sample_size=300
#' R1 = 3
#' R2 = 2
#' ratio = 0.5 #The ratio of zeroes in coefficients
#' # Set the true coefficients
#' zeroNum = round(ratio*p_size)
#' ind = sample(1:p_size,zeroNum)
#' beta_true = runif(p_size,0,R2)
#' beta_true[ind] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y=X%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' # Estimation
#' fit = GAGAs(X,y,alpha = 3,family="gaussian")
#' summary(fit)
#'
summary.GAGA <- function(object, ...){
  if(inherits(object, "GAGA", which = FALSE)==FALSE)stop("predict.GAGA need a GAGA object as its input")
  family = class(object)[2]
  if(family == 'multinomial'){
    cat("GAGAs algorithm was run for", object$itrNum, "iterations with alpha =",object$alpha,", and the obtained model has ",sum(object$beta!=0)," non-zero coefficients in total.
\n")
    cnames = colnames(object$beta)
    for(ii in 1:length(cnames)){
      cat(cnames[ii],"has",sum(object$beta[,ii]!=0)," non-zero coefficients:\n" )
      print(object$beta[object$beta[,ii]!=0,ii])
    }
    cat('\n')

  }else{
    cat("GAGAs algorithm was run for", object$itrNum, "iterations with alpha =",object$alpha,", and the obtained model has ",sum(object$beta!=0)," non-zero coefficients:
\n")
    print(object$beta[object$beta!=0])
    if(utils::hasName(object,'sigma'))cat("\n The estimated noise standard deviation is ", object$sigma)
    cat('\n')
  }
}
