#' GAGAs: A package for fiting a GLM with GAGA algorithm
#'
#' Fits linear, logistic and multinomial, poisson, and Cox regression models via Global Adaptive Generative Adjustment algorithm.
#'
#' @section Mypackage functions:
#' GAGAs
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @docType package
#' @name GAGAs
#' @useDynLib GAGAs, .registration=TRUE
NULL
#> NULL


#' Fit a GLM with the GAGA algorithm
#'
#' Fits linear, logistic and multinomial, poisson, and Cox regression models via the Global Adaptive Generative Adjustment algorithm.
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param y Response variable. Quantitative for \code{family="gaussian"}, or
#' \code{family="poisson"} (non-negative counts). For \code{family="binomial"} should be either a factor with two levels.
#' For \code{family="multinomial"} should be a one-hot matrix or a \code{nc>=2} level factor.
#' For \code{family="cox"} should be an n*2 matrix, one column should be named "time", indicating the survival time;
#' the other column must be named "status", and consists of 0 and 1, 0 indicates that the row of data is censored, 1 is opposite.
#' @param family Either a character string representing one of the built-in families,
#' \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"multinomial"} or \code{"cox"}.
#' @param alpha Hyperparameter. In general, alpha can be set to 1, 2 or 3.
#' for \code{family="gaussian"} and \code{family="cox"}, the suggested value for alpha is 2 or 3.
#' for \code{family="binomial"}, \code{family="poisson"} and \code{family="multinomial"}, the suggested value for alpha is 1 or 2.
#' but when the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
#' @param itrNum The number of iteration steps. In general, 20 steps are enough.
#' If the condition number of \code{X} is large, it is recommended to greatly increase the
#' number of iteration steps.
#' @param thresh Convergence threshold for beta Change, if \code{max(abs(beta-beta_old))<threshold}, return.
#' @param QR_flag It identifies whether to use QR decomposition to speed up the algorithm.
#' Currently only valid for linear models.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fdiag It identifies whether to use diag Approximation to speed up the algorithm.
#' @param frp Identifies if a method is preprocessed to reduce the number of parameters
#' @param subItrNum Maximum number of steps for subprocess iterations.
#'
#' @return Regression coefficients
#' @export GAGAs
#' @examples
#' # Gaussian
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
#' Eb = fit$beta
#' #Create testing data
#' X_t = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y_t=X_t%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' #Prediction
#' Ey = predict.GAGA(fit,newx=X_t)
#'
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n perr:", norm(Ey-y_t,type="2")/sqrt(sample_size))
#'
#' # binomial
#' set.seed(2022)
#' cat("\n")
#' cat("Test binomial GAGA\n")
#' p_size = 30
#' sample_size=600
#' test_size=1000
#' R1 = 1
#' R2 = 3
#' ratio = 0.5 #The ratio of zeroes in coefficients
#' #Set the true coefficients
#' zeroNum = round(ratio*p_size)
#' ind = sample(1:p_size,zeroNum)
#' beta_true = runif(p_size,R2*0.2,R2)
#' beta_true[ind] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' t = 1/(1+exp(-X%*%beta_true))
#' tmp = runif(sample_size,0,1)
#' y = rep(0,sample_size)
#' y[t>tmp] = 1
#' fit = GAGAs(X,y,family = "binomial", alpha = 1)
#' Eb = fit$beta
#' #Generate test samples
#' X_t = R1*matrix(rnorm(test_size * p_size), ncol = p_size)
#' X_t[1:test_size,1]=1
#' t = 1/(1+exp(-X_t%*%beta_true))
#' tmp = runif(test_size,0,1)
#' y_t = rep(0,test_size)
#' y_t[t>tmp] = 1
#' #Prediction
#' Ey = predict(fit,newx = X_t)
#' cat("\n--------------------")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n pacc:", cal.w.acc(as.character(Ey),as.character(y_t)))
#' cat("\n")
#' # multinomial
#' set.seed(2022)
#' cat("\n")
#' cat("Test multinomial GAGA\n")
#' p_size = 20
#' C = 3
#' classnames = c("C1","C2","C3","C4")
#' sample_size = 500
#' test_size = 1000
#' ratio = 0.5 #The ratio of zeroes in coefficients
#' Num = 10 # Total number of experiments
#' R1 = 1
#' R2 = 5
#' #Set the true coefficients
#' beta_true = matrix(rep(0,p_size*C),c(p_size,C))
#' zeroNum = round(ratio*p_size)
#' for(jj in 1:C){
#'   ind = sample(1:p_size,zeroNum)
#'   tmp = runif(p_size,0,R2)
#'   tmp[ind] = 0
#'   beta_true[,jj] = tmp
#' }
#' #Generate training samples
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' z = X%*%beta_true
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' tt = t(apply(t,1,cumsum))
#' tt = cbind(rep(0,sample_size),tt)
#' # y = matrix(rep(0,sample_size*(C+1)),c(sample_size,C+1))
#' y = rep(0,sample_size)
#' for(jj in 1:sample_size){
#'   tmp = runif(1,0,1)
#'   for(kk in 1:(C+1)){
#'     if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
#'       # y[jj,kk] = 1
#'       y[jj] = kk
#'       break
#'     }
#'   }
#' }
#' y = classnames[y]
#' fit = GAGAs(X, y,alpha=1,family = "multinomial")
#' Eb = fit$beta
#' #Prediction
#' #Generate test samples
#' X_t = R1*matrix(rnorm(test_size * p_size), ncol = p_size)
#' X_t[1:test_size,1]=1
#' z = X_t%*%beta_true
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' tt = t(apply(t,1,cumsum))
#' tt = cbind(rep(0,test_size),tt)
#' y_t = rep(0,test_size)
#' for(jj in 1:test_size){
#'   tmp = runif(1,0,1)
#'   for(kk in 1:(C+1)){
#'     if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
#'       y_t[jj] = kk
#'       break
#'     }
#'   }
#' }
#' y_t = classnames[y_t]
#' Ey = predict(fit,newx = X_t)
#' cat("\n--------------------")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n pacc:", cal.w.acc(as.character(Ey),as.character(y_t)))
#' cat("\n")
#'
#' # Poisson
#' set.seed(2022)
#' p_size = 30
#' sample_size=300
#' R1 = 1/sqrt(p_size)
#' R2 = 5
#' ratio = 0.5 #The ratio of zeroes in coefficients
#' # Set the true coefficients
#' zeroNum = round(ratio*p_size)
#' ind = sample(1:p_size,zeroNum)
#' beta_true = runif(p_size,0,R2)
#' beta_true[ind] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' y = rpois(sample_size,lambda = as.vector(exp(X%*%beta_true)))
#' y = as.vector(y)
#' # Estimate
#' fit = GAGAs(X,y,alpha = 2,family="poisson")
#' Eb = fit$beta
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#'
#' # cox
#' p_size = 50
#' sample_size = 500
#' test_size = 1000
#' R1 = 3
#' R2 = 1
#' ratio = 0.5 #The ratio of zeroes in coefficients
#' censoringRate = 0.25 #Proportion of censoring data in observation data
#' # Set the true coefficients
#' zeroNum = round(ratio*p_size)
#' ind = sample(1:p_size,zeroNum)
#' beta_true = runif(p_size,-R2,R2)
#' beta_true[ind] = 0
#' # Generate training samples
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' z = X%*%beta_true
#' u = runif(sample_size,0,1)
#' t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
#' cs = rep(0,sample_size)
#' csNum = round(censoringRate*sample_size)
#' ind = sample(1:sample_size,csNum)
#' cs[ind] = 1
#' t[ind] = runif(csNum,0,0.8)*t[ind]
#' y = cbind(t,1 - cs)
#' colnames(y) = c("time", "status")
#' #Estimation
#' fit = GAGAs(X,y,alpha=2,family="cox")
#' Eb = fit$beta
#'
#' #Generate testing samples
#' X_t = R1*matrix(rnorm(test_size * p_size), ncol = p_size)
#' z = X_t%*%beta_true
#' u = runif(test_size,0,1)
#' t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
#' cs = rep(0,test_size)
#' csNum = round(censoringRate*test_size)
#' ind = sample(1:test_size,csNum)
#' cs[ind] = 1
#' t[ind] = runif(csNum,0,0.8)*t[ind]
#' y_t = cbind(t,1 - cs)
#' colnames(y_t) = c("time", "status")
#' #Prediction
#' pred = predict(fit,newx=X_t)
#'
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n Cindex:", cal.cindex(pred,y_t))
GAGAs <- function(X,y,family=c("gaussian","binomial","poisson","multinomial","cox"),alpha=2,itrNum=100,thresh=1.e-3,QR_flag=FALSE,flag=TRUE,lamda_0=0.001,fdiag=TRUE,frp=TRUE,subItrNum = 20) {

  if(!is.character(family)){
    warning("Please check the input of family")
    family = "gaussian"
  }

  family=match.arg(family)
  beta=switch(family,
              "gaussian"=LM_GAGA(X,y,alpha=alpha,itrNum=itrNum,thresh=thresh,QR_flag=QR_flag,flag=flag,lamda_0=lamda_0,fix_sigma=FALSE,sigm2_0 = 1,fdiag=fdiag,frp=frp),
              "poisson"=poisson_GAGA(X,y,alpha=alpha,itrNum=itrNum,thresh=thresh,flag=flag,lamda_0=lamda_0,fdiag=fdiag,subItrNum),
              "binomial"=logistic_GAGA(X,y,alpha=alpha,itrNum=itrNum,thresh=thresh,flag=flag,lamda_0=lamda_0,fdiag=fdiag,subItrNum),
              "multinomial"=multinomial_GAGA(X,y,alpha=alpha,itrNum=itrNum,thresh=thresh,flag=flag,lamda_0=lamda_0,fdiag=fdiag,subItrNum),
              "cox"=cox_GAGA(X,y,alpha=alpha,itrNum=itrNum,thresh=thresh,flag=flag,lamda_0=lamda_0,fdiag=fdiag,subItrNum)
  )

  return(beta)
}
