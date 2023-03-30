#' Fit a multinomial model via the GAGA algorithm
#'
#' Fit a multinomial model the Global Adaptive Generative Adjustment algorithm
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param y a One-hot response matrix or a \code{nc>=2} level factor
#' @param alpha Hyperparameter. The suggested value for alpha is 1 or 2.
#' When the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
#' @param itrNum The number of iteration steps. In general, 20 steps are enough.
#' If the condition number of \code{X} is large, it is recommended to greatly increase the
#' number of iteration steps.
#' @param thresh Convergence threshold for beta Change, if \code{max(abs(beta-beta_old))<threshold}, return.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fdiag It identifies whether to use diag Approximation to speed up the algorithm.
#' @param subItrNum Maximum number of steps for subprocess iterations.
#'
#' @return Coefficient matrix with K-1 columns beta_1,...,beta_{K-1} where K is the class number.
#' For k=1,..,K-1, the probability
#' \deqn{Pr(G=k|x)=exp(x^T beta_k) /(1+sum_{k=1}^{K-1}exp(x^T beta_k))}.
#' For k=K, the probability \deqn{Pr(G=K|x)=1/(1+sum_{k=1}^{K-1}exp(x^T beta_k))}.
#' @export multinomial_GAGA
#'
#' @examples
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

multinomial_GAGA = function(X,y,alpha=1,itrNum=50,thresh=1.e-3,flag=TRUE,lamda_0=0.001,fdiag=TRUE, subItrNum = 20){

  exitflag = FALSE
  eps = 1.e-19
  N = nrow(X)
  P = ncol(X)

  # refers to "lognet.R" of glmnet
  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(P),sep="")

  if(!is.matrix(y)){
    y = as.factor(y)
    ntab = table(y)
    classnames=names(ntab)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ,drop=FALSE]
  }else{
    if(!(ncol(y)>=2&&min(y)==0&&max(y)==1))stop("y should be a one hot matrix or a nc>=2 level factor")
    classnames = colnames(y)
    if(is.null(classnames))classnames=paste("V",seq(ncol(y)),sep="")
  }

  K = ncol(y)
  C = K-1

  fit = list()
  fit$classnames = classnames
  class(fit) = c("GAGA","multinomial")

  tmpfit = cpp_multinomial_gaga(X, as.matrix(y), alpha, itrNum, thresh, flag, lamda_0, fdiag, subItrNum)

  fit$beta = tmpfit$beta
  rownames(fit$beta) = vnames
  colnames(fit$beta) = classnames[1:C]
  fit$alpha = alpha
  fit$itrNum = tmpfit$itrNum
  fit$fdiag = fdiag
  return(fit)
}

