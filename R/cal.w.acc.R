#' Calculate the weighted ACC of the classification, the inputs must be characters
#'
#' @param predictions predictions
#' @param truelabels true labels
#'
#' @return weighted ACC
#' @export cal.w.acc
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

cal.w.acc <- function(predictions, truelabels) {

  tmp = table(truelabels)
  # weights = tmp/sum(tmp)
  # print(weights)
  lvls <- levels(factor(truelabels))
  accs <- lapply(lvls, function(x) {
    idx <- which(truelabels == x)

    return(cal.acc(predictions[idx], truelabels[idx]))

  })
  #browser()
  acc <- mean(unlist(accs))

  return(acc)

}
