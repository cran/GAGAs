#' Calculate F1 score for classification, the inputs must be characters, and each of these elements must be either 'FALSE' or 'TRUE'.
#'
#' @param predictions predictions
#' @param truelabels true labels
#'
#' @return F1 score
#' @export cal.F1Score
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
#' cat("\n F1 score:", cal.F1Score(as.character(Eb!=0),as.character(beta_true!=0)))


cal.F1Score <- function(predictions, truelabels) {

  lvls <- levels(factor(truelabels))
  if (length(lvls)!=2) {
    stop("cal.F1Score error: The Number of factors must be 2!")
  }
  if(lvls[1]!='FALSE' || lvls[2]!='TRUE'){
    stop("cal.F1Score error: The factor must be a word of FALSE or TRUE !")
  }

  pos_true <- which(truelabels=='TRUE')
  pos2_true <- which(truelabels=='FALSE')
  pos <- which(predictions=='TRUE')
  pos2 <- which(predictions=='FALSE')

  TP=length(intersect(pos,pos_true))
  FP=length(intersect(pos,pos2_true))
  FN=length(intersect(pos2,pos_true))
  return(2*TP/(2*TP+FP+FN))
}
