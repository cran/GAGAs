//[[Rcpp::depends(RcppEigen)]]

#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "helper_function.h"

using namespace Rcpp;
using namespace std;
using namespace Eigen;

double Func_u_multinomial(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {
  double eps = 2.2204e-16;
  int C = u.cols();
  Eigen::MatrixXd Xu = X * u;
  double v1 = -(y.leftCols(C).array() * Xu.array()).sum();
  double v2 = (Xu.array().exp().rowwise().sum() + 1).log().sum();
  double v3 = 0.5* (b.array()*(u.array().pow(2))).sum();
  return  v1 + v2 + v3;
}

double Func_lambda_multinomial(double lambda, std::vector<Eigen::MatrixXd const*> const &Plist) {
  if (Plist.size() != 5) {
    //std::cerr << "Func_lambda_logistic need 5 input parameters!" << std::endl;
    return 10000;
  }
  Eigen::MatrixXd const*u = Plist[0];
  Eigen::MatrixXd const*X = Plist[1];
  Eigen::MatrixXd const*y = Plist[2];
  Eigen::MatrixXd const*b = Plist[3];
  Eigen::MatrixXd const*d = Plist[4];

  return Func_u_multinomial((*u).array() + lambda*(*d).array(), *X, *y, *b);
}

Eigen::MatrixXd Dfu_multinomial(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                                Eigen::MatrixXd const &y, Eigen::MatrixXd const &b) {

  Eigen::MatrixXd Xu = X * u;
  int C = u.cols();
  Eigen::MatrixXd tmp1, tmp2;
  tmp1 = Xu.array().exp();
  tmp2 = 1 / (tmp1.array().rowwise().sum() + 1);
  return X.transpose() * (-y.leftCols(C) + (tmp1.array().colwise()*tmp2.col(0).array()).matrix()) + (u.array()*b.array()).matrix();
}

Eigen::MatrixXd getDDfu_multinomial(Eigen::MatrixXd const &u, Eigen::MatrixXd const &X,
                                    Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, bool fdiag) {

  Eigen::MatrixXd Xu = X * u;
  int C = u.cols();
  int P = u.rows();
  Eigen::MatrixXd Xt = X.transpose();
  Eigen::MatrixXd invgg = Eigen::MatrixXd::Zero(P*C,P*C);
  Eigen::MatrixXd gg = Eigen::MatrixXd::Zero(P*C,1);
  Eigen::MatrixXd tmp2 = 1 + Xu.array().exp().rowwise().sum();
  Eigen::MatrixXd tmp22 = tmp2.array().pow(2);
  Eigen::MatrixXd tmp1, tmp3, tmp4, tmp5, tmp6;
  for (int ii = 0; ii < C; ii++) {
    int indexii = ii*P;
    for (int jj = 0; jj < C; jj++) {
      int indexjj = jj*P;
      if (ii == jj) {
        tmp1 = Xu.col(jj).array().exp();
        tmp4 = tmp1.array() * (tmp2.array() - tmp1.array()) / tmp22.array();
        if (fdiag == false) {
          tmp5 = Xt * (X.array().colwise() * tmp4.array().col(0)).matrix();
          tmp5.diagonal() += b.col(jj);
          invgg.block(indexjj, indexii, P, P) = tmp5;
        }
        else {
          tmp6 = (X.array().pow(2).colwise() * tmp4.array().col(0)).colwise().sum();
          gg.block(indexjj, 0, P, 1) = tmp6.transpose() + b.col(jj);
        }
      }
      else {
        if (fdiag == false) {
          tmp1 = (Xu.col(jj) + Xu.col(ii)).array().exp();
          tmp3 = Xt * (X.array().colwise()*(tmp1.array() / tmp22.array()).col(0)).matrix();
          invgg.block(indexjj, indexii, P, P) = tmp3;
          invgg.block(indexii, indexjj, P, P) = tmp3;
        }
      }
    }//for(jj in 0:C-1)
  }//for(ii in 0:C-1)
  if (fdiag == true) {
    gg = 1 / (gg.array());
    return Eigen::MatrixXd(gg.asDiagonal());
  }
  else {
    return invgg.inverse();
  }
}

Eigen::MatrixXd getEb_multinomial(Eigen::MatrixXd const &X, Eigen::MatrixXd const &y, Eigen::MatrixXd const &b, Eigen::MatrixXd const &beta, Eigen::MatrixXd &D0, int maxItr, bool fdiag) {
  int N = X.rows();
  int P = X.cols();
  int C = beta.cols();
  Eigen::MatrixXd D = D0;
  Eigen::MatrixXd u = beta;
  Eigen::MatrixXd g = Dfu_multinomial(u, X, y, b);
  Eigen::MatrixXd d;
  std::vector<Eigen::MatrixXd const*> Plist(5);
  Plist[0] = &u;
  Plist[1] = &X;
  Plist[2] = &y;
  Plist[3] = &b;
  Plist[4] = &d;

  for (int ii = 0; ii < maxItr; ii++) {
    g.resize(P*C, 1);
    d = -D * g;
    d.resize(P, C);
    double LL = myfmin(0, 2, Func_lambda_multinomial, 20, 1e-19, Plist);

    u += LL*d;
    D = getDDfu_multinomial(u, X, y, b, fdiag);

    g = Dfu_multinomial(u, X, y, b);
    if (sqrt(g.array().pow(2).sum()) / P < 1e-16) {
      break;
    }
  }
  D0 = D;
  return u;
}

//' Fit a multinomial model via the GAGA algorithm using cpp
//'
//' Fit a multinomial model the Global Adaptive Generative Adjustment algorithm
//' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
//' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
//' @param y a One-hot response matrix or a \code{nc>=2} level factor
//' @param s_alpha Hyperparameter. The suggested value for alpha is 1 or 2.
//' When the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
//' @param s_itrNum The number of iteration steps. In general, 20 steps are enough.
//' If the condition number of \code{X} is large, it is recommended to greatly increase the
//' number of iteration steps.
//' @param s_thresh Convergence threshold for beta Change, if \code{max(abs(beta-beta_old))<threshold}, return.
//' @param s_flag It identifies whether to make model selection. The default is \code{TRUE}.
//' @param s_lamda_0 The initial value of the regularization parameter for ridge regression.
//' The running result of the algorithm is not sensitive to this value.
//' @param s_fdiag It identifies whether to use diag Approximation to speed up the algorithm.
//'
//' @return Coefficient matrix with K-1 columns beta_1,...,beta_{K-1} where K is the class number.
//' For k=1,..,K-1, the probability
//' \deqn{Pr(G=k|x)=exp(x^T beta_k) /(1+sum_{k=1}^{K-1}exp(x^T beta_k))}.
//' For k=K, the probability \deqn{Pr(G=K|x)=1/(1+sum_{k=1}^{K-1}exp(x^T beta_k))}.
//'
// [[Rcpp::export]]
Rcpp::List cpp_multinomial_gaga(Eigen::MatrixXd X, Eigen::MatrixXd y, SEXP s_alpha, SEXP s_itrNum, SEXP s_thresh,
                                     SEXP s_flag, SEXP s_lamda_0, SEXP s_fdiag) {

  double alpha = Rcpp::as<double>(s_alpha);
  int itrNum = Rcpp::as<int>(s_itrNum);
  double thresh = Rcpp::as<double>(s_thresh);
  bool flag = Rcpp::as<bool>(s_flag);
  double lamda_0 = Rcpp::as<double>(s_lamda_0);
  bool fdiag = Rcpp::as<bool>(s_fdiag);

  bool exitflag = false;
  double eps = 1.e-19;
  int N = X.rows();
  int P = X.cols();

  int K = y.cols();
  int C = K - 1;

  Eigen::MatrixXd b, b_old, db, beta, beta_old, cov_beta, E_pow_beta, cov0;
  b = Eigen::MatrixXd::Ones(P, C).array()*lamda_0;
  b_old = b;

  int index = 1;
  for (index = 1; index <= itrNum; index++) {
    if (index == itrNum || exitflag) {
      db = b - b_old;
      b = b / alpha;
    }

    if (index == 1) {
      beta = Eigen::MatrixXd::Zero(P, C);
      cov_beta = getDDfu_multinomial(beta, X, y, b, fdiag);
    }
    int maxItr = 20;
    beta = getEb_multinomial(X, y, b, beta, cov_beta, maxItr, fdiag);
    Eigen::MatrixXd beta2 = beta.array().pow(2);
    beta2.resize(P*C,1);
    E_pow_beta = cov_beta.diagonal().array() + beta2.array();
    b = alpha / E_pow_beta.array();
    b.resize(P, C);

    if (flag && (index == itrNum || exitflag)) {
      cov0 = getDDfu_multinomial(beta, X, y, Eigen::MatrixXd::Zero(P,C), fdiag);
      Eigen::MatrixXd diagcov0 = cov0.diagonal();
      for (int k = 0; k < diagcov0.size(); k++) {
        if (E_pow_beta(k) < diagcov0(k)) beta(k) = 0;
      }
      break;
    }
    else {
      b_old = b;
    }

    if (index == 1) {
      beta_old = beta;
    }
    else {
      if ((beta - beta_old).array().abs().maxCoeff() < thresh)exitflag = true;
      beta_old = beta;
    }
  }
  return Rcpp::List::create(Rcpp::Named("itrNum") = index,
                            Rcpp::Named("beta") = beta);

}
