#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;
using namespace std;


// WEIGHT FUNCTIONS FOR GRADIENT

//' Linear weights
// [[Rcpp::export]]
mat lin_weights(mat Xc, mat theta) {
  return Xc * theta;
}


//' Linear weights
// [[Rcpp::export]]
mat lin_weights2(mat eta) {
  return eta;
}


// [[Rcpp::export]]
wptr make_lin_weights() {
  return wptr(new weightPtr(lin_weights));
}

// [[Rcpp::export]]
wptr2 make_lin_weights2() {
  return wptr2(new weightPtr2(lin_weights2));
}



//' normalized logit weights, numerically stable
// [[Rcpp::export]]
mat softmax_weights(mat Xc, mat theta) {
  mat eta = Xc * theta;
  double m = arma::max(arma::max(eta));

  return arma::exp(eta-m) / accu(exp(eta-m));
}



//' normalized logit weights, numerically stable
// [[Rcpp::export]]
mat softmax_weights2(mat eta) {
  double m = arma::max(arma::max(eta,0));
  eta -= m;
  eta = arma::exp(eta);
  return normalise(eta, 1, 0);
}


// [[Rcpp::export]]
wptr make_softmax_weights() {
  return wptr(new weightPtr(softmax_weights));
}


// [[Rcpp::export]]
wptr2 make_softmax_weights2() {
  return wptr2(new weightPtr2(softmax_weights2));
}

