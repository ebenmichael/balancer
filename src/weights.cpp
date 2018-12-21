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


//' Linear weights
// [[Rcpp::export]]
mat lin_weights_ipw(mat Xc, mat theta, mat q) {
  return q % (Xc * theta -1);
}


// [[Rcpp::export]]
wptr make_lin_weights() {
  return wptr(new weightPtr(lin_weights));
}

// [[Rcpp::export]]
wptr2 make_lin_weights2() {
  return wptr2(new weightPtr2(lin_weights2));
}

// [[Rcpp::export]]
wptripw make_lin_weights_ipw() {
  return wptripw(new weightPtrIPW(lin_weights_ipw));
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


//' normalized logit weights, numerically stable
// [[Rcpp::export]]
mat softmax_weights_ipw(mat Xc, mat theta, mat q) {
  mat eta = Xc * theta;
  double m = arma::max(arma::max(eta));

  return q % arma::exp(eta-m) / accu(q % exp(eta-m));
}


// [[Rcpp::export]]
wptr make_softmax_weights() {
  return wptr(new weightPtr(softmax_weights));
}


// [[Rcpp::export]]
wptr2 make_softmax_weights2() {
  return wptr2(new weightPtr2(softmax_weights2));
}


// [[Rcpp::export]]
wptripw make_softmax_weights_ipw() {
  return wptripw(new weightPtrIPW(softmax_weights_ipw));
}




//' un-normalized logit weights
// [[Rcpp::export]]
mat exp_weights(mat Xc, mat theta) {
  mat eta = Xc * theta;
  return arma::exp(eta-1);
}



//' un-normalized logit weights
// [[Rcpp::export]]
mat exp_weights2(mat eta) {
  return exp(eta - 1);
}


// [[Rcpp::export]]
wptr make_exp_weights() {
  return wptr(new weightPtr(exp_weights));
}


// [[Rcpp::export]]
wptr2 make_exp_weights2() {
  return wptr2(new weightPtr2(exp_weights2));
}




//' un-normalized logit weights
// [[Rcpp::export]]
mat exp_weights_ipw(mat Xc, mat theta, mat q) {
  mat eta = Xc * theta;

  return q % arma::exp(eta-1);
}



// [[Rcpp::export]]
wptripw make_exp_weights_ipw() {
  return wptripw(new weightPtrIPW(exp_weights_ipw));
}


//' Linear weights
// [[Rcpp::export]]
mat pos_lin_weights(mat Xc, mat theta) {
  mat eta = Xc * theta;
  return eta % (eta > 0);
}


//' Linear weights
// [[Rcpp::export]]
mat pos_lin_weights2(mat eta) {
  return eta % (eta > 0);
}


// [[Rcpp::export]]
wptr make_pos_lin_weights() {
  return wptr(new weightPtr(pos_lin_weights));
}

// [[Rcpp::export]]
wptr2 make_pos_lin_weights2() {
  return wptr2(new weightPtr2(pos_lin_weights2));
}



//' Linear weights
// [[Rcpp::export]]
mat pos_lin_weights_ipw(mat Xc, mat theta, mat q) {
  mat eta = Xc * theta;
  return q % (eta-1) % ((eta-1) > 0);
}


// [[Rcpp::export]]
wptripw make_pos_lin_weights_ipw() {
  return wptripw(new weightPtrIPW(pos_lin_weights_ipw));
}
