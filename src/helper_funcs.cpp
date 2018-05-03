
#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;


// GENERIC LOSS AND GRADIENT FUNCTIONS

//' Generic balancing loss function
// [[Rcpp::export]]
double balancing_loss(vec theta, List opts) {

  // control data
  mat Xc = as<mat>(opts["Xc"]);
  mat Xt = as<mat>(opts["Xt"]);
  
  // weighted loss portion
  wlossPtr wloss_func = *as<wlossptr>(opts["wloss"]);
  double wloss = wloss_func(Xc, theta);

  //combine with treated values
  return wloss - trace(Xt.t() * theta);
    
  
}


// [[Rcpp::export]]
lptr make_balancing_loss() {
  return lptr(new lossPtr(balancing_loss));
}



//' Generic balancing loss gradient
vec balancing_grad(vec theta, List opts) {

  // control data
  mat Xc = as<mat>(opts["Xc"]);
  mat Xt = as<mat>(opts["Xt"]);
  // weights function
  weightPtr weight_func = *as<wptr>(opts["weight_func"]);

  // get weights
  mat weights = weight_func(Xc, theta);

  //combine to get gradient
  return Xc.t() * weights - Xt;
    
}



// [[Rcpp::export]]
gptr make_balancing_grad() {
  return gptr(new gradPtr(balancing_grad));
}


// WEIGHTED PORITON OF LOSS FUNCTION

//' Loss function for linear link
//[Rcpp::export]]
double var_loss(mat X, vec theta) {
  return norm(X * theta, "fro");
}

//[Rcpp::export]]
wlossptr make_var_loss() {
  return wlossptr(new wlossPtr(var_loss));
}


// WEIGHT FUNCTIONS FOR GRADIENT

//' Linear weights
vec lin_weights(mat Xc, vec theta) {
  return Xc * theta;
}



//[Rcpp::export]]
wptr make_lin_weights() {
  return wptr(new weightPtr(lin_weights));
}


// PROX FUNCTIONS

// [[Rcpp::export]]
arma::vec no_prox(arma::vec theta, double t) {

  return theta;
     
}

// [[Rcpp::export]]
pptr make_no_prox() {
  return pptr(new proxPtr(no_prox));
}

// [[Rcpp::export]]
vec prox_l1(vec x, double lam) {
  return (x - lam) * (x > lam) + (x + lam) * (x < -lam);
}


// [[Rcpp::export]]
pptr make_prox_l1() {
  return pptr(new proxPtr(prox_l1));
}
