
#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
#include <string>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;
using namespace std;


// GENERIC GRADIENT FUNCTIONS


//' Generic balancing loss gradient
// [[Rcpp::export]]
mat balancing_grad(mat theta, List opts) {

  // control data
  mat Xc = as<mat>(opts["Xc"]);
  mat Xt = as<mat>(opts["Xt"]);
  // weights function
  // wptr weight_func = as<wptr>(opts["weight_func"]);

  // fullWeightPtr full_w_func = *as<fwptr>(opts["weight_type"]);

  // // get weights
  // mat weights = full_w_func(Xc, theta, weight_func, opts);

  

  mat grad;
  bool doridge = 1;  
  //compute weights in different ways
  if(as<string>(opts["weight_type"]) == "base") {
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    grad = Xc.t() * weight_func(Xc, theta);
    
  } else if(as<string>(opts["weight_type"])=="subgroup") {
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    // subgroup indicators
    vec z = as<vec>(opts["z"]);
    vec uni_z = unique(z);
    List z_ind = as<List>(opts["z_ind"]);
    
    // initialize gradient as zero
    grad = zeros(Xc.n_cols, uni_z.size());

    mat Xz;
    // iterate over subgroups and compute gradient
    for(int i=0; i < uni_z.size(); i++) {
      int k = uni_z(i);
      // Rcout << k<< "\n";
      // get indices for subgroups, passed in through opts
   
      uvec idxs = as<uvec>(z_ind[k]);
      uvec col(1);
      col.fill(i);
      // get gradient for subgroup
      Xz = Xc.rows(idxs);
      grad.col(i) = Xz.t() * weight_func(Xz, theta.col(i));

    }
    
  } else if(as<string>(opts["weight_type"])=="missing"){
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    // treatment assignment
    // here Xc is the matrix of covariates for units with outcomes observed
    vec trt = as<vec>(opts["trt"]);

    mat grad = zeros(Xc.n_cols, 2);
    
    // two sets of gradientss
    // ctrl -> trt gradients
    uvec idx_ctrl = as<uvec>(opts["idx_ctrl"]);
    mat Xz = Xc.rows(idx_ctrl);
    grad.col(0) = Xz.t() * weight_func(Xz, theta.col(0));

    // observed trt -> trt weights
    uvec idx_trt = as<uvec>(opts["idx_trt"]);
    Xz = Xc.rows(idx_trt);
    grad.col(1) = Xz.t() * weight_func(Xz, theta.col(1));

  } else if(as<string>(opts["weight_type"])=="hte"){

    weightPtr2 weight_func = *as<wptr2>(opts["weight_func"]);
    // iterate over the columns of theta
    int m = theta.n_cols;
    // special case for linear weights
    if(as<bool>(opts["linear"])) {
        grad = Xc.t() * Xc * theta;
    } else {
      mat eta = Xc * theta;
      grad = Xc.t() * weight_func(eta);
    }
  } else {
    throw runtime_error("weight_type must be one of 'base', 'subgroup'");
  }

  // include (generalized) ridge penalty
  if(as<bool>(opts["ridge"])) {
    if(as<bool>(opts["hasQ"])) {
      mat Q = as<mat>(opts["Q"]);
      grad += as<double>(opts["hyper"]) * theta * Q * Q.t();
    }
    else {
      grad += as<double>(opts["hyper"]) * theta;
    }
  }
  //combine to get gradient
  return grad - Xt;
    
}



// [[Rcpp::export]]
gptr make_balancing_grad() {
  return gptr(new gradPtr(balancing_grad));
}


//' Generic balancing loss gradient with kernels
// [[Rcpp::export]]
mat balancing_grad_kern(mat theta, List opts) {

  // control data
  mat Xc = as<mat>(opts["Xc"]);
  mat Xt = as<mat>(opts["Xt"]);

  // Only linear weights for now

  mat grad;

  mat Q = as<mat>(opts["Q"]);

  //compute weights in different ways
  if(as<string>(opts["weight_type"]) == "base") {
    // weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    grad = (Xc + Q) * theta;
    
  } else if(as<string>(opts["weight_type"])=="subgroup") {
    // weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    // subgroup indicators
    vec z = as<vec>(opts["z"]);
    vec uni_z = unique(z);
    List z_ind = as<List>(opts["z_ind"]);
    
    // initialize gradient as zero
    grad = zeros(Xc.n_cols, uni_z.size());

    mat Xz;
    // iterate over subgroups and compute gradient
    for(int i=0; i < uni_z.size(); i++) {
      int k = uni_z(i);
      // Rcout << k<< "\n";
      // get indices for subgroups, passed in through opts
   
      uvec idxs = as<uvec>(z_ind[k]);
      uvec col(1);
      col.fill(i);
      // get gradient for subgroup
      Xz = Xc.rows(idxs);
      grad.col(i) = Xz * (Xz + Q) * theta.col(i);

    }
    
  } else if(as<string>(opts["weight_type"])=="missing"){
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    // treatment assignment
    // here Xc is the matrix of covariates for units with outcomes observed
    vec trt = as<vec>(opts["trt"]);

    mat grad = zeros(Xc.n_cols, 2);
    
    // two sets of gradientss
    // ctrl -> trt gradients
    uvec idx_ctrl = as<uvec>(opts["idx_ctrl"]);
    mat Xz = Xc.rows(idx_ctrl);
    grad.col(0) = Xz * (Xz + Q) * theta.col(0);

    // observed trt -> trt weights
    uvec idx_trt = as<uvec>(opts["idx_trt"]);
    Xz = Xc.rows(idx_trt);
    grad.col(1) = Xz * (Xz + Q) * theta.col(1);

  } else if(as<string>(opts["weight_type"])=="hte"){

    weightPtr2 weight_func = *as<wptr2>(opts["weight_func"]);
    // iterate over the columns of theta
    int m = theta.n_cols;
    // special case for linear weights
    grad = Xc * (Xc + Q) * theta;
  } else {
    throw runtime_error("weight_type must be one of 'base', 'subgroup'");
  }


  //combine to get gradient
  return grad - Xt;
    
}



// [[Rcpp::export]]
gptr make_balancing_grad_kern() {
  return gptr(new gradPtr(balancing_grad_kern));
}


