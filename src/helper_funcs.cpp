
#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
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

  

  mat weights;
  //compute weights in different ways
  if(as<string>(opts["weight_type"]) == "base") {
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    weights = weight_func(Xc, theta);
    
  } else if(as<string>(opts["weight_type"])=="subgroup") {
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    // subgroup indicators
    vec z = as<vec>(opts["z"]);
    vec uni_z = unique(z);
    List z_ind = as<List>(opts["z_ind"]);
    
    // initialize gradient as zero
    mat grad = zeros(Xc.n_cols, uni_z.size());

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
    return grad - Xt;
    
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

    return grad - Xt;

  } else if(as<string>(opts["weight_type"])=="hte"){

    weightPtr2 weight_func = *as<wptr2>(opts["weight_func"]);
    // iterate over the columns of theta
    int m = theta.n_cols;

    // special case for linear weights
    if(as<bool>(opts["linear"])) {
      return (Xc.t() * Xc) * theta - Xt;
    } else {
      mat eta = Xc * theta;
      weights = weight_func(eta);
    }
  } else {
    throw runtime_error("weight_type must be one of 'base', 'subgroup'");
  }
  // Rcout << weights.n_rows << ", " << weights.n_cols << "\n";
  //combine to get gradient
  return Xc.t() * weights - Xt;
    
}



// [[Rcpp::export]]
gptr make_balancing_grad() {
  return gptr(new gradPtr(balancing_grad));
}


