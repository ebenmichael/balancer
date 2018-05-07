
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
    
    // initialize weights as zero
    weights = zeros(Xc.n_rows, uni_z.size());
    // Rcout << weights.n_rows << ", " << weights.n_cols << "\n";
    // iterate over subgroups and set weights
    for(int i=0; i < uni_z.size(); i++) {
      int k = uni_z(i);
      // Rcout << k<< "\n";
      // get indices for subgroups
      uvec idxs = find(z == k);
      uvec col = find(uni_z == k);
      // Rcout << col << "\n---\n";
      // Rcout << idxs << "\n---\n";
      // get weights for subgroup
      mat tmp_w = weight_func(Xc.rows(idxs),theta.col(i));
      // Rcout << tmp_w.n_rows << ", " << tmp_w.n_cols << "\n";
      weights.submat(idxs, col) = tmp_w;
      // Rcout << weights.n_rows << ", " << weights.n_cols << "\n--\n";
    }
    // Rcout << weights.n_rows << ", " << weights.n_cols << "\n--\n";
  } else if(as<string>(opts["weight_type"])=="missing"){
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    // treatment assignment
    // here Xc is the matrix of covariates for units with outcomes observed
    vec trt = as<vec>(opts["trt"]);
    
    // two sets of weights
    // ctrl -> trt weights
    uvec idx_ctrl = find(trt == 0);
    vec weights1 = zeros(Xc.n_rows);
    weights1.elem(idx_ctrl) = weight_func(Xc.rows(idx_ctrl), theta.col(0));

    // observed trt -> trt weights
    uvec idx_trt = find(trt == 1);
    vec weights2 = zeros(Xc.n_rows);
    weights2.elem(idx_trt) = weight_func(Xc.rows(idx_trt), theta.col(1));

    // combine weights
    weights = join_horiz(weights1, weights2);
                                             
    
    

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


