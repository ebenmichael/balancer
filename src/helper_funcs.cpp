
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
      // grad += as<double>(opts["hyper"]) * theta;
        grad += as<double>(opts["hyper"])  * theta;
    }
  }
  
  //combine to get gradient
  return grad - Xt;
    
}


// [[Rcpp::export]]
gptr make_balancing_grad() {
  return gptr(new gradPtr(balancing_grad));
}




//' Generic balancing loss gradient
// [[Rcpp::export]]
mat multilevel_grad(mat theta, List opts) {


  
  // control data
  mat Xc = as<mat>(opts["Xc"]);
  mat Xt = as<mat>(opts["Xt"]);

  int n_groups = as<int>(opts["n_groups"]);
  int dim = as<int>(opts["dim"]);

  // group level predictors
  int grp_cov_dim = as<int>(opts["grp_cov_dim"]);

  // initialize gradient as zero
  mat grad = zeros(Xc.n_cols, n_groups + 1);

  if(as<string>(opts["weight_type"])=="subgroup") {
    
    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    
    // subgroup indicators
    vec z = as<vec>(opts["z"]);
    vec uni_z = unique(z);
    List z_ind = as<List>(opts["z_ind"]);


    mat Xz;
    // iterate over subgroups and compute gradient
    vec weights = vec(size(Xc)[0]);
    for(int i=0; i < uni_z.size(); i++) {
      int k = uni_z(i);


      // get indices for subgroups, passed in through opts
   
      uvec idxs = as<uvec>(z_ind[k]);
      uvec col(1);
      col.fill(i);
      // get gradient for subgroup
      Xz = Xc.rows(idxs);
      weights.rows(idxs) = weight_func(Xz, theta.col(0) + theta.col(i+1));
      grad.col(i+1) = Xz.t() * weight_func(Xz, theta.col(0) + theta.col(i+1));

      // global parameters
      // grad.col(0) += grad.col(i+1);
    }
    
    grad.col(0) = Xc.t() * weights;
    
  } else if(as<string>(opts["weight_type"])=="multimatch") {

    weightPtr weight_func = *as<wptr>(opts["weight_func"]);
    
    // iterate over subgroups and compute gradient
    mat weights = mat(size(Xc)[0], n_groups);
    
    for(int i=0; i < n_groups; i++) {

      weights.col(i) = weight_func(Xc, theta.col(0) + theta.col(i+1));

      grad.col(i+1) = Xc.t() * weights.col(i);

      // global parameters
      // grad.col(0) += grad.col(i+1);
    }
    
    grad.col(0) = Xc.t() * sum(weights, 1);

    
  }

  if(as<bool>(opts["ridge"])) {
    // multi level structure with ridge penalties (Normal hyper priors)
    // global intercept
    vec global_int = as<vec>(opts["global_int"]);

    // group-level intercepts
    vec group_int = vec(n_groups);
    group_int.fill(as<double>(opts["group_int"]));

    // global parameters
    mat global_param = mat(dim + grp_cov_dim, 1);
    global_param.fill(as<double>(opts["global_param"]));


    // group parameters
    mat group_param = mat(dim + grp_cov_dim, n_groups);
    group_param.fill(as<double>(opts["group_param"]));
    mat hyper = join_cols(join_cols(global_int, group_int).t(), join_rows(global_param, group_param));

    grad += hyper % theta;
  }


  //combine to get gradient
  grad -= Xt;
  
  // zero out group parameters for group level predictors
  if(grp_cov_dim > 0) {
    grad.submat(dim+1, 1, dim + grp_cov_dim, n_groups) = zeros(grp_cov_dim, n_groups); 
  }

  return grad;
                 
    
  
    
}



// [[Rcpp::export]]
gptr make_multilevel_grad() {
  return gptr(new gradPtr(multilevel_grad));
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


