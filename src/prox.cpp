
#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

// PROX FUNCTIONS


mat no_prox(mat theta, double t, List opts) {

  return theta;
     
}

// [[Rcpp::export]]
pptr make_no_prox() {
  return pptr(new proxPtr(no_prox));
}

// L1 PROX
mat prox_l1(mat x, double lam, List opts) {
  lam = lam * as<double>(opts["lam"]);
  return (x - lam) % (x > lam) + (x + lam) % (x < -lam);
}


// [[Rcpp::export]]
pptr make_prox_l1() {
  return pptr(new proxPtr(prox_l1));
}


// GROUP L1 PROX
mat prox_l1_grp(mat x, double lam, List opts) {
  lam = lam * as<double>(opts["lam"]);
  double rownorm;
  for(int i=0; i < x.n_rows; i++) {
    rownorm = norm(x.row(i), 2);
    x.row(i) = (x.row(i) - x.row(i) / rownorm * lam) * (rownorm > lam);
  }
  return x;
}

// [[Rcpp::export]]
pptr make_prox_l1_grp() {
  return pptr(new proxPtr(prox_l1_grp));
}



// Nuclear norm PROX
mat prox_nuc(mat x, double lam, List opts) {
  mat U; vec s; mat V;
  // SVD then threshold singular values
  svd_econ(U, s, V, x);
  // Rcout << U.n_rows << ", " << U.n_cols << "\n";
  // Rcout << V.n_rows << ", " << V.n_cols << "\n";
  // Rcout << s.size() << "\n";
  int d = x.n_rows;
  int m = x.n_cols;

  // threshold singular values
  lam = lam * as<double>(opts["lam"]);
  s = (s - lam) % (s > lam) + (s + lam) % (s < -lam);
  // mat smat;
  // if(d >= m) {
  //   //smat = join_vert(diagmat(s), zeros<mat>(d - m, m));
  //   x = U.cols(0, m-1) * (diagmat(s) * V.t());
  // } else {
  //   // smat = join_horiz(diagmat(s), zeros<mat>(d, m - d));
    
  //   x = (U * diagmat(s)) * V.cols(0, d-1).t();
  // }
  // return x;

  return U * diagmat(s) * V.t();
}

// [[Rcpp::export]]
pptr make_prox_nuc() {
  return pptr(new proxPtr(prox_nuc));
}
