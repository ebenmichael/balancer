
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
