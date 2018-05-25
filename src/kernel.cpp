#include <RcppArmadillo.h>
#include <numeric>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;
using namespace std;



//' Polynomial kernel
//' @param x,y Input vectors
//' @param d Degree of polynomial
// [[Rcpp::export]]
double poly_kernel(mat x, mat y, double d) {
  return pow(dot(x, y), d);
}



//' Radial Basis Function kernel
//' @param x,y Input vectors
//' @param sig
// [[Rcpp::export]]
double rbf_kernel(mat x, mat y, double sig) {
  return exp(-norm(x - y, 2) / ( 2 * pow(sig, 2)));
}



//' Internal function to compute a gram matrix between X and Y
//' @param X, Y, matrices
//' @param kernel Kernel function
//' @param param Hyper parameter for the kernel
mat compute_gram_(mat X, mat Y, kptr kernel, double param) {
  double nx = X.n_rows;
  double ny = Y.n_rows;
  mat gram = zeros(nx, ny);
  kernelPtr kern_func = *kernel;
  for(int i = 0; i < nx; i++) {
    for(int j = 0; j < ny; j++) {
      gram(i, j) = kern_func(X.row(i), Y.row(j), param);
    }
  }
  return gram;
}



//' Compute a gram matrix between X and Y
//' for a given kernel
//' @param X, Y, matrices
//' @param kernel Name of kernel
//' @param param Hyper parameter for the kernel
// [[Rcpp::export]]
mat compute_gram(mat X, mat Y, std::string kernel, double param) {

  if(kernel == "rbf") {
    kptr kern = kptr(new kernelPtr(rbf_kernel));
    return compute_gram_(X, Y, kern, param);
  } else if(kernel == "poly") {
    kptr kern = kptr(new kernelPtr(poly_kernel));
    return compute_gram_(X, Y, kern, param);
  } else {
    throw invalid_argument("kernel must be one of 'rbf', 'poly'");
  }

  
}
