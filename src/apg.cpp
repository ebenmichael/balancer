
#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;


// vec tmp(vec x) {
//   return x;
// }


//' Accelerated Proximal Gradient method
//' @export
// [[Rcpp::export]]
vec apg2(Function f,
                   Function grad_f,
                   Function prox_h,
                   int dim, int max_it, double eps, double beta) {

  // initialize at zero
  vec x(dim);
  vec y = vec(x);
  vec oldx;
  vec gtx;
  vec grad;
  double fx;
  double fgtx;
  double improve;
  double diff;
  double t = 1;
  bool backcond;
  int j;
  for(int i = 1; i <= max_it; i++) {
    
    oldx = vec(x);

    fx = as<double>(f(y));
    grad = as<vec>(grad_f(x));
    gtx = 1 / t * (y - as<vec>(prox_h(y - t * grad, t)));    

    // backtracking line search
    backcond = as<double>(f(y - gtx)) <= fx - t * dot(grad, gtx) + t / 2 * dot(gtx, gtx);

    while(!backcond) {
      t = beta * t;

      gtx = 1 / t * (y -
                     as<vec>(prox_h(y -
                                    t * grad, t)));
      fgtx = as<double>(f(y - t * gtx));
      improve = fx - t * dot(grad , gtx) +
        t / 2 * dot(gtx, gtx);
      backcond =  fgtx <= improve;

    }

    x = y - t * gtx;
    y = x + (i - 1) / (i + 2) * (x - oldx);

    diff = norm(gtx, 2);
    Rcout << diff << "\n";
    Rcout << (diff <= eps) << "\n";
    j = i;
    if(diff == 0) {
      diff = 2 * eps;
    }
    if(diff <= eps) {
      break;
    }
    
  }
  Rcout << eps << "\n";
  Rcout << j << "\n";
  return(x);
}


// vec tmp(XPtr<funcPtr> f, const vec& x) {
//   funcPtr fun = *f;
//   return(fun(x));
// }


//' Accelerated proximal gradient method
//'
//' @param loss_ptr Pointer to loss function
//' @param grad_ptr Pointer to gradient function
//' @param prox_ptr Pointer to prox function
//' @param loss_opts List of options for loss (input data, tuning params, etc.)
//' @param dim Dimension
//' @param max_it Maximum number of iterations
//' @param eps Convergence tolerance
//' @param beta Backtracking line search parameter
//'
//' @return Optimal value
//' @export
// [[Rcpp::export]]
vec apg(lptr loss_ptr,
        gptr grad_ptr,
        pptr prox_ptr,
        List loss_opts,
        int dim, int max_it,
        double eps, double beta) {
  // grab the functions from pointers
  lossPtr f = *loss_ptr;
  gradPtr grad_f = *grad_ptr;
  proxPtr prox_h = *prox_ptr;
  
  // initialize at zero
  vec x(dim, fill::zeros);
  // accelerated
  vec y = vec(x);
  vec oldx;
  vec gtx;
  vec grad;
  double fx;
  double fgtx;
  double improve;
  double diff;
  double t = 1;
  bool backcond;
  for(int i = 1; i <= max_it; i++) {
    
    oldx = vec(x);
    fx = f(y, loss_opts);
    grad = grad_f(x, loss_opts);
    gtx = 1 / t * (y - prox_h(y - t * grad, t));
    

    // backtracking line search
    backcond = f(y - gtx, loss_opts) <= fx - t * dot(grad, gtx) + t / 2 * dot(gtx, gtx);

    while(!backcond) {
      t = beta * t;

      gtx = 1 / t * (y -
                     prox_h(y -
                            t * grad, t));
      fgtx = f(y - t * gtx, loss_opts);
      improve = fx - t * dot(grad , gtx) +
        t / 2 * dot(gtx, gtx);
      backcond =  fgtx <= improve;

    }

    x = y - t * gtx;
    y = x + (i - 1) / (i + 2) * (x - oldx);

    diff = norm(gtx, 2);

    if(diff <= eps) {
      break;
    }
  }

  return(x);
}
