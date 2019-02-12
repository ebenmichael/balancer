

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <numeric>
#include "balancer_types.h"

using namespace arma;
using namespace Rcpp;


//' Accelerated proximal gradient method
//'
//' @param grad_ptr Pointer to gradient function
//' @param prox_ptr Pointer to prox function
//' @param loss_opts List of options for loss (input data, tuning params, etc.)
//' @param prox_opts List of options for prox (regularization parameter)
//' @param x Initial value
//' @param max_it Maximum number of iterations
//' @param eps Convergence tolerance
//' @param beta Backtracking line search parameter
//' @param verbose How much information to print
//'
//' @return Optimal value
//' @export
// [[Rcpp::export]]
mat apg(gptr grad_ptr,
        pptr prox_ptr,
        List loss_opts,
        List prox_opts,
        mat x,
        int max_it,
        double eps, double alpha,
        double beta, bool accel, bool verbose) {

  // don't do anything if the initial x is infinite
  if(x.has_inf()) {
    return x;
  }
  
  
  int dim1 = x.n_rows;
  int dim2 = x.n_cols;
  // grab the functions from pointers
  gradPtr grad_f = *grad_ptr;
  proxPtr prox_h = *prox_ptr;
  

  mat y = x;
  // accelerated
  double theta = 1;
  mat oldx;
  mat oldy;
  mat gtx;
  mat grad;
  mat oldg;
  double fx;
  double fgtx;
  double improve;
  double diff;
  double t = 1.0;
  double oldt;
  double t_hat;
  bool backcond;
  int j;

  // if fail more than once, stop
  bool fail = 0;
  
  // step size initialization
  grad = grad_f(y, loss_opts);
  // Rcout << grad << "\n";
  t = 1 / sqrt(accu(pow(grad,2)));
  mat x_hat = x - t * grad;
  mat g_hat = grad_f(x_hat, loss_opts);
  double num = accu((x - x_hat) % (grad - g_hat));
  double denom = accu(pow(grad - g_hat,2));
  t = fabs(num / denom);
  // Rcout << num << ", " << denom << ", " << t << "\n";
  for(int i = 1; i <= max_it; i++) {
    // Rcout << i << "\n";
    // Rcout << t << "\n";
    oldx = mat(x);
    oldy = mat(y);

    x = prox_h(y - t * grad, t, prox_opts);

    // Rcout << accu(pow(y - x,2)) << "\n";
    // stopping criterion
    if(accu(pow(y - x,2)) < eps) {
      break;
    }

    if(accel) {
      theta = 2 / (1 + sqrt(1 + 4/(pow(theta,2))));
    } else {
      theta = 1;
    }

    // restart
    if(dot(grad, (x - oldx)) > 0) {
      x = oldx;
      y = x;
      theta = 1;
    }
    

    y = x + (1-theta) * (x - oldx);

    oldg = mat(grad);
    grad = grad_f(y, loss_opts);
    
    t_hat = 0.5 * accu(pow(y - oldy, 2)) /
      fabs(accu((y - oldy) % (oldg - grad)));
    double maxval = (t_hat > beta * t) ? t_hat : beta * t;
    t = (alpha * t < maxval) ? alpha * t : maxval;

    if((i % 100) == 0) {
      if(verbose) {
        Rcout << i << "\n";
        Rcout << t << "\n";
      }
      Rcpp::checkUserInterrupt();
    }

    // if x has any non-finite element, restart from zero
    // if it happens again, quit
    // TODO: figure out what is up with this and fix this hack
    // This probably has to do with the step size initialization
    if(!x.is_finite() & !fail) {
      Rcpp::warning("One or more parameters are non-finite. Restarting from zero");
      x.zeros();
      y = x;
      fail = 1;

      // restart step size
      grad = grad_f(y, loss_opts);
      t = 1 / sqrt(accu(pow(grad,2)));
      mat x_hat = x - t * grad;
      mat g_hat = grad_f(x_hat, loss_opts);
      double num = accu((x - x_hat) % (grad - g_hat));
      double denom = accu(pow(grad - g_hat,2));
      t = fabs(num / denom);
      
    } else if(!x.is_finite() & fail) {
      Rcpp::warning("One or more parameters are non-finite for the second time. Passing NaN.");
      break;
    }
    
  }

  return(x);
}



//' Accelerated proximal gradient method
//'
//' @param grad_ptr Pointer to gradient function
//' @param prox_ptr Pointer to prox function
//' @param loss_opts List of options for loss (input data, tuning params, etc.)
//' @param prox_opts List of options for prox (regularization parameter)
//' @param lams Vector of hyperparameters to use
//' @param x Initial value
//' @param max_it Maximum number of iterations
//' @param eps Convergence tolerance
//' @param beta Backtracking line search parameter
//' @param verbose How much information to print
//'
//' @return Optimal value
//' @export
// [[Rcpp::export]]
List apg_warmstart(gptr grad_ptr,
                  pptr prox_ptr,
                  List loss_opts,
                  List prox_opts,
                  vec lams,
                  mat x,
                  int max_it,
                  double eps, double alpha,
                  double beta, bool accel, bool verbose) {

  // keep different fitted values in a list
  int nlam = lams.n_elem;
  List output(nlam);
  // iterate over values of lambda
  for(int j = 0; j < nlam; j++) {
    prox_opts["lam"] = lams[j];
    // use previous optimizer as warm start
    x = apg(grad_ptr, prox_ptr, loss_opts, prox_opts,
            x, max_it, eps, alpha, beta, accel, verbose);

    output[j] = x;

    // replace large numbers with infinity
    if(arma::any(arma::abs(arma::vectorise(x)) > 1e23)) {
      x.fill(datum::inf);
    }
  }

  return(output);
}
