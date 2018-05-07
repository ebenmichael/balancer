
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
         int dim, int max_it, double eps, double beta, bool accel) {

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

  // step size initialization
  grad = as<vec>(grad_f(y));
  t = 1 / sqrt(accu(pow(grad,2)));
  Rcout << t << "\n";
  vec x_hat = x - t * grad;
  vec g_hat = as<vec>(grad_f(x_hat));
  t = abs(accu( (x - x_hat) % (grad - g_hat)) / accu(pow(grad - g_hat,2)));
  Rcout << abs(accu( (x - x_hat) % (grad - g_hat))) << "\n";
  Rcout << accu(pow(grad - g_hat,2)) << "\n";
  Rcout << t << "\n---\n\n\n";

  if(accel) {
    for(int i = 1; i <= max_it; i++) {
      
      oldx = vec(x);
      
      fx = as<double>(f(y));
      grad = as<vec>(grad_f(y));
      gtx = 1 / t * (y - as<vec>(prox_h(y - t * grad, t)));    
      // Rcout << y(0) << "\n";
      // Rcout << fx << "\n";
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
      // Rcout << diff << "\n";
      // Rcout << (diff <= eps) << "\n";
      // Rcout << "---------\n\n";
      j = i;
      if(diff == 0) {
        diff = 2 * eps;
      }
      if(diff <= eps) {
        ;
      }
      if((i % 50) == 0) {
        Rcout << t << "\n---\n\n";
      }
    
    }
  } else {
    for(int i = 1; i <= max_it; i++) {
      
      oldx = vec(x);
      
      // fx = as<double>(f(x));
      grad = as<vec>(grad_f(x));
      

      // x = as<vec>(prox_h(oldx - t * grad, t));
      
      gtx = 1 / t * (x - as<vec>(prox_h(x - t * grad, t)));    
      // Rcout << x(0) << "\n";
      // Rcout << fx << "\n";
      // backtracking line search
      backcond = as<double>(f(x - gtx)) <= fx - t * dot(grad, gtx) + t / 2 * dot(gtx, gtx);
      
      while(!backcond) {
        t = beta * t;
        
        gtx = 1 / t * (x -
                       as<vec>(prox_h(x -
                                      t * grad, t)));
        fgtx = as<double>(f(x - t * gtx));
        improve = fx - t * dot(grad , gtx) +
        t / 2 * dot(gtx, gtx);
        backcond =  fgtx <= improve;
        
      }
      
      x = x - t * gtx;
      //y = x + (i - 1) / (i + 2) * (x - oldx);

      diff = norm(gtx, 2);
      // Rcout << diff << "\n";
      // Rcout << (diff <= eps) << "\n";
      // Rcout << "---------\n\n";
      j = i;
      if(diff == 0) {
        diff = 2 * eps;
      }
      if(diff <= eps) {
        //break;
        ;
      }
      if((i % 50) == 0) {
        Rcout << t << "\n---\n\n";
      }
      
      
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
//' @param grad_ptr Pointer to gradient function
//' @param prox_ptr Pointer to prox function
//' @param loss_opts List of options for loss (input data, tuning params, etc.)
//' @param prox_opts List of options for prox (regularization parameter)
//' @param x Initial value
//' @param max_it Maximum number of iterations
//' @param eps Convergence tolerance
//' @param beta Backtracking line search parameter
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
        double beta, bool accel) {

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

  // step size initialization
  grad = grad_f(y, loss_opts);
  t = 1 / sqrt(accu(pow(grad,2)));
  mat x_hat = x - t * grad;
  mat g_hat = grad_f(x_hat, loss_opts);
  double num = accu((x - x_hat) % (grad - g_hat));
  double denom = accu(pow(grad - g_hat,2));
  t = fabs(num / denom);

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
      // Rcout << i << "\n";
      // Rcout << t << "\n";
      Rcpp::checkUserInterrupt();
    }
  }

  return(x);
}


//' Accelerated Proximal Gradient method
//' @export
// [[Rcpp::export]]
vec apg3(Function grad_f,
         Function prox_h,
         int dim, int max_it, double eps,
         double beta, bool accel, double alpha) {

  // initialize at zero
  vec x= zeros<vec>(dim);
  vec y = vec(x);
  double theta = 1;
  vec oldx;
  vec oldy;
  vec gtx;
  vec grad;
  vec oldg;
  double fx;
  double fgtx;
  double improve;
  double diff;
  double t = 1;
  double oldt;
  double t_hat;
  bool backcond;
  int j;

  // step size initialization
  grad = as<vec>(grad_f(y));
  t = 1 / sqrt(accu(pow(grad,2)));
  // Rcout << t << "\n";
  vec x_hat = x - t * grad;
  vec g_hat = as<vec>(grad_f(x_hat));
  t = abs(accu( (x - x_hat) % (grad - g_hat)) / accu(pow(grad - g_hat,2)));
  // Rcout << abs(accu( (x - x_hat) % (grad - g_hat))) << "\n";
  // Rcout << accu(pow(grad - g_hat,2)) << "\n";
  // Rcout << t << "\n---\n\n\n";
  
  for(int i = 1; i <= max_it; i++) {

    // prox update
    oldx = vec(x);
    oldy = vec(y);

    grad = as<vec>(grad_f(y));

    x = as<vec>(prox_h(y - t * grad, t));

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
    
    // Rcout << "---\n";
    y = x + (1-theta) * (x - oldx);
    // Rcout << y.n_rows << ", "<< y.n_cols << "\n";
    oldg = vec(grad);
    grad = as<vec>(prox_h(y - t * grad, t));
    // Rcout << grad.n_rows << ", "<< grad.n_cols << "\n";
    
     t_hat = 0.5 * accu(pow(y - oldy, 2)) / abs(accu((y - oldy) % (oldg - grad)));
    //t_hat = 0.5 * abs(accu((y - oldy) % (oldg - grad))) / accu(pow(grad - oldg, 2));

    double maxval = (t_hat > beta * t) ? t_hat : beta * t;
    t = (alpha * t < maxval) ? alpha * t : maxval;
    // if((i % 50) == 0) {
    //   // Rcout << abs(accu((y - oldy) % (oldg - grad))) << "\n";
    //   // Rcout << accu(pow(grad - oldg, 2)) << "\n---\n";
    //   Rcout << t_hat << "\n";
    //   Rcout << maxval << "\n";
    //   Rcout << t << "\n---\n\n";
    // }
      
  }

  return(x);
}
