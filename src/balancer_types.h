#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace arma;
using namespace Rcpp;


typedef double (*lossPtr)(arma::vec x, List opts);
typedef arma::vec (*gradPtr)(arma::vec x, List opts);
typedef arma::vec (*proxPtr)(arma::vec x, double t);
typedef vec (*weightPtr)(mat X, vec theta);
typedef double (*wlossPtr)(mat X, vec theta);
typedef double (*regPtr)(vec theta);

typedef XPtr<lossPtr> lptr;
typedef XPtr<gradPtr> gptr;
typedef XPtr<proxPtr> pptr;
typedef XPtr<weightPtr> wptr;
typedef XPtr<wlossPtr> wlossptr;
typedef XPtr<regPtr> regptr;
