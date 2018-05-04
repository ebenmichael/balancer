#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace arma;
using namespace Rcpp;


typedef double (*lossPtr)(arma::vec x, List opts);
typedef mat (*gradPtr)(mat x, List opts);
typedef mat (*proxPtr)(mat x, double t, List opts);
typedef mat (*weightPtr)(mat X, mat theta);
typedef double (*wlossPtr)(mat X, vec theta);
typedef double (*regPtr)(vec theta);

typedef XPtr<lossPtr> lptr;
typedef XPtr<gradPtr> gptr;
typedef XPtr<proxPtr> pptr;
typedef XPtr<weightPtr> wptr;
typedef XPtr<wlossPtr> wlossptr;
typedef XPtr<regPtr> regptr;


typedef mat (*fullWeightPtr)(mat Xc, mat theta, wptr weight_func, List opts);
typedef XPtr<fullWeightPtr> fwptr;
