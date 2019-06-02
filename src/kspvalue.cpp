#include <Rcpp.h>
using namespace Rcpp;

//' Compute p-value of KS.
//'
//' This is not in use.
//'
// [[Rcpp::export]]

double ks_pvalue(double d, int n, int imax = 100000, double epsi = 1.0e-12) {
  const double pi = 4.0 * std::atan(1.0);
  double x = std::sqrt(n) * d;
  double tmp = std::sqrt(2*pi) / x;
  double sum = 0.0;
  for (int i=1; i<=imax; i++) {
    double tmp = exp(-(2*i-1)*(2*i-1)*pi*pi/(8*x*x));
    sum += tmp;
    if (tmp < epsi) {
      break;
    }
  }
  return 1.0 - tmp * sum;
}
