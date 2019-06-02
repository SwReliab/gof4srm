#include <Rcpp.h>
using namespace Rcpp;

//' Compute P(D->d-) for the generalized KS statistic
//'
//' This is an experimental function instead of compute_Pdminus.
//' This is not used yet.
//'
//' @param ctime A sequence represents time slots (bins)
//' @param count A sequence indicates the number of samples falls int a bin
//' @param dplus A value of d+
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A value of the probability
// // [[Rcpp::export]]

List compute_Pdminus_arma(NumericVector ctime, IntegerVector count,
                          double dminus, Function cdf, Function solve) {
  int n = sum(count);
  int jmax = n * (1 - dminus);

  // step 1
  NumericVector c(jmax+1);
  NumericVector cp = cdf(ctime);
  for (int j=0; j<=jmax; j++) {
    double p = dminus + (double) j / n;
    for (int i=0; i<ctime.length(); i++) {
      if (cp[i] >= p) {
        c[j] = 1 - cp[i];
        break;
      }
    }
  }

  // step 2
  NumericMatrix B(jmax+1,jmax+1);
  NumericVector one(jmax+1);
  for (int j=0; j<=jmax; j++) {
    one(j) = 1;
    double x = 1;
    for (int i=j; i<=jmax; i++) {
      B(i,j) = x;
      x *= c(j) * (i+1) / (i-j+1);
    }
  }
  NumericVector b = solve(Named("a")=B, Named("b")=one);

  // step 3
  double prob = 0.0;
  double binom = 1.0;
  for (int j=0; j<=jmax; j++) {
    prob += binom * std::pow(c(j), n-j) * b(j);
    binom *= (double) (n-j) / (j+1);
  }

  return List::create(Named("B")=B, Named("prob")=prob);
}

//' Compute P(D+>d+) for the generalized KS statistic
//'
//' This is an experimental function instead of compute_Pdplus.
//' This is not used yet.
//'
//' @param ctime A sequence represents time slots (bins)
//' @param count A sequence indicates the number of samples falls int a bin
//' @param dplus A value of d+
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A value of the probability
// // [[Rcpp::export]]

List compute_Pdplus_arma(NumericVector ctime, IntegerVector count,
                         double dplus, Function cdf, Function solve) {
  int n = sum(count);
  int jmax = n * (1 - dplus);

  // step 1
  // step 1
  NumericVector c(jmax+1);
  NumericVector cp = cdf(ctime);
  double cp0 = 0;
  for (int j=0; j<=jmax; j++) {
    double p = 1 - dplus - (double) j / n;
    for (int i=0; i<ctime.length(); i++) {
      if (cp[i] >= p) {
        c[j] = cp0;
        break;
      }
      cp0 = cp[i];
    }
  }

  // step 2
  NumericMatrix B(jmax+1,jmax+1);
  NumericVector one(jmax+1);
  for (int j=0; j<=jmax; j++) {
    one(j) = 1;
    double x = 1;
    for (int i=j; i<=jmax; i++) {
      B(i,j) = x;
      x *= c(j) * (i+1) / (i-j+1);
    }
  }
  NumericVector b = solve(Named("a")=B, Named("b")=one);

  // step 3
  double prob = 0.0;
  double binom = 1.0;
  for (int j=0; j<=jmax; j++) {
    prob += binom * std::pow(c(j), n-j) * b(j);
    binom *= (double) (n-j) / (j+1);
  }

  return List::create(Named("B")=B, Named("b")=b, Named("prob")=prob);
}
