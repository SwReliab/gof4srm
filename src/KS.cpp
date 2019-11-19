#include <Rcpp.h>
using namespace Rcpp;

//' KS distance for IID sampes drawn from a continuous distribution.
//'
//' @param sample A numeric vector of IID samples from a continuous distribution
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A list with components;
//' \item{dplus}{A value of statistic D+.}
//' \item{dminus}{A value of statistic D-.}
//' \item{d}{A value of statistic max(D-,D+).}
//' @export
// [[Rcpp::export]]

List KSdistance_point(NumericVector sample, Function cdf) {
  int n = sample.length();
  sample.sort();
  NumericVector fv = cdf(sample);
  // t = 0
  double ft = 0.0;
  double et0 = 0.0;
  double et1 = 0.0;
  double dplus = 0.0;
  double dminus = 0.0;
  for (int i=0; i<n; i++) {
    ft = fv[i];
    et0 = et1;
    et1 = (i+1.0) / n;
    if (et1 - ft > dplus) {
      dplus = et1 - ft;
    }
    if (ft - et0 > dminus) {
      dminus = ft - et0;
    }
  }
  return List::create(Named("dplus")=dplus, Named("dminus")=dminus,
                      Named("d")=std::max(dplus,dminus));
}

//' KS distance for IID sampes drawn from a continuous distribution.
//'
//' @param msample A numeric matrix of IID samples from a continuous distribution (m series)
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A list with components;
//' \item{dplus}{A value of statistic D+.}
//' \item{dminus}{A value of statistic D-.}
//' \item{d}{A value of statistic max(D-,D+).}
//' @export
// [[Rcpp::export]]

List KSdistance_pointM(NumericMatrix msample, Function cdf) {
  int n = msample.nrow();
  int m = msample.ncol();
  NumericVector dplusv(m);
  NumericVector dminusv(m);
  NumericVector dv(m);
  for (int k=0; k<m; k++) {
    NumericVector sample = msample(_,k);
    sample.sort();
    NumericVector fv = cdf(sample);
    // t = 0
    double ft = 0.0;
    double et0 = 0.0;
    double et1 = 0.0;
    double dplus = 0.0;
    double dminus = 0.0;
    for (int i=0; i<n; i++) {
      ft = fv[i];
      et0 = et1;
      et1 = (i+1.0) / n;
      if (et1 - ft > dplus) {
        dplus = et1 - ft;
      }
      if (ft - et0 > dminus) {
        dminus = ft - et0;
      }
    }
    dplusv(k) = dplus;
    dminusv(k) = dminus;
    dv(k) = std::max(dplus,dminus);
  }
  return List::create(Named("dplus")=dplusv, Named("dminus")=dminusv,
                      Named("d")=dv);
}

//' KS distance for grouped samples drawn from a continuous distribution.
//'
//' @param ctime A sequence represents time slots (bins)
//' @param count A sequence indicates the number of samples falls int a bin
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A list with components;
//' \item{dplus}{A value of statistic D+.}
//' \item{dminus}{A value of statistic D-.}
//' \item{d}{A value of statistic max(D-,D+).}
// [[Rcpp::export]]

List KSdistance_group(NumericVector ctime, IntegerVector count, Function cdf) {
  int n = sum(count);
  int m = ctime.length();
  NumericVector fv = cdf(ctime);
  // t = 0
  int cumsum = 0;
  double ft = 0.0;
  double dplus = 0.0;
  double dminus = 0.0;
  for (int i=0; i<m; i++) {
    cumsum += count[i];
    ft = fv[i];
    double tmp = (double) cumsum / n;
    if (tmp - ft > dplus) {
      dplus = tmp - ft;
    }
    if (ft - tmp > dminus) {
      dminus = ft - tmp;
    }
  }
  return List::create(Named("dplus")=dplus, Named("dminus")=dminus,
                      Named("d")=std::max(dplus,dminus));
}

//' KS distance for grouped samples drawn from a continuous distribution.
//' (experimental: We do not know this is mathematically correct or not.)
//'
//' @param ctime A sequence represents time slots (bins)
//' @param count A sequence indicates the number of samples falls int a bin
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A list with components;
//' \item{dplus}{A value of statistic D+.}
//' \item{dminus}{A value of statistic D-.}
//' \item{d}{A value of statistic max(D-,D+).}
//' @export
// [[Rcpp::export]]

List KSdistance_group2(NumericVector ctime, IntegerVector count, Function cdf) {
  int n = sum(count);
  int m = ctime.length();
  NumericVector fv = cdf(ctime);
  // t = 0
  int cumsum = 0;
  double ft = 0.0;
  double et0 = 0.0;
  double et1 = 0.0;
  double dplus = 0.0;
  double dminus = 0.0;
  for (int i=0; i<m; i++) {
    cumsum += count[i];
    ft = fv[i];
    et0 = et1;
    et1 = (double) cumsum / n;
    if (et1 - ft > dplus) {
      dplus = et1 - ft;
    }
    if (ft - et0 > dminus) {
      dminus = ft - et0;
    }
  }
  return List::create(Named("dplus")=dplus, Named("dminus")=dminus,
                      Named("d")=std::max(dplus,dminus));
}

//' Compute P(C->d-) for the generalized KS statistic
//'
//' @param ctime A sequence represents time slots (bins)
//' @param count A sequence indicates the number of samples falls int a bin
//' @param dminus A value of d-
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A value of the probability
//' @export
// [[Rcpp::export]]

double compute_Pdminus(NumericVector ctime, IntegerVector count,
                       double dminus, Function cdf) {
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
  NumericVector b(jmax+1);
  b(0) = 1;
  for (int k=1; k<=jmax; k++) {
    double prob = 0.0;
    double binom = 1.0;
    for (int j=0; j<=k-1; j++) {
      prob += binom * std::pow(c(j), k-j) * b(j);
      binom *= (double) (k-j) / (j+1);
    }
    b(k) = 1 - prob;
  }

  // step 3
  double prob = 0.0;
  double binom = 1.0;
  for (int j=0; j<=jmax; j++) {
    prob += binom * std::pow(c(j), n-j) * b(j);
    binom *= (double) (n-j) / (j+1);
  }

  return prob;
}

//' Compute P(D+>d+) for the generalized KS statistic
//'
//' @param ctime A sequence represents time slots (bins)
//' @param count A sequence indicates the number of samples falls int a bin
//' @param dplus A value of d+
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A value of the probability
//' @export
// [[Rcpp::export]]

double compute_Pdplus(NumericVector ctime, IntegerVector count,
                       double dplus, Function cdf) {
  int n = sum(count);
  int jmax = n * (1 - dplus);

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
  NumericVector b(jmax+1);
  b(0) = 1;
  for (int k=1; k<=jmax; k++) {
    double prob = 0.0;
    double binom = 1.0;
    for (int j=0; j<=k-1; j++) {
      prob += binom * std::pow(c(j), k-j) * b(j);
      binom *= (double) (k-j) / (j+1);
    }
    b(k) = 1 - prob;
  }

  // step 3
  double prob = 0.0;
  double binom = 1.0;
  for (int j=0; j<=jmax; j++) {
    prob += binom * std::pow(c(j), n-j) * b(j);
    binom *= (double) (n-j) / (j+1);
  }

  return prob;
}

//' KS distance for grouped samples drawn from a continuous distribution.
//'
//' @param ctime A sequence represents time slots (bins)
//' @param size An integer for the number of total samples.
//' @param sample A set of sampels drawn from the multinomial distribution.
//' @param cdf A function of CDF. It is allowd to get a closure.
//' @return A list with components;
//' \item{dplus}{A vector of statistic D+.}
//' \item{dminus}{A vector of statistic D-.}
//' \item{d}{A vector of statistic max(D-,D+).}
//' @export
// [[Rcpp::export]]

List KSdistance_groupM(NumericVector ctime, int size, IntegerMatrix sample, Function cdf) {
  int n = size;
  int m = ctime.length();
  int N = sample.ncol();
  NumericVector fv = cdf(ctime);
  NumericVector dplusv(N);
  NumericVector dminusv(N);
  NumericVector dv(N);
  for (int k=0; k<N; k++) {
    int cumsum = 0;
    double ft = 0.0;
    double dplus = 0.0;
    double dminus = 0.0;
    for (int i=0; i<m; i++) {
      cumsum += sample(i,k);
      ft = fv[i];
      double tmp = (double) cumsum / n;
      if (tmp - ft > dplus) {
        dplus = tmp - ft;
      }
      if (ft - tmp > dminus) {
        dminus = ft - tmp;
      }
    }
    dplusv(k) = dplus;
    dminusv(k) = dminus;
    dv(k) = std::max(dplus, dminus);
  }
  return List::create(Named("dplus")=dplusv, Named("dminus")=dminusv,
                      Named("d")=dv);
}
