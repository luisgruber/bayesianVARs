#include <Rcpp.h>
#include <GIGrvg.h>
using namespace Rcpp;

//' @export
//[[Rcpp::export]]
NumericMatrix my_gig(int n, NumericVector lambda, NumericVector chi, NumericVector psi) {

  NumericVector mtmp = NumericVector::create(lambda.size(), chi.size(), psi.size());
  int m = Rcpp::max(mtmp);
  NumericVector lambda1 =  rep_len(lambda, m);
  NumericVector chi1 =  rep_len(chi, m);
  NumericVector psi1 =  rep_len(psi, m);

  NumericMatrix out(n,m);

  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");

  for(int i = 0; i<n; ++i) {
    for(int j = 0; j<m; ++j){

      out(i,j) = as<double>(fun(1, lambda1[j], chi1[j], psi1[j]));

    }
  }

  return out;

}
