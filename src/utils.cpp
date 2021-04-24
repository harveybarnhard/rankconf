#include <Rcpp.h>
using namespace Rcpp;

// Goal here is to update the reject matrix in place based on the difference
// matrix and the critical value. The output will be the number of new
// rejections.

//' Update the symmetric rejection matrix in place based on what values of
//' the difference matrix are greater than the provided critical value
//' @param rejmat The rejection matrix
//' @param diffmat The difference matrix
//' @param c A constant
//' @export
// [[Rcpp::export]]
int rejupdate(LogicalMatrix &rejmat,
              NumericMatrix const &diffmat,
              double const &c) {
  int nr = rejmat.nrow();
  int out = 0;
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < i; j++) {
      if(diffmat(i,j) > c & rejmat(i,j)==0) {
        rejmat(i,j) = 1;
        rejmat(j,i) = 1;
        out += 2;
      }
    }
  }
  return out;
}

//' Update sampling indicator matrix
//' @param rejmat The rejection matrix
//' @param indmat The sample indicator matrix
//' @param diffmat The difference matrix
//' @param c The constant for checking already-rejected values
//' @param k The k in k-FWER
//' @export
// [[Rcpp::export]]
int indupdate(LogicalMatrix const &rejmat,
              LogicalMatrix &indmat,
              NumericMatrix const &diffmat,
              double const c,
              int const k) {
  int nr = rejmat.nrow();
  int out = 0;
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < i; j++) {
      if((rejmat(i,j)==1 & diffmat(i,j) <= c & k  > 1) | (rejmat(i,j)==0)) {
        indmat(i,j) = 1;
        indmat(j,i) = 1;
      } else {
        indmat(i,j) = 0;
        indmat(j,i) = 0;
      }
    }
  }
  return out;
}
