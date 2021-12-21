#include <Rcpp.h>
#include <queue>
#include <utility>
using namespace Rcpp;
using namespace std;

//' Update the symmetric rejection matrix in place based on what values of
//' the difference matrix are greater than the provided critical value.
//' @param rejmat The rejection matrix
//' @param diffmat The difference matrix
//' @param c A constant
// [[Rcpp::export]]
int rejupdate(LogicalMatrix &rejmat,
              NumericMatrix const &diffmat,
              double const c) {
  int nr = rejmat.nrow();
  int out = 0;

  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < i; j++) {
      if((diffmat(i,j) > c) & !rejmat(i,j)) {
        rejmat(i,j) = 1;
        rejmat(j,i) = 1;
        out += 2;
      }
    }
  }
  return out;
}

//' Update sampling indicator matrix, in place modifications
//' @param rejmat The rejection matrix
//' @param sigmat The sample indicator matrix
//' @param diffmat The difference matrix
//' @param c The constant for checking already-rejected values
//' @param numind Size of the sample
//' @param k The k in k-FWER
// [[Rcpp::export]]
int sigupdate(LogicalMatrix const &rejmat,
               NumericMatrix &sigmat,
               NumericMatrix const &diffmat,
               double const c,
               int &numind,
               int const k) {
  int nr = rejmat.nrow();
  if(k > 1) {
    for(int i = 0; i < nr; i++) {
      for(int j = 0; j < i; j++) {
        if((sigmat(i,j) < 0) & (diffmat(i,j) <= c)) {
          sigmat(i,j) = -sigmat(i,j);
          sigmat(j,i) = -sigmat(i,j);
          numind += 2;
        } else if((sigmat(i,j) > 0) & rejmat(i,j)){
          sigmat(i,j) = -sigmat(i,j);
          sigmat(j,i) = -sigmat(i,j);
          numind -= 2;
        }
      }
    }
  } else if(k==1) {
    for(int i = 0; i < nr; i++) {
      for(int j = 0; j < i; j++) {
        if((sigmat(i,j) > 0) & rejmat(i,j)) {
          sigmat(i,j) = -sigmat(i,j);
          sigmat(j,i) = -sigmat(i,j);
          numind -= 2;
        }
      }
    }
  }
  return numind;
}

bool compfun (double a, double b) { return (a > b); }

//' Returns the kth largest value by sorting in place
//' @param x Numeric vector.
//' @param k Specifies the kth largest value.
// [[Rcpp::export]]
double kmax(NumericVector &x, int const k) {
  std::nth_element(x.begin(), x.begin() + k - 1, x.end(), compfun);
  return x(k - 1);
}
