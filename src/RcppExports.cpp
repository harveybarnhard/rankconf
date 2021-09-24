// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rejupdate
int rejupdate(LogicalMatrix& rejmat, NumericMatrix const& diffmat, double const c);
RcppExport SEXP _rankconf_rejupdate(SEXP rejmatSEXP, SEXP diffmatSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix& >::type rejmat(rejmatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix const& >::type diffmat(diffmatSEXP);
    Rcpp::traits::input_parameter< double const >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(rejupdate(rejmat, diffmat, c));
    return rcpp_result_gen;
END_RCPP
}
// sigupdate
int sigupdate(LogicalMatrix const& rejmat, NumericMatrix& sigmat, NumericMatrix const& diffmat, double const c, int& numind, int const k);
RcppExport SEXP _rankconf_sigupdate(SEXP rejmatSEXP, SEXP sigmatSEXP, SEXP diffmatSEXP, SEXP cSEXP, SEXP numindSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix const& >::type rejmat(rejmatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type sigmat(sigmatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix const& >::type diffmat(diffmatSEXP);
    Rcpp::traits::input_parameter< double const >::type c(cSEXP);
    Rcpp::traits::input_parameter< int& >::type numind(numindSEXP);
    Rcpp::traits::input_parameter< int const >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(sigupdate(rejmat, sigmat, diffmat, c, numind, k));
    return rcpp_result_gen;
END_RCPP
}
// kmax
double kmax(NumericVector& x, int const k);
RcppExport SEXP _rankconf_kmax(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int const >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kmax(x, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rankconf_rejupdate", (DL_FUNC) &_rankconf_rejupdate, 3},
    {"_rankconf_sigupdate", (DL_FUNC) &_rankconf_sigupdate, 6},
    {"_rankconf_kmax", (DL_FUNC) &_rankconf_kmax, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rankconf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}