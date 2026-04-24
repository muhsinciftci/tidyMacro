#include "fGetBands.h"
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

static double lin_interp(const std::vector<double>& sv, int B, double pct) {
    double raw = (pct / 100.0) * (B - 1);
    int lo = static_cast<int>(raw);
    double frac = raw - lo;
    if (lo + 1 >= B) return sv[lo];
    return sv[lo] * (1.0 - frac) + sv[lo + 1] * frac;
}

FGetBandsResult fGetBands_cpp(const arma::cube& bootirf, double prc) {
    int N = static_cast<int>(bootirf.n_rows);
    int H = static_cast<int>(bootirf.n_cols);
    int B = static_cast<int>(bootirf.n_slices);

    double up_pct  = 50.0 + prc / 2.0;
    double low_pct = 50.0 - prc / 2.0;

    arma::mat upper(N, H, arma::fill::none);
    arma::mat lower(N, H, arma::fill::none);
    arma::mat med(N, H, arma::fill::none);

    for (int i = 0; i < N; ++i) {
        for (int h = 0; h < H; ++h) {
            std::vector<double> sv(B);
            for (int b = 0; b < B; ++b) sv[b] = bootirf(i, h, b);
            std::sort(sv.begin(), sv.end());
            upper(i, h) = lin_interp(sv, B, up_pct);
            lower(i, h) = lin_interp(sv, B, low_pct);
            med(i, h)   = lin_interp(sv, B, 50.0);
        }
    }

    FGetBandsResult result;
    result.upper  = upper;
    result.lower  = lower;
    result.median = med;
    return result;
}

// Sort once, extract all five quantiles — avoids double sort when called for two bands
FGetBands2Result fGetBands2_cpp(const arma::cube& bootirf, double prc, double prc2) {
    int N = static_cast<int>(bootirf.n_rows);
    int H = static_cast<int>(bootirf.n_cols);
    int B = static_cast<int>(bootirf.n_slices);

    double up1  = 50.0 + prc  / 2.0;
    double lo1  = 50.0 - prc  / 2.0;
    double up2  = 50.0 + prc2 / 2.0;
    double lo2  = 50.0 - prc2 / 2.0;

    FGetBands2Result res;
    res.upper  = arma::mat(N, H, arma::fill::none);
    res.lower  = arma::mat(N, H, arma::fill::none);
    res.upper2 = arma::mat(N, H, arma::fill::none);
    res.lower2 = arma::mat(N, H, arma::fill::none);
    res.median = arma::mat(N, H, arma::fill::none);

    for (int i = 0; i < N; ++i) {
        for (int h = 0; h < H; ++h) {
            std::vector<double> sv(B);
            for (int b = 0; b < B; ++b) sv[b] = bootirf(i, h, b);
            std::sort(sv.begin(), sv.end());
            res.upper(i, h)  = lin_interp(sv, B, up1);
            res.lower(i, h)  = lin_interp(sv, B, lo1);
            res.upper2(i, h) = lin_interp(sv, B, up2);
            res.lower2(i, h) = lin_interp(sv, B, lo2);
            res.median(i, h) = lin_interp(sv, B, 50.0);
        }
    }
    return res;
}

//' Compute Confidence Bands from Bootstrap IRF Cube
//'
//' Computes upper, lower, and median bands from a 3-D bootstrap IRF array
//' along the third (bootstrap) dimension using linear interpolation of order
//' statistics.
//'
//' @param bootirf N x (hor+1) x nboot cube of bootstrapped IRFs.
//' @param prc Confidence level (e.g. 68 for 68\% band). Upper and lower
//'   quantiles are \eqn{(50 + prc/2)}\% and \eqn{(50 - prc/2)}\%.
//'   Default 68.
//'
//' @return List with three N x (hor+1) matrices: upper, lower, median.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fGetBands(const arma::cube& bootirf, double prc = 68.0) {
    FGetBandsResult r = fGetBands_cpp(bootirf, prc);
    return Rcpp::List::create(
        Rcpp::Named("upper")  = r.upper,
        Rcpp::Named("lower")  = r.lower,
        Rcpp::Named("median") = r.median
    );
}
