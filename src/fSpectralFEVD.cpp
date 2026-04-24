#include "fSpectralFEVD.h"
#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec fSpectralFEVD_cpp(const arma::cube& D, const arma::mat& irf_s,
                              const arma::vec& band, int J, bool fourier) {
    const int K = static_cast<int>(D.n_rows);
    const int H = static_cast<int>(D.n_slices);
    const int N = static_cast<int>(D.n_cols);

    // Build frequency grid
    // fourier = true:  pi/J, 2*pi/J, ..., pi  (J points, then filter to band)
    // fourier = false: linspace(2*pi/high, 2*pi/low, J) (all points averaged)
    arma::vec grid;
    double omega_l = 0.0, omega_u = 0.0;
    if (fourier) {
        const double step = arma::datum::pi / J;
        grid    = arma::regspace(step, step, arma::datum::pi);
        omega_l = 2.0 * arma::datum::pi / band(1);
        omega_u = 2.0 * arma::datum::pi / band(0);
    } else {
        grid = arma::linspace(2.0 * arma::datum::pi / band(1),
                              2.0 * arma::datum::pi / band(0), J);
    }

    arma::vec totvar(K, arma::fill::zeros);
    arma::vec shockvar(K, arma::fill::zeros);
    int count = 0;

    // Pre-allocate cos/sin arrays and real/imag work buffers (no per-j allocation)
    arma::vec c_arr(H), s_arr(H);
    arma::vec temp_re(N), temp_im(N);

    for (int j = 0; j < J; ++j) {
        const double omega = grid(j);
        if (fourier && (omega < omega_l || omega > omega_u)) continue;

        // Recursive trig: avoids H calls to std::exp per frequency
        // cos(omega*h) and sin(omega*h) via angle-addition recurrence
        const double c_step = std::cos(omega);
        const double s_step = std::sin(omega);
        c_arr(0) = 1.0;  s_arr(0) = 0.0;
        for (int h = 1; h < H; ++h) {
            c_arr(h) = c_arr(h-1) * c_step - s_arr(h-1) * s_step;
            s_arr(h) = s_arr(h-1) * c_step + c_arr(h-1) * s_step;
        }

        for (int k = 0; k < K; ++k) {
            // Accumulate real/imaginary parts of:
            //   sum_h D[k,:,h] * exp(-i*omega*h)   (exp(-i*x) = cos(x) - i*sin(x))
            // Avoids conv_to<cx_vec> — work entirely in real arithmetic
            temp_re.zeros();
            temp_im.zeros();
            double irf_re = 0.0, irf_im = 0.0;

            for (int h = 0; h < H; ++h) {
                const double ch =  c_arr(h);   // cos(omega*h)
                const double sh = -s_arr(h);   // -sin(omega*h)  [imag part of exp(-i*omega*h)]
                temp_re += ch * D.slice(h).row(k).t();
                temp_im += sh * D.slice(h).row(k).t();
                irf_re  += ch * irf_s(k, h);
                irf_im  += sh * irf_s(k, h);
            }

            totvar(k)   += arma::dot(temp_re, temp_re) + arma::dot(temp_im, temp_im);
            shockvar(k) += irf_re * irf_re + irf_im * irf_im;
        }
        ++count;
    }

    if (!fourier && count > 0) {
        totvar  /= count;
        shockvar /= count;
    }

    return shockvar / totvar;
}

//' Spectral Forecast Error Variance Decomposition
//'
//' Computes the fraction of spectral variance attributable to a single
//' structural shock within a given frequency band.
//'
//' @param D Cholesky IRF cube of dimension \eqn{N \times N \times H}
//'   (output of \code{fcholeskyIRF}).
//' @param irf_s Matrix of dimension \eqn{N \times H}: IRF of each variable to
//'   the structural shock of interest (unit-variance normalisation).
//' @param band Numeric vector of length 2: \code{c(low_period, high_period)}
//'   in the same time units as the data (e.g. months).
//'   Example: \code{c(18, 96)} for the business-cycle band.
//' @param J Integer. Number of points on the frequency grid.
//' @param fourier Logical. If \code{TRUE} (default), uses Fourier frequencies
//'   \eqn{(\pi/J,\, 2\pi/J,\, \ldots,\, \pi)} and sums over those falling
//'   inside \code{band}. If \code{FALSE}, places \code{J} points uniformly
//'   across the band via linspace and averages over all of them.
//'
//' @return Numeric vector of length \eqn{N}: share of spectral variance
//'   explained by the shock in the specified band.
//'
//' @seealso \code{\link{fcholeskyIRF}}, \code{\link{fPolyConvolve}}
//'
//' @export
// [[Rcpp::export]]
arma::vec fSpectralFEVD(const arma::cube& D, const arma::mat& irf_s,
                         const arma::vec& band, int J, bool fourier = true) {
    return fSpectralFEVD_cpp(D, irf_s, band, J, fourier);
}
