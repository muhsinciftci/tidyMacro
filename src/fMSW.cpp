// [[Rcpp::depends(RcppArmadillo)]]

#include "fMSW.h"
#include "flagmakerMatrix.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>

// =========================================================================
// Internal helpers
// =========================================================================

// Newey-West HAC covariance estimator (STATA convention, same as NW_hac_STATA.m)
// vars : T1 x d matrix (each row is an observation)
// lags : number of Newey-West lags (0 = heteroskedasticity-only correction)
// returns d x d long-run variance matrix
static arma::mat nw_hac_cpp(const arma::mat& vars, int lags) {
    int T = (int)vars.n_rows;
    arma::mat S = (vars.t() * vars) / T;
    for (int l = 1; l <= lags; ++l) {
        arma::mat Gl = (vars.head_rows(T - l).t() * vars.tail_rows(T - l)) / T;
        double w = 1.0 - (double)l / (lags + 1);
        S += w * (Gl + Gl.t());
    }
    return S;
}

// MA representation (identical recursion to MARep.m)
// AL   : n x n*p matrix [A1 | A2 | ... | Ap]
// p    : lag order
// hori : maximum horizon
// returns n x n*hori where block h (1-indexed, h=1..hori) = C_h
static arma::mat ma_rep_cpp(const arma::mat& AL, int p, int hori) {
    int n = (int)AL.n_rows;
    std::vector<arma::mat> C(hori + 1);
    C[0] = arma::eye(n, n);
    for (int h = 1; h <= hori; ++h) {
        C[h].zeros(n, n);
        for (int j = 1; j <= std::min(h, p); ++j) {
            C[h] += AL.cols((j - 1) * n, j * n - 1) * C[h - j];
        }
    }
    arma::mat result(n, n * hori);
    for (int h = 1; h <= hori; ++h) {
        result.cols((h - 1) * n, h * n - 1) = C[h];
    }
    return result;
}

// G matrices: derivative of vec(C_h) wrt vec(AL) (free VAR parameters)
// G_h = sum_{k=0}^{h-1} kron( (Alut^k * J')' , C_{h-1-k} )
// where Alut = companion matrix (np x np), J = [I_n | 0] (n x np)
// G  : n^2 x n^2*p x (hor+1) cube  (G.slice(0) = zeros, G.slice(h) for h>=1)
// Gcum: same shape, cumulative sum over slices
static void g_matrices_cpp(const arma::mat& AL,
                            const std::vector<arma::mat>& C_vec,  // C_0..C_hor
                            int p, int hor, int n,
                            arma::cube& G, arma::cube& Gcum) {

    // Companion matrix Alut (np x np)
    arma::mat Alut(n * p, n * p, arma::fill::zeros);
    Alut.rows(0, n - 1) = AL;
    if (p > 1) {
        Alut.submat(n, 0, n * p - 1, n * (p - 1) - 1) =
            arma::eye(n * (p - 1), n * (p - 1));
    }

    // J = [I_n | 0_{n,(p-1)n}]  (n x np)
    arma::mat J = arma::join_rows(arma::eye(n, n), arma::zeros(n, (p - 1) * n));

    // Precompute AJ_list[k] = (Alut^k * J')' = J * (Alut^k)'  (n x np) for k=0..hor-1
    // Note: (Alut^k * J')' = (J Alut'^k)  [transpose of np x n]
    std::vector<arma::mat> AJt(hor);   // AJt[k] = (Alut^k * J')' = n x np
    arma::mat Apower = arma::eye(n * p, n * p);
    for (int k = 0; k < hor; ++k) {
        AJt[k] = (Apower * J.t()).t();  // J * (Alut^k)' → n x np
        if (k < hor - 1) Apower = Apower * Alut;
    }

    int n2 = n * n;
    int n2p = n * n * p;

    G.zeros(n2, n2p, hor + 1);

    for (int h = 1; h <= hor; ++h) {
        arma::mat Gh(n2, n2p, arma::fill::zeros);
        // G_h = sum_{k=0}^{h-1} kron(AJt[k], C_{h-1-k})
        for (int k = 0; k < h; ++k) {
            Gh += arma::kron(AJt[k], C_vec[h - 1 - k]);
        }
        G.slice(h) = Gh;
    }

    // Cumulative G (manual cumsum along slice dimension)
    Gcum.zeros(n2, n2p, hor + 1);
    for (int h = 1; h <= hor; ++h) {
        Gcum.slice(h) = Gcum.slice(h - 1) + G.slice(h);
    }
}

// Vech selector matrix V such that vech(Sigma) = V * vec(Sigma)
// V: n*(n+1)/2 x n^2
static arma::mat vech_selector(int n) {
    int vech_len = n * (n + 1) / 2;
    arma::mat V(vech_len, n * n, arma::fill::zeros);
    arma::mat In = arma::eye(n, n);
    int row = 0;
    for (int i = 0; i < n; ++i) {
        arma::rowvec ei = In.row(i);
        arma::mat Isub = In.rows(i, n - 1);             // (n-i) x n
        arma::mat block = arma::kron(ei, Isub);          // (n-i) x n^2
        V.rows(row, row + (n - i) - 1) = block;
        row += n - i;
    }
    return V;
}

// =========================================================================
// CovAhat_Sigmahat_Gamma equivalent
// Returns WHatall, WHat, and stores them in output args
// WHatall: asymptotic variance of [vec(Ahat)', vech(Sigmahat)', vec(Gammahat)']'
// WHat   : asymptotic variance of [vec(Ahat)', vec(Gammahat)']'
// =========================================================================
static void cov_ahat_sigma_gamma(
        int p,
        const arma::mat& X,      // T1 x (m + n*p) design matrix
        const arma::mat& Z,      // T1 x k instrument
        const arma::mat& eta,    // n x T1 residuals (transposed)
        int NWlags,
        arma::mat& WHatall,
        arma::mat& WHat) {

    int n   = (int)eta.n_rows;
    int T1  = (int)eta.n_cols;
    int k   = (int)Z.n_cols;
    int m   = (int)X.n_cols - n * p;   // deterministic regressors (intercept count)

    int T2aux = (int)X.n_cols + n + k;  // m + n*p + n + k

    // Build matagg: T2aux x T1, each column = [x_t; eta_t; z_t]
    // x_t = X[t,:]', eta_t = eta[:,t], z_t = Z[t,:]'
    arma::mat matagg(T2aux, T1);
    for (int t = 0; t < T1; ++t) {
        matagg.col(t) = arma::join_cols(
            arma::join_cols(X.row(t).t(), eta.col(t)),
            Z.row(t).t());
    }

    // AuxHAC2: T1 x (n * T2aux)
    // Row t = vec(eta_t * matagg_t')' = vectorised outer product (row-vector)
    arma::mat AuxHAC2(T1, n * T2aux);
    for (int t = 0; t < T1; ++t) {
        arma::mat outer = eta.col(t) * matagg.col(t).t();  // n x T2aux
        AuxHAC2.row(t) = arma::vectorise(outer).t();
    }

    // Demean
    arma::rowvec mean_row = arma::mean(AuxHAC2, 0);
    AuxHAC2.each_row() -= mean_row;

    // NW HAC
    arma::mat WhatAss1 = nw_hac_cpp(AuxHAC2, NWlags);  // (n*T2aux) x (n*T2aux)

    // Vech selector V  (n*(n+1)/2 x n^2)
    arma::mat V = vech_selector(n);

    // Q1 = X'X / T1  ((m+np) x (m+np))
    arma::mat Q1 = (X.t() * X) / T1;
    arma::mat Q1inv = arma::inv_sympd(Q1);

    // Q2 = Z'X / T1  (k x (m+np))
    arma::mat Q2 = (Z.t() * X) / T1;

    // Shat: output rows = n^2*p + n*(n+1)/2 + n*k
    //       input cols  = n*m + n^2*p + n^2 + n*k = n*T2aux
    int n2p      = n * n * p;
    int vech_len = n * (n + 1) / 2;
    int nk       = n * k;
    int ncols    = n * T2aux;   // = n*m + n^2*p + n^2 + n*k
    int nrows    = n2p + vech_len + nk;

    arma::mat Shat(nrows, ncols, arma::fill::zeros);

    // --- Row block 1: vec(Ahat), rows 0..n2p-1 ---
    // Selector [0_{n*p,m} | I_{n*p}] is (n*p x (m+n*p)), rows m..m+np-1 of Q1inv
    arma::mat sel_lag = Q1inv.rows(m, m + n * p - 1);  // np x (m+np)
    // kron(sel_lag, I_n) is n^2*p x n*(m+np)
    arma::mat kron_A = arma::kron(sel_lag, arma::eye(n, n));
    // Place in cols 0..n*(m+np)-1
    Shat.submat(0, 0, n2p - 1, n * (m + n * p) - 1) = kron_A;

    // --- Row block 2: vech(Sigmahat), rows n2p..n2p+vech_len-1 ---
    // V acts on eta block: cols n*(m+np)..n*(m+np)+n^2-1
    int eta_col_start = n * ((int)X.n_cols);  // = n*(m+np)
    Shat.submat(n2p, eta_col_start, n2p + vech_len - 1, eta_col_start + n * n - 1) = V;

    // --- Row block 3: vec(Gammahat), rows n2p+vech_len..nrows-1 ---
    // -kron(Q2*Q1inv, I_n) acts on cols 0..n*(m+np)-1
    arma::mat kron_G = arma::kron(Q2 * Q1inv, arma::eye(n, n));  // nk x n*(m+np)
    Shat.submat(n2p + vech_len, 0, nrows - 1, n * (m + n * p) - 1) = -kron_G;
    // I_{nk} acts on Z block: cols n*(m+np)+n^2..ncols-1
    int z_col_start = eta_col_start + n * n;
    Shat.submat(n2p + vech_len, z_col_start, nrows - 1, ncols - 1) =
        arma::eye(nk, nk);

    // WHatall = Shat * WhatAss1 * Shat'
    WHatall = Shat * WhatAss1 * Shat.t();

    // WHat = extract rows/cols for [vec(A), vec(Gamma)] (drop vech(Sigma) block)
    // Rows/cols: 0..n2p-1 and n2p+vech_len..nrows-1
    arma::uvec idx(n2p + nk);
    for (int i = 0; i < n2p; ++i) idx(i) = i;
    for (int i = 0; i < nk;  ++i) idx(n2p + i) = n2p + vech_len + i;

    WHat = WHatall(idx, idx);
}

// =========================================================================
// fMSW_cpp  — main Montiel-Olea, Stock, Watson (2021) inference
// =========================================================================
MSWResult fMSW_cpp(const arma::mat& AL,
                   const arma::mat& Sigma_sub,
                   const arma::mat& eta,
                   const arma::mat& X,
                   const arma::mat& Z,
                   int p, int hor,
                   int nvar,
                   double scale,
                   double confidence,
                   int NWlags) {

    int n  = (int)AL.n_rows;
    int T1 = (int)eta.n_cols;
    int k  = (int)Z.n_cols;

    // Critical value (two-sided chi^2 with 1 df)
    // norminv(1-(1-conf)/2)^2 = qnorm(0.5+conf/2)^2
    // Use the approximation / exact value via R's qnorm via Rcpp
    double alpha_tail = (1.0 - confidence) / 2.0;
    // Beasley-Springer-Moro approximation for the normal quantile
    auto qnorm_approx = [](double p_val) -> double {
        // Rational approximation for tail quantiles
        const double a[] = {-3.969683028665376e+01,  2.209460984245205e+02,
                            -2.759285104469687e+02,  1.383577518672690e+02,
                            -3.066479806614716e+01,  2.506628277459239e+00};
        const double b[] = {-5.447609879822406e+01,  1.615858368580409e+02,
                            -1.556989798598866e+02,  6.680131188771972e+01,
                            -1.328068155288572e+01};
        const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                            -2.400758277161838e+00, -2.549732539343734e+00,
                             4.374664141464968e+00,  2.938163982698783e+00};
        const double d[] = { 7.784695709041462e-03,  3.224671290700398e-01,
                              2.445134137142996e+00,  3.754408661907416e+00};
        double q, r;
        if (p_val < 0.02425) {
            q = std::sqrt(-2.0 * std::log(p_val));
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                   ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
        } else if (p_val <= 0.97575) {
            q = p_val - 0.5;
            r = q * q;
            return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                   (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
        } else {
            q = std::sqrt(-2.0 * std::log(1.0 - p_val));
            return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
        }
    };
    double z_crit  = qnorm_approx(1.0 - alpha_tail);
    double critval = z_crit * z_crit;   // chi^2(1) critical value

    // Gamma = eta * Z / T1   (n x k)
    arma::mat Gamma = (eta * Z) / T1;

    // Compute WHat and WHatall
    arma::mat WHatall, WHat;
    cov_ahat_sigma_gamma(p, X, Z, eta, NWlags, WHatall, WHat);

    int n2p  = n * n * p;
    int nk   = n * k;

    // Sub-blocks of WHat (ordered as [vec(Ahat); vec(Gammahat)])
    arma::mat W1  = WHat.submat(0,    0,    n2p - 1,        n2p - 1);        // n^2*p x n^2*p
    arma::mat W12 = WHat.submat(0,    n2p,  n2p - 1,        n2p + nk - 1);   // n^2*p x n*k
    arma::mat W2  = WHat.submat(n2p,  n2p,  n2p + nk - 1,  n2p + nk - 1);   // n*k x n*k

    // Waldstat and (standard) F-stat for nvar-th element of Gamma (1-based)
    int nv0 = nvar - 1;   // 0-based index
    double Gamma_nvar = Gamma(nv0, 0);
    double Waldstat = ((double)T1 * Gamma_nvar * Gamma_nvar) /
                       WHat(n2p + nv0, n2p + nv0);
    // Standard F-stat (not HAC-robust): T1 * Gamma_nvar^2 / var_basic
    // Use Sigma_sub diagonal element as basic variance estimate
    double Fstat = Waldstat;  // report Waldstat as F (HAC robust)

    // MA representation: C_vec[h] = C_h for h=0..hor
    arma::mat MA_mat = ma_rep_cpp(AL, p, hor);  // n x n*hor
    std::vector<arma::mat> C_vec(hor + 1);
    C_vec[0] = arma::eye(n, n);
    for (int h = 1; h <= hor; ++h) {
        C_vec[h] = MA_mat.cols((h - 1) * n, h * n - 1);
    }

    // G and Gcum cubes: n^2 x n^2*p x (hor+1)
    arma::cube G, Gcum;
    g_matrices_cpp(AL, C_vec, p, hor, n, G, Gcum);

    // Allocate output matrices: n x (hor+1)
    arma::mat IRF(n, hor + 1, arma::fill::zeros);
    arma::mat IRFstderr(n, hor + 1, arma::fill::zeros);
    arma::mat Dmlbound(n, hor + 1, arma::fill::zeros);
    arma::mat Dmubound(n, hor + 1, arma::fill::zeros);
    arma::mat MSWlbound(n, hor + 1, arma::fill::zeros);
    arma::mat MSWubound(n, hor + 1, arma::fill::zeros);

    // ahat is the same for all (j, h)
    double ahat = (double)T1 * Gamma_nvar * Gamma_nvar - critval * W2(nv0, nv0);

    // Unit vectors
    arma::mat In = arma::eye(n, n);

    for (int j = 0; j < n; ++j) {
        arma::vec e_j = In.col(j);               // n x 1
        arma::vec e_nv = In.col(nv0);             // n x 1 (normalization variable)

        for (int h = 0; h <= hor; ++h) {
            const arma::mat& Ch = C_vec[h];       // n x n

            // 1 x n^2 row vector
            arma::rowvec kron_Gt_ej = arma::kron(Gamma.t(), e_j.t());

            // Plugin IRF: scale * e_j' * C_h * Gamma / Gamma_nvar
            double ej_Ch_Gamma = arma::as_scalar(e_j.t() * Ch * Gamma);
            double lambda_jh   = scale * ej_Ch_Gamma / Gamma_nvar;
            IRF(j, h) = lambda_jh;

            // Useful quantities
            arma::rowvec kGh   = kron_Gt_ej * G.slice(h);    // 1 x n^2*p
            arma::rowvec ej_Ch = e_j.t() * Ch;                // 1 x n

            arma::vec W2_col  = W2.col(nv0);                  // n*k x 1
            arma::vec W12_col = W12.col(nv0);                  // n^2*p x 1

            // Delta-method variance
            arma::vec d1 = (kGh * scale).t();                 // n^2*p x 1
            arma::vec d2 = (scale * ej_Ch - lambda_jh * e_nv.t()).t(); // n x 1
            arma::vec d  = arma::join_cols(d1, d2);            // (n^2*p + n) x 1
            double DmVar = arma::as_scalar(d.t() * WHat * d);

            double sqrt_T1    = std::sqrt((double)T1);
            double abs_Gnv    = std::abs(Gamma_nvar);
            IRFstderr(j, h)   = std::sqrt(DmVar) / (sqrt_T1 * abs_Gnv);
            double half_width = std::sqrt(critval / T1) * std::sqrt(DmVar) / abs_Gnv;
            Dmlbound(j, h)    = lambda_jh - half_width;
            Dmubound(j, h)    = lambda_jh + half_width;

            // MSW (Anderson-Rubin) bounds
            // bhat
            double b1 = -2.0 * T1 * scale * ej_Ch_Gamma * Gamma_nvar;
            double b2 =  2.0 * critval * scale * arma::as_scalar(kGh * W12_col);
            double b3 =  2.0 * critval * scale * arma::as_scalar(ej_Ch * W2_col);
            double bhat = b1 + b2 + b3;

            // chat
            double c1 = std::pow(sqrt_T1 * scale * ej_Ch_Gamma, 2.0);
            double c2 = critval * scale * scale *
                        arma::as_scalar(kGh * W1 * kGh.t());
            double c3 = 2.0 * critval * scale * scale *
                        arma::as_scalar(kGh * W12 * ej_Ch.t());
            double c4 = critval * scale * scale *
                        arma::as_scalar(ej_Ch * W2 * ej_Ch.t());
            double chat = c1 - c2 - c3 - c4;

            double Delta = bhat * bhat - 4.0 * ahat * chat;

            double lb, ub;
            if (ahat > 0 && Delta > 0) {
                lb = (-bhat - std::sqrt(Delta)) / (2.0 * ahat);
                ub = (-bhat + std::sqrt(Delta)) / (2.0 * ahat);
            } else if (ahat < 0 && Delta > 0) {
                lb = (-bhat + std::sqrt(Delta)) / (2.0 * ahat);
                ub = (-bhat - std::sqrt(Delta)) / (2.0 * ahat);
            } else if (ahat > 0 && Delta <= 0) {
                lb = arma::datum::nan;
                ub = arma::datum::nan;
            } else {
                lb = -arma::datum::inf;
                ub =  arma::datum::inf;
            }
            MSWlbound(j, h) = lb;
            MSWubound(j, h) = ub;
        }
    }

    // Force normalization variable at h=0 to exactly scale
    MSWlbound(nv0, 0) = scale;
    MSWubound(nv0, 0) = scale;

    MSWResult result;
    result.IRF       = IRF;
    result.IRFstderr = IRFstderr;
    result.Dmlbound  = Dmlbound;
    result.Dmubound  = Dmubound;
    result.MSWlbound = MSWlbound;
    result.MSWubound = MSWubound;
    result.Waldstat  = Waldstat;
    result.Fstat     = Fstat;
    return result;
}

// =========================================================================
// R wrapper
// =========================================================================

//' Weak-IV Robust Impulse Response Inference (Montiel-Olea, Stock, Watson 2021)
//'
//' Computes both plug-in (delta method) and Anderson-Rubin weak-IV robust
//' confidence sets for IV-identified structural IRFs following
//' Montiel-Olea, Stock and Watson (2021).
//'
//' @param var_result List from \code{fVAR}.
//' @param Z Instrument matrix (proxy sample, no missing values).
//' @param finaldata T x N matrix of original endogenous variables (full sample).
//' @param adjustu Integer vector \code{c(start, end)} (1-based) selecting the
//'   proxy-sample rows of the VAR residuals.
//' @param hor Forecast horizon. Default 48.
//' @param nvar 1-based index of the normalization variable. Default 1.
//' @param scale Scale of the structural shock. Default 1.
//' @param confidence Confidence level (0-1). Default 0.9.
//' @param NWlags Number of Newey-West lags for HAC correction: 0 = White
//'   heteroskedasticity-consistent only (no autocorrelation), k > 0 = Newey-West
//'   with k lags. Must be >= 0. Default 0.
//'
//' @return A list with elements:
//'   \item{IRF}{n x (hor+1) plug-in IRF point estimates.}
//'   \item{IRFstderr}{n x (hor+1) delta-method standard errors.}
//'   \item{Dmlbound}{n x (hor+1) delta-method lower confidence bounds.}
//'   \item{Dmubound}{n x (hor+1) delta-method upper confidence bounds.}
//'   \item{MSWlbound}{n x (hor+1) Anderson-Rubin lower bounds.}
//'   \item{MSWubound}{n x (hor+1) Anderson-Rubin upper bounds.}
//'   \item{Waldstat}{Wald statistic for instrument relevance (HAC-robust F-stat).}
//'   \item{Fstat}{Same as Waldstat (reported for compatibility).}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fMSW(const Rcpp::List& var_result,
                const arma::mat& Z,
                const arma::mat& finaldata,
                const arma::ivec& adjustu,
                int hor        = 48,
                int nvar       = 1,
                double scale   = 1.0,
                double confidence = 0.9,
                int NWlags     = 0) {

    if (NWlags < 0)
        Rcpp::stop("NWlags must be >= 0. Use 0 for White HC-only, or a positive integer for Newey-West with that many lags.");

    // Extract VAR components
    arma::mat beta      = Rcpp::as<arma::mat>(var_result["beta"]);
    arma::mat residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
    int p               = Rcpp::as<int>(var_result["p"]);
    int c               = Rcpp::as<int>(var_result["c"]);

    int N = (int)finaldata.n_cols;
    int T = (int)finaldata.n_rows;

    // Proxy-sample residuals (0-based)
    int u_start = adjustu(0) - 1;
    int u_end   = adjustu(1) - 1;
    arma::mat u = residuals.rows(u_start, u_end);  // T1 x N
    int T1 = (int)u.n_rows;

    // AL: n x n*p (strip constant row from beta, then transpose)
    // beta is (c + N*p) x N; first c rows are deterministic
    arma::mat AL = beta.rows(c, (int)beta.n_rows - 1).t();  // N x N*p

    // Sigma_sub = cov(u) with denominator T1-1 (standard sample covariance)
    arma::mat u_dm = u.each_row() - arma::mean(u, 0);
    arma::mat Sigma_sub = (u_dm.t() * u_dm) / (T1 - 1);

    // eta = u' (N x T1)
    arma::mat eta = u.t();

    // Build full design matrix for VAR sample, then restrict to proxy sample
    // X_full: (T-p) x (c + N*p)
    arma::mat lags_full = flagmakerMatrix(finaldata, p);   // (T-p) x N*p
    arma::mat X_full;
    if (c == 1) {
        X_full = arma::join_rows(arma::ones<arma::mat>(T - p, 1), lags_full);
    } else {
        X_full = lags_full;
    }
    // Restrict to proxy sample (same rows as residuals[u_start:u_end])
    arma::mat X = X_full.rows(u_start, u_end);  // T1 x (c + N*p)

    MSWResult res = fMSW_cpp(AL, Sigma_sub, eta, X, Z,
                             p, hor, nvar, scale, confidence, NWlags);

    return Rcpp::List::create(
        Rcpp::Named("IRF")       = res.IRF,
        Rcpp::Named("IRFstderr") = res.IRFstderr,
        Rcpp::Named("Dmlbound")  = res.Dmlbound,
        Rcpp::Named("Dmubound")  = res.Dmubound,
        Rcpp::Named("MSWlbound") = res.MSWlbound,
        Rcpp::Named("MSWubound") = res.MSWubound,
        Rcpp::Named("Waldstat")  = res.Waldstat,
        Rcpp::Named("Fstat")     = res.Fstat
    );
}
