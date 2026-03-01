#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Print SVAR Identification Steps
//'
//' Prints a numbered guide to the proxy-SVAR identification and
//' invertibility-testing procedure to the R console.
//'
//' @return Called for its side-effect of printing to the console.
//'   Returns \code{NULL} invisibly.
//'
//' @details
//' Step 1: Estimate the VAR(p) and extract the residuals.
//'
//' Step 2: Regress the instrument on the current residuals plus f leads;
//' save the fitted value.
//'
//' Step 3: F-test whether the lead coefficients are jointly zero. Rejection
//' implies the shock is non-invertible.
//'
//' Step 4: If invertible, apply the standard Stock-Watson IV procedure to
//' obtain IRFs and unit-variance shocks.
//'
//' Step 5: If non-invertible, test recoverability via the Ljung-Box test on
//' the fitted instrument series.
//'
//' Step 6: If recoverable (null of no autocorrelation not rejected), repeat
//' the IV procedure adding residual leads to the first-stage regression.
//'
//' Step 7: If recoverability is rejected, SVAR analysis with this instrument
//' is not feasible.
//'
//' @examples
//' fSVAR_steps()
//'
//' @export
// [[Rcpp::export]]
void fSVAR_steps() {
    Rcpp::Rcout << "Step 1: Estimate the VAR(p) and extract the residuals." << std::endl;
    Rcpp::Rcout << "Step 2: Regress the instrument on current residuals plus f leads. Save the fitted value." << std::endl;
    Rcpp::Rcout << "Step 3: F-test that lead coefficients are jointly zero. Rejection implies non-invertibility." << std::endl;
    Rcpp::Rcout << "Step 4: If invertible, apply the Stock-Watson IV procedure for IRFs and unit-variance shocks." << std::endl;
    Rcpp::Rcout << "Step 5: If non-invertible, test recoverability via the Ljung-Box test on the fitted instrument." << std::endl;
    Rcpp::Rcout << "Step 6: If recoverable, redo the IV procedure adding residual leads to the first-stage regression." << std::endl;
    Rcpp::Rcout << "Step 7: If recoverability is rejected, SVAR analysis with this instrument is not feasible." << std::endl;
}
