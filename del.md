🎉 tidyMacro is growing fast. After the initial release, the second version of the package comes with huge updates.

---

📐 **Identification via Long-Run Restrictions**
  1. Impulse response functions
  2. Bias-corrected impulse response functions
  3. Variance decomposition
  4. Replications:
     1. Galí, J. (1999). "Technology, Employment, and the Business Cycle: Do Technology Shocks Explain Aggregate Fluctuations?" American Economic Review, 89(1), 249–271
     2. Beaudry, P. & Portier, F. (2014). "News-Driven Business Cycles: Insights and Challenges." Journal of Economic Literature, 52(4), 993–1074

🎯 **Identification via External Instruments (Extended)**
  1. Historical decompositions
  2. Instrument Strength: First-stage F-statistics
  3. Recovering the structural shocks
     1. Unit normalization
     2. One SD normalization
  4. Weak IV Robust IRF following Montiel Olea, Stock & Watson (2021), "Inference in Structural Vector Autoregressions Identified with an External Instrument," Journal of Econometrics, 225(1), 74–87
     1. Delta method
     2. Anderson-Rubin
  5. External instrument SVAR analysis for noninvertible shocks following Forni, Gambetti & Ricco (2022), "External Instrument SVAR Analysis for Noninvertible Shocks"
     1. Invertibility test
     2. Recoverability test
     3. IRF, HD, FEVD, and relative IRFs

🔗 **Identification via Internal Instruments**
  1. Instrument enters VAR as first variable; identification via recursive ordering
  2. All short-run tools apply
  3. Replications:
     1. Känzig, D. R. (2021). "The Macroeconomic Effects of Oil Supply News: Evidence from OPEC Announcements." American Economic Review, 111(4), 1092–1125

⚡ **Identification via Heteroskedasticity**
  1. Impulse response functions with MBB bootstrap bands following Rigobon (2003), "Identification Through Heteroskedasticity," Review of Economics and Statistics, 85(4), 777–792
  2. Replications:
     1. Känzig, D. R. (2021). "The Macroeconomic Effects of Oil Supply News: Evidence from OPEC Announcements." American Economic Review, 111(4), 1092–1125

🚀 **Performance & Visuals**
  1. Bootstraps are now significantly faster
  2. Refined IRF charts — 68% and 90% confidence bands are now the default
  3. tidyMacro has its own website

---

📦 GitHub: github.com/muhsinciftci/tidyMacro
🌐 Website: https://tidymacro.netlify.app/
